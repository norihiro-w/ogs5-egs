/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_std.h"

#include <Eigen/Eigen>

#include "rf_mmp_new.h"
#include "rf_pcs.h"

namespace FiniteElement
{

void CFiniteElementStd::AssembleTHResidual()
{
	const unsigned c_dim = dim;
	const int c_nnodes = nnodes;
	const double dt = pcs->Tim->time_step_length;
	const double theta = pcs->m_num->ls_theta;
	const bool hasGravity =
		(coordinate_system) % 10 == 2 && FluidProp->CheckGravityCalculation();
	const double g_const = hasGravity ? gravity_constant : .0;
	const bool useSUPG = (pcs->m_num->ele_supg_method > 0);
	const bool useLumpedMass = (pcs->m_num->ele_mass_lumping > 0);
	const bool isTransient =
		(this->pcs->tim_type == FiniteElement::TIM_TRANSIENT);
	const bool isMatrixFlowInactive = pcs->deactivateMatrixFlow;
	const bool isMatrixElement = (MeshElement->GetDimension() == pcs->m_msh->GetMaxElementDim());

	Matrix t_transform_tensor(3, 3);
	if (dim > MediaProp->geo_dimension)
	{
		if (MeshElement->getTransformTensor() == NULL)
		{
			std::cout << "***Error: Geometric dimension in MMP is not "
						 "consistent with element."
					  << "\n";
			exit(0);
		}
		t_transform_tensor.LimitSize(dim, dim);
		MeshElement->getTransformTensor()->GetTranspose(t_transform_tensor);
	}

	Eigen::VectorXd nodal_p0(nnodes), nodal_p1(nnodes);
	Eigen::VectorXd nodal_T0(nnodes), nodal_T1(nnodes);
	for (int i = 0; i < nnodes; i++)
	{
		nodal_p0[i] = pcs->GetNodeValue(nodes[i], idxp0);
		nodal_p1[i] = pcs->GetNodeValue(nodes[i], idxp1);
		nodal_T0[i] = pcs->GetNodeValue(nodes[i], idxT0);
		nodal_T1[i] = pcs->GetNodeValue(nodes[i], idxT1);
	}

	// Calculate matrices
	const int offset_p = 0;
	const int offset_T = c_nnodes;
	(*RHS) = 0.0;          // Residual

//#define TH_DEBUG
#ifdef TH_DEBUG
	// for debugging
	Mass->LimitSize(nnodes, nnodes);
	Laplace->LimitSize(nnodes, nnodes);
	Advection->LimitSize(nnodes, nnodes);
	(*Mass) = .0;
	(*Laplace) = .0;
	(*Advection) = .0;
#endif
	double* tmp_k_tensor = MediaProp->PermeabilityTensor(Index);
	if (c_dim > MediaProp->geo_dimension)
	{
		Matrix local_tensor(c_dim, c_dim), temp_tensor(c_dim, c_dim),
			global_tensor(c_dim, c_dim);
		const unsigned c_ele_dim = ele_dim;
		for (size_t i = 0; i < c_ele_dim; i++)
			for (size_t j = 0; j < c_ele_dim; j++)
				local_tensor(i, j) = tmp_k_tensor[j + i * c_ele_dim];
		// cout << "K':" << endl; local_tensor.Write();
		local_tensor.multi(t_transform_tensor, temp_tensor);
		for (size_t i = 0; i < c_dim; i++)
			for (size_t j = 0; j < c_dim; j++)
				for (size_t k = 0; k < c_dim; k++)
					global_tensor(i, j) +=
						(*MeshElement->getTransformTensor())(i, k) *
						temp_tensor(k, j);
		// cout << "K:" << endl; global_tensor.Write();
		for (size_t i = 0; i < c_dim; i++)
			for (size_t j = 0; j < c_dim; j++)
				tmp_k_tensor[c_dim * i + j] = global_tensor(i, j);
	}
	double const* const k_tensor = tmp_k_tensor;
	const double Ss =
		MediaProp->StorageFunction(Index, unit, theta);
	double dummy[3] = {};
	double const* const lambda_tensor =
		MediaProp->HeatDispersionTensorNew(0, dummy);
	Eigen::MatrixXd k(c_dim, c_dim);
	for (unsigned i=0; i<c_dim; i++)
		for (unsigned j=0; j<c_dim; j++)
			k(i,j) = k_tensor[i*c_dim + j];
	Eigen::MatrixXd lambda(c_dim, c_dim);
	for (unsigned i=0; i<c_dim; i++)
		for (unsigned j=0; j<c_dim; j++)
			lambda(i,j) = lambda_tensor[i*c_dim + j];
	Eigen::Vector3d vec_g;
	vec_g << 0, 0, -g_const;

	//======================================================================
	Eigen::VectorXd r_p(c_nnodes), r_T(c_nnodes);
	r_p.setZero();
	r_T.setZero();

	//======================================================================
	// Mass lumping
	if (isTransient && useLumpedMass)
	{
		const double vol =
			MeshElement->GetVolume() * MeshElement->GetFluxArea();
		// Center of the reference element
		SetCenterGP();
		ComputeShapefct(1);
		double const* const c_shapefct = shapefct;
		Eigen::RowVectorXd N(c_nnodes);
		N.setZero();
		for (int i=0; i<c_nnodes; i++)
			N(i) = c_shapefct[i];
		const double gp_p1 = N * nodal_p1;
		const double gp_T1 = N * nodal_T1;
		double var[3] = {};
		var[0] = gp_p1;
		var[1] = gp_T1;
		const double rhocp = MediaProp->HeatCapacity(Index, theta, this, var);
		const double fkt = vol / (double)nnodes;
		r_T += rhocp * (nodal_T1 - nodal_T0)/dt * fkt;
	}

	//======================================================================
	// Loop over Gauss points
	int gp_r, gp_s, gp_t;
	double var[3] = {};
	Eigen::RowVectorXd N(c_nnodes);
	Eigen::RowVectorXd W_T(c_nnodes), W_SUPG(c_nnodes);
	Eigen::MatrixXd dN(c_dim, c_nnodes);
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		const double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeShapefct(1);
		ComputeGradShapefct(1);
		double const* const c_shapefct = shapefct;
		double const* const c_dshapefct = dshapefct;
		N.setZero();
		dN.setZero();
		for (int i=0; i<c_nnodes; i++)
			N(i) = c_shapefct[i];
		for (unsigned i = 0; i < c_dim; i++)
			for (int j = 0; j < c_nnodes; j++)
				dN(i,j) = c_dshapefct[i * c_nnodes + j];
		W_T = N;

		//---------------------------------------------------------
		//  Get state variables
		//---------------------------------------------------------
		const double gp_p0 = N * nodal_p0;
		const double gp_T0 = N * nodal_T0;
		const double gp_p1 = N * nodal_p1;
		const double gp_T1 = N * nodal_T1;
		Eigen::VectorXd grad_p1 = dN * nodal_p1;
		Eigen::VectorXd grad_T1 = dN * nodal_T1;

		//---------------------------------------------------------
		//  Get material properties
		//---------------------------------------------------------
		var[0] = gp_p1;
		var[1] = gp_T1;
		// Fluid properties
		const double rho_w = FluidProp->Density(var);
		const double vis = FluidProp->Viscosity(var);
		const double cp_w = FluidProp->SpecificHeatCapacity(var);
		// Medium properties
		const double rhocp = MediaProp->HeatCapacity(Index, theta, this, var);

		//---------------------------------------------------------
		//  Set velocity
		//---------------------------------------------------------
		Eigen::VectorXd grad_h1 = grad_p1;
		if (hasGravity)
			grad_h1 -= rho_w * vec_g;
		Eigen::VectorXd q = - k / vis * grad_h1;

		//---------------------------------------------------------
		//  SUPG coefficients
		//---------------------------------------------------------
		if (useSUPG)
		{
			W_SUPG.setZero();
			double supg_tau = 0;
			CalcSUPGWeightingFunction(&q[0], gp, supg_tau, &W_SUPG[0]);
			W_T += supg_tau * W_SUPG;
		}

		//---------------------------------------------------------
		//  Assemble Liquid flow equation
		//  original: N^T*S*N*dp/dt + dN^T*k/mu*dN*p+dN^T*rho*g*z = 0
		//---------------------------------------------------------
		// Rp += [1/dt*N^T*Ss*N + theta*dN^T*k/mu*dN]*p1 -
		// [1/dt*N^T*Ss*N - (1-theta)*dN^T*k/mu*dN]*p0 +
		// dN^T*k/mu*rho*g*z
		if (isTransient)
		{
			// 1/dt*N^T*Ss*(p1-p0)
			r_p.noalias() += N.transpose() * Ss * (gp_p1 - gp_p0)/dt * fkt;
		}
		// - dN^T*vel
		r_p.noalias() += - fkt * dN.transpose() * q;

		//---------------------------------------------------------
		//  Assemble Heat transport equation
		//---------------------------------------------------------

		// Rt += [1/dt*N^T*Cp*N + theta*(dN^T*lambda*dN+N^T*Cp_w*dN)]*T1 -
		// [1/dt*N^T*Cp*N - (1-theta)*(dN^T*lambda*dN+N^T*Cp_w*dN)]*T0
		if (isTransient && !useLumpedMass)
		{
			r_T.noalias() += W_T.transpose() * fkt * rhocp * (gp_T1 - gp_T0)/dt;
		}
		r_T.noalias() += fkt * dN.transpose() * lambda * grad_T1;
		r_T.noalias() += fkt * W_T.transpose() * rho_w * cp_w * q.transpose() * grad_T1;
	}

	if (isMatrixElement && isMatrixFlowInactive)
		r_p.setZero();

	for (int i = 0; i < c_nnodes; i++)
	{
		(*RHS)(offset_p + i) = r_p[i];
		(*RHS)(offset_T + i) = r_T[i];
	}

	if (pcs->scaleEQS)
	{
		for (int ii = 0; ii < 2; ii++)
		{
			const double scale_eqs = pcs->vec_scale_eqs[ii];
			for (int i = 0; i < c_nnodes; i++)
				(*RHS)(ii* c_nnodes + i) *= scale_eqs;
		}
	}

	// RHS should be - residual
	(*RHS) *= -1.;

#ifdef NEW_EQS
	for (size_t ii = 0; ii < pcs->GetPrimaryVNumber(); ii++)
		for (long i = 0; i < c_nnodes; i++)
			eqs_rhs[NodeShift[ii] + eqs_number[i]] += (*RHS)(i + ii * c_nnodes);
#endif
}

void CFiniteElementStd::AssembleTHJacobian()
{
	const unsigned c_dim = dim;
	const int c_nnodes = nnodes;
	const double dt = pcs->Tim->time_step_length;
	const double theta = pcs->m_num->ls_theta;
	const bool hasGravity =
		(coordinate_system) % 10 == 2 && FluidProp->CheckGravityCalculation();
	const double g_const = hasGravity ? gravity_constant : .0;
	const bool useSUPG = (pcs->m_num->ele_supg_method > 0);
	const bool useLumpedMass = (pcs->m_num->ele_mass_lumping > 0);
	const bool isTransient =
		(this->pcs->tim_type == FiniteElement::TIM_TRANSIENT);
	const bool isMatrixFlowInactive = pcs->deactivateMatrixFlow;
	const bool isMatrixElement = (MeshElement->GetDimension() == pcs->m_msh->GetMaxElementDim());

	Matrix t_transform_tensor(3, 3);
	if (dim > MediaProp->geo_dimension)
	{
		if (MeshElement->getTransformTensor() == NULL)
		{
			std::cout << "***Error: Geometric dimension in MMP is not "
						 "consistent with element."
					  << "\n";
			exit(0);
		}
		t_transform_tensor.LimitSize(dim, dim);
		MeshElement->getTransformTensor()->GetTranspose(t_transform_tensor);
	}

	Eigen::VectorXd nodal_p0(nnodes), nodal_p1(nnodes);
	Eigen::VectorXd nodal_T0(nnodes), nodal_T1(nnodes);
	for (int i = 0; i < nnodes; i++)
	{
		nodal_p0[i] = pcs->GetNodeValue(nodes[i], idxp0);
		nodal_p1[i] = pcs->GetNodeValue(nodes[i], idxp1);
		nodal_T0[i] = pcs->GetNodeValue(nodes[i], idxT0);
		nodal_T1[i] = pcs->GetNodeValue(nodes[i], idxT1);
	}

	// Calculate matrices
	const int offset_p = 0;
	const int offset_T = c_nnodes;
	(*StiffMatrix) = 0.0;  // Jacobian

//#define TH_DEBUG
#ifdef TH_DEBUG
	// for debugging
	Mass->LimitSize(nnodes, nnodes);
	Laplace->LimitSize(nnodes, nnodes);
	Advection->LimitSize(nnodes, nnodes);
	(*Mass) = .0;
	(*Laplace) = .0;
	(*Advection) = .0;
#endif
	double* tmp_k_tensor = MediaProp->PermeabilityTensor(Index);
	if (c_dim > MediaProp->geo_dimension)
	{
		Matrix local_tensor(c_dim, c_dim), temp_tensor(c_dim, c_dim),
			global_tensor(c_dim, c_dim);
		const unsigned c_ele_dim = ele_dim;
		for (size_t i = 0; i < c_ele_dim; i++)
			for (size_t j = 0; j < c_ele_dim; j++)
				local_tensor(i, j) = tmp_k_tensor[j + i * c_ele_dim];
		// cout << "K':" << endl; local_tensor.Write();
		local_tensor.multi(t_transform_tensor, temp_tensor);
		for (size_t i = 0; i < c_dim; i++)
			for (size_t j = 0; j < c_dim; j++)
				for (size_t k = 0; k < c_dim; k++)
					global_tensor(i, j) +=
						(*MeshElement->getTransformTensor())(i, k) *
						temp_tensor(k, j);
		// cout << "K:" << endl; global_tensor.Write();
		for (size_t i = 0; i < c_dim; i++)
			for (size_t j = 0; j < c_dim; j++)
				tmp_k_tensor[c_dim * i + j] = global_tensor(i, j);
	}
	double const* const k_tensor = tmp_k_tensor;
	const double Ss = MediaProp->StorageFunction(Index, unit, theta);
	double dummy[3] = {};
	double const* const lambda_tensor =
		MediaProp->HeatDispersionTensorNew(0, dummy);
	Eigen::MatrixXd k(c_dim, c_dim);
	for (unsigned i=0; i<c_dim; i++)
		for (unsigned j=0; j<c_dim; j++)
			k(i,j) = k_tensor[i*c_dim + j];
	Eigen::MatrixXd lambda(c_dim, c_dim);
	for (unsigned i=0; i<c_dim; i++)
		for (unsigned j=0; j<c_dim; j++)
			lambda(i,j) = lambda_tensor[i*c_dim + j];
	Eigen::Vector3d vec_g;
	vec_g << 0, 0, -g_const;

	//======================================================================
	Eigen::MatrixXd J_pp(c_nnodes, c_nnodes), J_pT(c_nnodes, c_nnodes);
	Eigen::MatrixXd J_Tp(c_nnodes, c_nnodes), J_TT(c_nnodes, c_nnodes);
	J_pp.setZero();
	J_pT.setZero();
	J_Tp.setZero();
	J_TT.setZero();

	//======================================================================
	// Mass lumping
	if (isTransient && useLumpedMass)
	{
		const double vol =
			MeshElement->GetVolume() * MeshElement->GetFluxArea();
		// Center of the reference element
		SetCenterGP();
		ComputeShapefct(1);
		double const* const c_shapefct = shapefct;
		Eigen::RowVectorXd N(c_nnodes);
		N.setZero();
		for (int i=0; i<c_nnodes; i++)
			N(i) = c_shapefct[i];
		const double gp_p1 = N * nodal_p1;
		const double gp_T1 = N * nodal_T1;
		double var[3] = {};
		var[0] = gp_p1;
		var[1] = gp_T1;
		const double rhocp = MediaProp->HeatCapacity(Index, theta, this, var);
		const double drho_w_dp = FluidProp->drhodP(var);
		const double drho_w_dT = FluidProp->drhodT(var);
		const double cp_w = FluidProp->SpecificHeatCapacity(var);
		const double porosity = MediaProp->Porosity(Index, theta);
		const double drhocp_dp = porosity * cp_w * drho_w_dp;
		const double drhocp_dT = porosity * cp_w * drho_w_dT;
		const double fkt = vol / (double)nnodes;
		J_TT.diagonal().setConstant(1/dt * rhocp * fkt);
		J_TT += fkt * drhocp_dT * (nodal_T1 - nodal_T0)/dt * N;
		J_Tp += fkt * drhocp_dp * (nodal_T1 - nodal_T0)/dt * N;
	}

	//======================================================================
	// Loop over Gauss points
	int gp_r, gp_s, gp_t;
	double var[3] = {};
	Eigen::RowVectorXd N(c_nnodes);
	Eigen::MatrixXd dN(c_dim, c_nnodes);
	Eigen::RowVectorXd W_T(c_nnodes), W_SUPG(c_nnodes);
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		const double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeShapefct(1);
		ComputeGradShapefct(1);
		double const* const c_shapefct = shapefct;
		double const* const c_dshapefct = dshapefct;
		N.setZero();
		dN.setZero();
		for (int i=0; i<c_nnodes; i++)
			N(i) = c_shapefct[i];
		for (unsigned i = 0; i < c_dim; i++)
			for (int j = 0; j < c_nnodes; j++)
				dN(i,j) = c_dshapefct[i * c_nnodes + j];
		W_T = N;

		//---------------------------------------------------------
		//  Get state variables
		//---------------------------------------------------------
		// const double gp_p0 = N * nodal_p0;
		const double gp_T0 = N * nodal_T0;
		const double gp_p1 = N * nodal_p1;
		const double gp_T1 = N * nodal_T1;
		//const double gp_dp = gp_p1 - gp_p0;
		const double gp_dT = gp_T1 - gp_T0;
		//Eigen::VectorXd grad_p0 = dN * nodal_p0;
		Eigen::VectorXd grad_p1 = dN * nodal_p1;
		//Eigen::VectorXd grad_T0 = dN * nodal_T0;
		Eigen::VectorXd grad_T1 = dN * nodal_T1;
		if (pcs->m_num->nls_jacobian_level == 1)
			grad_T1.setZero();

		//---------------------------------------------------------
		//  Get material properties
		//---------------------------------------------------------
		var[0] = gp_p1;
		var[1] = gp_T1;
		// Fluid properties
		const double rho_w = FluidProp->Density(var);
		const double drho_w_dp = FluidProp->drhodP(var);
		const double drho_w_dT = FluidProp->drhodT(var);
		const double vis = FluidProp->Viscosity(var);
		const double dvis_dp = FluidProp->dViscositydP(var);
		const double dvis_dT = FluidProp->dViscositydT(var);
		const double cp_w = FluidProp->SpecificHeatCapacity(var);
		// Medium properties
		const double rhocp = MediaProp->HeatCapacity(Index, theta, this, var);
		const double porosity = MediaProp->Porosity(Index, theta);
		const double drhocp_dp = porosity * cp_w * drho_w_dp; //TODO d(cp)/dp
		const double drhocp_dT = porosity * cp_w * drho_w_dT;

		//---------------------------------------------------------
		//  Set velocity
		//---------------------------------------------------------
		Eigen::VectorXd grad_h1 = grad_p1;
		if (hasGravity)
			grad_h1 -= rho_w * vec_g;
		Eigen::VectorXd q = - k / vis * grad_h1;

		//-----------------------------------------
		// Derivatives of flow velocity and heat flux
		//-----------------------------------------
		Eigen::MatrixXd dq_dp = - k / vis * dN;
		if (hasGravity && drho_w_dp != 0.0)
			dq_dp.noalias() += - k / vis * drho_w_dp * vec_g * N;
		if (dvis_dp != 0.0)
			dq_dp.noalias() += - dvis_dp / vis * q * N;
		Eigen::MatrixXd dq_dT = - dvis_dT * q * N;
		if (hasGravity && drho_w_dT != 0.0)
			dq_dT.noalias() += - k / vis * drho_w_dT * vec_g * N;

		Eigen::MatrixXd djDiff_dT = - lambda * dN;
		Eigen::MatrixXd djAdv_dT = rhocp * q.transpose() * dN;
		if (dq_dT.size() > 0)
			djAdv_dT.noalias() += rhocp * grad_T1.transpose() * dq_dT;
		if (drho_w_dT != .0)
			djAdv_dT.noalias() += drhocp_dT * q.transpose() * grad_T1 * N;
		Eigen::MatrixXd djAdv_dp = rhocp * grad_T1.transpose() * dq_dp;
		if (drho_w_dp != .0)
			djAdv_dp.noalias() += drhocp_dp * q.transpose() * grad_T1 * N;

		//---------------------------------------------------------
		//  SUPG coefficients
		//---------------------------------------------------------
		if (useSUPG)
		{
			W_SUPG.setZero();
			double supg_tau = 0;
			CalcSUPGWeightingFunction(&q[0], gp, supg_tau, &W_SUPG[0]);
			W_T += supg_tau * W_SUPG;
		}

		//---------------------------------------------------------
		//  Assemble Jacobian
		//---------------------------------------------------------

		if (isTransient)
			J_pp.noalias() += 1/dt * fkt * N.transpose() * Ss * N;
		J_pp.noalias() += - fkt * dN.transpose() * dq_dp;
		if (dq_dT.size() > 0)
			J_pT.noalias() += - fkt * dN.transpose() * dq_dT;

		if (isTransient && !useLumpedMass)
			J_TT.noalias() += 1/dt * W_T.transpose() * (rhocp + gp_dT * drhocp_dT) * N * fkt;
		J_TT.noalias() += - dN.transpose() * djDiff_dT * fkt;
		J_TT.noalias() += W_T.transpose() * djAdv_dT * fkt;
		if (isTransient && !useLumpedMass)
			J_Tp.noalias() += 1/dt * W_T.transpose() * (gp_dT * drhocp_dp) * N * fkt;
		J_Tp.noalias() += W_T.transpose() * djAdv_dp * fkt;
	}

	if (isMatrixElement && isMatrixFlowInactive)
	{
		J_pp.setZero();
		J_pT.setZero();
		J_Tp.setZero();
	}

	for (int i=0; i<c_nnodes; i++)
	{
		for (int j=0; j<c_nnodes; j++)
		{
			(*StiffMatrix)(offset_p + i, offset_p + j) = J_pp(i, j);
			(*StiffMatrix)(offset_p + i, offset_T + j) = J_pT(i, j);
			(*StiffMatrix)(offset_T + i, offset_p + j) = J_Tp(i, j);
			(*StiffMatrix)(offset_T + i, offset_T + j) = J_TT(i, j);
		}
	}

	if (pcs->scaleUnknowns || pcs->scaleEQS)
	{
		// scale p by 1e6
		const double scale_p =
			1. / pcs->vec_scale_dofs[0] * pcs->vec_scale_eqs[0];
		for (long i = 0; i < c_nnodes; i++)
			for (long j = 0; j < c_nnodes; j++)
			{
				(*StiffMatrix)(i, j) *= scale_p;
				(*StiffMatrix)(i, j + c_nnodes) *= scale_p;
				//				(*StiffMatrix)(i+c_nnodes,j) *= scale_p;
			}
		// scale T
		const double scale_T =
			1. / pcs->vec_scale_dofs[1] * pcs->vec_scale_eqs[1];
		for (long i = 0; i < c_nnodes; i++)
			for (long j = 0; j < c_nnodes; j++)
			{
				(*StiffMatrix)(i + c_nnodes, j) *= scale_T;
				(*StiffMatrix)(i + c_nnodes, j + c_nnodes) *= scale_T;
				//				(*StiffMatrix)(i,j+c_nnodes) *= scale_T;
				//				(*StiffMatrix)(i+c_nnodes,j+c_nnodes) *=
				// scale_T;
			}
	}
}

}  // end namespace

