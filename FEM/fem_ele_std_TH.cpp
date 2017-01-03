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
	// Loop over Gauss points
	int gp_r, gp_s, gp_t;
	double var[3] = {};
	Eigen::RowVectorXd N(c_nnodes);
	Eigen::MatrixXd dN(c_dim, c_nnodes);
	Eigen::VectorXd r_p(c_nnodes), r_T(c_nnodes);
	r_p.setZero();
	r_T.setZero();
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
		//  Assemble Liquid flow equation
		//  original: N^T*S*N*dp/dt + dN^T*k/mu*dN*p+dN^T*rho*g*z = 0
		//---------------------------------------------------------
		// Rp += [1/dt*N^T*Ss*N + theta*dN^T*k/mu*dN]*p1 -
		// [1/dt*N^T*Ss*N - (1-theta)*dN^T*k/mu*dN]*p0 +
		// dN^T*k/mu*rho*g*z
		if (isTransient)
		{
			// 1/dt*N^T*Ss*(p1-p0)
			r_p += N.transpose() * Ss * (gp_p1 - gp_p0)/dt * fkt;
		}
		// - dN^T*vel
		r_p += - fkt * dN.transpose() * q;
		//---------------------------------------------------------
		//  Assemble Heat transport equation
		//---------------------------------------------------------

		// Rt += [1/dt*N^T*Cp*N + theta*(dN^T*lambda*dN+N^T*Cp_w*dN)]*T1 -
		// [1/dt*N^T*Cp*N - (1-theta)*(dN^T*lambda*dN+N^T*Cp_w*dN)]*T0
		if (isTransient)
		{
			r_T += N.transpose() * fkt * rhocp * (gp_T1 - gp_T0)/dt;
		}
		r_T += fkt * dN.transpose() * lambda * grad_T1;
		r_T += fkt * N.transpose() * rho_w * cp_w * q.transpose() * grad_T1;
	}

	if (isMatrixElement && isMatrixFlowInactive)
		r_p.setZero();

	for (int i = 0; i < c_nnodes; i++)
	{
		(*RHS)(offset_p + i) = r_p[i];
		(*RHS)(offset_T + i) = r_T[i];
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
	// Loop over Gauss points
	int gp_r, gp_s, gp_t;
	double var[3] = {};
	Eigen::RowVectorXd N(c_nnodes);
	Eigen::MatrixXd dN(c_dim, c_nnodes);
	Eigen::MatrixXd J_pp(c_nnodes, c_nnodes), J_pT(c_nnodes, c_nnodes);
	Eigen::MatrixXd J_Tp(c_nnodes, c_nnodes), J_TT(c_nnodes, c_nnodes);
	J_pp.setZero();
	J_pT.setZero();
	J_Tp.setZero();
	J_TT.setZero();
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
		//  Assemble Jacobian
		//---------------------------------------------------------

		if (isTransient)
			J_pp += 1/dt * fkt * N.transpose() * Ss * N;
		J_pp += - fkt * dN.transpose() * dq_dp;
		if (dq_dT.size() > 0)
			J_pT.noalias() += - fkt * dN.transpose() * dq_dT;

		if (isTransient)
			J_TT.noalias() += 1/dt * N.transpose() * (rhocp + gp_dT * drhocp_dT) * N * fkt;
		J_TT.noalias() += - dN.transpose() * djDiff_dT * fkt;
		J_TT.noalias() += N.transpose() * djAdv_dT * fkt;
		if (isTransient)
			J_Tp.noalias() += 1/dt * N.transpose() * (gp_dT * drhocp_dp) * N * fkt;
		J_Tp.noalias() += N.transpose() * djAdv_dp * fkt;
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
}

#if 0
void CFiniteElementStd::AssembleTHEquation(bool updateA, bool updateRHS)
{
	const unsigned c_dim = dim;
	const int c_nnodes = nnodes;
	const double pcs_time_step = pcs->Tim->time_step_length;
	const double dt_inverse = 1.0 / pcs_time_step;
	const double theta = pcs->m_num->ls_theta;
	const bool hasGravity =
	    (coordinate_system) % 10 == 2 && FluidProp->CheckGravityCalculation();
	const double g_const = hasGravity ? gravity_constant : .0;
	const bool useSUPG = (pcs->m_num->ele_supg_method > 0);
	const bool useLumpedMass = (pcs->m_num->ele_mass_lumping > 0);
	const bool isTransient =
	    (this->pcs->tim_type == FiniteElement::TIM_TRANSIENT);
	//	ElementValue* gp_ele = ele_gp_value[Index];
	static Matrix t_transform_tensor(3, 3);
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

	for (int i = 0; i < nnodes; i++)
	{
		NodalVal_p0[i] = pcs->GetNodeValue(nodes[i], idxp0);
		NodalVal_p1[i] = pcs->GetNodeValue(nodes[i], idxp1);
		NodalVal_T0[i] = pcs->GetNodeValue(nodes[i], idxT0);
		NodalVal_T1[i] = pcs->GetNodeValue(nodes[i], idxT1);
	}

	// Calculate matrices
	const int offset_p = 0;
	const int offset_T = c_nnodes;
	(*StiffMatrix) = 0.0;  // Jacobian
	(*RHS) = 0.0;          // Residual
                           //	if (myrank==2 && updateA) {
                           //		StiffMatrix->Write();
                           //	}

#define TH_DEBUG
#ifdef TH_DEBUG
	// for debugging
	Mass->LimitSize(nnodes, nnodes);
	Laplace->LimitSize(nnodes, nnodes);
	Advection->LimitSize(nnodes, nnodes);
	(*Mass) = .0;
	(*Laplace) = .0;
	(*Advection) = .0;
#endif
#define ELE_CONST_K
#ifdef ELE_CONST_K
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
	    MediaProp->StorageFunction(Index, unit, pcs->m_num->ls_theta);
	double dummy[3] = {};
	double const* const lambda_tensor =
	    MediaProp->HeatDispersionTensorNew(0, dummy);
#endif

#define USE_LMASS1
#ifdef USE_LMASS1
	if (isTransient && useLumpedMass)
	{
		const double vol =
		    MeshElement->GetVolume() * MeshElement->GetFluxArea();
		// Center of the reference element
		SetCenterGP();
		ComputeShapefct(1);
		const double gp_p1 = interpolate(NodalVal_p1);
		const double gp_T1 = interpolate(NodalVal_T1);
		double var[3] = {};
		var[0] = gp_p1;
		var[1] = gp_T1;
		const double rhocp =
		    MediaProp->HeatCapacity(Index, pcs->m_num->ls_theta, this, var);
		const double drho_w_dp = FluidProp->drhodP(var);
		const double drho_w_dT = FluidProp->drhodT(var);
		const double cp_w = FluidProp->SpecificHeatCapacity(var);
		const double porosity =
		    MediaProp->Porosity(Index, pcs->m_num->ls_theta);
		const double drhocp_dp = porosity * cp_w * drho_w_dp;
		const double drhocp_dT = porosity * cp_w * drho_w_dT;
		const double fkt = vol / (double)nnodes;
		const double M_tt_coeff1 = fkt * dt_inverse * rhocp;
		const double M_tt_coeff2 = fkt * dt_inverse * drhocp_dT;
		const double M_tp_coeff2 = fkt * dt_inverse * drhocp_dp;
		if (updateA)
			for (int i = 0; i < nnodes; i++)
			{
				(*StiffMatrix)(offset_T + i, offset_T + i) += M_tt_coeff1;
				for (int j = 0; j < nnodes; j++)
				{
					(*StiffMatrix)(offset_T + i, offset_T + j) +=
					    M_tt_coeff2 * (NodalVal_T1[i] - NodalVal_T0[i]) *
					    shapefct[j];
					(*StiffMatrix)(offset_T + i, offset_p + j) +=
					    M_tp_coeff2 * (NodalVal_T1[i] - NodalVal_T0[i]) *
					    shapefct[j];
				}
			}
		if (updateRHS)
			for (int i = 0; i < nnodes; i++)
			{
				(*RHS)(offset_T + i) +=
				    M_tt_coeff1 * (NodalVal_T1[i] - NodalVal_T0[i]);
			}
	}
#endif

	//======================================================================
	// Loop over Gauss points
	int gp_r, gp_s, gp_t;
	double var[3] = {}, vel[3] = {};
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

		//---------------------------------------------------------
		//  Get state variables
		//---------------------------------------------------------
		const double gp_p0 = interpolate(NodalVal_p0);
		const double gp_T0 = interpolate(NodalVal_T0);
		const double gp_p1 = interpolate(NodalVal_p1);
		const double gp_T1 = interpolate(NodalVal_T1);
		//		if (myrank==2 && gp==0 && updateA)
		//		ScreenMessage2("%d: p0=%.3e, p1=%.3e, T0=%.3e, T1=%.3e\n",
		//Index,
		// gp_p0, gp_p1, gp_T0, gp_T1);
		var[0] = gp_p1;
		var[1] = gp_T1;
		double grad_p0[3] = {}, grad_T0[3] = {};
		double grad_p1[3] = {}, grad_T1[3] = {};
		double grad_T1_r[3] = {};
		//		static double max_gradT1 = .0;
		for (unsigned i = 0; i < c_dim; i++)
		{
			for (int j = 0; j < c_nnodes; j++)
			{
				grad_p0[i] += c_dshapefct[i * nnodes + j] * NodalVal_p0[j];
				grad_p1[i] += c_dshapefct[i * nnodes + j] * NodalVal_p1[j];
				grad_T0[i] += c_dshapefct[i * nnodes + j] * NodalVal_T0[j];
				grad_T1_r[i] += c_dshapefct[i * nnodes + j] * NodalVal_T1[j];
			}
#define TH_INCLUDE_GRAD_T
#ifdef TH_INCLUDE_GRAD_T
			if (pcs->m_num->nls_jacobian_level == 0) grad_T1[i] = grad_T1_r[i];
#ifdef USE_PETSC
//			if (fabs(grad_T1[i]) > fabs(max_gradT1) ) {
//				max_gradT1 = grad_T1[i];
////				ScreenMessage2("-> max grad T updated = %g\n", max_gradT1);
//			}
#endif
#endif
		}

//---------------------------------------------------------
//  Get material properties
//---------------------------------------------------------
#ifdef CONST_HEAT
		{
#endif
			// Fluid properties
			const double rho_w = FluidProp->Density(var);
			const double drho_w_dp = FluidProp->drhodP(var);
			const double drho_w_dT = FluidProp->drhodT(var);
			const double vis = FluidProp->Viscosity(var);
			const double dvis_dp = FluidProp->dViscositydP(var);
			const double dvis_dT = FluidProp->dViscositydT(var);
			const double cp_w = FluidProp->SpecificHeatCapacity(var);
// Medium properties
#ifndef ELE_CONST_K
			const double Ss =
			    MediaProp->StorageFunction(Index, unit, pcs->m_num->ls_theta);
			double* k_tensor = MediaProp->PermeabilityTensor(Index);
			if (dim > MediaProp->geo_dimension)
			{
				Matrix local_tensor(dim, dim), temp_tensor(dim, dim),
				    global_tensor(dim, dim);
				for (size_t i = 0; i < ele_dim; i++)
					for (size_t j = 0; j < ele_dim; j++)
						local_tensor(i, j) = k_tensor[j + i * ele_dim];
				// cout << "K':" << endl; local_tensor.Write();
				local_tensor.multi(t_transform_tensor, temp_tensor);
				for (size_t i = 0; i < dim; i++)
					for (size_t j = 0; j < dim; j++)
						for (size_t k = 0; k < dim; k++)
							global_tensor(i, j) +=
							    (*MeshElement->transform_tensor)(i, k) *
							    temp_tensor(k, j);
				// cout << "K:" << endl; global_tensor.Write();
				for (size_t i = 0; i < dim; i++)
					for (size_t j = 0; j < dim; j++)
						k_tensor[dim * i + j] = global_tensor(i, j);
			}
			double const* const lambda_tensor =
			    MediaProp->HeatDispersionTensorNew(gp, var);
#endif
			const double rhocp =
			    MediaProp->HeatCapacity(Index, pcs->m_num->ls_theta, this, var);
			const double porosity =
			    MediaProp->Porosity(Index, pcs->m_num->ls_theta);
			const double drhocp_dp =
			    porosity * cp_w *
			    drho_w_dp;  //  MediaProp->dHeatCapacitydP(Index,pcs->m_num->ls_theta,this,var);
			const double drhocp_dT =
			    porosity * cp_w *
			    drho_w_dT;  // MediaProp->dHeatCapacitydT(Index,pcs->m_num->ls_theta,this,var);

//---------------------------------------------------------
//  Set velocity
//---------------------------------------------------------
#if 0
		vel[0] = gp_ele->Velocity(0, gp);
		vel[1] = gp_ele->Velocity(1, gp);
		vel[2] = gp_ele->Velocity(2, gp);
		if (myrank==0 && gp==0 && updateRHS)
		ScreenMessage2("%d: vx=%.3e, vy=%.3e, vz=%.3e\n", Index, vel[0], vel[1], vel[2]);
#else

#if 1
		for (unsigned i = 0; i < c_dim; i++)
		{
			vel[i] = .0;
			for (unsigned j = 0; j < c_dim; j++)
				vel[i] -= k_tensor[c_dim * i + j] / vis * grad_p1[j];
			vel[i] -= k_tensor[c_dim * i + c_dim - 1] / vis * rho_w * g_const;
		}
#endif
//		if (gp==0 && Index<2)
//		ScreenMessage2("%d: vx=%.3e, vy=%.3e, vz=%.3e\n", Index, vel[0], vel[1],
// vel[2]);
//		if (gp==0 && Index<2)
//		ScreenMessage2("%d: gx=%.3e, gy=%.3e, gz=%.3e\n", Index,
// gp_ele->Velocity(0, gp), gp_ele->Velocity(1, gp), gp_ele->Velocity(2, gp));
//		if (c_dim > MediaProp->geo_dimension) {
//			if (gp==0 && Index%100==0)
//				ScreenMessage2("%d: vx=%.3e, vy=%.3e, vz=%.3e\n", Index, vel[0],
// vel[1], vel[2]);
//			if (gp==0 && Index%100==0)
//				ScreenMessage2("%d: gx=%.3e, gy=%.3e, gz=%.3e\n", Index,
// gp_ele->Velocity(0, gp), gp_ele->Velocity(1, gp), gp_ele->Velocity(2, gp));
//		}
#endif

//---------------------------------------------------------
//  Assemble Liquid flow equation
//  original: N^T*S*N*dp/dt + dN^T*k/mu*dN*p+dN^T*rho*g*z = 0
//---------------------------------------------------------
#define ENABLE_PT
#define ENABLE_TP
#define ENABLE_TP_MASS

#if 1
			if (updateA)
			{
				const double Mpp_coeff =
				    isTransient ? dt_inverse * fkt * Ss : 0;
				const double Lpp_coeff1 = theta * fkt / vis;
				// const double Lpp_coeff2 = theta * fkt/(vis*vis);
				const double Lpp_coeff3 = theta * fkt / vis * dvis_dp;
				const double LpT_coeff3 = theta * fkt / vis * dvis_dT;
				// const double Gpp_coeff1 = - theta * fkt / (vis*vis) * dvis_dp
				// * rho_w * g_const;
				const double Gpp_coeff2 =
				    theta * fkt / vis * drho_w_dp * g_const;
				// const double GpT_coeff3 = - theta * fkt / (vis*vis) * dvis_dT
				// * rho_w * g_const;
				const double GpT_coeff4 =
				    theta * fkt / vis * drho_w_dT * g_const;

				for (int i = 0; i < c_nnodes; i++)
				{
					// symmetric
					for (int j = i; j < c_nnodes; j++)
					{
						double v_pp = .0;
						// Jpp = 1/dt*N^T*Ss*N
						v_pp = Mpp_coeff * c_shapefct[i] * c_shapefct[j];
						for (unsigned k = 0; k < c_dim; k++)
						{
							for (unsigned l = 0; l < c_dim; l++)
							{
								// Jpp = theta*dN^T*k/mu*dN
								v_pp += Lpp_coeff1 *
								        c_dshapefct[k * c_nnodes + i] *
								        k_tensor[c_dim * k + l] *
								        c_dshapefct[l * c_nnodes + j];
							}
						}

						(*StiffMatrix)(offset_p + i, offset_p + j) += v_pp;
						if (i != j)
						{  // make it symmetric
							(*StiffMatrix)(offset_p + j, offset_p + i) += v_pp;
						}
					}
					// unsymmetric
					for (int j = 0; j < c_nnodes; j++)
					{
						double v_pp = .0;
						double v_pT = .0;
						for (unsigned k = 0; k < c_dim; k++)
						{
							const double dN_N =
							    c_dshapefct[k * c_nnodes + i] * c_shapefct[j];
							const double dN_kz_N =
							    dN_N * k_tensor[c_dim * k + c_dim - 1];
							const double dN_v_N =
							    c_dshapefct[k * c_nnodes + i] * vel[k] *
							    c_shapefct[j];
							// Jpp += theta/mu*dmu_dp*dN^T*q*N
							v_pp += Lpp_coeff3 * dN_v_N;
							// JpT += theta/mu*dmu_dT*dN^T*q*N
							v_pT += LpT_coeff3 * dN_v_N;
							// Jpp += theta*dN^T*k/mu*drho_w_dp*g*N
							v_pp += Gpp_coeff2 * dN_kz_N;
							// Jpt += theta*dN^T*k/mu*drho_w_dT*g*N
							v_pT += GpT_coeff4 * dN_kz_N;
						}
						(*StiffMatrix)(offset_p + i, offset_p + j) += v_pp;
#ifdef ENABLE_PT
						(*StiffMatrix)(offset_p + i, offset_T + j) += v_pT;
#endif
					}
				}
			}  // update A

			if (updateRHS)
			{
				// Rp += [1/dt*N^T*Ss*N + theta*dN^T*k/mu*dN]*p1 -
				// [1/dt*N^T*Ss*N - (1-theta)*dN^T*k/mu*dN]*p0 +
				// dN^T*k/mu*rho*g*z
				const double Rp_coeff1 =
				    isTransient ? dt_inverse * fkt * Ss * (gp_p1 - gp_p0) : 0;
				for (int i = 0; i < c_nnodes; i++)
				{
					double v_p = 0;
					// 1/dt*N^T*Ss*(p1-p0)
					v_p += Rp_coeff1 * c_shapefct[i];
					for (unsigned k = 0; k < c_dim; k++)
					{
						// - dN^T*vel
						v_p += -fkt * c_dshapefct[k * c_nnodes + i] * vel[k];
					}
					(*RHS)(offset_p + i) += v_p;
				}
			}  // update RHS
               //---------------------------------------------------------
               //  Assemble Heat transport equation
//---------------------------------------------------------
#ifdef CONST_HEAT
		}
		double varT[3] = {1750000., 333.15, .0};
		const double rho_w = FluidProp->Density(varT);
		const double drho_w_dp = 0;  // FluidProp->drhodP(varT);
		const double drho_w_dT = 0;  // FluidProp->drhodT(varT);
		const double vis = FluidProp->Viscosity(varT);
		const double dvis_dp = 0;  // FluidProp->dViscositydP(varT);
		const double dvis_dT = 0;  // FluidProp->dViscositydT(varT);
		const double cp_w = FluidProp->SpecificHeatCapacity(varT);
		const double rhocp =
		    MediaProp->HeatCapacity(Index, pcs->m_num->ls_theta, this, varT);
		const double drhocp_dp =
		    0;  // MediaProp->dHeatCapacitydP(Index,pcs->m_num->ls_theta,this,varT);
		const double drhocp_dT =
		    0;  // MediaProp->dHeatCapacitydT(Index,pcs->m_num->ls_theta,this,varT);
		const int offset_p = 0;
		const int offset_T = c_nnodes;
#endif

		double supg_tau = 0;
		if (useSUPG) CalcSUPGWeightingFunction(vel, gp, supg_tau, weight_func);
		double const* const c_weight_func = weight_func;

		if (updateA)
		{
			// Jtt = 1/dt*N^T*rho*cp*N
			// Jtt += 1/dt*N^T*drhocp_dT*T1*N
			// Jtp += 1/dt*N^T*drhocp_dp*T1*N
			//		const double M_tt_coeff1 = isTransient ? dt_inverse * fkt *
			// rhocp : 0;
			const double M_tt_coeff1 =
			    isTransient
			        ? dt_inverse * fkt * (rhocp + drhocp_dT * (gp_T1 - gp_T0))
			        : 0;
			const double M_tp_coeff1 =
			    isTransient ? dt_inverse * fkt * drhocp_dp * (gp_T1 - gp_T0)
			                : 0;
			if (useLumpedMass)
			{
#ifdef USE_LMASS2
				const double R_t_coeff1 =
				    isTransient ? dt_inverse * fkt * rhocp : 0;
				//			double vol = MeshElement->GetVolume() *
				// MeshElement->GetFluxArea();
				//			// Center of the reference element
				//			SetCenterGP();
				//			double factor = CalCoefMass();
				//			factor *= vol / (double)nnodes;
				//			for (int i = 0; i < nnodes; i++)
				//				(*StiffMatrix)(offset_T+i,offset_T+i) =  factor;

				for (int i = 0; i < c_nnodes; i++)
				{
					double NN = .0;
					for (int j = 0; j < c_nnodes; j++)
					{
						NN += c_shapefct[i] * c_shapefct[j];
					}
					(*StiffMatrix)(offset_T + i, offset_T + i) +=
					    M_tt_coeff1 * NN;
#ifdef ENABLE_TP_MASS1
					for (int k = 0; k < c_nnodes; k++)
					{
						//					(*StiffMatrix)(offset_T+i, offset_T+k)
						//+=
						// dt_inverse * fkt * drhocp_dT * NodalVal_T1[i] * NN *
						// c_shapefct[k];
						(*StiffMatrix)(offset_T + i, offset_p + k) +=
						    dt_inverse * fkt * drhocp_dp * NodalVal_T1[i] * NN *
						    c_shapefct[k];
					}
#endif
					(*RHS)(offset_T + i) +=
					    R_t_coeff1 * NN * (NodalVal_T1[i] - NodalVal_T0[i]);
#ifdef TH_DEBUG
					(*Mass)(i, i) += 1. / dt_inverse * R_t_coeff1 * NN;
#endif
				}
#endif
			}
			else
			{
				for (int i = 0; i < c_nnodes; i++)
					for (int j = 0; j < c_nnodes; j++)
					{
						const double NN = c_shapefct[i] * c_shapefct[j];
						(*StiffMatrix)(offset_T + i, offset_T + j) +=
						    M_tt_coeff1 * NN;
						(*StiffMatrix)(offset_T + i, offset_p + j) +=
						    M_tp_coeff1 * NN;
#ifdef TH_DEBUG
						(*Mass)(i, j) += 1. / dt_inverse * M_tt_coeff1 * NN;
#endif
					}
				if (useSUPG)
				{
					// tau*({v}[dN])^T*[N]
					for (int i = 0; i < c_nnodes; i++)
						for (int j = 0; j < c_nnodes; j++)
						{
							const double NN = c_weight_func[i] * c_shapefct[j];
							(*StiffMatrix)(offset_T + i, offset_T + j) +=
							    supg_tau * M_tt_coeff1 * NN;
							(*StiffMatrix)(offset_T + i, offset_p + j) +=
							    supg_tau * M_tp_coeff1 * NN;
#ifdef TH_DEBUG
							(*Mass)(i, j) +=
								supg_tau * 1. / dt_inverse * M_tt_coeff1 * NN;
#endif
						}
				}
			}

			// L = theta*dN^T*lambda*dN*T
			// Jtt += theta*dN^T*lambda*dN
			const double theta_fkt = theta * fkt;
			const double theta_fkt_Cw = theta_fkt * rho_w * cp_w;
			const double Att_coeff2 = theta_fkt * drho_w_dT * cp_w;
			const double Att_coeff3 = theta_fkt * drho_w_dp * cp_w;
			//		const double Att_coeff4 = theta_fkt_Cw / (vis*vis) *
			// dvis_dT;
			//		const double Att_coeff5 = theta_fkt_Cw / (vis*vis) *
			// dvis_dp;
			//		const double Att_coeff6 = - theta_fkt_Cw / vis;
			//		double G_tt_coeff1 = theta_fkt_Cw / vis
			//*(1./vis*dvis_dT*rho_w
			//- drho_w_dT) * g_const;
			//		double G_tp_coeff1 = theta_fkt_Cw / vis
			//*(1./vis*dvis_dp*rho_w
			//- drho_w_dp) * g_const;
			for (int i = 0; i < c_nnodes; i++)
			{
				// symmetric
				for (int j = i; j < c_nnodes; j++)
				{
					double v_TT = 0;
					for (unsigned k = 0; k < c_dim; k++)
					{
						for (unsigned l = 0; l < c_dim; l++)
						{
							v_TT += theta_fkt * c_dshapefct[k * c_nnodes + i] *
							        lambda_tensor[c_dim * k + l] *
							        c_dshapefct[l * c_nnodes + j];
#ifdef TH_DEBUG
							(*Laplace)(i, j) += fkt *
							                    lambda_tensor[dim * k + l] *
							                    dshapefct[k * nnodes + i] *
							                    dshapefct[l * nnodes + j];
							if (i != j)
								(*Laplace)(j, i) += fkt *
													lambda_tensor[dim * k + l] *
													dshapefct[k * nnodes + i] *
													dshapefct[l * nnodes + j];
#endif
						}
					}
					(*StiffMatrix)(offset_T + i, offset_T + j) += v_TT;
					if (i != j)
						(*StiffMatrix)(offset_T + j, offset_T + i) += v_TT;
				}
				// unsymmetric
				for (int j = 0; j < c_nnodes; j++)
				{
					double v_TT = 0;
					double v_Tp = 0;
					const double NN = c_shapefct[i] * c_shapefct[j];
					const double WN =
					    useSUPG ? c_weight_func[i] * c_shapefct[j] : 0;

					for (unsigned k = 0; k < c_dim; k++)
					{
						const double NdN =
						    c_shapefct[i] * c_dshapefct[k * c_nnodes + j];
						// Jtt += theta*N^T*rho_w*cp_w*q^T*dN
						v_TT += theta_fkt_Cw * vel[k] * NdN;
						// Jtt += theta*N^T*drho_w_dT*cp_w*q^T*grad(T)*N
						v_TT += Att_coeff2 * vel[k] * grad_T1[k] * NN;
						// Jtp += theta*N^T*drho_w_dp*cp_w*q^T*grad(T)*N
						v_Tp += Att_coeff3 * vel[k] * grad_T1[k] * NN;
#ifdef TH_DEBUG
						(*Advection)(i, j) += fkt * shapefct[i] * rho_w * cp_w *
						                      vel[k] *
						                      dshapefct[k * nnodes + j];
#endif
						const double N_gT_kz_N =
						    grad_T1[k] * k_tensor[c_dim * k + c_dim - 1] * NN;
						// Jtt += -
						// theta*N^T*rho_w*cp_w*grad(T)*k/mu*g*drho_dT*N
						v_TT += -theta_fkt_Cw / vis * drho_w_dT * g_const *
						        N_gT_kz_N;
						// Jtt += - theta*N^T*rho_w*cp_w*grad(T)*1/mu*q*dmu_dT*N
						v_TT += -theta_fkt_Cw * grad_T1[k] / vis * vel[k] *
						        dvis_dT * NN;

						// Jtp += -
						// theta*N^T*rho_w*cp_w*grad(T)*k/mu*g*drho_dp*N
						v_Tp += -theta_fkt_Cw / vis * drho_w_dp * g_const *
						        N_gT_kz_N;
						// Jtp += - theta*N^T*rho_w*cp_w*grad(T)/mu*q*dmu_dp*N
						v_Tp += -theta_fkt_Cw * grad_T1[k] / vis * vel[k] *
						        dvis_dp * NN;
						for (unsigned l = 0; l < c_dim; l++)
						{
							// Jtp += - theta*N^T*rho_w*cp_w*grad(T)*k/mu*dN
							v_Tp += -theta_fkt_Cw * grad_T1[k] *
							        k_tensor[c_dim * k + l] / vis *
							        c_shapefct[i] *
							        c_dshapefct[l * c_nnodes + j];
						}

						if (useSUPG)
						{
							// Calculate mat_factor*tau*({v}[dN])^T*({v}[dN])
							const double WdN = c_weight_func[i] *
							                   c_dshapefct[k * c_nnodes + j];
							v_TT += supg_tau * theta_fkt_Cw * vel[k] * WdN;
							v_TT += supg_tau * Att_coeff2 * vel[k] *
							        grad_T1[k] * WN;
							v_Tp += supg_tau * Att_coeff3 * vel[k] *
							        grad_p1[k] * WN;
#ifdef TH_DEBUG
							(*Advection)(i, j) +=
							    supg_tau * fkt * rho_w * cp_w * vel[k] * WdN;
#endif
							const double W_gT_kz_N =
							    grad_T1[k] * k_tensor[c_dim * k + c_dim - 1] *
							    WN;
							// Jtt += -
							// theta*N^T*rho_w*cp_w*grad(T)*k/mu*g*drho_dT*N
							v_TT += -theta_fkt_Cw / vis * drho_w_dT * g_const *
							        W_gT_kz_N;
							// Jtt += -
							// theta*N^T*rho_w*cp_w*grad(T)*1/mu*q*dmu_dT*N
							v_TT += -theta_fkt_Cw * grad_T1[k] / vis * vel[k] *
							        dvis_dT * WN;

							// Jtp += -
							// theta*N^T*rho_w*cp_w*grad(T)*k/mu*g*drho_dp*N
							v_Tp += -theta_fkt_Cw / vis * drho_w_dp * g_const *
							        W_gT_kz_N;
							// Jtp += -
							// theta*N^T*rho_w*cp_w*grad(T)/mu*q*dmu_dp*N
							v_Tp += -theta_fkt_Cw * grad_T1[k] / vis * vel[k] *
							        dvis_dp * WN;
							for (unsigned l = 0; l < c_dim; l++)
							{
								// Jtp += - theta*N^T*rho_w*cp_w*grad(T)*k/mu*dN
								v_Tp += -theta_fkt_Cw * grad_T1[k] *
								        k_tensor[c_dim * k + l] / vis *
								        c_weight_func[i] *
								        c_dshapefct[l * c_nnodes + j];
							}
						}
					}
#define ADV
#ifdef ADV
					(*StiffMatrix)(offset_T + i, offset_T + j) += v_TT;
#ifdef ENABLE_TP
					(*StiffMatrix)(offset_T + i, offset_p + j) += v_Tp;
#endif
#endif
				}
			}
		}  // update A

		if (updateRHS)
		{
			// Rt += [1/dt*N^T*Cp*N + theta*(dN^T*lambda*dN+N^T*Cp_w*dN)]*T1 -
			// [1/dt*N^T*Cp*N - (1-theta)*(dN^T*lambda*dN+N^T*Cp_w*dN)]*T0
			if (isTransient)
			{
				if (useLumpedMass)
				{
#ifdef USE_LMASS2
					for (int i = 0; i < c_nnodes; i++)
					{
						double row_sum = .0;
						for (int j = 0; j < c_nnodes; j++)
							row_sum += dt_inverse * fkt * rhocp *
							           c_shapefct[i] * c_shapefct[j];
						(*RHS)(offset_T + i) +=
						    row_sum * (NodalVal_T1[i] - NodalVal_T0[i]);
					}
#endif
				}
				else
				{
					const double coeff =
					    dt_inverse * fkt * rhocp * (gp_T1 - gp_T0);
					for (int i = 0; i < c_nnodes; i++)
						(*RHS)(offset_T + i) += c_shapefct[i] * coeff;
					if (useSUPG)
					{
						// tau*({v}[dN])^T*[N]
						for (int i = 0; i < c_nnodes; i++)
							(*RHS)(offset_T + i) +=
							    supg_tau * c_weight_func[i] * coeff;
					}
				}
			}
			for (int i = 0; i < c_nnodes; i++)
			{
				double v_T = 0;
				for (unsigned k = 0; k < c_dim; k++)
				{
					for (unsigned l = 0; l < c_dim; l++)
					{
						// dN^T*lambda*(theta*gradT1+(1-theta)*gradT0)
						v_T +=
						    fkt * c_dshapefct[k * c_nnodes + i] *
						    lambda_tensor[c_dim * k + l] *
						    (theta * grad_T1_r[l] + (1. - theta) * grad_T0[l]);
					}
#ifdef ADV
					// N^T*Cp_w*(theta*gradT1+(1-theta)*gradT0)
					const double cp_v_gradT =
					    fkt * rho_w * cp_w * vel[k] *
					    (theta * grad_T1_r[k] + (1. - theta) * grad_T0[k]);
					v_T += c_shapefct[i] * cp_v_gradT;
					if (useSUPG)
					{
						// Calculate mat_factor*tau*({v}[dN])^T*({v}[dN])
						v_T += supg_tau * c_weight_func[i] * cp_v_gradT;
					}
#endif
				}
				(*RHS)(offset_T + i) += v_T;
			}
		}  // update RHS
#else

#endif
	}

	// RHS should be - residual
	if (updateRHS) (*RHS) *= -1.;

#if !defined(USE_PETSC)
	for (size_t ii = 0; ii < pcs->GetPrimaryVNumber(); ii++)
		for (long i = 0; i < c_nnodes; i++)
			eqs_rhs[NodeShift[ii] + eqs_number[i]] += (*RHS)(i + ii * c_nnodes);
#endif

	if (updateA && (pcs->scaleUnknowns || pcs->scaleEQS))
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
	if (updateRHS && pcs->scaleEQS)
	{
		for (int ii = 0; ii < 2; ii++)
		{
			const double scale_eqs = pcs->vec_scale_eqs[ii];
			for (int i = 0; i < c_nnodes; i++)
				(*RHS)(ii* c_nnodes + i) *= scale_eqs;
		}
	}
}
#endif

}  // end namespace

