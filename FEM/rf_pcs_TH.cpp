/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_pcs_TH.h"

#ifdef USE_PETSC
#include <petscksp.h>
#include <petsctime.h>
#endif

#include "display.h"
#include "StringTools.h"
#if defined(NEW_EQS)
#include "equation_class.h"
#elif defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#include "fem_ele_std.h"

CRFProcessTH::CRFProcessTH() : CRFProcess(), error_k0(0.0)
{
}

CRFProcessTH::~CRFProcessTH()
{
}

void CRFProcessTH::Initialization()
{
	int Axisymm = 1;                           // ani-axisymmetry
	if (m_msh->isAxisymmetry()) Axisymm = -1;  // Axisymmetry is true
	fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
	fem->SetGaussPointNumber(m_num->ele_gauss_points);

	if (vec_scale_dofs.empty()) vec_scale_dofs.resize(2, 1.);
	if (vec_scale_eqs.empty()) vec_scale_eqs.resize(2, 1.);
	if (scaleUnknowns)
	{
		std::stringstream ss;
		ss << "-> scale DOF:";
		for (unsigned i = 0; i < vec_scale_dofs.size(); i++)
			ss << "[" << i << "] " << vec_scale_dofs[i] << " ";
		ScreenMessage("%s\n", ss.str().c_str());
	}
	if (scaleEQS)
	{
		std::stringstream ss;
		ss << "-> scale EQS:";
		for (unsigned i = 0; i < vec_scale_eqs.size(); i++)
			ss << "[" << i << "] " << vec_scale_eqs[i] << " ";
		ScreenMessage("%s\n", ss.str().c_str());
	}

#if 0
	// calculate initial velocity
	ScreenMessage("-> compute velocity from initial pressure\n");
	CalIntegrationPointValue();
#endif
}

double CRFProcessTH::Execute(int loop_process_number)
{
	ScreenMessage("\n================================================\n");
	ScreenMessage("->Process %d : %s\n", loop_process_number,
	              convertProcessTypeToString(getProcessType()).c_str());
	ScreenMessage("================================================\n");

	clock_t dm_time = -clock();

	m_msh->SwitchOnQuadraticNodes(false);
	if (hasAnyProcessDeactivatedSubdomains || Deactivated_SubDomain.size() > 0)
		CheckMarkedElement();

#if defined(USE_PETSC)
	// set Dirichlet BC to nodal values
	ScreenMessage("-> set initial guess\n");
	IncorporateBoundaryConditions(-1, false, false, true);
	copyNodalValuesToVec(eqs_new->total_x);
	// VecView(eqs_new->total_x, PETSC_VIEWER_STDOUT_WORLD);

	if (m_num->petsc_use_snes)
	{
		bool error = ExecuteNonlinearWithPETsc();
		dm_time += clock();
		ScreenMessage("CPU time elapsed in this process: %gs\n",
		              (double)dm_time / CLOCKS_PER_SEC);
		ScreenMessage("------------------------------------------------\n");
		return error;
	}

	eqs_x = eqs_new->GetGlobalSolution();
#endif
#if defined(NEW_EQS)  // WW
//
	eqs_new->ConfigNumerics(m_num);  // 27.11.2007 WW
//
#endif

	// Begin Newton-Raphson steps
	double Error = 1.0;
	//	double Error1 = 0.0;
	//	double ErrorU = 1.0;
	//	double ErrorU1 = 0.0;
	double InitialNorm = 0.0;
	//	double InitialNormDx = 0.0;
	//	double InitialNormU = 0.0;
	double NormDx = std::numeric_limits<double>::max();
#ifdef USE_PETSC
	static double rp0 = .0, rT0 = 0;
	static double rp0_L2 = .0, rT0_L2 = 0;
	double dp_max = std::numeric_limits<double>::max(),
	       dT_max = std::numeric_limits<double>::max();
	double dp_L2 = std::numeric_limits<double>::max(),
	       dT_L2 = std::numeric_limits<double>::max();
	double p_max = std::numeric_limits<double>::max(),
	       T_max = std::numeric_limits<double>::max();
	double p_L2 = std::numeric_limits<double>::max(),
	       T_L2 = std::numeric_limits<double>::max();
	const double tol_dp = m_num->nls_error_tolerance[1];
	const double tol_dT = m_num->nls_error_tolerance[2];
#endif

	const double newton_tol = m_num->nls_error_tolerance[0];
	const int n_max_iterations = m_num->nls_max_iterations;

	iter_nlin = 0;
	bool converged = false;
	while (iter_nlin < n_max_iterations)
	{
		iter_nlin++;

//----------------------------------------------------------------------
// Solve
//----------------------------------------------------------------------
// Refresh solver
#if defined(NEW_EQS)
		eqs_new->Initialize();  // 27.11.2007 WW
#elif defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
		eqs_new->Initialize();
#endif

		ScreenMessage("-> Assembling equation system...\n");
		GlobalAssembly();

//
#if defined(NEW_EQS)
		const double NormR = eqs_new->ComputeNormRHS();
//		const double NormR = eqs_new->NormRHS();
#elif defined(USE_PETSC)
		const double NormR = eqs_new->GetVecNormRHS();
#endif
#if defined(USE_PETSC)
		double rp_max = std::numeric_limits<double>::max(),
		       rT_max = std::numeric_limits<double>::max();
		double rp_L2 = std::numeric_limits<double>::max(),
		       rT_L2 = std::numeric_limits<double>::max();
		Vec sub_x;
		VecGetSubVector(eqs_new->b, eqs_new->vec_isg[0], &sub_x);
		VecNorm(sub_x, NORM_2, &rp_L2);
		VecNormBegin(sub_x, NORM_INFINITY, &rp_max);
		VecNormEnd(sub_x, NORM_INFINITY, &rp_max);
		VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[0], &sub_x);
		VecGetSubVector(eqs_new->b, eqs_new->vec_isg[1], &sub_x);
		VecNorm(sub_x, NORM_2, &rT_L2);
		VecNormBegin(sub_x, NORM_INFINITY, &rT_max);
		VecNormEnd(sub_x, NORM_INFINITY, &rT_max);
		VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[1], &sub_x);
#endif
		if (iter_nlin == 1 && this->first_coupling_iteration)
		{
			InitialNorm = NormR;
#ifdef USE_PETSC
			rp0 = rp_max;
			rT0 = rT_max;
			rp0_L2 = rp_L2;
			rT0_L2 = rT_L2;
#endif
			static bool firstime = true;
			if (firstime)
				firstime = false;
		}
#ifdef USE_PETSC
		Error = std::max(rp_max / rp0, rT_max / rT0);  // NormR / InitialNorm;
		const double Error_L2 = std::max(rp_L2 / rp0_L2, rT_L2 / rT0_L2);
		const double dx_i = std::max(dp_max / p_max, dT_max / T_max);
		const double dx_L2 = std::max(dp_L2 / p_L2, dT_L2 / T_L2);
		ScreenMessage(
		    "-> Newton-Raphson: r_i=%.3e, r_2=%.3e, dx_i=%.3e, dx_2=%.3e "
		    "(tol=%g) %d/%d \n",
		    Error, Error_L2, dx_i, dx_L2, newton_tol, iter_nlin - 1,
		    n_max_iterations);
		ScreenMessage("|r0|=%.3e, |r|=%.3e, |r|/|r0|=%.3e, |dx|=%.3e\n",
		              InitialNorm, NormR, NormR / InitialNorm, NormDx);
		ScreenMessage(
		    "|rp|_i=%.3e, |rT|_i=%.3e, |rp/r0|_i=%.3e, |rT/r0|_i=%.3e\n",
		    rp_max, rT_max, rp_max / rp0, rT_max / rT0);
		ScreenMessage(
		    "|rp|_2=%.3e, |rT|_2=%.3e, |rp/r0|_2=%.3e, |rT/r0|_2=%.3e\n", rp_L2,
		    rT_L2, rp_L2 / rp0_L2, rT_L2 / rT0_L2);
		ScreenMessage(
		    "|dp|_i=%.3e, |dT|_i=%.3e, |dp/p|_i=%.3e, |dT/T|_i=%.3e "
		    "(tol.p=%.1e,T=%.1e)\n",
		    dp_max, dT_max, dp_max / p_max, dT_max / T_max, tol_dp, tol_dT);
		ScreenMessage(
		    "|dp|_2=%.3e, |dT|_2=%.3e, |dp/p|_2=%.3e, |dT/T|_2=%.3e "
		    "(tol.p=%.1e,T=%.1e)\n",
		    dp_L2, dT_L2, dp_L2 / p_L2, dT_L2 / T_L2, tol_dp, tol_dT);
		ScreenMessage("------------------------------------------------\n");
		ScreenMessage("-> Nonlinear iteration: %d/%d, |r|=%.3e, |r|/|r0|=%.3e\n", iter_nlin - 1,
		              n_max_iterations, NormR, NormR / InitialNorm);
		ScreenMessage("------------------------------------------------\n");
#else
		Error = NormR / InitialNorm;
#if 1
		double r_dof[2] = {};
		for (size_t ii = 0; ii < this->GetPrimaryVNumber(); ii++) {
			const size_t nnodes = this->m_msh->GetNodesNumber(false);
			const size_t shift = nnodes * ii;
			for (size_t i = 0; i < nnodes; i++)
				r_dof[ii] += std::pow(eqs_new->getRHS()[shift + i], 2);
			r_dof[ii] = std::sqrt(r_dof[ii]);
		}
		ScreenMessage("------------------------------------------------\n");
		ScreenMessage("-> Nonlinear iteration: %d/%d, |r|=%.3e, |r|/|r0|=%.3e, |rp|=%.3e, |rT|=%.3e\n", iter_nlin - 1,
		              n_max_iterations, NormR, NormR / InitialNorm, r_dof[0], r_dof[1]);
		ScreenMessage("------------------------------------------------\n");
#else
		ScreenMessage("------------------------------------------------\n");
		ScreenMessage("-> Nonlinear iteration: %d/%d, |r|=%.3e, |r|/|r0|=%.3e\n", iter_nlin - 1,
		              n_max_iterations, NormR, NormR / InitialNorm);
		ScreenMessage("------------------------------------------------\n");
#endif
#endif


		if (Error < newton_tol)
		{
			ScreenMessage("-> Newton-Raphson converged\n");
			converged = true;
			break;
		}

		ScreenMessage("-> Calling linear solver...\n");
// Linear solver
#if defined(NEW_EQS)
		bool compress_eqs =
			(type / 10 == 4 || this->Deactivated_SubDomain.size() > 0);
		iter_lin = eqs_new->Solver(this->m_num, compress_eqs);  // NW
#elif defined(USE_PETSC)
		//		if (write_leqs) {
		//			std::string fname = FileName + "_" +
		// convertProcessTypeToString(this->getProcessType()) +
		//"_leqs_assembly.txt";
		//			eqs_new->EQSV_Viewer(fname);
		//		}
		iter_lin = eqs_new->Solver();
		//		if (write_leqs) {
		//			std::string fname = FileName + "_" +
		// convertProcessTypeToString(this->getProcessType()) +
		//"_leqs_solution.txt";
		//			eqs_new->EQSV_Viewer(fname);
		//		}
		//		if (iter_nlin==1) {
		//			std::string fname = FileName + "_" +
		// convertProcessTypeToString(this->getProcessType()) +
		//"_leqs_residual.txt";
		//			eqs_new->Residual_Viewer(fname);
		//		}
		if (!m_num->petsc_split_fields)
		{
			eqs_new->MappingSolution();
			eqs_x = eqs_new->GetGlobalSolution();
		}
		else
		{
			VecAXPY(eqs_new->total_x, 1.0, eqs_new->x);
			//			VecView(eqs_new->total_x, PETSC_VIEWER_STDOUT_SELF);
		}
#endif
		if (iter_lin < 0)
		{
			accepted = false;
			Tim->last_dt_accepted = false;
			break;
		}
		iter_lin_max = std::max(iter_lin_max, iter_lin);

//----------------------------------------------------------------------
// Check convergence
//----------------------------------------------------------------------
#if defined(USE_PETSC)
		if (m_num->petsc_split_fields)
		{
			Vec sub_x;
			// dx
			VecGetSubVector(eqs_new->x, eqs_new->vec_isg[0], &sub_x);
			VecNorm(sub_x, NORM_2, &dp_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &dp_max);
			VecNormEnd(sub_x, NORM_INFINITY, &dp_max);
			VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[0], &sub_x);
			VecGetSubVector(eqs_new->x, eqs_new->vec_isg[1], &sub_x);
			VecNorm(sub_x, NORM_2, &dT_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &dT_max);
			VecNormEnd(sub_x, NORM_INFINITY, &dT_max);
			VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[1], &sub_x);
			// total x
			VecGetSubVector(eqs_new->total_x, eqs_new->vec_isg[0], &sub_x);
			VecNorm(sub_x, NORM_2, &p_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &p_max);
			VecNormEnd(sub_x, NORM_INFINITY, &p_max);
			VecRestoreSubVector(eqs_new->total_x, eqs_new->vec_isg[0], &sub_x);
			VecGetSubVector(eqs_new->total_x, eqs_new->vec_isg[1], &sub_x);
			VecNorm(sub_x, NORM_2, &T_L2);
			VecNormBegin(sub_x, NORM_INFINITY, &T_max);
			VecNormEnd(sub_x, NORM_INFINITY, &T_max);
			VecRestoreSubVector(eqs_new->total_x, eqs_new->vec_isg[1], &sub_x);
			// relative error
		}
#endif

// Norm of dx
#if defined(NEW_EQS)
		NormDx = eqs_new->NormX();
#elif defined(USE_PETSC)
		NormDx = eqs_new->GetVecNormX();
#endif
		ScreenMessage("-> |dx|=%.3e\n", NormDx);

// Check the convergence
//		Error1 = Error;
//		ErrorU1 = ErrorU;
//		if(iter_nlin == 1 && this->first_coupling_iteration)
//		{
//			InitialNormDx = NormDx;
//			static bool firstime = true;
//			if (firstime) {
////				InitialNormU = NormDx;
//				firstime = false;
//			}
//		}

#if 0
		ErrorU = NormDx / InitialNormDx;
		if(NormR < newton_tol && Error > NormR)
			Error = NormR;
		//           if(Norm<TolNorm)  Error = 0.01*Tolerance_global_Newton;
		if((NormDx / InitialNormU) <= newton_tol)
			Error = NormDx / InitialNormU;
		if(ErrorU < Error)
			Error = ErrorU;
#endif
		// JT: Store the process and coupling errors
		pcs_num_dof_errors = 1;
		if (iter_nlin == 1)
		{
			pcs_absolute_error[0] = NormDx;
			pcs_relative_error[0] = pcs_absolute_error[0] / newton_tol;
			cpl_max_relative_error = pcs_relative_error[0];
			cpl_num_dof_errors = 1;
		}
		else
		{
			pcs_absolute_error[0] = Error;
			pcs_relative_error[0] = Error / newton_tol;
		}
//
// Screan printing:
//		ScreenMessage("-> update solutions\n");
#if 0
		if(Error > 100.0 && iter_nlin > 1)
		{
			ScreenMessage ("\n  Attention: Newton-Raphson step is diverged. Programme halt!\n");
			accepted = false;
			Tim->last_dt_accepted = false;
			return -1;
		}
#endif
#ifdef USE_PETSC
		if (std::max(rp_max / rp0, rT_max / rT0) < newton_tol ||
		    (dp_max < tol_dp && dT_max < tol_dT))
		{
			ScreenMessage("-> Newton-Raphson converged\n");
			converged = true;
			break;
		}
#endif
		//		if(InitialNorm < 10 * newton_tol
		//			|| NormR < 0.001 * InitialNorm
		//			|| Error <= newton_tol)
		//			break;

		// x^k1 = x^k + dx
		UpdateIterativeStep(1.0);

		//		ScreenMessage("-> update velocity\n");
		//		CalIntegrationPointValue();
	}  // Newton-Raphson iteration

	// x^k1 = x^k + dx
	UpdateIterativeStep(1.0);

	iter_nlin_max = std::max(iter_nlin_max, iter_nlin);

	if (!converged && m_num->nls_max_iterations > 1) accepted = false;

	//
	dm_time += clock();
#if defined(USE_MPI) || defined(USE_PETSC)
	if (myrank == 0)
	{
#endif
		std::cout << "CPU time elapsed in this process: "
		          << (double)dm_time / CLOCKS_PER_SEC << "s"
		          << "\n";
		std::cout << "------------------------------------------------"
		          << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
	}
#endif
// Recovery the old solution.  Temp --> u_n	for flow proccess
//	RecoverSolution();
//
#ifdef NEW_EQS  // WW
	// Also allocate temporary memory for linear solver. WW
	eqs_new->Clean();
#endif

	// For coupling control
	Error = 0.0;
	for (size_t n = 0; n < m_msh->GetNodesNumber(false); n++)
	{
		for (int l = 0; l < pcs_number_of_primary_nvals; l++)
		{
			double NormU = GetNodeValue(n, 2 * l) - GetNodeValue(n, 2 * l + 1);
			Error += NormU * NormU;
		}
	}
	double sqrt_norm = sqrt(Error);

	if (this->first_coupling_iteration)
		Error = fabs(sqrt_norm - error_k0) / error_k0;
	else
		Error = sqrt_norm;
	error_k0 = sqrt_norm;

	if (!accepted || Tim->isDynamicTimeFailureSuggested(this))
	{
		accepted = false;
		Tim->last_dt_accepted = false;
	}

	return Error;
}

/**************************************************************************
   Update solution in Newton-Raphson procedure
**************************************************************************/
void CRFProcessTH::UpdateIterativeStep(const double damp)
{
//	long shift = 0;
#if defined(NEW_EQS)
	const double* eqs_x = eqs_new->getX();
#endif

	// x^k1 = x^k + dx
	const long number_of_nodes = num_nodes_p_var[0];
#if defined(USE_PETSC)
	if (m_num->petsc_split_fields)
	{
		Vec sub_x;
		// std::vector<double> array_x(m_msh->getNumNodesGlobal());
		double* array_x = eqs_new->global_x0;
		for (int i = 0; i < 2; i++)
		{
			VecGetSubVector(eqs_new->x, eqs_new->vec_isg[i], &sub_x);
			eqs_new->getGlobalVectorArray(sub_x, array_x);
			const double inv_scaling = 1. / vec_scale_dofs[i];
			for (long j = 0; j < number_of_nodes; j++)
			{
				int k = m_msh->Eqs2Global_NodeIndex[j];
				SetNodeValue(j, p_var_index[i],
				             GetNodeValue(j, p_var_index[i]) +
				                 array_x[k] * damp * inv_scaling);
			}
			VecRestoreSubVector(eqs_new->x, eqs_new->vec_isg[i], &sub_x);
		}
	}
	else
	{
#endif
		for (int i = 0; i < 2; i++)
		{
			const double inv_scaling = 1. / vec_scale_dofs[i];
			for (long j = 0; j < number_of_nodes; j++)
			{
#if defined(USE_PETSC)
				double dx = eqs_x[m_msh->Eqs2Global_NodeIndex[j] *
				                      pcs_number_of_primary_nvals +
				                  i];
#else
			double dx =
			    eqs_x[j + number_of_nodes * i] * damp * vec_scale_dofs[i];
//			std::cout << GetNodeValue(j, ColIndex) << " ";
#endif
				SetNodeValue(
				    j, p_var_index[i],
				    GetNodeValue(j, p_var_index[i]) + dx * damp * inv_scaling);
			}
		}
#if defined(USE_PETSC)
	}
#endif
}


#ifdef USE_PETSC
void CRFProcessTH::setSolver(petsc_group::PETScLinearSolver* petsc_solver)
{
	eqs_new = petsc_solver;

	char str_precon_type[128];
	PetscBool set;
	PetscOptionsGetString(NULL, "-pc_type", str_precon_type,
	                      sizeof(str_precon_type), &set);
	if (set)
	{
		std::string str(str_precon_type);
		if (str.compare("fieldsplit") == 0)
			this->m_num->petsc_split_fields = true;
	}
	if (this->m_num->petsc_split_fields)
	{
		ScreenMessage("-> prepare field splits in PETSc\n");
		PetscErrorCode ierr;
		int dof = this->GetPrimaryVNumber();
		eqs_new->vec_subA.resize(dof * dof);
		eqs_new->vec_isg.resize(dof);
		eqs_new->vec_subRHS.resize(dof);
		// setup matrix
		int n_nodes = this->m_msh->getNumNodesGlobal();
		for (int i = 0; i < dof * dof; i++)
		{
			ierr = MatCreate(PETSC_COMM_WORLD, &eqs_new->vec_subA[i]);
			CHKERRCONTINUE(ierr);
			std::string str =
			    "a" + number2str(i / dof) + number2str(i % dof) + "_";
			ierr = MatSetOptionsPrefix(eqs_new->vec_subA[i], str.c_str());
			CHKERRCONTINUE(ierr);
			ierr = MatSetSizes(eqs_new->vec_subA[i], PETSC_DECIDE, PETSC_DECIDE,
			                   n_nodes, n_nodes);
			CHKERRCONTINUE(ierr);
			ierr = MatSetType(eqs_new->vec_subA[i], MATMPIAIJ);
			CHKERRCONTINUE(ierr);
			ierr = MatSetFromOptions(eqs_new->vec_subA[i]);
			CHKERRCONTINUE(ierr);
			ierr = MatMPIAIJSetPreallocation(eqs_new->vec_subA[i],
			                                 eqs_new->d_nz, PETSC_NULL,
			                                 eqs_new->o_nz, PETSC_NULL);
			CHKERRCONTINUE(ierr);
			MatSetOption(eqs_new->vec_subA[i], MAT_NEW_NONZERO_ALLOCATION_ERR,
			             PETSC_FALSE);
			MatSetOption(eqs_new->vec_subA[i], MAT_KEEP_NONZERO_PATTERN,
			             PETSC_TRUE);  // for MatZeroRows()
			                           //			int i_start, i_end;
			//			MatGetOwnershipRange(eqs_new->vec_subA[i],&i_start,&i_end);
			//			ScreenMessage2("-> Sub A[%d]: start=%d, end=%d\n", i,
			// i_start, i_end);
		}
		ierr = MatDestroy(&eqs_new->A);
		CHKERRCONTINUE(ierr);
		ierr = VecDestroy(&eqs_new->b);
		CHKERRCONTINUE(ierr);
		ierr = VecDestroy(&eqs_new->x);
		CHKERRCONTINUE(ierr);
		ierr = MatCreateNest(PETSC_COMM_WORLD, dof, NULL, dof, NULL,
		                     &eqs_new->vec_subA[0], &eqs_new->A);
		CHKERRCONTINUE(ierr);
		ierr = MatGetVecs(eqs_new->A, &eqs_new->b, &eqs_new->x);
		CHKERRCONTINUE(ierr);
		ierr = VecDuplicate(eqs_new->x, &eqs_new->total_x);
		CHKERRCONTINUE(ierr);
		// setup index sets
		ierr = MatNestGetISs(eqs_new->A, &eqs_new->vec_isg[0], NULL);
		CHKERRCONTINUE(ierr);
#if 0
		for (int i=0; i<dof; i++) {
			int is_global_size, is_local_size;
			ISGetSize(eqs_new->vec_isg[i], &is_global_size);
			ISGetLocalSize(eqs_new->vec_isg[i], &is_local_size);
			ScreenMessage2("-> IS[%d]: global size=%d, local size=%d\n", i, is_global_size, is_local_size);
		}
#endif
	}
	else
	{
		ScreenMessage("-> do not use field splits in PETSc\n");
		int ierr = VecDuplicate(eqs_new->x, &eqs_new->total_x);
		CHKERRCONTINUE(ierr);
		int dof = this->GetPrimaryVNumber();
		eqs_new->vec_isg.resize(dof);
		int n_local_row, rstart;
		VecGetLocalSize(eqs_new->b, &n_local_row);
		VecGetOwnershipRange(eqs_new->b, &rstart, NULL);
		// p
		if (rstart % 2 == 0)
			ISCreateStride(PETSC_COMM_WORLD, std::ceil(n_local_row * 0.5),
			               rstart, 2, &eqs_new->vec_isg[0]);
		else
			ISCreateStride(PETSC_COMM_WORLD, std::floor(n_local_row * 0.5),
			               rstart + 1, 2, &eqs_new->vec_isg[0]);
		// T
		if (rstart % 2 == 0)
			ISCreateStride(PETSC_COMM_WORLD, std::floor(n_local_row * 0.5),
			               rstart + 1, 2, &eqs_new->vec_isg[1]);
		else
			ISCreateStride(PETSC_COMM_WORLD, std::ceil(n_local_row * 0.5),
			               rstart, 2, &eqs_new->vec_isg[1]);
	}
	//	ISView(eqs_new->vec_isg[0],PETSC_VIEWER_STDOUT_WORLD);
	//	ISView(eqs_new->vec_isg[1],PETSC_VIEWER_STDOUT_WORLD);

	const int n_local_nodes = this->m_msh->GetNodesNumber(false);
	std::vector<int> idx_from(n_local_nodes), idx_to(n_local_nodes);
	for (long j = 0; j < n_local_nodes; j++)
	{
		idx_from[j] = m_msh->Eqs2Global_NodeIndex[j];
		idx_to[j] = j;
	}
	ISCreateGeneral(PETSC_COMM_SELF, n_local_nodes, &idx_from[0],
	                PETSC_COPY_VALUES, &eqs_new->is_global_node_id);
	ISCreateGeneral(PETSC_COMM_SELF, n_local_nodes, &idx_to[0],
	                PETSC_COPY_VALUES, &eqs_new->is_local_node_id);
	//		ISView(eqs_new->is_global_node_id,
	// PETSC_VIEWER_STDOUT_SELF);CHKERRCONTINUE(ierr);
	//		ISView(eqs_new->is_local_node_id,
	// PETSC_VIEWER_STDOUT_SELF);CHKERRCONTINUE(ierr);

	PetscBool usePetscSNES = PETSC_FALSE;
	PetscOptionsHasName(NULL, "-use_snes", &usePetscSNES);
	if (usePetscSNES)
	{
		ScreenMessage("-> setup SNES\n");
		this->m_num->petsc_use_snes = true;
		eqs_new->ConfigWithNonlinear(
		    m_num->ls_error_tolerance, m_num->ls_max_iterations,
		    m_num->getLinearSolverName(), m_num->getPreconditionerName(),
		    m_num->ls_extra_arg, "");
		//	SNESGSSetTolerances(eqs_new->snes, PETSC_DECIDE,
		// m_num->nls_error_tolerance[0], PETSC_DECIDE,
		// m_num->nls_max_iterations);
		SNESSetFunction(eqs_new->snes, eqs_new->b, FormFunctionTH, (void*)this);
#if (PETSC_VERSION_NUMBER >= 3050)
		SNESSetJacobian(eqs_new->snes, eqs_new->A, eqs_new->A, FormJacobianTH,
		                (void*)this);
#else
		SNESSetJacobian(eqs_new->snes, eqs_new->A, eqs_new->A, FormJacobianTH,
		                (void*)this);
#endif
		SNESSetFromOptions(eqs_new->snes);
	}
	else
	{
		eqs_new->Config(m_num->ls_error_tolerance, m_num->ls_max_iterations,
		                m_num->getLinearSolverName(),
		                m_num->getPreconditionerName(), m_num->ls_extra_arg, "");
	}
}

void CRFProcessTH::copyVecToNodalValues(Vec x)
{
	const long number_of_nodes = num_nodes_p_var[0];
	Vec sub_x_local;
	VecCreateSeq(PETSC_COMM_SELF, number_of_nodes, &sub_x_local);

	Vec sub_x;
	for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
		VecGetSubVector(x, eqs_new->vec_isg[ii], &sub_x);
		VecScatter scatter;
		VecScatterCreate(sub_x, eqs_new->is_global_node_id, sub_x_local,
		                 eqs_new->is_local_node_id, &scatter);
		VecScatterBegin(scatter, sub_x, sub_x_local, INSERT_VALUES,
		                SCATTER_FORWARD);
		VecScatterEnd(scatter, sub_x, sub_x_local, INSERT_VALUES,
		              SCATTER_FORWARD);
		// VecView(sub_x_local, PETSC_VIEWER_STDOUT_SELF);
		const double* array_x;
		VecGetArrayRead(sub_x_local, &array_x);
		const int p_var_id_ii = p_var_index[ii];
		const double inv_scaling = 1. / vec_scale_dofs[ii];
		for (long j = 0; j < number_of_nodes; j++)
			SetNodeValue(j, p_var_id_ii, array_x[j] * inv_scaling);
		VecRestoreArrayRead(x, &array_x);
		VecRestoreSubVector(x, eqs_new->vec_isg[ii], &sub_x);
		VecScatterDestroy(&scatter);
	}

	VecDestroy(&sub_x_local);
}

void CRFProcessTH::copyNodalValuesToVec(Vec x)
{
	//	ScreenMessage2("-> copyNodalValuesToVec\n");
	const long number_of_nodes = num_nodes_p_var[0];
	if (vec_pos.empty())
	{
		vec_pos.resize(number_of_nodes);
		for (long j = 0; j < number_of_nodes; j++)
			vec_pos[j] = m_msh->Eqs2Global_NodeIndex[j];
	}

	Vec sub_x;
	double* vec_values = eqs_new->global_x1;
	for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
		const double scaling = vec_scale_dofs[ii];
		for (long j = 0; j < number_of_nodes; j++)
		{
			vec_values[j] = scaling * GetNodeValue(j, p_var_index[ii]);
		}
		VecGetSubVector(x, eqs_new->vec_isg[ii], &sub_x);
		int subx_size, rstart, rend;
		VecGetSize(sub_x, &subx_size);
		VecGetOwnershipRange(sub_x, &rstart, &rend);
		// ScreenMessage2("-> sub vec: gn=%d, ln=%d, rs=%d\n", subx_size,
		// rend-rstart, rstart);
		VecSetValues(sub_x, number_of_nodes, &vec_pos[0], &vec_values[0],
		             INSERT_VALUES);
		VecAssemblyBegin(sub_x);
		VecAssemblyEnd(sub_x);
		if (!this->m_num->petsc_split_fields)
		{
			IS is;
			ISCreateStride(PETSC_COMM_WORLD, rend - rstart, rstart, 1, &is);
			VecScatter scatter;
			VecScatterCreate(sub_x, is, x, eqs_new->vec_isg[ii], &scatter);
			VecScatterBegin(scatter, sub_x, x, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterEnd(scatter, sub_x, x, INSERT_VALUES, SCATTER_FORWARD);
			VecScatterDestroy(&scatter);
			ISDestroy(&is);
		}
		VecRestoreSubVector(x, eqs_new->vec_isg[ii], &sub_x);
		// VecView(sub_x, PETSC_VIEWER_STDOUT_SELF);
		// VecView(x, PETSC_VIEWER_STDOUT_SELF);
	}
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);
}

double CRFProcessTH::ExecuteNonlinearWithPETsc()
{
	PetscLogDouble v1, v2;
	PetscPrintf(PETSC_COMM_WORLD,
	            "------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "*** PETSc nonlinear solver\n");
	// VecView(eqs_new->total_x, PETSC_VIEWER_STDOUT_SELF);
	PetscTime(&v1);
	SNESSolve(eqs_new->snes, NULL, eqs_new->total_x);
	PetscTime(&v2);
	eqs_new->time_elapsed = v2 - v1;
	SNESConvergedReason reason;
	SNESGetConvergedReason(eqs_new->snes, &reason);
	SNESGetIterationNumber(eqs_new->snes, &iter_nlin);
	iter_nlin_max = std::max(iter_nlin_max, iter_nlin);
	PetscPrintf(PETSC_COMM_WORLD,
	            "------------------------------------------------\n");
	// VecView(eqs_new->total_x, PETSC_VIEWER_STDOUT_SELF);
	if (reason > 0)
	{
		copyVecToNodalValues(eqs_new->total_x);
	}

	if (reason <= 0 || Tim->isDynamicTimeFailureSuggested(this))
	{
		accepted = false;
		Tim->last_dt_accepted = false;
	}

	return 0;  // TODO error
}

PetscErrorCode FormFunctionTH(SNES /*snes*/, Vec x, Vec f, void* ctx)
{
	bool debugOutput = false;
	if (debugOutput) ScreenMessage("-> form a residual function\n");
	// VecView(x, PETSC_VIEWER_STDOUT_WORLD); //PETSC_VIEWER_STDOUT_SELF
	CRFProcessTH* pcs = (CRFProcessTH*)ctx;
	MeshLib::CFEMesh* m_msh = pcs->m_msh;
	petsc_group::PETScLinearSolver* eqs_new = pcs->eqs_new;
	CFiniteElementStd* fem = pcs->GetAssembler();

	Vec b_org = eqs_new->b;
	if (eqs_new->b != f) eqs_new->b = f;

	// update nodal values
	if (debugOutput) ScreenMessage("-> update nodal values\n");
	pcs->copyVecToNodalValues(x);
	//	pcs->IncorporateBoundaryConditions(-1, false, false, true); // set bc
	// node values

	//	if (debugOutput)
	//	ScreenMessage("-> update integration point values\n");
	//	pcs->CalIntegrationPointValue();

	const size_t dn = m_msh->ele_vector.size() / 10;
	const bool print_progress = (dn >= 100) && debugOutput;
	if (print_progress)
		ScreenMessage("start local assembly for %d elements...\n",
		              m_msh->ele_vector.size());

	VecSet(f, 0.0);
	//	VecAssemblyBegin(f);
	//	VecAssemblyEnd(f);
	for (size_t i = 0; i < eqs_new->vec_subRHS.size(); i++)
	{
		VecGetSubVector(f, eqs_new->vec_isg[i], &eqs_new->vec_subRHS[i]);
	}

	const size_t n_ele = m_msh->ele_vector.size();
	for (size_t i = 0; i < n_ele; i++)
	{
		if (print_progress && (i + 1) % dn == 0) ScreenMessage("* ");
		// ScreenMessage("%d \%\n", ((i+1)*100/n_eles));
		MeshLib::CElem* elem = m_msh->ele_vector[i];
		// Marked for use //WX: modified for coupled excavation
		if (elem->GetMark())
		{
			elem->SetOrder(false);
			fem->ConfigElement(elem, false);
			fem->Assembly(false, true);
		}
	}
	//	if (debugOutput)
	if (print_progress) ScreenMessage("done\n");

	for (size_t i = 0; i < eqs_new->vec_subRHS.size(); i++)
	{
		VecAssemblyBegin(eqs_new->vec_subRHS[i]);
		VecAssemblyEnd(eqs_new->vec_subRHS[i]);
		VecRestoreSubVector(f, eqs_new->vec_isg[i], &eqs_new->vec_subRHS[i]);
	}
	VecAssemblyBegin(f);
	VecAssemblyEnd(f);

	if (debugOutput)
		ScreenMessage("-> impose Neumann BC and source/sink terms\n");
	pcs->IncorporateSourceTerms();
	VecAssemblyBegin(f);
	VecAssemblyEnd(f);

	if (debugOutput) ScreenMessage("-> impose Dirichlet BC\n");
	pcs->IncorporateBoundaryConditions(-1, false, true,
	                                   true);  // set bc residual = 0
	                                           //
	VecAssemblyBegin(f);
	VecAssemblyEnd(f);

	// assembly set -r
	VecScale(f, -1.);
//	VecAssemblyBegin(f);
//	VecAssemblyEnd(f);

// VecView(f, PETSC_VIEWER_STDOUT_WORLD); //PETSC_VIEWER_STDOUT_SELF

#if 0
	int offset = 0;
	if (eqs_new->vec_subA.size()>0) {
		int dummy;
		MatGetLocalSize(eqs_new->vec_subA[0], &offset, &dummy);
	}
	unsigned n = 5;
	std::vector<int> idx_from(n*2);
	for (unsigned i=0; i<n; i++)
		idx_from[i] = 1000*i;
	for (unsigned i=n; i<n*2; i++)
		idx_from[i] = offset + idx_from[i-n];
	std::vector<int> idx_to(idx_from.size());
	for (unsigned i=0; i<idx_from.size(); i++)
		idx_to[i] = i;
	IS is_from, is_to;
	Vec sub_x_local;
	VecCreateSeq(PETSC_COMM_SELF, idx_from.size(), &sub_x_local);

	ISCreateGeneral(PETSC_COMM_SELF, idx_from.size(), &idx_from[0], PETSC_COPY_VALUES, &is_from);
	ISCreateGeneral(PETSC_COMM_SELF, idx_to.size(), &idx_to[0], PETSC_COPY_VALUES, &is_to);
	VecScatter scatter;
	VecScatterCreate(x, is_from, sub_x_local, is_to, &scatter);
	VecScatterBegin(scatter, x, sub_x_local, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, x, sub_x_local, INSERT_VALUES, SCATTER_FORWARD);
	const double *array_x;
	VecGetArrayRead(sub_x_local, &array_x);
	std::vector<double> vec_xi(idx_from.size());
	for (unsigned i=0; i<vec_xi.size(); i++)
		vec_xi[i] = array_x[i];
//	double x34 = array_x[0];
	VecRestoreArrayRead(x, &array_x);
	VecScatterBegin(scatter, f, sub_x_local, INSERT_VALUES, SCATTER_FORWARD);
	VecScatterEnd(scatter, f, sub_x_local, INSERT_VALUES, SCATTER_FORWARD);
	VecGetArrayRead(sub_x_local, &array_x);
	std::vector<double> vec_fi(idx_from.size());
	for (unsigned i=0; i<vec_fi.size(); i++)
		vec_fi[i] = array_x[i];
//	double f34 = array_x[0];
	VecRestoreArrayRead(x, &array_x);
	VecScatterDestroy(&scatter);
	ISDestroy(&is_from);
	ISDestroy(&is_to);
	VecDestroy(&sub_x_local);

	for (unsigned i=0; i<idx_from.size(); i++)
		ScreenMessage("-> row %d: x=%g, r=%g\n", idx_from[i], vec_xi[i], vec_fi[i]);
//	ScreenMessage("-> row 34: x=%g, r=%g\n", x34, f34);
//	VecView(x, PETSC_VIEWER_STDOUT_WORLD);
//	ScreenMessage("-> evaluated r\n");
//	VecView(f, PETSC_VIEWER_STDOUT_WORLD);

	double r_max, r_min;
	VecMax(f, NULL, &r_max);
	VecMin(f, NULL, &r_min);
	ScreenMessage("-> r: max=%g, min=%g\n", r_max, r_min);
#endif

	eqs_new->b = b_org;

#if 0
	// check if given x is properly set to node values
	Vec x2;
	VecDuplicate(x, &x2);
	pcs->copyNodalValuesToVec(x2);
	VecAXPY(x2, -1, x);
	double diff_norm;
	VecNorm(x2, NORM_2, &diff_norm);
	ScreenMessage("-> ||x1-x2||=%g\n", diff_norm);
#endif

	return 0;
}

#if (PETSC_VERSION_NUMBER >= 3050)
PetscErrorCode FormJacobianTH(SNES snes, Vec x, Mat jac_, Mat B_, void* ctx)
{
	(void)snes;
	Mat* jac = &jac_;
	Mat* B = &B_;
#else
PetscErrorCode FormJacobianTH(SNES /*snes*/, Vec x, Mat* jac, Mat* B,
                              MatStructure* flag, void* ctx)
{
#endif
	bool debugOutput = false;
	if (debugOutput) ScreenMessage("-> form a Jacobian matrix\n");
	CRFProcessTH* pcs = (CRFProcessTH*)ctx;
	petsc_group::PETScLinearSolver* eqs_new = pcs->eqs_new;
	MeshLib::CFEMesh* m_msh = pcs->m_msh;
	CFiniteElementStd* fem = pcs->GetAssembler();

	Mat* mat = jac;
	mat = B;
	eqs_new->B = eqs_new->A;
	eqs_new->A = *B;

	if (debugOutput) ScreenMessage("-> update nodal values\n");
	pcs->copyVecToNodalValues(x);
	//	pcs->IncorporateBoundaryConditions(-1, false, false, true); // set bc
	// node values
	//	if (debugOutput)
	//	ScreenMessage("-> update integration point values\n");
	//	pcs->CalIntegrationPointValue();

	const std::vector<Mat> vec_subA_org(eqs_new->vec_subA);
	if (eqs_new->vec_subA.size() > 0)
	{
		Mat** mats;
		MatNestGetSubMats(*mat, NULL, NULL, &mats);
		eqs_new->vec_subA[0] = mats[0][0];
		eqs_new->vec_subA[1] = mats[0][1];
		eqs_new->vec_subA[2] = mats[1][0];
		eqs_new->vec_subA[3] = mats[1][1];
		for (unsigned i = 0; i < eqs_new->vec_subA.size(); i++)
		{
			MatZeroEntries(eqs_new->vec_subA[i]);
		}
	}

	MatZeroEntries(*mat);

	const size_t dn = m_msh->ele_vector.size() / 10;
	const bool print_progress = (dn >= 100) && debugOutput;
	// if (debugOutput)
	if (print_progress)
		ScreenMessage("-> start local assembly for %d elements...\n",
		              m_msh->ele_vector.size());

	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		if (print_progress && (i + 1) % dn == 0) ScreenMessage("* ");
		// ScreenMessage("%d \%\n", ((i+1)*100/n_eles));
		MeshLib::CElem* elem = m_msh->ele_vector[i];
		// Marked for use //WX: modified for coupled excavation
		if (elem->GetMark())
		{
			elem->SetOrder(false);
			fem->ConfigElement(elem, false);
			fem->Assembly(true, false);
		}
	}
	//	if (debugOutput)
	if (print_progress) ScreenMessage("done\n");

	for (size_t i = 0; i < eqs_new->vec_subA.size(); i++)
	{
		MatAssemblyBegin(eqs_new->vec_subA[i], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(eqs_new->vec_subA[i], MAT_FINAL_ASSEMBLY);
	}
	MatAssemblyBegin(*mat, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*mat, MAT_FINAL_ASSEMBLY);

	if (debugOutput) ScreenMessage("-> impose Dirichlet BC\n");
	pcs->IncorporateBoundaryConditions(-1, true, false);

	//	MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);
	//	MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(*jac, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*jac, MAT_FINAL_ASSEMBLY);
	MatAssemblyBegin(*B, MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(*B, MAT_FINAL_ASSEMBLY);

#if (PETSC_VERSION_NUMBER < 3050)
	*flag = SAME_NONZERO_PATTERN;
#endif

	//	for (size_t i=0; i<eqs_new->vec_subA.size(); i++)
	//		MatView(eqs_new->vec_subA[i], PETSC_VIEWER_STDOUT_WORLD);

	eqs_new->A = eqs_new->B;
	for (unsigned i = 0; i < vec_subA_org.size(); i++)
		eqs_new->vec_subA[i] = vec_subA_org[i];

	return 0;
}

#endif
