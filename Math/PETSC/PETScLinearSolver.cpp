/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PETScLinearSolver.h"

#include <iostream>
#include <list>

#include <petsc-private/petscimpl.h>
#include <petscversion.h>
#include <petsctime.h>

#include "display.h"
#include "StringTools.h"

namespace petsc_group
{
PETScLinearSolver::PETScLinearSolver(const int size)
    : A(NULL),
      B(NULL),
      b(NULL),
      x(NULL),
      snes(NULL),
      lsolver(NULL),
      prec(NULL),
      total_x(NULL),
      global_x0(NULL),
      global_x1(NULL),
      global_buff(NULL)
{
	i_start = i_end = 0;
	ltolerance = 1.e-10;
	m_size = size;
	time_elapsed = 0.0;
	d_nz = 10;
	o_nz = 10;
	//nz = 10;
	m_size_loc = PETSC_DECIDE;
	mpi_size = 0;
	rank = 0;
	is_global_node_id = NULL;
	is_local_node_id = NULL;
}

PETScLinearSolver::~PETScLinearSolver()
{
	VecDestroy(&b);
	VecDestroy(&x);
	VecDestroy(&total_x);
	for (size_t i = 0; i < vec_subA.size(); i++)
		MatDestroy(&vec_subA[i]);
	MatDestroy(&A);
	if (snes)
		SNESDestroy(&snes);
	else if (lsolver)
		KSPDestroy(&lsolver);
	// if(prec) PCDestroy(&prec);

	if (vec_subA.empty())
	{
		for (size_t i = 0; i < vec_isg.size(); i++)
			ISDestroy(&vec_isg[i]);
	}
	ISDestroy(&is_global_node_id);
	ISDestroy(&is_local_node_id);

	if (global_x0) delete[] global_x0;
	if (global_x1) delete[] global_x1;
	if (global_buff) delete[] global_buff;

	PetscPrintf(PETSC_COMM_WORLD, "\n>>Number of Unknows: %d\n", m_size);
	PetscPrintf(PETSC_COMM_WORLD, ">>Elapsed time in linear solver: %fs\n",
	            time_elapsed);
}

void PETScLinearSolver::Init(SparseIndex* sparse_index)
{
	if (sparse_index)
	{
		d_nz = sparse_index->d_nz;
		o_nz = sparse_index->d_nz;
		//nz = sparse_index->nz;
		m_size_loc = sparse_index->m_size_loc;
	}

	CreateMatrix();

	global_x0 = new PetscScalar[m_size];
	global_x1 = new PetscScalar[m_size];
	global_buff = new PetscScalar[m_size];
}

void PETScLinearSolver::Config(const PetscReal tol, const PetscInt maxits,
                               const KSPType lsol, const PCType prec_type,
                               const std::string& misc_setting, const std::string& prefix)
{
	KSPCreate(PETSC_COMM_WORLD, &lsolver);
	ConfigLinear(tol, maxits, lsol, prec_type, misc_setting, prefix);
	KSPSetFromOptions(lsolver);
}

void PETScLinearSolver::ConfigWithNonlinear(const PetscReal tol,
                                            const PetscInt maxits,
                                            const KSPType lsol,
                                            const PCType prec_type,
                                            const std::string& misc_setting, const std::string& prefix)
{
	SNESCreate(PETSC_COMM_WORLD, &snes);
	SNESGetKSP(snes, &lsolver);
	ConfigLinear(tol, maxits, lsol, prec_type, misc_setting, prefix);
}

void PETScLinearSolver::ConfigLinear(const PetscReal tol, const PetscInt maxits,
                                     const KSPType lsol, const PCType prec_type,
                                     const std::string& misc_setting, const std::string& prefix)
{
	if (lsolver == NULL) return;

	ltolerance = tol;
	sol_type = lsol;
	pc_type = prec_type;

#if (PETSC_VERSION_NUMBER >= 3050)
	KSPSetOperators(lsolver, A, A);
#elif(PETSC_VERSION_NUMBER > 3040)
	KSPSetOperators(lsolver, A, A, SAME_NONZERO_PATTERN);
#else
	KSPSetOperators(lsolver, A, A, DIFFERENT_NONZERO_PATTERN);
#endif

	KSPSetType(lsolver, lsol);
	KSPSetTolerances(lsolver, ltolerance, PETSC_DEFAULT, PETSC_DEFAULT, maxits);

	KSPGetPC(lsolver, &prec);
	if (vec_subA.size() > 0)
	{
		PCSetType(prec, PCFIELDSPLIT);
		for (size_t i = 0; i < vec_isg.size(); i++)
		{
			PCFieldSplitSetIS(prec, number2str(i).c_str(), vec_isg[i]);
		}
	}
	else
	{
		PCSetType(prec, prec_type);  //  PCJACOBI); //PCNONE);
	}
	if (!prefix.empty())
	{
		KSPSetOptionsPrefix(lsolver, prefix.c_str());
		PCSetOptionsPrefix(prec, prefix.c_str());
	}

#if 1
	if (vec_subA.empty())
	{
		if (!misc_setting.empty())
		{
			PetscPrintf(PETSC_COMM_WORLD, "-> additional PETSc arguments:\n");
			std::list<std::string> lst = splitString(misc_setting, ' ');
			for (std::list<std::string>::iterator itr = lst.begin();
			     itr != lst.end();
			     ++itr)
			{
				// key-value or only key
				std::string& str1 = *itr;
				std::string val = "";
				if (str1.length() < 2 || str1.find('-') == std::string::npos) continue;
				++itr;
				if (itr != lst.end())
				{
					std::string& str2 = *itr;
					if (str2[0] != '-')
						val = str2;
					else
						--itr;
				}
				std::string key = "-" + prefix + str1.substr(1);
				vec_para.push_back(std::make_pair(key, val));
				PetscPrintf(PETSC_COMM_WORLD, "\t %s = %s\n", key.c_str(),
				            val.c_str());
				if (itr == lst.end())
					break;
			}
		}
		for (std::vector<Para>::iterator itr = vec_para.begin();
		     itr != vec_para.end();
		     ++itr)
		{
			PetscOptionsSetValue(itr->first.c_str(), itr->second.c_str());
		}
		char* copts;
		PetscOptionsGetAll(&copts);
		PetscPrintf(PETSC_COMM_WORLD, "-> PETSc options = %s\n", copts);
		PetscFree(copts);
	}
#endif
#if 0
   // reset options
   for (std::vector<Para>::iterator itr=vec_para.begin(); itr!=vec_para.end(); ++itr) {
	   PetscOptionsSetValue(itr->first.c_str(),"");
   }
#endif

}

void PETScLinearSolver::CreateMatrix()
{
	PetscErrorCode ierr;
	MatCreate(PETSC_COMM_WORLD, &A);

	ierr = MatSetSizes(A, m_size_loc, m_size_loc, PETSC_DECIDE, PETSC_DECIDE);
	//ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n_global_rows, n_global_cols);
	PetscInt M, N;
	MatGetSize(A, &M, &N);
	ScreenMessage("-> matrix: global nrows=%d, ncols=%d\n", M, N);
	PetscInt m, n;
	MatGetLocalSize(A, &m, &n);
	ScreenMessage2("-> matrix: local nrows=%d, ncols=%d\n", m, n);
	CHKERRCONTINUE(ierr);

	MatSetType(A, MATMPIAIJ);
	MatSetFromOptions(A);
#if 1
	ScreenMessage2("-> preallocate PETSc matrix with d_nz=%d and o_nz=%d\n", d_nz, o_nz);
	MatMPIAIJSetPreallocation(A, d_nz, PETSC_NULL, o_nz, PETSC_NULL);
	MatSeqAIJSetPreallocation(A, d_nz, PETSC_NULL);
	MatSetOption(A, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_TRUE);
	MatSetOption(A, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);

	MatInfo info;
	MatGetInfo(A, MAT_LOCAL, &info);
	ScreenMessage2("-> MatInfo: number of nonzeros=%g, allocated memory=%g MB\n", info.nz_allocated, info.memory/(1024*1024));
#else
	ScreenMessage("-> do not preallocate PETSc\n");
	MatSetUp(A);
#endif
	MatSetOption(A, MAT_KEEP_NONZERO_PATTERN, PETSC_TRUE);  // for MatZeroRows()
	MatGetOwnershipRange(A, &i_start, &i_end);
	ScreenMessage2d("-> PETSc linear solver range: start=%d, end=%d\n", i_start,
	                i_end);

	//  std::cout<<"sub_a  "<<i_start<<";   sub_d "<<i_end<<"\n";

	// PetscErrorCode ierr;  // returned value from PETSc functions
#if (PETSC_VERSION_NUMBER >= 3060)
	MatCreateVecs(A, &x, &b);
#else
	MatGetVecs(A, &x, &b);
#endif
	VecGetSize(x, &M);
	ScreenMessage("-> x: global nrows=%d\n", M);
	VecGetLocalSize(x, &m);
	ScreenMessage2("-> x: local nrows=%d\n", m);
	VecGetSize(b, &M);
	ScreenMessage("-> b: global nrows=%d\n", M);
	VecGetLocalSize(b, &m);
	ScreenMessage2("-> b: local nrows=%d\n", m);
//	VecCreate(PETSC_COMM_WORLD, &b);
//	////VecCreateMPI(PETSC_COMM_WORLD,m_size_loc, m, &b);
//	// VecSetSizes(b, m_size_loc, m);
//	VecSetSizes(b, PETSC_DECIDE, m);
//	VecSetFromOptions(b);
//	VecSetUp(b);  // kg44 for PETSC 3.3
//	VecDuplicate(b, &x);

	// VecGetOwnershipRange(b, &i_start,&i_end);

}

void PETScLinearSolver::getLocalRowColumnSizes(int* m, int* n)
{
	MatGetLocalSize(A, m, n);
}
void PETScLinearSolver::getOwnerRange(int* start_r, int* end_r)
{
	*start_r = i_start;
	*end_r = i_end;
}

void PETScLinearSolver::CheckIfMatrixIsSame(const std::string& filename)
{
	ScreenMessage(
	    "-> Check if the assembled matrix is the same as the one in a file: "
	    "%s\n",
	    filename.c_str());
	PetscErrorCode ierr;
	PetscViewer fd;
	ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD, filename.c_str(),
	                             FILE_MODE_READ, &fd);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);

	PETSc_Mat Y;
	ierr = MatCreate(PETSC_COMM_WORLD, &Y);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatSetFromOptions(Y);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	ierr = MatLoad(Y, fd);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);

	ierr = MatAXPY(Y, -1., A, SAME_NONZERO_PATTERN);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);
	PetscReal diff = .0;
	ierr = MatNorm(Y, NORM_1, &diff);
	CHKERRABORT(PETSC_COMM_WORLD, ierr);

	ScreenMessage("\t||A_assembled - A_file|| = %e\n", diff);
}

int PETScLinearSolver::Solver()
{
// TEST
#define TEST_MEM_PETSC
#ifdef TEST_MEM_PETSC
	PetscLogDouble mem1, mem2;
	PetscMemoryGetCurrentUsage(&mem1);
	PetscPrintf(PETSC_COMM_WORLD, "-> memory usage: %f MB\n", mem1/(1024*1024));
#endif

	/*
	//TEST
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x.txt", &viewer);
	PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB);
	PetscObjectSetName((PetscObject)x,"Solution");
	VecView(x, viewer);
	*/

	int its;
	PetscLogDouble v1, v2;
	KSPConvergedReason reason;

// #define PETSC34
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER >= 3040)
	PetscTime(&v1);
#else
	PetscGetTime(&v1);
#endif

	PetscPrintf(PETSC_COMM_WORLD,
	            "------------------------------------------------\n");
	PetscPrintf(PETSC_COMM_WORLD, "*** PETSc linear solver\n");
#if (PETSC_VERSION_NUMBER >= 3050)
	KSPSetOperators(lsolver, A, A);
#elif(PETSC_VERSION_NUMBER >= 3040)
	KSPSetOperators(lsolver, A, A, SAME_NONZERO_PATTERN);
#else
	KSPSetOperators(lsolver, A, A, DIFFERENT_NONZERO_PATTERN);
#endif
	KSPSolve(lsolver, b, x);

	const char* slv_type;
	const char* prc_type;
	KSPGetType(lsolver, &slv_type);
	PCGetType(prec, &prc_type);
	KSPGetConvergedReason(lsolver, &reason);  // CHKERRQ(ierr);
	KSPGetIterationNumber(lsolver, &its);     // CHKERRQ(ierr);
	PetscReal rtol = .0, abstol = .0, dtol = .0;
	PetscInt maxits = 0;
	KSPGetTolerances(lsolver, &rtol, &abstol, &dtol, &maxits);
	PetscReal r_norm = .0;
	KSPGetResidualNorm(lsolver, &r_norm);
//	PetscReal b_norm = .0;
//	VecNorm(b, NORM_2, &b_norm);
//	PetscReal error_r = r_norm / b_norm;

	PetscPrintf(PETSC_COMM_WORLD, "solver    : %s\n", slv_type);
	PetscPrintf(PETSC_COMM_WORLD, "precon    : %s\n", prc_type);
	PetscPrintf(PETSC_COMM_WORLD, "iteration : %d/%d\n", its, maxits);
//	PetscPrintf(PETSC_COMM_WORLD, "||b||     : %e\n", b_norm);
	PetscPrintf(PETSC_COMM_WORLD,
	            "residual  : ||r||=%e, rel. tol=%e\n", r_norm, rtol);
//	PetscPrintf(PETSC_COMM_WORLD,
//	            "residual  : ||r||=%e, ||r||/||b||=%e(approx.), rtol=%e\n",
//	            r_norm, error_r, rtol);
	if (reason >= 0)
	{
		PetscPrintf(PETSC_COMM_WORLD, "status    : Converged (reason=%d)\n",
		            reason);
	}
	else
	{
		if (reason == KSP_DIVERGED_INDEFINITE_PC)
		{
			PetscPrintf(PETSC_COMM_WORLD,
			            "status    : Diverged (indefinite precon)\n", reason);
			PetscPrintf(PETSC_COMM_WORLD,
			            "            Run the executable again but with "
			            "-pc_factor_shift_positive_definite option.\n");
		}
		else if (reason == KSP_DIVERGED_ITS)
		{
			PetscPrintf(PETSC_COMM_WORLD,
			            "status    : Diverged (max iteration)\n", reason);
		}
		else if (reason == KSP_DIVERGED_BREAKDOWN)
		{
			PetscPrintf(PETSC_COMM_WORLD, "status    : Diverged (breakdown)\n",
			            reason);
		}
		else
		{
			PetscPrintf(PETSC_COMM_WORLD, "status    : Diverged (reason=%d)\n",
			            reason);
		}
		//      PetscFinalize();
		//      exit(1);
		return -1;
	}

	PetscPrintf(PETSC_COMM_WORLD,
	            "------------------------------------------------\n");

// VecAssemblyBegin(x);
// VecAssemblyEnd(x);

// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER >= 3040)
	PetscTime(&v2);
#else
	PetscGetTime(&v2);
#endif

	time_elapsed += v2 - v1;

#define aTEST_OUT
#ifdef TEST_OUT
	// TEST
	PetscViewer viewer;
	PetscViewerASCIIOpen(PETSC_COMM_WORLD, "x2.txt", &viewer);
	PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);
	PetscObjectSetName((PetscObject)A, "Matrix");
	MatView(A, viewer);
	PetscObjectSetName((PetscObject)x, "Solution");
	VecView(x, viewer);
	PetscObjectSetName((PetscObject)b, "RHS");
	VecView(b, viewer);
	VecDestroy(&b);
	VecDestroy(&x);
	MatDestroy(&A);
	if (lsolver) KSPDestroy(&lsolver);
	// if(prec) PCDestroy(&prec);
	if (global_x0) delete[] global_x0;
	if (global_x1) delete[] global_x1;
	PetscFinalize();
	exit(0);
#endif

#ifdef TEST_MEM_PETSC
	PetscMemoryGetCurrentUsage(&mem2);
	PetscPrintf(PETSC_COMM_WORLD, "-> memory usage: %f MB\n", mem1/(1024*1024));
#endif
#undef TEST_MEM_PETSC
	return its;
}

void PETScLinearSolver::AssembleRHS_PETSc(bool assemble_subvec)
{
	if (assemble_subvec)
	{
		for (size_t i = 0; i < vec_subRHS.size(); i++)
			VecRestoreSubVector(b, vec_isg[i], &vec_subRHS[i]);
	}

	VecAssemblyBegin(b);
	VecAssemblyEnd(b);
}
void PETScLinearSolver::AssembleUnkowns_PETSc()
{
	VecAssemblyBegin(x);
	VecAssemblyEnd(x);
}
void PETScLinearSolver::AssembleMatrixPETSc(const MatAssemblyType type)
{
	for (size_t i = 0; i < vec_subA.size(); i++)
	{
		MatAssemblyBegin(vec_subA[i], MAT_FINAL_ASSEMBLY);
		MatAssemblyEnd(vec_subA[i], MAT_FINAL_ASSEMBLY);
	}

	MatAssemblyBegin(A, type);
	MatAssemblyEnd(A, type);
}

void PETScLinearSolver::getGlobalVectorArray(Vec& vec, PetscScalar* u1)
{
	PetscInt local_size;
	VecGetLocalSize(vec, &local_size);

	PetscScalar* xp = NULL;
	VecGetArray(vec, &xp);

	// number of elements to be sent for each rank
	int size_rank;
	MPI_Comm_size(PETSC_COMM_WORLD, &size_rank);
	std::vector<PetscInt> i_cnt(size_rank);

	MPI_Allgather(&local_size, 1, MPI_INT, &i_cnt[0], 1, MPI_INT,
	              PETSC_COMM_WORLD);

	// collect local array
	PetscInt offset = 0;
	// offset in the receive vector of the data from each rank
	std::vector<PetscInt> i_disp(size_rank);
	for (PetscInt i = 0; i < size_rank; i++)
	{
		i_disp[i] = offset;
		offset += i_cnt[i];
	}

	MPI_Allgatherv(xp, local_size, MPI_DOUBLE, u1, &i_cnt[0], &i_disp[0],
	               MPI_DOUBLE, PETSC_COMM_WORLD);

	VecRestoreArray(vec, &xp);
}

void PETScLinearSolver::UpdateSolutions(PetscScalar* u0, PetscScalar* u1)
{
#ifdef TEST_MEM_PETSC
	// TEST
	PetscLogDouble mem1, mem2;
	PetscMemoryGetCurrentUsage(&mem1);
#endif

	int i, j;
	PetscScalar* xp;

	int receivecount;
	PetscInt low, high, otherlow;
	MPI_Status status;
	PetscInt count;
	int tag = 9999;
	VecGetOwnershipRange(x, &low, &high);
	VecGetLocalSize(x, &count);

	VecGetArray(x, &xp);
	for (i = 0; i < count; i++)
		u1[i] = xp[i];

	// double *global_buff = new double[m_size];

	// Collect solution from processes.
	for (j = 0; j < count; j++)
		global_buff[low + j] = u1[j];
	for (i = 0; i < mpi_size; i++)
	{
		if (i != rank)
		{
			MPI_Sendrecv(&count, 1, MPI_INT, i, tag, &receivecount, 1, MPI_INT,
			             i, tag, PETSC_COMM_WORLD, &status);
			MPI_Sendrecv(&low, 1, MPI_INT, i, tag, &otherlow, 1, MPI_INT, i,
			             tag, PETSC_COMM_WORLD, &status);
			MPI_Sendrecv(u1, count, MPI_DOUBLE, i, tag, u0, receivecount,
			             MPI_DOUBLE, i, tag, PETSC_COMM_WORLD, &status);
			for (j = 0; j < receivecount; j++)
				global_buff[otherlow + j] = u0[j];
		}
	}

	// MPI_Barrier(PETSC_COMM_WORLD);
	// Copy the collected solution to the array for the previous solution
	for (i = 0; i < m_size; i++)
	{
		u1[i] = global_buff[i];
		u0[i] = global_buff[i];
	}

	// delete [] global_buff;

	VecRestoreArray(x, &xp);

// TEST
#ifdef TEST_MEM_PETSC
	PetscMemoryGetCurrentUsage(&mem2);
	PetscPrintf(
	    PETSC_COMM_WORLD,
	    "### Memory usage by Updating. Before :%f After:%f Increase:%d\n", mem1,
	    mem2, (int)(mem2 - mem1));
#endif
}

void PETScLinearSolver::MappingSolution()
{
	UpdateSolutions(global_x0, global_x1);
}

int PETScLinearSolver::GetLocalSolution(PetscScalar* x_l)
{
	PetscInt count;
	VecGetLocalSize(x, &count);

	VecGetArray(x, &x_l);

	return count;
}

int PETScLinearSolver::GetLocalRHS(PetscScalar* rhs_l)
{
	PetscInt count;
	VecGetLocalSize(b, &count);

	VecGetArray(b, &rhs_l);

	return count;
}

double* PETScLinearSolver::GetGlobalSolution() const
{
	return global_x1;
}

/*!
  Get values of the specified elements from a global vector

  @param v_type - Indicator for vector: 0: x; 1: rhs
  @param ni 	- number of elements to get
  @param ix 	- indices where to get them from (in global 1d numbering)
*/
void PETScLinearSolver::GetVecValues(const int v_type, PetscInt ni,
                                     const PetscInt ix[], PetscScalar y[]) const
{
	if (v_type == 0)
		VecGetValues(x, ni, ix, y);
	else
		VecGetValues(b, ni, ix, y);
}

/*!
    Get norm of RHS
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i|
    06.2012. WW
*/
PetscReal PETScLinearSolver::GetVecNormRHS(NormType nmtype)
{
	PetscReal norm = 0.;
	VecNorm(b, nmtype, &norm);
	return norm;
}
/*!
    Get norm of x
    @param nmtype  - norm type
                     NORM_1 denotes sum_i |x_i|
                     NORM_2 denotes sqrt(sum_i (x_i)^2)
                     NORM_INFINITY denotes max_i |x_i|
    06.2012. WW
*/
PetscReal PETScLinearSolver::GetVecNormX(NormType nmtype)
{
	PetscReal norm = 0.;
	VecNorm(x, nmtype, &norm);
	return norm;
}

void PETScLinearSolver::RestoreLocalSolutionArray(PetscScalar* x_l)
{
	VecRestoreArray(x, &x_l);
}
void PETScLinearSolver::RestoreLocalRHSArray(PetscScalar* rhs_l)
{
	VecRestoreArray(b, &rhs_l);
}

void PETScLinearSolver::set_bVectorEntry(const int i, const double value)
{
	VecSetValues(b, 1, &i, &value, INSERT_VALUES);
}
void PETScLinearSolver::set_xVectorEntry(const int i, const double value)
{
	VecSetValues(x, 1, &i, &value, INSERT_VALUES);
}

void PETScLinearSolver::setArrayValues(int arr_idx, PetscInt ni,
                                       const PetscInt ix[],
                                       const PetscScalar y[], InsertMode iora)
{
	if (arr_idx == 0)
		VecSetValues(x, ni, ix, y, iora);
	else if (arr_idx == 1)
		VecSetValues(b, ni, ix, y, iora);
}

void PETScLinearSolver::add_bVectorEntry(const int i, const double value,
                                         InsertMode mode)
{
	VecSetValue(b, i, value, mode);
}
void PETScLinearSolver::add_xVectorEntry(const int i, const double value,
                                         InsertMode mode)
{
	VecSetValue(x, i, value, mode);
}

void PETScLinearSolver::Initialize()
{
	VecSet(b, 0.0);
	VecSet(x, 0.0);
	for (unsigned i = 0; i < vec_subA.size(); i++)
	{
		MatZeroEntries(vec_subA[i]);
	}
	MatZeroEntries(A);
#if 0
   for (unsigned i=0; i<vec_subRHS.size(); i++) {
      VecGetSubVector(b, vec_isg[i], &vec_subRHS[i]);
      int gstart, gend, start;
      PetscBool contiguous;
      VecGetOwnershipRange(vec_subRHS[i],&gstart,&gend);
      ISContiguousLocal(vec_isg[i],gstart,gend,&start,&contiguous);
      ScreenMessage2("-> sub Vec[%d]: start=%d, end=%d, contiguous=%s\n", i, gstart, gend, contiguous ? "true" : "false");
   }
#endif
}

void PETScLinearSolver::addMatrixEntry(const int i, const int j,
                                       const double value)
{
	MatSetValue(A, i, j, value, ADD_VALUES);
}

void PETScLinearSolver::addMatrixEntries(const int m, const int idxm[],
                                         const int n, const int idxn[],
                                         const PetscScalar v[])
{
	int err = MatSetValues(A, m, idxm, n, idxn, v, ADD_VALUES);
	CHKERRABORT(PETSC_COMM_WORLD, err);
}

void PETScLinearSolver::zeroRows_in_Matrix(const int nrows,
                                           const PetscInt* rows)
{
	PetscScalar one = 1.0;
	// Each process indicates only rows it owns that are to be zeroed
	// MatSetOption(A, MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE);
	if (nrows > 0)
		MatZeroRows(A, nrows, rows, one, PETSC_NULL, PETSC_NULL);
	else
		MatZeroRows(A, 0, PETSC_NULL, one, PETSC_NULL, PETSC_NULL);
}

void PETScLinearSolver::EQSV_Viewer(const std::string& file_name, bool ascii)
{
	if (ascii)
	{
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
		PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

		AssembleRHS_PETSc();
		AssembleUnkowns_PETSc();
		AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);

//	  PetscBool flg = PETSC_FALSE;
//	  PetscOptionsGetBool(((PetscObject) A)->prefix, "-mat_ascii_output_large",
//&flg,NULL);
//	  ScreenMessage2("A: mat_ascii_output_large found: %s\n", flg==PETSC_TRUE ?
//"true" : "false");

// PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
#if 0
	  PetscObjectSetName((PetscObject)A,"Stiffness_matrix");
	  MatView(A,viewer);
#if 1
	  for (size_t i=0; i<vec_subA.size(); i++) {
		  std::string name = "SubMatrix" + number2str(i);
//		  PetscOptionsGetBool(((PetscObject) vec_subA[i])->prefix, "-mat_ascii_output_large", &flg,NULL);
//		  ScreenMessage2("subA[%d]: mat_ascii_output_large found: %s\n", i, flg==PETSC_TRUE ? "true" : "false");
		  PetscObjectSetName((PetscObject)vec_subA[i],name.c_str());
		  MatView(vec_subA[i],viewer);
	  }
#endif
#endif
		for (size_t i = 0; i < vec_subRHS.size(); i++)
		{
			std::string name = "SubVector" + number2str(i);
			VecGetSubVector(b, vec_isg[i], &vec_subRHS[i]);
			PetscObjectSetName((PetscObject)vec_subRHS[i], name.c_str());
			VecView(vec_subRHS[i], viewer);
			VecRestoreSubVector(b, vec_isg[i], &vec_subRHS[i]);
		}
		for (size_t i = 0; i < vec_subRHS.size(); i++)
		{
			std::string name = "SubVectorX" + number2str(i);
			VecGetSubVector(x, vec_isg[i], &vec_subRHS[i]);
			PetscObjectSetName((PetscObject)vec_subRHS[i], name.c_str());
			VecView(vec_subRHS[i], viewer);
			VecRestoreSubVector(b, vec_isg[i], &vec_subRHS[i]);
		}
		if (vec_subRHS.empty())
		{
			PetscObjectSetName((PetscObject)b, "RHS");
			PetscObjectSetName((PetscObject)x, "Solution");
			VecView(b, viewer);
			VecView(x, viewer);
		}

//#define  EXIT_TEST
#ifdef EXIT_TEST
		VecDestroy(&b);
		VecDestroy(&x);
		MatDestroy(&A);
		if (lsolver) KSPDestroy(&lsolver);
		// if(prec) PCDestroy(&prec);
		if (global_x0) delete[] global_x0;
		if (global_x1) delete[] global_x1;
		PetscFinalize();
		exit(0);
#endif
	}
	else
	{
		AssembleRHS_PETSc();
		AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);

		{
			PetscViewer viewer;
			std::string fnameA = file_name + "_eqs_A.dat";
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, fnameA.c_str(),
			                      FILE_MODE_WRITE, &viewer);
			MatView(A, viewer);
			PetscViewerDestroy(&viewer);
		}

		{
			PetscViewer viewer;
			std::string fnameRHS = file_name + "_eqs_rhs.dat";
			PetscViewerBinaryOpen(PETSC_COMM_WORLD, fnameRHS.c_str(),
			                      FILE_MODE_WRITE, &viewer);
			VecView(b, viewer);
			PetscViewerDestroy(&viewer);
		}
	}
}

void PETScLinearSolver::Residual_Viewer(const std::string& file_name,
                                        bool ascii)
{
	if (ascii)
	{
		PetscViewer viewer;
		PetscViewerASCIIOpen(PETSC_COMM_WORLD, file_name.c_str(), &viewer);
		PetscViewerPushFormat(viewer, PETSC_VIEWER_ASCII_MATLAB);

		AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);

		// residual vector
		Vec Br, v, w;
		MatGetVecs(A, &w, &v);
		KSPBuildResidual(lsolver, v, w, &Br);

		// PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_VTK);
		PetscObjectSetName((PetscObject)Br, "Residual");
		VecView(Br, viewer);
	}
}

}  // end of namespace
