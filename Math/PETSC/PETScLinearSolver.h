/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*!
   \brief Declaration of class PETScLinearSolver

   11.2011. WW

*/
#ifndef PETSC_LSOLVER_INC
#define PETSC_LSOLVER_INC

#include <string>
#include <vector>
// kg needed for memcpy in petsc libs
#include <stdio.h>
#include <string.h>

#include "petscmat.h"
#include "petscksp.h"
#include "petscsnes.h"

#ifdef USEPETSC34
#include "petsctime.h"
#endif

typedef Mat PETSc_Mat;
typedef Vec PETSc_Vec;

namespace petsc_group
{

struct SparseIndex
{
	int d_nz = 0;
	int o_nz = 0;
	//int nz = 0;
	int m_size_loc = 0;
};


class PETScLinearSolver
{
public:
	PETScLinearSolver(const int size);
	~PETScLinearSolver();

	void Config(const PetscReal tol, const PetscInt maxits, const KSPType lsol,
	            const PCType prec_type, const std::string& misc_setting, const std::string& prefix);
	void ConfigWithNonlinear(const PetscReal tol, const PetscInt maxits,
	                         const KSPType lsol, const PCType prec_type,
	                         const std::string& misc_setting, const std::string& prefix);
	void ConfigLinear(const PetscReal tol, const PetscInt maxits,
	                  const KSPType lsol, const PCType prec_type,
	                  const std::string& misc_setting, const std::string& prefix);

	void Init(SparseIndex* sparse_index);

	void CheckIfMatrixIsSame(const std::string& filename);
	int Solver();
	void AssembleRHS_PETSc(bool assemble_subvec = false);
	void AssembleUnkowns_PETSc();
	void AssembleMatrixPETSc(
	    const MatAssemblyType type = MAT_FINAL_ASSEMBLY);  // MAT_FLUSH_ASSEMBLY

	void getGlobalVectorArray(Vec& vec, PetscScalar* u1);
	void UpdateSolutions(PetscScalar* u0, PetscScalar* u1);
	void MappingSolution();
	int GetLocalSolution(PetscScalar* x_l);
	int GetLocalRHS(PetscScalar* rhs_l);
	double* GetGlobalSolution() const;
	void GetVecValues(const int v_type, PetscInt ni, const PetscInt ix[],
	                  PetscScalar y[]) const;
	PetscReal GetVecNormRHS(NormType nmtype = NORM_2);
	PetscReal GetVecNormX(NormType nmtype = NORM_2);

	void RestoreLocalSolutionArray(PetscScalar* x_l);
	void RestoreLocalRHSArray(PetscScalar* rhs_l);
	void getLocalRowColumnSizes(int* m, int* n);
	void getOwnerRange(int* start_r, int* end_r);

	int Size() const { return m_size; }

	void set_xVectorEntry(const int i, const double value);
	void set_bVectorEntry(const int i, const double value);
	void setArrayValues(int arr_idx, PetscInt ni, const PetscInt ix[],
	                    const PetscScalar y[], InsertMode iora = ADD_VALUES);

	void add_xVectorEntry(const int i, const double value, InsertMode mode);
	void add_bVectorEntry(const int i, const double value, InsertMode mode);
	void addMatrixEntry(const int i, const int j, const double value);
	void addMatrixEntries(const int m, const int idxm[], const int n,
	                      const int idxn[], const PetscScalar v[]);

	void Initialize();

	void zeroRows_in_Matrix(const int nrow, const PetscInt* rows);
	void zeroMatrix() { MatZeroEntries(A); }

	void set_rank_size(const int m_rank, const int size)
	{
		mpi_size = size;
		rank = m_rank;
	}

	PetscInt getStartRow() const { return i_start; }
	PetscInt getEndRow() const { return i_end; }

	PetscInt getMPI_Size() const { return mpi_size; }
	PetscInt getMPI_Rank() const { return rank; }

	void EQSV_Viewer(const std::string& file_name, bool ascii = true);
	void Residual_Viewer(const std::string& file_name, bool ascii = true);

public:
	//  private:
	PETSc_Mat A;
	PETSc_Mat B;
	PETSc_Vec b;
	PETSc_Vec x;
	SNES snes;
	KSP lsolver;
	PC prec;
	std::vector<Mat> vec_subA; /* the four blocks */
	std::vector<IS> vec_isg;   /* index sets of split "0" and "1" */
	std::vector<Vec> vec_subRHS;
	PETSc_Vec total_x;

	PetscInt i_start;
	PetscInt i_end;

	PetscScalar* global_x0;
	PetscScalar* global_x1;
	PetscScalar* global_buff;

	// Slover and preconditioner names, only for log
	std::string sol_type;
	std::string pc_type;

	PetscLogDouble time_elapsed;

	PetscInt m_size;
	PetscInt m_size_loc;
	float ltolerance;
	// Number of nonzeros per row in DIAGONAL portion of
	// local submatrix (same value is used for all local rows)
	PetscInt d_nz;
	// Number of nonzeros per row in the OFF-DIAGONAL portion of
	// local submatrix (same value is used for all local rows).
	PetscInt o_nz;
	// Number of nonzeros per row (same for all rows)
	//PetscInt nz;

	int mpi_size;
	int rank;

	typedef std::pair<std::string, std::string> Para;
	std::vector<Para> vec_para;

	void CreateMatrix();

	IS is_global_node_id, is_local_node_id;
};

// extern std::vector<PETScLinearSolver*> EQS_Vector;
}  // end namespace
#endif
