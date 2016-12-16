/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "equation_class.h"

#include <algorithm>
#include <cfloat>
#include <iomanip>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef LIS
#ifndef LIS_INT
#define LIS_INT int
#endif
#endif

#if defined(_WIN32) || defined(_WIN64)
#define pardiso_ PARDISO
#else
#ifdef OGS_USE_LONG
#define PARDISO pardiso_64
#else
#define PARDISO pardiso_
#endif
#endif

#ifdef MKL
#ifdef _WIN32
/* PARDISO prototype. */
extern "C" int PARDISOINIT(void*, int*, int*, int*, double*, int*);
extern "C" int PARDISO(void*, int*, int*, int*, int*, int*, double*, int*, int*,
                       int*, int*, int*, int*, double*, double*, int*, double*);
#else
#include "mkl.h"

/* PARDISO prototype. */
//#define PARDISO pardiso_
extern int omp_get_max_threads();
extern int PARDISO(int*, int*, int*, int*, int*, int*, double*, int*, int*,
                   int*, int*, int*, int*, double*, double*, int*);
#endif
#endif

#include "Configure.h"
#include "display.h"
#include "makros.h"

#include "matrix_class.h"
#include "rf_num_new.h"
#include "rf_pcs.h"

std::vector<Math_Group::Linear_EQS*> EQS_Vector;
using namespace std;

//
namespace Math_Group
{
/**************************************************************************
   Task: Linear equation::Constructor
   Programing:
   10/2007 WW/
**************************************************************************/
Linear_EQS::Linear_EQS(const SparseTable& sparse_table,
					   const long dof, bool /*messg*/)
{
	long i;
	A = new CSparseMatrix(sparse_table, dof);
	size_A = A->Dim();

	x = new double[size_A];
	b = new double[size_A];
	//
	for (i = 0; i < size_A; i++)
	{
		x[i] = 0.;
		b[i] = 0.;
	}
}

/**************************************************************************
   Task: Linear equation::Destructor
   Programing:
   10/2007 WW/
**************************************************************************/
Linear_EQS::~Linear_EQS()
{
	if (A) delete A;
	if (x) delete[] x;
	if (b) delete[] b;
	//
	A = NULL;
	x = NULL;
	b = NULL;
}
/**************************************************************************
   Task: Linear equation::
   Programing:
   10/2007 WW/
**************************************************************************/
void Linear_EQS::ConfigNumerics(CNumerics* /*m_num*/, const long /*n*/)
{
}
/**************************************************************************
   Task: Linear equation::Alocate memory for solver
   Programing:
   11/2007 WW/
**************************************************************************/
void Linear_EQS::Initialize()
{
	if (A) (*A) = 0.;
	for (long i = 0; i < size_A; i++)
		b[i] = 0.;
}
/**************************************************************************
   Task: Linear equation::Alocate memory for solver
   Programing:
   11/2007 WW/
**************************************************************************/
void Linear_EQS::Clean()
{
}
/**************************************************************************
   Task: Linear equation::Write
   Programing:
   11/2007 WW/
**************************************************************************/
void Linear_EQS::Write(std::ostream& os)
{
	A->Write(os);
	//
	os << " b ( RHS): " << endl;
	os.width(10);
	os.precision(6);
	//
	for (long i = 0; i < size_A; i++)
		os << setw(10) << i << " " << setw(15) << b[i] << endl;

	os << " x : " << endl;
	for (long i = 0; i < size_A; i++)
		os << setw(10) << i << " " << setw(15) << x[i] << endl;
}
/**************************************************************************
   Task: Linear equation::Write
   Programing:
   07/2010 NW
**************************************************************************/
void Linear_EQS::WriteRHS(ostream& os)
{
	os.width(10);
	os.precision(6);
	//
	for (long i = 0; i < A->Dim(); i++)
		os << setw(15) << b[i] << endl;
}
//**************************************************************************
/*!
    \brief Write the equation into a binary file

     Programing:
     03/2011 WW
 */
//**************************************************************************
void Linear_EQS::Write_BIN(ostream& os)
{
	if ((A->GetStorageType() != CRS) || (!A)) return;

	A->Write_BIN(os);
	os.write((char*)b, A->Dim() * sizeof(double));
}

/**************************************************************************
   Task: Linear equation::Write
   Programing:
   07/2010 NW
**************************************************************************/
void Linear_EQS::WriteX(ostream& os)
{
	os.width(10);
	os.precision(6);
	//
	for (long i = 0; i < A->Dim(); i++)
		os << setw(15) << x[i] << endl;
}

/**************************************************************************
   Task: Linear equation::Solver
   Programing:

   PARDISO openmp-paralle direct solver: 805

   LIS matrix solver options
   CG -i {cg|1}
   BiCG -i {bicg|2}
   CGS -i {cgs|3}
   BiCGSTAB -i {bicgstab|4}
   BiCGSTAB(l) -i {bicgstabl|5} -ell [2] Value for l
   GPBiCG -i {gpbicg|6}
   TFQMR -i {tfqmr|7}
   Orthomin(m) -i {orthomin|8} -restart [40] Value for Restart m
   GMRES(m) -i {gmres|9} -restart [40] Value for Restart m
   Jacobi -i {jacobi|10}
   Gauss-Seidel -i {gs|11}
   SOR -i {sor|12} -omega [1.9] Value for Relaxation Coefficient  (0 <  < 2)
   BiCGSafe -i {bicgsafe|13}
   CR -i {cr|14}
   BiCR -i {bicr|15}
   CRS -i {crs|16}
   BiCRSTAB -i {bicrstab|17}
   GPBiCR -i {gpbicr|18}
   BiCRSafe -i {bicrsafe|19}
   FGMRES(m) -i {fgmres|20} -restart [40] Value for Restart m
   IDR(s) -i {idrs|21} -restart [40] Value for Restart s

   Preconditioner Option Auxiliary Option
   None -p {none|0}
   Jacobi -p {jacobi|1}
   ILU(k) -p {ilu|2} -ilu_fill [0] Fill level k
   SSOR -p {ssor|3} -ssor_w [1.0] Relaxation Coefficient  (0 <  < 2)
   Hybrid -p {hybrid|4} -hybrid_i [sor] Iterative method
   -hybrid_maxiter [25] Maximum number of iterations
   -hybrid_tol [1.0e-3] Convergence criteria
   -hybrid_w [1.5] Relaxation Coefficient  for
   the SOR method (0 <  < 2)
   -hybrid_ell [2] Value for l of the BiCGSTAB(l) method
   -hybrid_restart [40] Restart values for GMRES and Orthomin
   I+S -p {is|5} -is_alpha [1.0] Parameter ?for preconditioner
   of a I + ?(m) type
   -is_m [3] Parameter m for preconditioner
   of a I + ?(m) type
   SAINV -p {sainv|6} -sainv_drop [0.05] Drop criteria
   SA-AMG -p {saamg|7} -saamg_unsym [false] Selection of asymmetric version
   Crout ILU -p {iluc|8} -iluc_drop [0.05] Drop criteria
   -iluc_rate [5.0] Ratio of Maximum fill-in
   ILUT -p {ilut|9} -ilut_drop [0.05] Drop criteria
   -ilut_rate [5.0] Ratio of Maximum fill-in
   additive Schwarz -adds true -adds_iter [1] Number of iterations

   02/2008 PCH OpenMP parallelization by LIS
   03/2009 PCH Solver type and precondition options added for .num file
**************************************************************************/
#if defined(USE_MPI)
#else  // if defined(USE_MPI)

#if defined(LIS) || defined(MKL)
Linear_EQS::IndexType Linear_EQS::searcgNonZeroEntries(
    IndexType nrows, IndexType* ptr, double* value,
    std::vector<IndexType>& vec_nz_rows, std::vector<IndexType>& vec_z_rows)
{
	vec_nz_rows.reserve(nrows);
	IndexType n_nz_entries = 0;
	for (LIS_INT i = 0; i < nrows; i++)
	{
		IndexType const j_row_begin = ptr[i];
		IndexType const j_row_end = ptr[i + 1];
		IndexType const old_n_nz_entries = n_nz_entries;
		for (IndexType j = j_row_begin; j < j_row_end; j++)
		{
			if (value[j] == .0) continue;
			n_nz_entries++;
		}
		if (n_nz_entries > old_n_nz_entries)
			vec_nz_rows.push_back(i);
		else
			vec_z_rows.push_back(i);
	}
	return n_nz_entries;
}

void Linear_EQS::compressCRS(const IndexType* org_ptr,
                             const IndexType* org_col_idx,
                             const double* org_value,
                             const std::vector<IndexType>& vec_nz_rows,
                             const std::vector<IndexType>& vec_z_rows,
                             IndexType n_nz_entries, IndexType*& new_ptr,
                             IndexType*& new_col_index, double*& new_value)
{
	IndexType n_new_rows = vec_nz_rows.size();
	// lis_matrix_malloc_crs(n_new_rows, n_nz_entries, &new_ptr, &new_col_index,
	// &new_value);
	new_value = new double[n_nz_entries];
	new_ptr = (IndexType*)malloc((vec_nz_rows.size() + 1) * sizeof(IndexType));
	new_col_index = (IndexType*)malloc((n_nz_entries) * sizeof(IndexType));

	LIS_INT nnz_counter = 0;
	for (LIS_INT i = 0; i < n_new_rows; i++)
	{
		const LIS_INT old_i = vec_nz_rows[i];
		const LIS_INT j_row_begin = org_ptr[old_i];
		const LIS_INT j_row_end = org_ptr[old_i + 1];

		new_ptr[i] = nnz_counter;

		for (LIS_INT j = j_row_begin; j < j_row_end; j++)
		{
			if (org_value[j] == .0) continue;
			// count how many columns are excluded before current
			LIS_INT offset_col = 0;
			for (size_t k = 0; k < vec_z_rows.size(); k++)
			{
				if (vec_z_rows[k] < org_col_idx[j]) offset_col++;
			}
			new_col_index[nnz_counter] = org_col_idx[j] - offset_col;
			new_value[nnz_counter] = org_value[j];
			nnz_counter++;
		}
	}
	new_ptr[vec_nz_rows.size()] = nnz_counter;
	// assert(n_nz_entries==nnz_counter);
}
#endif

#ifdef MKL
void Linear_EQS::solveWithPARDISO(CNumerics* num, bool compress_if_possible)
{
	ScreenMessage2(
	    "------------------------------------------------------------------\n");
	ScreenMessage2("*** PARDISO solver computation\n");

	// Prepare CRS data
	_INTEGER_t nonzero = A->nnz();
	_INTEGER_t nrows = A->Size() * A->Dof();
	double* value = new double[nonzero];
	A->GetCRSValue(value);
	_INTEGER_t* ptr = A->ptr;
	_INTEGER_t* index = A->col_idx;

	double* tmp_x = x;
	double* tmp_b = b;

	bool is_compressed = false;
	std::vector<IndexType> vec_nz_rows;
	if (compress_if_possible)
	{
		// check non-zero rows, non-zero entries
		ScreenMessage2("-> Check non-zero entries\n");
		std::vector<_INTEGER_t> vec_z_rows;
		IndexType n_nz_entries =
		    searcgNonZeroEntries(nrows, ptr, value, vec_nz_rows, vec_z_rows);

		if (vec_nz_rows.size() != (std::size_t)nrows)
		{
			ScreenMessage2(
			    "-> Found %d empty rows. Delete them from total %d rows... \n",
			    nrows - vec_nz_rows.size(), nrows);
			is_compressed = true;
			double* new_value;
			_INTEGER_t* new_ptr, *new_index;
			compressCRS(ptr, index, value, vec_nz_rows, vec_z_rows,
			            n_nz_entries, new_ptr, new_index, new_value);

			// update
			nrows = (int)vec_nz_rows.size();
			nonzero = n_nz_entries;
			ptr = new_ptr;
			index = new_index;
			delete[] value;
			value = new_value;

			// compress also x and RHS
			tmp_x = new double[nrows];
			tmp_b = new double[nrows];
#pragma omp parallel for
			for (int i = 0; i < nrows; ++i)
			{
				tmp_x[i] = x[vec_nz_rows[i]];
				tmp_b[i] = b[vec_nz_rows[i]];
			}
		}
	}

//#define PARDISO_INDEX_FROM_ONE
#ifdef PARDISO_INDEX_FROM_ONE  // in case PARDISO requires index start from one
	const bool is_index_from_one = true;
#else
	const bool is_index_from_one = false;
#endif
	if (is_index_from_one)
	{
		ptr = (_INTEGER_t*)malloc((nrows + 1) * sizeof(_INTEGER_t));
		index = (_INTEGER_t*)malloc((nonzero) * sizeof(_INTEGER_t));
		// Reindexing ptr according to Fortran-based PARDISO
		_INTEGER_t i = 0;
		for (i = 0; i < nrows; ++i)
			ptr[i] = A->ptr[i] + 1;
		// ptr needs one more storage
		ptr[i] = A->ptr[i] + 1;
		// Reindexing index according to Fortran-based PARDISO
		// and zonzero of Matrix A
		for (_INTEGER_t i = 0; i < nonzero; ++i)
			index[i] = A->col_idx[i] + 1;
	}

#if 0
	{
		std::ofstream ofs("pardiso.txt");
		ofs << "ptr=\n";
		for(i = 0; i < nrows+1; ++i)
			ofs << ptr[i] << " ";
		ofs << "\ncol_idx=\n";
		for(i = 0; i < nonzero; ++i)
			ofs << index[i] << " ";
		ofs << "\nval=\n";
		for(i = 0; i < nonzero; ++i)
			ofs << value[i] << " ";
		ofs << "\n";
		ofs.close();
	}
#endif

	_INTEGER_t mtype = 11; /* Real unsymmetric matrix */
	if (num->ls_storage_method == 102)
		mtype = 1;       // Real and structurally symmetric
	_INTEGER_t nrhs = 1; /* Number of right hand sides. */
	/* Internal solver memory pointer pt, */
	/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
	/* or void *pt[64] should be OK on both architectures */
	void* pt[64];
	/* Pardiso control parameters.*/
	_INTEGER_t iparm[64];
	_INTEGER_t maxfct, mnum, phase, error, msglvl;

	/* Auxiliary variables.*/
	double ddum;     /* Double dummy */
	_INTEGER_t idum; /* Integer dummy. */

#ifdef _WIN32
	double dparm[64];
	_INTEGER_t solver;
	// Check the license and initialize the solver
	{
		// static bool done = false;
		// if (!done) {
		PARDISOINIT(pt, &mtype, &solver, iparm, dparm, &error);
		if (error != 0)
		{
			if (error == -10) printf("->No license file found \n");
			if (error == -11) printf("->License is expired \n");
			if (error == -12) printf("->Wrong username or hostname \n");
			exit(1);
		}
		else
			printf("->PARDISO license check was successful ... \n");

		//  done = true;
		//}
	}
#endif

	/* --------------------------------------------------------------------*/
	/* .. Setup Pardiso control parameters.*/
	/* --------------------------------------------------------------------*/
	for (int i = 0; i < 64; i++)
		iparm[i] = 0;
	iparm[0] = 1; /* No solver default */
	iparm[1] = 2; /* Fill-in reordering from METIS */
/* Numbers of processors, value of MKL_NUM_THREADS */
#ifdef _WIN32
	iparm[2] = omp_get_max_threads();
#else
	iparm[2] = mkl_get_max_threads();
#endif
	iparm[3] = 0;   /* No iterative-direct algorithm */
	iparm[4] = 0;   /* No user fill-in reducing permutation */
	iparm[5] = 0;   /* Write solution into x */
	iparm[6] = 0;   /* Not in use */
	iparm[7] = 2;   /* Max numbers of iterative refinement steps */
	iparm[8] = 0;   /* Not in use */
	iparm[9] = 13;  /* Perturb the pivot elements with 1E-13 */
	iparm[10] = 1;  /* Use nonsymmetric permutation and scaling MPS */
	iparm[11] = 0;  /* Not in use */
	iparm[12] = 1;  /* Use (non-) symmetric weighted matching  */
	iparm[13] = 0;  /* Output: Number of perturbed pivots */
	iparm[14] = 0;  /* Not in use */
	iparm[15] = 0;  /* Not in use */
	iparm[16] = 0;  /* Not in use */
	iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	iparm[18] = -1; /* Output: Mflops for LU factorization */
	iparm[19] = 0;  /* Output: Numbers of CG Iterations */
	iparm[34] = 1;  /* Input: Zero-based indexing */
	iparm[59] = 1;  /* PARDISO mode - in-core (0) or out-core (2) */
	maxfct = 1;     /* Maximum number of numerical factorizations. */
	mnum = 1;       /* Which factorization to use. */
	msglvl = 0;     /* Print statistical information in file */
	if (nrows > 1e6) msglvl = 1;  // output log for large problems
	error = 0;                    /* Initialize error flag */

	/* --------------------------------------------------------------------*/
	/* .. Initialize the internal solver memory pointer. This is only */
	/* necessary for the FIRST call of the PARDISO solver. */
	/* --------------------------------------------------------------------*/
	for (int i = 0; i < 64; i++)
		pt[i] = 0;

	/* --------------------------------------------------------------------*/
	/* .. Reordering and Symbolic Factorization. This step also allocates */
	/* all memory that is necessary for the factorization. */
	/* --------------------------------------------------------------------*/
	ScreenMessage("-> Executing PARDISO\n");
	phase = 11;
#ifdef _WIN32
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
#else
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
#endif

	if (msglvl == 1)
	{
		printf("< Memory usage >\n");
		printf("             Peak memory on symbolic factorization = %d kb\n",
		       iparm[14]);
		printf(
		    "             Permanent memory on symbolic factorization = %d kb\n",
		    iparm[15]);
		printf(
		    "             Size of factors/Peak memory on numerical "
		    "factorization and solution = %d\n",
		    iparm[16]);
		printf("             Total peak memory = %d kb\n",
		       std::max(iparm[14], iparm[15] + iparm[16]));
	}
	if (error != 0)
	{
		printf("\nERROR during symbolic factorization: %d", error);
		exit(1);
	}

	/* --------------------------------------------------------------------*/
	/* .. Numerical factorization.*/
	/* --------------------------------------------------------------------*/
	phase = 22;
#ifdef _WIN32
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
#else
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
#endif
	if (msglvl == 1)
	{
		printf("< Memory usage >\n");
		printf("             Peak memory on symbolic factorization = %d kb\n",
		       iparm[14]);
		printf(
		    "             Permanent memory on symbolic factorization = %d kb\n",
		    iparm[15]);
		printf(
		    "             Size of factors/Peak memory on numerical "
		    "factorization and solution = %d\n",
		    iparm[16]);
		printf("             Total peak memory = %d kb\n",
		       std::max(iparm[14], iparm[15] + iparm[16]));
	}
	if (error != 0)
	{
		printf("\nERROR during numerical factorization: %d", error);
		exit(2);
	}

	/* --------------------------------------------------------------------*/
	/* .. Back substitution and iterative refinement. */
	/* --------------------------------------------------------------------*/
	phase = 33;
	iparm[7] = 2; /* Max numbers of iterative refinement steps. */

/* Set right hand side to one. */

#ifdef _WIN32
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, b, x, &error, dparm);
#else
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, tmp_b, tmp_x, &error);
#endif
	if (msglvl == 1)
	{
		printf("< Memory usage >\n");
		printf("             Peak memory on symbolic factorization = %d kb\n",
		       iparm[14]);
		printf(
		    "             Permanent memory on symbolic factorization = %d kb\n",
		    iparm[15]);
		printf(
		    "             Size of factors/Peak memory on numerical "
		    "factorization and solution = %d\n",
		    iparm[16]);
		printf("             Total peak memory = %d kb\n",
		       std::max(iparm[14], iparm[15] + iparm[16]));
	}
	if (error != 0)
	{
		printf("\nERROR during solution: %d", error);
		exit(3);
	}

	phase = -1; /* Release internal memory. */
#ifdef _WIN32
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error, dparm);
#else
	PARDISO(pt, &maxfct, &mnum, &mtype, &phase, &nrows, value, ptr, index,
	        &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
#endif

	if (is_compressed)
	{
#pragma omp parallel for
		for (int i = 0; i < nrows; ++i)
		{
			x[vec_nz_rows[i]] = tmp_x[i];
		}
	}

	// Releasing the local memory
	delete[] value;
	if (is_compressed || is_index_from_one)
	{
		free(ptr);
		free(index);
	}
	if (is_compressed)
	{
		delete tmp_x;
		delete tmp_b;
	}

	//		MKL_FreeBuffers();
	ScreenMessage2("-> Finished PARDISO computation\n");
	ScreenMessage2(
	    "------------------------------------------------------------------\n");
}
#endif

#ifdef LIS
int Linear_EQS::solveWithLIS(CNumerics* m_num, bool compress)
{
	ScreenMessage2(
	    "------------------------------------------------------------------\n");
	ScreenMessage2("*** LIS solver computation\n");

	// Prepare CRS data
	LIS_INT nrows = A->Size() * A->Dof();
	LIS_INT nonzero = A->nnz();
	// ScreenMessage2("-> copying CRS data with dim=%ld and nnz=%ld\n", nrows,
	// nonzero);
	LIS_REAL* value = new LIS_REAL[nonzero];
	A->GetCRSValue(value);
	LIS_INT* ptr = A->ptr;
	LIS_INT* col_idx = A->col_idx;

	bool is_compressed = false;
	std::vector<LIS_INT> vec_nz_rows;
	if (compress)
	{
		// check non-zero rows, non-zero entries
		ScreenMessage2("-> Check non-zero entries\n");
		std::vector<LIS_INT> vec_z_rows;
		IndexType n_nz_entries =
		    searcgNonZeroEntries(nrows, ptr, value, vec_nz_rows, vec_z_rows);

		const LIS_INT n_new_rows = vec_nz_rows.size();
		if (n_new_rows < nrows)
		{
			ScreenMessage2(
			    "-> Found %d empty rows. Delete them from total %d rows\n",
			    nrows - n_new_rows, nrows);
			is_compressed = true;
			LIS_INT* new_ptr, *new_col_index;
			LIS_REAL* new_value;
			compressCRS(ptr, col_idx, value, vec_nz_rows, vec_z_rows,
			            n_nz_entries, new_ptr, new_col_index, new_value);

			nrows = n_new_rows;
			nonzero = n_nz_entries;
			ptr = new_ptr;
			col_idx = new_col_index;
			delete[] value;
			value = new_value;
		}
	}

	// Creating a matrix.
	int ierr = lis_matrix_create(0, &AA);
	CHKERR(ierr);
#ifndef OGS_USE_LONG
	ierr = lis_matrix_set_type(AA, LIS_MATRIX_CRS);
	CHKERR(ierr);
#else
	ierr = lis_matrix_set_type(AA, LIS_MATRIX_CSR);
	CHKERR(ierr);
#endif
	ierr = lis_matrix_set_size(AA, 0, nrows);
	CHKERR(ierr);

	// Matrix solver and Precondition can be handled better way.
	char solver_options[MAX_ZEILE], tol_option[MAX_ZEILE];
	sprintf(solver_options,
	        "-i %d -p %d %s",
	        m_num->ls_method,
	        m_num->ls_precond,
	        m_num->ls_extra_arg.c_str());
	// tolerance and other setting parameters are same
	// NW add max iteration counts
	sprintf(tol_option,
	        "-tol %e -maxiter %d",
	        m_num->ls_error_tolerance,
	        m_num->ls_max_iterations);

#ifndef OGS_USE_LONG
	ierr = lis_matrix_set_crs(nonzero, ptr, col_idx, value, AA);
	CHKERR(ierr);
//		ierr = lis_matrix_set_crs(nonzero,A->ptr,A->col_idx, value,AA);
#else
	ierr = lis_matrix_set_csr(nonzero, ptr, col_idx, value, AA);
	CHKERR(ierr);
//		ierr = lis_matrix_set_csr(nonzero,A->ptr,A->col_idx, value,AA);
#endif
	ierr = lis_matrix_assemble(AA);
	CHKERR(ierr);

	// Assemble the vector, b, x
	ierr = lis_vector_duplicate(AA, &bb);
	CHKERR(ierr);
	ierr = lis_vector_duplicate(AA, &xx);
	CHKERR(ierr);

	if (!is_compressed)
	{
#pragma omp parallel for
		for (int i = 0; i < nrows; ++i)
		{
			ierr = lis_vector_set_value(LIS_INS_VALUE, i, x[i], xx);
			ierr = lis_vector_set_value(LIS_INS_VALUE, i, b[i], bb);
		}
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i < nrows; ++i)
		{
			ierr =
			    lis_vector_set_value(LIS_INS_VALUE, i, x[vec_nz_rows[i]], xx);
			ierr =
			    lis_vector_set_value(LIS_INS_VALUE, i, b[vec_nz_rows[i]], bb);
		}
	}

	// lis_output(AA, bb, xx, LIS_FMT_MM, "leqs.txt");

	// Create solver
	ierr = lis_solver_create(&solver);
	CHKERR(ierr);

	ierr = lis_solver_set_option(solver_options, solver);
	ierr = lis_solver_set_option(tol_option, solver);
	ierr = lis_solver_set_option("-print mem", solver);
	ierr = lis_solver_set_optionC(solver);
	ScreenMessage2("-> Execute Lis\n");
	ierr = lis_solve(AA, bb, xx, solver);
	CHKERR(ierr);
	ScreenMessage2("-> done\n");
	int status = 0;
	lis_solver_get_status(solver, &status);
	printf("status   : %d\n", status);
	int iter = 0;
	ierr = lis_solver_get_iters(solver, &iter);
	// NW
	printf("iteration: %d/%d\n", iter, m_num->ls_max_iterations);
	double resid = 0.0;
	ierr = lis_solver_get_residualnorm(solver, &resid);
	printf("residuals: %e\n", resid);
	//	lis_vector_print(xx);
	//	lis_vector_print(bb);

	// Update the solution (answer) into the x vector
	if (!is_compressed)
	{
#pragma omp parallel for
		for (int i = 0; i < nrows; ++i)
			lis_vector_get_value(xx, i, &(x[i]));
	}
	else
	{
#pragma omp parallel for
		for (int i = 0; i < nrows; ++i)
		{
			lis_vector_get_value(xx, i, &(x[vec_nz_rows[i]]));
		}
	}

	// Clear memory
	if (is_compressed)
	{
#if 0
	    lis_matrix_destroy(AA);
#endif
	}
	else
	{
		delete[] value;
	}
	//	lis_matrix_destroy(AA);
	lis_vector_destroy(bb);
	lis_vector_destroy(xx);
	lis_solver_destroy(solver);
	ScreenMessage2(
	    "------------------------------------------------------------------\n");

	return (status == 0) ? iter : -1;
}
#endif

#if defined(LIS) || defined(MKL)
int Linear_EQS::Solver(CNumerics* num, bool compress)
{
	CNumerics* m_num = (num == NULL) ? num_vector[0] : num;
#define ENABLE_COMPRESS_EQS
#ifndef ENABLE_COMPRESS_EQS
	compress = false;
#endif

	int iter = 0;
#ifdef _OPENMP
	// omp_set_num_threads (1);
	ScreenMessage2("-> Use OpenMP with %d threads\n", omp_get_max_threads());
#endif
#ifdef OGS_USE_LONG
	ScreenMessage2("-> 64bit integer is used in PARDISO\n");
#endif

	if (m_num->ls_method == 805)  // Then, PARDISO parallel direct solver
	{
#ifdef MKL
		solveWithPARDISO(num, compress);
#endif
	}
	else  // LIS parallel solver
	{
#ifdef LIS
		iter = solveWithLIS(num, compress);
#endif
	}

	return iter;
}
#endif
#endif

/**************************************************************************
   Task: Linear equation::SetKnownXi
      Configure equation system when one entry of the vector of
      unknown is given
   Programing:
   10/2007 WW/
**************************************************************************/
void Linear_EQS::SetKnownX_i(const long i, const double x_i)
{
	A->Diagonize(i, x_i, b);
}

/*\!
 ********************************************************************
   Dot production of two vectors
   Programm:
   10/2007 WW
   12/2007 WW  Parallel
 ********************************************************************/
double Linear_EQS::dot(const double* xx, const double* yy)
{
	long i;
	double val = 0.;
	for (i = 0; i < size_A; i++)
		val += xx[i] * yy[i];
	return val;
}
/*\!
 ********************************************************************
   Dot production of two vectors
   Programm:
   01/2008 WW
 ********************************************************************/
double Linear_EQS::NormX()
{
	return sqrt(dot(x, x));
}


//------------------------------------------------------------------------
}  // namespace
