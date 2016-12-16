/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   Task: Linear equation
   Programing:
   11/2007 WW/
**************************************************************************/
// There is a name conflict between stdio.h and the MPI C++ binding
// with respect to the names SEEK_SET, SEEK_CUR, and SEEK_END.  MPI
// wants these in the MPI namespace, but stdio.h will #define these
// to integer values.  #undef'ing these can cause obscure problems
// with other include files (such as iostream), so we instead use
// #error to indicate a fatal error.  Users can either #undef
// the names before including mpi.h or include mpi.h *before* stdio.h
// or iostream.
#if defined(USE_MPI)
#include "par_ddc.h"
#include <mpi.h>
#endif

#include "equation_class.h"

#include <cfloat>
#include <iomanip>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

// NEW_EQS To be removed
#ifdef NEW_EQS  // 1.11.2007 WW

#ifdef LIS  // 07.02.2008 PCH
#include "lis.h"
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
#include "makros.h"
#include "display.h"

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
                       const long dof, bool messg)
    : message(messg)
{
	long i;
	A = new CSparseMatrix(sparse_table, dof);
	size_A = A->Dim();
	size_global = 0;

	prec_M = NULL;

#if defined(USE_MPI)
	x = NULL;
	cpu_time = 0.0;
	// WW
	border_buffer0 = NULL;
	border_buffer1 = NULL;
#else
	x = new double[size_A];
#endif
	b = new double[size_A];
	//
	for (i = 0; i < size_A; i++)
	{
#ifndef USE_MPI
		x[i] = 0.;
#endif
		b[i] = 0.;
	}
	iter = 0;
	bNorm = 1.0;
	error = 1.0e10;
}
#if defined(USE_MPI)
/**************************************************************************
   Task: Linear equation::Constructor
   Programing:
   12/2007 WW/
**************************************************************************/
Linear_EQS::Linear_EQS(const long size)
{
	long i;
	A = NULL;
	b = NULL;
	size_global = size;
	x = new double[size];
	//
	for (i = 0; i < size; i++)
		x[i] = 0.;
	iter = 0;
	bNorm = 1.0;
	error = 1.0e10;
}
#endif
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

	/// GMRES. 30.06.2010. WW
	if (solver_type == 13) H.ReleaseMemory();
}
/**************************************************************************
   Task: Linear equation::
   Programing:
   10/2007 WW/
**************************************************************************/
void Linear_EQS::ConfigNumerics(CNumerics* m_num, const long n)
{
	(void)n;
	int i, nbuffer = 0;  // Number of temperary float arrays
	precond_type = m_num->ls_precond;
	solver_type = m_num->ls_method;
	switch (solver_type)
	{
		case 1:
			solver_name = "Gauss";
			break;
		case 2:
			solver_name = "BiCGSTab";
			nbuffer = 8;
			break;
		case 3:
			solver_name = "BiCG";
			nbuffer = 8;  // 20.10.2010. WW
			break;
		case 4:
			solver_name = "QMRCGStab";
			break;
		case 5:
			solver_name = "CG";
			nbuffer = 3;
			break;
		case 6:
			solver_name = "CGNR";
			break;
		case 7:
			solver_name = "CGS";
			nbuffer = 9;
			break;
		case 8:
			solver_name = "Richardson";
			break;
		case 9:
			solver_name = "JOR";
			break;
		case 10:
			solver_name = "SOR";
			break;
		case 11:
			solver_name = "AMG1R5";
			break;
		case 12:
			solver_name = "UMF";
			break;
		case 13:  // 06.2010. WW
			solver_name = "GMRES";
			m_gmres = m_num->Get_m();
			for (i = 0; i < 4; i++)
			{
				double* new_array = new double[m_gmres + 1];
				f_buffer.push_back(new_array);
			}
			H.resize(m_gmres + 1, m_gmres + 1);
			nbuffer = m_gmres + 4;
			break;
	}
	// Buffer
	/*
	   #if defined(USE_MPI)
	   long size = A->Dim()+A->Dof()*dom->BSize();
	   #else
	 */
	//#endif
	for (i = 0; i < nbuffer; i++)
	{
		double* new_array = new double[size_A];
		f_buffer.push_back(new_array);
	}
#if defined(USE_MPI)
	if (!x) x = new double[size_A];
#if defined(NEW_BREDUCE)
	// For concatenate border entries
	double* catb = NULL;
	f_buffer.push_back(catb);
#endif
	// Buffer for a global array
	double* x_g_buffer = new double[n];
	f_buffer.push_back(x_g_buffer);
	// Buffer for border array
	// Will be removed if topo is ready
	border_buffer0 = new double[A->Dof() * dom->BSize()];
	border_buffer1 = new double[A->Dof() * dom->BSize()];
//
#endif
	//---------------------------------------------
	switch (precond_type)
	{
		case 1:
			precond_name = "Jacobi";
#if defined(USE_MPI)
			prec_M = new double[size_A];
#endif
			break;
		case 100:
			// precond_name = "ILU"; break;
			// If ILU is ready, remove follows
			// ----------------------------------------------
			precond_name = "ILU not available. Use Jacobi";
			precond_type = 1;
#if defined(USE_MPI)
			prec_M = new double[size_A];
#endif
			// ----------------------------------------------
			break;
		default:
			precond_name = "No preconditioner";
			break;
	}
	//
	//
	max_iter = m_num->ls_max_iterations;
	tol = m_num->ls_error_tolerance;
	//
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
	error = 1.0e10;
}
/**************************************************************************
   Task: Linear equation::Alocate memory for solver
   Programing:
   11/2007 WW/
**************************************************************************/
void Linear_EQS::Clean()
{
#if defined(USE_MPI)
	double cpu_time_local = -MPI_Wtime();
#endif
	for (int i = 0; i < (int)f_buffer.size(); i++)
	{
		if (f_buffer[i]) delete[] f_buffer[i];
		f_buffer[i] = NULL;
	}
	f_buffer.clear();
	if (prec_M) delete[] prec_M;
	prec_M = NULL;
#if defined(USE_MPI)
	//
	if (border_buffer0) delete[] border_buffer0;
	border_buffer0 = NULL;
	if (border_buffer1) delete[] border_buffer1;
	border_buffer1 = NULL;
	//
	cpu_time_local += MPI_Wtime();
	cpu_time += cpu_time_local;
#endif
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
int Linear_EQS::Solver(double* xg, const long n)
{
	//
	double cpu_time_local = -MPI_Wtime();
	iter = 0;
	ComputePreconditioner();
	size_global = n;
	switch (solver_type)
	{
		case 2:
			iter = BiCGStab(xg, n);
			break;
		case 3:
			iter = BiCG(xg, n);
			break;
		case 5:
			iter = CG(xg, n);
			break;
		case 7:
			iter = CGS(xg, n);
			break;
	}
	cpu_time_local += MPI_Wtime();
	cpu_time += cpu_time_local;
	return iter;
}

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

	// get RHS norm
	this->bNorm = Norm(b);

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
#else  // ifdef LIS
int Linear_EQS::Solver()
{
	//
	iter = 0;
	ComputePreconditioner();
	switch (solver_type)
	{
		case 1:
			return Gauss();
		case 2:
			iter = BiCGStab();
			return iter;  // kg44 only to make sure here is iter returned
		case 3:
			return BiCG();
		case 4:
			return QMRCGStab();
		case 5:
			return CG();
		case 6:
			return CGNR();
		case 7:
			return CGS();
		case 8:
			return Richardson();
		case 9:
			return JOR();
		case 10:
			return SOR();
		case 11:
			return AMG1R5();
		case 12:
			return UMF();
		case 13:
			return GMRES();
			break;
	}
	return -1;
}
#endif
#endif
// Preconditioners
/**************************************************************************
   Task: Preconditioners
   Programing:
   10/2007 WW
**************************************************************************/
void Linear_EQS::ComputePreconditioner()
{
	switch (precond_type)
	{
		case 1:
#if defined(USE_MPI)
			ComputePreconditioner_Jacobi();
#endif
			return;
		case 100:
			ComputePreconditioner_ILU();
			return;
		default:
			return;
	}
}
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
/**************************************************************************
   Task: Linear equation::Preconditioner
   Programing:
   08/2007 WW/
**************************************************************************/
void Linear_EQS::Precond(double* vec_s, double* vec_r)
{
	bool pre = true;
	switch (precond_type)
	{
		case 1:
#if defined(USE_MPI)
			Precond_Jacobi(vec_s, vec_r);
#else
				A->Precond_Jacobi(vec_s, vec_r);
#endif
			break;
		case 100:
			pre = false;  // A->Precond_ILU(vec_s, vec_r);
			break;
		default:
			pre = false;  // A->Precond_ILU(vec_s, vec_r);
			break;
	}
	if (!pre)
		for (long i = 0; i < size_A; i++)
			vec_r[i] = vec_s[i];
}
/**************************************************************************
   Task: Linear equation:: M^T x
   Transpose of preconditioner times a vector
   Programing:
   02/2010 WW/
**************************************************************************/
void Linear_EQS::TransPrecond(double* vec_s, double* vec_r)
{
	Precond(vec_s, vec_r);
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
#if defined(USE_MPI)
	double val_i = dom->Dot_Interior(xx, yy);
	val_i += dom->Dot_Border_Vec(xx, yy);
	//
	MPI_Allreduce(&val_i, &val, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
	for (i = 0; i < size_A; i++)
		val += xx[i] * yy[i];
#endif
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
#if defined(USE_MPI)
	return sqrt(dot(x, x, size_global));
#else
	return sqrt(dot(x, x));
#endif
}
//
#if defined(USE_MPI)
/*\!
 ********************************************************************
   Dot production of two vectors
   Programm:
   12/2007 WW
 ********************************************************************/
double Linear_EQS::dot(const double* xx, const double* yy, const long n)
{
	long i;
	double val = 0.;
	for (i = 0; i < n; i++)
		val += xx[i] * yy[i];
	return val;
}
/*\!
 ********************************************************************
   Dot production of two vectors
   Programm:
   12/2007 WW
 ********************************************************************/
inline void Linear_EQS::MatrixMulitVec(double* xx, double* yy)
{
	long i;  //, size = A->Dim();
	//
	A->multiVec(xx, yy);
#if defined(NEW_BREDUCE)
	dom->ReduceBorderV(yy);
#else
	dom->Local2Border(yy, border_buffer0);
	MPI_Allreduce(border_buffer0, border_buffer1, A->Dof() * dom->BSize(),
	              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dom->Border2Local(border_buffer1, yy);
#endif
}
/*!
 ********************************************************************
   Dot production of two vectors
   Programm:
   12/2007 WW
   02/2010 WW Revise
 ********************************************************************/
inline void Linear_EQS::TransMatrixMulitVec(double* xx, double* yy)
{
	//
	A->Trans_MultiVec(xx, yy);
#if defined(NEW_BREDUCE)
	dom->ReduceBorderV(yy);
#else
	dom->Local2Border(yy, border_buffer0);
	MPI_Allreduce(border_buffer0, border_buffer1, A->Dof() * dom->BSize(),
	              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dom->Border2Local(border_buffer1, yy);
#endif
}
#endif
/*!
 ********************************************************************
   ConvergeTest
   Programm:
   09/2007 WW
 ********************************************************************/
void Linear_EQS::Message()
{
#ifdef USE_MPI
	if (myrank > 0) return;
#endif
	if (!message) return;
	cout.width(10);
	cout.precision(3);
	cout.setf(ios::scientific);
	//
	// system("color 0B");
	cout << "      ------------------------------------------------\n";
	cout << "      Linear solver " << solver_name << " with " << precond_name
	     << ":\n";
	cout << "      Iterations |"
	     << " Max Iters |"
	     << " Norm of b |"
	     << " Error\n";
	cout << "      " << setw(11) << iter << "|" << setw(11) << max_iter << "|"
	     << setw(11) << bNorm << "|" << setw(11) << error << "\n";
	if (iter == max_iter)
		cout << "      WARNING: Maximum iterations reached !!! \n";
	cout << "      ------------------------------------------------\n";
	cout.flush();
}
/*\!
 ********************************************************************
   Check if the norm of b is samll enough for convengence.
   normb_new is given to bNorm;
   Programm:
   09/2007 WW
 ********************************************************************/
inline bool Linear_EQS::CheckNormRHS(const double normb_new)
{
	if (bNorm > 0.0)
		if ((normb_new / bNorm) < tol)
		{
			error = normb_new / bNorm;
			bNorm = normb_new;
			Message();
			return true;
		}
	bNorm = normb_new;
	if (bNorm < DBL_MIN)
	{
		error = 0.;
		Message();
		return true;
	}
	return false;
}
#ifndef USE_MPI
/**************************************************************************
   Task: Linear equation::CG
   Programing:
   11/2007 WW/
**************************************************************************/
int Linear_EQS::CG()
{
	long i, size;
	double rrM1;
	double* p, *r, *s;
	//
	size = A->Dim();
	p = f_buffer[0];
	r = f_buffer[1];
	s = f_buffer[2];
	//
	double bNorm_new = Norm(b);
	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//
	// r0 = b-Ax
	A->multiVec(x, s);
	for (i = 0; i < size; i++)
		r[i] = b[i] - s[i];
	//
	// Preconditioning: M^{-1}r
	Precond(r, s);
	for (i = 0; i < size; i++)
		p[i] = s[i];
	// Check the convergence
	if ((error = Norm(r) / bNorm) < tol)
	{
		Message();
		return 1;
	}
	//
	double rr = dot(r, s);
	//
	for (iter = 1; iter <= max_iter; ++iter)
	{
		A->multiVec(p, s);
		const double alpha = rr / dot(p, s);
		// Update
		for (i = 0; i < size; i++)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * s[i];
		}
		if ((error = Norm(r) / bNorm) < tol)
		{
			Message();
			return iter <= max_iter;
		}
		//
		Precond(r, s);
		//
		rrM1 = rr;
		rr = dot(s, r);
		const double beta = rr / rrM1;
		for (i = 0; i < size; i++)
			p[i] = s[i] + beta * p[i];
	}
	//
	Message();
	return iter <= max_iter;
}
/**************************************************************************
   Task: Linear equation::BiCG
   Programing:
   10/2010 WW/
**************************************************************************/
int Linear_EQS::BiCG()
{
	long i, size;
	double rho1, rho2 = 0., alpha, beta;
	double* z, *zt, *p, *pt, *q, *qt, *r, *rt;
	//
	size = A->Dim();
	z = f_buffer[0];
	zt = f_buffer[1];
	p = f_buffer[2];
	pt = f_buffer[3];
	q = f_buffer[4];
	qt = f_buffer[5];
	r = f_buffer[6];
	rt = f_buffer[7];
	//
	double bNorm_new = Norm(b);
	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//
	// r0 = b-Ax
	A->multiVec(x, rt);
	for (i = 0; i < size; i++)
	{
		r[i] = b[i] - rt[i];
		rt[i] = r[i];
	}
	//
	// Check the convergence
	if ((error = Norm(r) / bNorm) < tol)
	{
		Message();
		return 1;
	}
	//
	//
	for (iter = 1; iter <= max_iter; ++iter)
	{
		Precond(r, z);
		TransPrecond(rt, zt);
		rho1 = dot(z, rt);
		//
		if (fabs(rho1) < DBL_MIN)
		{
			Message();
			return iter <= max_iter;
		}
		//
		if (iter == 1)
			for (i = 0; i < size; i++)
			{
				p[i] = z[i];
				pt[i] = zt[i];
			}
		else
		{
			beta = rho1 / rho2;
			for (i = 0; i < size; i++)
			{
				p[i] = z[i] + beta * p[i];
				pt[i] = zt[i] + beta * pt[i];
			}
		}
		//
		A->multiVec(p, q);
		A->Trans_MultiVec(pt, qt);
		alpha = rho1 / dot(pt, q);
		//
		for (i = 0; i < size; i++)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * q[i];
			rt[i] -= alpha * qt[i];
		}
		//
		rho2 = rho1;
		if ((error = Norm(r) / bNorm) < tol)
		{
			Message();
			return iter <= max_iter;
		}
		//
	}
	//
	Message();
	return iter <= max_iter;
}

/*************************************************************************
   GeoSys-Function:
   Task: BiCGStab solver
   Programming:
   10/2007 WW
 **************************************************************************/
int Linear_EQS::BiCGStab()
{
	long i, size;
	double rho_0, rho_1, alpha, beta, omega, tt = 0., norm_r = 0.;
	double* r0, *r, *s, *s_h, *t, *v, *p, *p_h;
	//
	size = size_A;
	r0 = f_buffer[0];
	r = f_buffer[1];
	s = f_buffer[2];
	s_h = f_buffer[3];
	t = f_buffer[4];
	v = f_buffer[5];
	p = f_buffer[6];
	p_h = f_buffer[7];
	//
	rho_0 = alpha = omega = 1.0;
	//
	double bNorm_new = Norm(b);
	// Check if the norm of b is small enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
//
// Norm of M r
	A->multiVec(x, s);  // s as buffer
	for (i = 0; i < size; i++)
		r0[i] = b[i] - s[i];                    // r = b-Ax

	for (i = 0; i < size; i++)
	{
		r[i] = r0[i];
		v[i] = 0.;
		p[i] = 0.;
	}
	if ((error = Norm(r) / bNorm) < tol)
	{
		Message();
		return 0;
	}
	//
	for (iter = 1; iter <= max_iter; iter++)
	{
		rho_1 = dot(r0, r);
		if (fabs(rho_1) < DBL_MIN)  // DBL_EPSILON
		{
			Message();
			return 0;
		}
		if (iter == 1)
			for (i = 0; i < size; i++)
				p[i] = r[i];
		else
		{
			beta = (rho_1 / rho_0) * (alpha / omega);
			for (i = 0; i < size; i++)
				p[i] = r[i] + beta * (p[i] - omega * v[i]);
		}
		// Preconditioner
		Precond(p, p_h);
// A M^{-1}p-->v
			A->multiVec(p_h, v);
		//
		alpha = rho_1 / dot(r0, v);
		//
		for (i = 0; i < size; i++)
			s[i] = r[i] - alpha * v[i];
		if ((error = Norm(s) / bNorm) < tol)
		{
			for (i = 0; i < size; i++)
				x[i] += alpha * p_h[i];
			Message();
			return iter;
		}
		//  M^{-1}s,
		Precond(s, s_h);
// A* M^{-1}s
			A->multiVec(s_h, t);
		//
		tt = dot(t, t);
		if (tt > DBL_MIN)
			omega = dot(t, s) / tt;
		else
			omega = 1.0;
		// Update solution
		for (i = 0; i < size; i++)
		{
			x[i] += alpha * p_h[i] + omega * s_h[i];
			r[i] = s[i] - omega * t[i];
		}
		rho_0 = rho_1;
		//
		norm_r = Norm(r);
		if ((error = norm_r / bNorm) < tol)
		{
			Message();
			return iter;
		}
		if (fabs(omega) < DBL_MIN)
		{
			error = norm_r / bNorm;
			Message();
			return iter;
		}
	}
	//
	Message();
	//
	return iter;
}

/*************************************************************************
   GeoSys-Function:
   Task: CGS solver
   Programming:
   11/2007 WW
 **************************************************************************/
int Linear_EQS::CGS()
{
	long i, size;
	double rho_1, rho_2, alpha, beta;
	double* r0, *r, *p, *p_h, *q, *q_h, *v, *u, *u_h;
	//
	size = A->Dim();
	r0 = f_buffer[0];
	r = f_buffer[1];
	p = f_buffer[2];
	p_h = f_buffer[3];
	q = f_buffer[4];
	q_h = f_buffer[5];
	v = f_buffer[6];
	u = f_buffer[7];
	u_h = f_buffer[8];
	//
	rho_1 = rho_2 = 1.0;
	//
	double bNorm_new = Norm(b);
	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//
	A->multiVec(x, v);  // v as buffer
	for (i = 0; i < size; i++)
	{
		r0[i] = b[i] - v[i];  // r = b-Ax
		r[i] = r0[i];
		v[i] = 0.;
	}
	if ((error = Norm(r) / bNorm) < tol)
	{
		Message();
		return 0;
	}
	//
	for (iter = 1; iter <= max_iter; iter++)
	{
		rho_1 = dot(r0, r);
		if (fabs(rho_1) < DBL_MIN)  //  DBL_EPSILON
		{
			Message();
			return 0;
		}
		if (iter == 1)
			for (i = 0; i < size; i++)
				p[i] = u[i] = r[i];
		else
		{
			beta = rho_1 / rho_2;
			for (i = 0; i < size; i++)
			{
				u[i] = r[i] + beta * q[i];
				p[i] = u[i] + beta * (q[i] + beta * p[i]);
			}
		}
		// Preconditioner
		Precond(p, p_h);
		// A M^{-1}p-->v
		A->multiVec(p_h, v);
		//
		alpha = rho_1 / dot(r0, v);
		//
		for (i = 0; i < size; i++)
		{
			q[i] = u[i] - alpha * v[i];
			q_h[i] = u[i] + q[i];
		}
		// Preconditioner
		Precond(q_h, u_h);
		for (i = 0; i < size; i++)
			x[i] += alpha * u_h[i];
		//
		A->multiVec(u_h, q_h);
		//
		for (i = 0; i < size; i++)
			r[i] -= alpha * q_h[i];
		rho_2 = rho_1;
		if ((error = Norm(r) / bNorm) < tol)
		{
			Message();
			return iter <= max_iter;
		}
	}
	//
	Message();
	//
	return iter <= max_iter;
}
//
//------------------------------------------------------------------------
#define aGMRES
#ifdef aGMRES

//-----------------------------------------------------------------
/*!
     GMRES solver.

     by WW. 06.2010
 */
//-----------------------------------------------------------------
/// For GMRES
inline void Linear_EQS::Get_Plane_Rotation(double& dx, double& dy, double& cs,
                                           double& sn)
{
	if (dy == 0.0)
	{
		cs = 1.0;
		sn = 0.0;
	}
	else if (fabs(dy) > fabs(dx))
	{
		double temp = dx / dy;
		sn = 1.0 / sqrt(1.0 + temp * temp);
		cs = temp * sn;
	}
	else
	{
		double temp = dy / dx;
		cs = 1.0 / sqrt(1.0 + temp * temp);
		sn = temp * cs;
	}
}

/// For GMRES.
inline void Linear_EQS::Set_Plane_Rotation(double& dx, double& dy, double& cs,
                                           double& sn)
{
	double temp = cs * dx + sn * dy;
	dy = -sn * dx + cs * dy;
	dx = temp;
}

/// Update solution in GMRES
inline void Linear_EQS::Update(double* x, int k, Matrix& h, double* s)
{
	long i;
	long m, j;
	long size = 0;
	int v_idx0 = 7;

	if (A)
		size = A->Dim();
	else
		size = size_global;

	double* v_j;

	m = m_gmres;

	double* y = f_buffer[3];
	for (j = 0; j < m + 1; j++)
		y[j] = s[j];

	// Back solve
	for (i = k; i >= 0; i--)
	{
		y[i] /= h(i, i);
		for (j = i - 1; j >= 0; j--)
			y[j] -= h(j, i) * y[i];
	}

	for (j = 0; j <= k; j++)
	{
		v_j = f_buffer[v_idx0 + j];
		for (i = 0; i < size; i++)
			x[i] += v_j[i] * y[j];
	}
}

//#ifndef USE_MPI
/// GMRES solver. WW
int Linear_EQS::GMRES()
{
	long i, k, l, m;
	double normb, beta;
	double* s, *cs, *sn, *v, *w, *r, *t, *v_k;

	normb = Norm(b);
	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(normb)) return 0;

	//
	m = m_gmres;

	s = f_buffer[0];
	cs = f_buffer[1];
	sn = f_buffer[2];
	w = f_buffer[4];
	r = f_buffer[5];
	t = f_buffer[6];  // Buffer array

	int v_idx0 = 7;

	// Norm of Mb
	Precond(b, r);
	// Here Mb-->r
	normb = Norm(r);

// Norm of M r
	A->multiVec(x, w);  // Ax-->w
	for (l = 0; l < size_A; l++)
		r[l] = b[l] - w[l];  // r = b-Ax.

	Precond(r, w);  // Mr-->w
	beta = Norm(w);

	if (normb < DBL_MIN) normb = 1;

	// if ((error = Norm(r) / normb) <= tol)
	if ((error = beta / normb) <= tol)
	{
		Message();
		return 0;
	}

	iter = 1;
	while (iter <= max_iter)
	{
		v = f_buffer[v_idx0];
		for (l = 0; l < size_A; l++)
			v[l] = r[l] / beta;  //  r/beta
		for (l = 0; l < m + 1; l++)
			s[l] = 0.0;
		s[0] = beta;

		for (i = 0; i < m && iter <= max_iter; i++, iter++)
		{
			v = f_buffer[v_idx0 + i];
			A->multiVec(v, t);
			Precond(t, w);

			for (k = 0; k <= i; k++)
			{
				v_k = f_buffer[v_idx0 + k];
				H(k, i) = dot(w, v_k);

				for (l = 0; l < size_A; l++)
					w[l] -= H(k, i) * v_k[l];
			}
			H(i + 1, i) = Norm(w);
			v_k = f_buffer[v_idx0 + i + 1];
			for (l = 0; l < size_A; l++)
				v_k[l] = w[l] / H(i + 1, i);

			for (k = 0; k < i; k++)
				Set_Plane_Rotation(H(k, i), H(k + 1, i), cs[k], sn[k]);

			Get_Plane_Rotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
			Set_Plane_Rotation(H(i, i), H(i + 1, i), cs[i], sn[i]);
			Set_Plane_Rotation(s[i], s[i + 1], cs[i], sn[i]);

			if ((error = fabs(s[i + 1]) / normb) < tol)
			{
				Update(x, i, H, s);
				Message();
				return iter <= max_iter;
			}
		}

		Update(x, i - 1, H, s);
		A->multiVec(x, t);

		for (l = 0; l < size_A; l++)
			w[l] = b[l] - t[l];  // r = b-Ax.
		Precond(w, r);           // M*r

		beta = Norm(r);
		if ((error = beta / normb) < tol)
		{
			Message();
			return iter <= max_iter;
		}
	}

	Message();
	return iter <= max_iter;
}
//-----------------------------------------------------------------
//#endif // USE_MPI
#endif  // GMRES

#endif  // If not defined USE_MPI
#if defined(USE_MPI)
/*************************************************************************
   GeoSys-Function:
   Task: Parallel preconditioner
   Programming:
   12/2007 WW
 **************************************************************************/
void Linear_EQS::ComputePreconditioner_Jacobi()
{
	//
	A->DiagonalEntries(prec_M);
//
#if defined(NEW_BREDUCE)
	dom->ReduceBorderV(prec_M);
#else
	dom->Local2Border(prec_M, border_buffer0);  // to buffer
	MPI_Allreduce(border_buffer0, border_buffer1, A->Dof() * dom->BSize(),
	              MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
	dom->Border2Local(border_buffer1, prec_M);
}
/*************************************************************************
   GeoSys-Function:
   Task: Parallel preconditioner, inverse
   Programming:
   12/2007 WW
 **************************************************************************/
void Linear_EQS::Precond_Jacobi(const double* vec_s, double* vec_r)
{
	double val;
	for (long i = 0; i < A->Dim(); i++)
	{
		val = prec_M[i];
		if (val > DBL_EPSILON)
			vec_r[i] = vec_s[i] / val;
		else
			vec_r[i] = vec_s[i];
	}
}
#define TEST_MPII
/**************************************************************************
   Task: Linear equation::CG
   Programing:
   01/2008 WW/
**************************************************************************/
int Linear_EQS::CG(double* xg, const long n)
{
	long i, size, size_b, size_t;
	double rrM1;
	double* p, *r, *s;
	double* p_b, *r_b, *s_b;
	double* x_b;
	//
	size = A->Dim();
	size_b = dom->BSize() * A->Dof();
	// size_t = size + size_b; //[0, size_i): internal; [size_i, size_i+size_b):
	// border.
	//
	p = f_buffer[0];
	r = f_buffer[1];
	s = f_buffer[2];
	//

	//*** Norm b
	double bNorm_new;
	double buff_fl = dom->Dot_Interior(b, b);
	MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dom->Local2Border(b, border_buffer0);  // p_b s_b as buffer
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	// (rhs on border)
	buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
	bNorm_new = sqrt(buff_fl);
	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//*** r = b-Ax
	//    A*x
	dom->Global2Local(xg, x, n);
	A->multiVec(x, s);  // s as buffer
	for (i = 0; i < size; i++)
		r[i] = b[i] - s[i];  // r = b-Ax
	//   Collect border r
	dom->Local2Border(r, border_buffer0);  //
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	dom->Border2Local(border_buffer1, r);
	//
	// Preconditioning: M^{-1}r
	Precond(r, s);
	for (i = 0; i < size; i++)
		p[i] = s[i];
	// Check the convergence
	if ((error = Norm(r) / bNorm) < tol)
	{
		Message();
		return 1;
	}
	//
	double rr = dot(r, s);
	//
	for (iter = 1; iter <= max_iter; ++iter)
	{
		MatrixMulitVec(p, s);
		const double alpha = rr / dot(p, s);
		// Update
		for (i = 0; i < size; i++)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * s[i];
		}
		if ((error = Norm(r) / bNorm) < tol) break;
		// Preconditioner
		Precond(r, s);
		//
		rrM1 = rr;
		rr = dot(s, r);
		//
		const double beta = rr / rrM1;
		for (i = 0; i < size; i++)
			p[i] = s[i] + beta * p[i];
	}
	//
	// concancert internal x
	dom->CatInnerX(xg, x, n);
	//
	Message();
	return iter <= max_iter;
}

/*************************************************************************
   GeoSys-Function:
   Task: CGS solver
   Programming:
   02/2008 WW
 **************************************************************************/
int Linear_EQS::CGS(double* xg, const long n)
{
	long i, size, size_b, size_t;
	double rho_1, rho_2, alpha, beta;
	double* r0, *r, *p, *p_h, *q, *q_h, *v, *u, *u_h;
	double* r0_b, *r_b, *p_b, *p_h_b, *q_b, *q_h_b, *v_b, *u_b, *u_h_b;
	double* x_b;
	//
	size = A->Dim();
	size_b = dom->BSize() * A->Dof();
	// size_t = size + size_b; //[0, size_i): internal; [size_i, size_i+size_b):
	// border.
	r0 = f_buffer[0];
	r = f_buffer[1];
	p = f_buffer[2];
	p_h = f_buffer[3];
	q = f_buffer[4];
	q_h = f_buffer[5];
	v = f_buffer[6];
	u = f_buffer[7];
	u_h = f_buffer[8];
	//
	//
	rho_1 = rho_2 = 1.0;
	//
	//*** Norm b
	double bNorm_new;
	double buff_fl = dom->Dot_Interior(b, b);
	MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dom->Local2Border(b, border_buffer0);  //
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	//  (rhs on border)
	buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
	bNorm_new = sqrt(buff_fl);
	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//*** r = b-Ax
	//    A*x
	dom->Global2Local(xg, x, n);
	A->multiVec(x, v);  // v as buffer
	for (i = 0; i < size; i++)
		r0[i] = b[i] - v[i];  // r = b-Ax
	//   Collect border r
	dom->Local2Border(r0, border_buffer0);  //  buffer
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	dom->Border2Local(border_buffer1, r0);
	//
	for (i = 0; i < size; i++)
	{
		r[i] = r0[i];
		v[i] = 0.;
	}
	if ((error = Norm(r) / bNorm) < tol)
	{
		Message();
		return 0;
	}
	//
	for (iter = 1; iter <= max_iter; iter++)
	{
		rho_1 = dot(r0, r);
		if (fabs(rho_1) < DBL_MIN)  //  DBL_EPSILON
			break;
		//
		if (iter == 1)
			for (i = 0; i < size; i++)
				p[i] = u[i] = r[i];
		else
		{
			beta = rho_1 / rho_2;
			for (i = 0; i < size; i++)
			{
				u[i] = r[i] + beta * q[i];
				p[i] = u[i] + beta * (q[i] + beta * p[i]);
			}
		}
		// Preconditioner
		Precond(p, p_h);
		// A M^{-1}p-->v
		MatrixMulitVec(p_h, v);
		//
		alpha = rho_1 / dot(r0, v);
		//
		for (i = 0; i < size; i++)
		{
			q[i] = u[i] - alpha * v[i];
			q_h[i] = u[i] + q[i];
		}
		// Preconditioner
		Precond(q_h, u_h);
		for (i = 0; i < size; i++)
			x[i] += alpha * u_h[i];
		//
		MatrixMulitVec(u_h, q_h);
		//
		for (i = 0; i < size; i++)
			r[i] -= alpha * q_h[i];
		rho_2 = rho_1;
		if ((error = Norm(r) / bNorm) < tol) break;
	}
	// concancert internal x
	dom->CatInnerX(xg, x, n);
	//
	Message();
	//
	return iter <= max_iter;
}
#define TEST_MPII
/*************************************************************************
   GeoSys-Function:
   Task: Parallel BiCGStab solver
    xg:  Global solution
    n:  Size of x
   Programming:
   12/2007 WW
 **************************************************************************/
int Linear_EQS::BiCGStab(double* xg, const long n)
{
	long i, size, size_b;  //, size_t;
	double rho_0, rho_1, alpha, beta, omega, tt = 0., norm_v = 0.;
	double* r0, *r, *s, *s_h, *t, *v, *p, *p_h;
	double* r0_b, *r_b, *s_b, *s_h_b, *t_b, *v_b, *p_b, *p_h_b;
	double* x_b;
	//
	size = A->Dim();
	//
	size_b = dom->BSize() * A->Dof();
	// size_t = size + size_b; //[0, size_i): internal; [size_i, size_i+size_b):
	// border.
	r0 = f_buffer[0];
	r = f_buffer[1];
	s = f_buffer[2];
	s_h = f_buffer[3];
	t = f_buffer[4];
	v = f_buffer[5];
	p = f_buffer[6];
	p_h = f_buffer[7];
	//
	//
	rho_0 = alpha = omega = 1.0;

#ifdef TEST_MPI
	// TEST
	string test = "rank";
	char stro[64];
	sprintf(stro, "%d", myrank);
	string test1 = test + (string)stro + "Assemble.txt";
	ofstream Dum(test1.c_str(), ios::out);
	Dum.width(20);
	Dum.precision(15);
	Dum.setf(ios::scientific);

	Dum << "Time step: " << aktueller_zeitschritt << endl;
	Dum << "Norm b inner  " << dom->Dot_Interior(b, b) << endl;
//  if(A->Dof()==1)
//  Dum.close();
#endif

	//*** Norm b
	double bNorm_new;
	double buff_fl = dom->Dot_Interior(b, b);
	MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dom->Local2Border(b, border_buffer0);  // buffer
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	// (rhs on border)
	buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
	bNorm_new = sqrt(buff_fl);

	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//*** r = b-Ax
	//    A*x
	dom->Global2Local(xg, x, n);
	A->multiVec(x, s);  // s as buffer
	for (i = 0; i < size; i++)
		r0[i] = b[i] - s[i];  // r = b-Ax
	//   Collect border r
	dom->Local2Border(r0, border_buffer0);  // buffer
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	dom->Border2Local(border_buffer1, r0);

#ifdef TEST_MPI
	// TEST
	// if(A->Dof()>1)
	{
		Dum << " |r_0|= " << Norm(r0) << endl;
		Dum << " |b_0|= " << bNorm << endl;
		Dum << " x  " << endl;
		for (i = 0; i < n; i++)
			Dum << xg[i] << endl;
		/*
		   Dum<<" b  "<<endl;
		   for(i=0; i<size; i++)
		   Dum<<b[i]<<endl;
		   Dum<<"inter: r  Ax  "<<endl;
		   for(i=0; i<size; i++)
		   Dum<<r0[i]<<endl;
		   Dum<<"border: r  Ax  "<<endl;
		   for(i=0; i<size_b; i++)
		   Dum<<r0_b[i]<<endl;

		   Dum<<"inter: x  "<<endl;
		   for(i=0; i<size; i++)
		   Dum<<x[i]<<endl;

		   Dum<<"border: x  "<<endl;
		   for(i=0; i<size_b; i++)
		   Dum<<x_b[i]<<endl;
		 */
	}
#endif

	// Initial. [0, size_i): internal; [size_i, size_i+size_b): border.
	for (i = 0; i < size; i++)
	{
		r[i] = r0[i];
		v[i] = 0.;
		p[i] = 0.;
	}
	if ((error = Norm(r) / bNorm) < tol)
	{
		if (myrank == 0)  // Make screen output only by processor 0
			Message();
		return 0;
	}

	/*
	   //TEST
	   //if(A->Dof()>1)
	   {
	    Dum<<" bNorm  "<<bNorm <<"  Norm(r) "<<  Norm(r) <<endl;

	   }
	 */

	//
	for (iter = 1; iter <= max_iter; iter++)
	{
		rho_1 = dot(r0, r);
		if (fabs(rho_1) < DBL_MIN) break;

#ifdef TEST_MPI
		// TEST
		//  if(A->Dof()>1)
		Dum << " rho_1  " << rho_1 << endl;
#endif

		if (iter == 1)
			// p[0, size_i): internal; p[size_i, size_i+size_b)-->p_b: border.
			for (i = 0; i < size; i++)
				p[i] = r[i];
		else
		{
			beta = (rho_1 / rho_0) * (alpha / omega);
			// [0, size_i): internal; [size_i, size_i+size_b): border.
			for (i = 0; i < size; i++)
				p[i] = r[i] + beta * (p[i] - omega * v[i]);
		}
		// Preconditioner
		Precond(p, p_h);
		// A M^{-1}p-->v
		MatrixMulitVec(p_h, v);
		//
		alpha = rho_1 / dot(r0, v);

#ifdef TEST_MPI
		// TEST
		//  if(A->Dof()>1)
		{
			// TEST
			Dum << "  alpha  " << alpha << " dot(r0, v) " << dot(r0, v)
			    << " dot(r0, r) " << dot(r0, r) << endl;

			Dum << "\n r0, r,  v   " << endl;
			//  for(i=0; i<size_t; i++)
			//   Dum<<r0[i]<<"       "<<r[i]<<"       "<<v[i]<<endl;
		}
#endif

		//
		for (i = 0; i < size; i++)
			s[i] = r[i] - alpha * v[i];
		norm_v = sqrt(dot(s, s));
		if ((error = norm_v / bNorm) < tol)
		{
			for (i = 0; i < size; i++)
				x[i] += alpha * p_h[i];
			break;
		}

#ifdef TEST_MPI
/*
   //TEST
   if(A->Dof()>1)
   {
   //TEST
   Dum<<"\n  norm_v/bNorm  "<<error <<endl;
   for(i=0; i<size_t; i++)
   Dum<<s[i]<<endl;
   exit(0);
   }
 */
#endif

		//  M^{-1}s,
		Precond(s, s_h);
		// A* M^{-1}s
		MatrixMulitVec(s_h, t);
		//
		tt = dot(t, t);

#ifdef TEST_MPI
		// TEST
		// TEST
		//  if(A->Dof()>1)
		Dum << "  tt  " << tt << endl;
#endif

		if (tt > DBL_MIN)
			omega = dot(t, s) / tt;
		else
			omega = 1.0;
		// Update solution
		for (i = 0; i < size; i++)
		{
			x[i] += alpha * p_h[i] + omega * s_h[i];
			r[i] = s[i] - omega * t[i];
		}
		rho_0 = rho_1;
		//
		norm_v = sqrt(dot(r, r));

#ifdef TEST_MPI
		// TEST
		// TEST
		// if(A->Dof()>1)
		{
			Dum << " sqrt(dot(r,r))  " << norm_v << endl;
			// exit(0);
		}
#endif

		if ((error = norm_v / bNorm) < tol) break;
		if (fabs(omega) < DBL_MIN)
		{
			error = norm_v / bNorm;
			break;
		}

#ifdef TEST_MPI
// TEST
// Dum.close();
// MPI_Finalize();
// exit(0);
#endif
	}
	//
	// concancert internal x
	dom->CatInnerX(xg, x, n);
// Form the local to global
// dom->Border2Global(x_b, xg, n);

#ifdef TEST_MPI
	// if(A->Dof()>1)
	{
		Dum << " x " << endl;
		for (i = 0; i < n; i++)
			Dum << xg[i] << endl;
		Dum.close();
		// exit(0);
	}
#endif

	//
	Message();
	//
	//  return iter <= max_iter;
	return iter;
}

//#define  CG_test
/*************************************************************************
   GeoSys-Function:
   Task: Parallel BiCG solver
    xg:  Global solution
    n:  Size of x
   Programming:
   08/2008 WW
 **************************************************************************/
int Linear_EQS::BiCG(double* xg, const long n)
{
	long i, size, size_b;
	double rho1, rho2, alpha, beta;
	double* z, *zt, *p, *pt, *q, *qt, *r, *rt;
	//
	size = A->Dim();
	//
	size_b = dom->BSize() * A->Dof();
	z = f_buffer[0];
	zt = f_buffer[1];
	p = f_buffer[2];
	pt = f_buffer[3];
	q = f_buffer[4];
	qt = f_buffer[5];
	r = f_buffer[6];
	rt = f_buffer[7];
	//
	//*** Norm b
	double bNorm_new;
	double buff_fl = dom->Dot_Interior(b, b);
	MPI_Allreduce(&buff_fl, &bNorm_new, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dom->Local2Border(b, border_buffer0);  //
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	//  (rhs on border)
	buff_fl = bNorm_new + dot(border_buffer1, border_buffer1, size_b);
	bNorm_new = sqrt(buff_fl);
	// Check if the norm of b is samll enough for convengence

	// Check if the norm of b is samll enough for convengence
	if (CheckNormRHS(bNorm_new)) return 0;
	//*** r = b-Ax
	//    A*x
	dom->Global2Local(xg, x, n);
	A->multiVec(x, rt);  // rt as buffer
	for (i = 0; i < size; i++)
		r[i] = b[i] - rt[i];  // r = b-Ax
	//   Collect border r
	dom->Local2Border(r, border_buffer0);  //  buffer
	MPI_Allreduce(border_buffer0, border_buffer1, size_b, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	dom->Border2Local(border_buffer1, r);
	//
	// Initial.
	for (i = 0; i < size; i++)
	{
		r[i] = b[i] - rt[i];
		rt[i] = r[i];
	}
	if ((error = Norm(r) / bNorm) < tol)
	{
		if (myrank == 0)  // Make screen output only by processor 0
			Message();
		return 0;
	}

#ifdef CG_test
	// TEST
	string test = "rank";
	char stro[64];
	sprintf(stro, "%d", myrank);
	string test1 = test + (string)stro + "_Assemble.txt";
	ofstream Dum(test1.c_str(), ios::out);
	Dum.width(20);
	Dum.precision(15);
	Dum.setf(ios::scientific);

	Dum << " Norm(r) " << Norm(r) << " bNorm " << bNorm << endl;
#endif

	//
	for (iter = 1; iter <= max_iter; iter++)
	{
		Precond(r, z);
		TransPrecond(rt, zt);
		rho1 = dot(z, rt);

#ifdef CG_test
		Dum << " rho1 " << rho1 << endl;
#endif

		//
		if (fabs(rho1) < DBL_MIN)
		{
			Message();
			break;
		}
		//
		if (iter == 1)
			for (i = 0; i < size; i++)
			{
				p[i] = z[i];
				pt[i] = zt[i];
			}
		else
		{
			beta = rho1 / rho2;
			for (i = 0; i < size; i++)
			{
				p[i] = z[i] + beta * p[i];
				pt[i] = zt[i] + beta * pt[i];
			}
		}
		MatrixMulitVec(p, q);
		TransMatrixMulitVec(pt, qt);
		alpha = rho1 / dot(pt, q);

#ifdef CG_test
		Dum << " alpha " << alpha << endl;
#endif

		for (i = 0; i < size; i++)
		{
			x[i] += alpha * p[i];
			r[i] -= alpha * q[i];
			rt[i] -= alpha * qt[i];
		}
		//
		rho2 = rho1;
		if ((error = Norm(r) / bNorm) < tol) break;

#ifdef CG_test
		Dum << " error = Norm(r)/bNorm " << error << endl;
#endif

		//
	}
	//
	// concancert internal x
	dom->CatInnerX(xg, x, n);

	//
	Message();
#ifdef CG_test
	exit(0);
#endif

	//
	return iter <= max_iter;
}
#endif
//------------------------------------------------------------------------
}  // namespace
#endif  // if defined(NEW_EQS)
