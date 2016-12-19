/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef eqs_class_INC
#define eqs_class_INC

#include <cmath>
#include <iostream>
#include <vector>

#ifdef LIS
#include <lis.h>
#endif

#include "sparse_matrix.h"

namespace Math_Group
{
using namespace std;

struct SparseTable;
//
class Linear_EQS
{
#ifndef OGS_USE_LONG
	typedef int IndexType;
#else
	typedef long long int IndexType;
#endif
public:
	Linear_EQS(const SparseTable& sparse_table, const long dof,
	           bool messg = true);
	~Linear_EQS();

	void ConfigNumerics(int ls_precond, int ls_method, int ls_max_iterations, double ls_error_tolerance, int storage_type, std::string const& extra_arg);
	int Solver(bool compress = false);
	//
	void Initialize();
	void Clean();

	void SetDOF(const int dof_n)
	{
		A->SetDOF(dof_n);
	}
	void SetKnownX_i(const long i, const double x_i);
	double X(const long i) const { return x[i]; }
	const double* getX() const { return x; }
	double* getX() { return x; }
	const double* getRHS() const { return b; }
	double* getRHS() { return b; }
	const CSparseMatrix* getA() const { return A; }
	CSparseMatrix* getA() { return A; }
	double RHS(const long i) const { return b[i]; }
	double NormX();
	double ComputeNormRHS() { return Norm(b); }
	// Write
	void Write(std::ostream& os = std::cout);
	void WriteRHS(std::ostream& os = std::cout);
	void WriteX(std::ostream& os = std::cout);
	void Write_BIN(ostream& os);

private:
	CSparseMatrix* A;
	double* b;
	double* x;
//
#ifdef LIS
	// lis solver interface starts here
	LIS_MATRIX AA;
	LIS_VECTOR bb, xx;
	LIS_SOLVER solver;
#endif

	// Controls
	int precond_type;
	int solver_type;
	int max_iter;
	double tol;
	long size_A;
	int storage_type;
	std::string extra_arg;

	// Operators
	double dot(const double* xx, const double* yy);
	inline double Norm(const double* xx) { return sqrt(dot(xx, xx)); }
#ifdef MKL
	void solveWithPARDISO(bool compress);
#endif
#ifdef LIS
	int solveWithLIS(bool compress);
#endif
	IndexType searcgNonZeroEntries(IndexType nrows, IndexType* ptr,
	                               double* value,
	                               std::vector<IndexType>& vec_nz_rows,
	                               std::vector<IndexType>& vec_z_rows);
	void compressCRS(const IndexType* org_ptr, const IndexType* org_col_idx,
	                 const double* org_value,
	                 const std::vector<IndexType>& vec_nz_rows,
	                 const std::vector<IndexType>& vec_z_rows,
	                 IndexType n_nz_entries, IndexType*& new_ptr,
	                 IndexType*& new_col_index, double*& new_value);

};
}

using Math_Group::Linear_EQS;
extern std::vector<Math_Group::Linear_EQS*> EQS_Vector;
extern std::vector<Math_Group::SparseTable*> SparseTable_Vector;

#endif
