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

#ifdef USE_PARALUTION
#include <paralution.hpp>
#endif

#include "Configure.h"
#include "sparse_table.h"
#include "sparse_matrix.h"

class CNumerics;

namespace Math_Group
{
class SparseTable;

class Linear_EQS
{
#ifndef OGS_USE_LONG
	typedef int IndexType;
#else
	typedef long long int IndexType;
#endif
public:
	Linear_EQS(const SparseTable& sparse_table, const long dof);
	~Linear_EQS();

	void ConfigNumerics(CNumerics* m_num, const long n = 0);

	int Solver(CNumerics* num = NULL, bool compress = false);

	//
	void Initialize();
	void Clean();
	//
	// Access to the members
	void SetDOF(const int dof_n)  // For different processes with different DOF
	                              // of OPDE. _new. 02/2010. WW
	{
		A->SetDOF(dof_n);
	}
	void SetKnownX_i(const long i, const double x_i);
	double X(const long i) const { return x[i]; }
	const double* getX() const { return x; }
	const double* getRHS() const { return b; }
	double RHS(const long i) const { return b[i]; }
	double NormX();
	double ComputeNormRHS() { return Norm(b); }
	double NormRHS() { return bNorm; }

	// Write
	void Write(std::ostream& os = std::cout);
	void WriteRHS(std::ostream& os = std::cout);
	void WriteX(std::ostream& os = std::cout);
	void Write_BIN(std::ostream& os);

private:
	CSparseMatrix* A;
	double* b;
	double* x;

#ifdef LIS
	LIS_MATRIX AA;
	LIS_VECTOR bb, xx;
	LIS_SOLVER solver;
#endif
#ifdef USE_PARALUTION
	paralution::LocalMatrix<double> AA;
	paralution::LocalVector<double> bb;
	paralution::LocalVector<double> xx;
#endif

	// Controls
	int precond_type;
	int solver_type;
	int iter, max_iter;
	double tol, bNorm, error;
	long size_A;

	// Operators
	double dot(const double* xx, const double* yy);
	inline double Norm(const double* x_) { return sqrt(dot(x_, x_)); }
	inline bool CheckNormRHS(const double normb_new);

#ifdef MKL
	void solveWithPARDISO(CNumerics* num, bool compress);
#endif
#ifdef LIS
	int solveWithLIS(CNumerics* num, bool compress);
#endif
#ifdef USE_PARALUTION
	int solveWithParalution(CNumerics* num, bool compress);
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

#endif // eqs_class_INC
