/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef sparse_matrix_INC
#define sparse_matrix_INC

#include <iostream>
#include "sparse_table.h"

namespace Math_Group
{

class CSparseMatrix
{
public:
#ifndef OGS_USE_LONG
	typedef int IndexType;
#else
	typedef long long int IndexType;
#endif

	CSparseMatrix(const SparseTable& sparse_table, const int dof);
	~CSparseMatrix();
	// Preconditioner
	void Precond_Jacobi(double* vec_s, double* vec_r);
	// TEMP
	void Precond_ILU(double* /*vec_s*/, double* /*vec_r*/) {}
	// Operator
	void operator=(const double a);
	void operator*=(const double a);
	void operator+=(const double a);
	void operator=(const CSparseMatrix& m);
	void operator+=(const CSparseMatrix& m);
	void operator-=(const CSparseMatrix& m);
	// Vector pass through augment and bring results back.
	void multiVec(double* vec_s, double* vec_r);
	void Trans_MultiVec(double* vec_s, double* vec_r);
	void Diagonize(const long idiag, const double b_given, double* b);
	//
	// Access to members
	double& operator()(const long i, const long j = 0) const;
	//
	StorageType GetStorageType() const { return storage_type; }  // 05.2011. WW
	long Dim() const { return DOF * rows; }
	int Dof() const { return DOF; }
	void SetDOF(const int dof_n)  //_new. 02/2010. WW
	{
		DOF = dof_n;
	}
	long Size() const { return rows; }
#if defined(LIS) || defined(MKL) || defined(USE_PARALUTION)
	// These two pointers are in need for Compressed Row Storage
	IndexType nnz() const  // PCH
	{
		return DOF * DOF * size_entry_column;
	}
	IndexType* ptr;
	IndexType* col_idx;
	IndexType* entry_index;
	int GetCRSValue(double* value);
#endif
	// Print
	void Write(std::ostream& os = std::cout);
	void Write_BIN(std::ostream& os);
private:
	// Data
	double* entry;
	mutable double zero_e;
	/// 0. 03.2011. WW
	StorageType storage_type;
	//
	bool symmetry;
	// Topology mapping from data array to matrix. All are only pointers to the
	// correpinding members in SparseTable, and no memory are allocated for them
	long* entry_column;
	long* num_column_entries;     // number of entries of each columns in sparse
	                              // table
	long* row_index_mapping_n2o;  // Row index of sparse table to row index of
	                              // matrix
	long* row_index_mapping_o2n;  // Inverse of last
	long* diag_entry;
	long size_entry_column;
	long max_columns;
	long rows;
	//
	int DOF;
};
// Since the pointer to member funtions gives lower performance

} // Math_Group

#endif // sparse_table_INC
