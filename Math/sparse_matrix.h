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
	// Operator
	void operator=(const double a);
	void operator*=(const double a);
	void operator+=(const double a);
	void operator=(const CSparseMatrix& m);
	void operator+=(const CSparseMatrix& m);
	void operator-=(const CSparseMatrix& m);
	// Vector pass through augment and bring results back.
	void Diagonize(const long idiag, const double b_given, double* b);
	//
	// Access to members
	double& operator()(const long i, const long j = 0) const;
	//
	long Dim() const { return DOF * rows; }
	int Dof() const { return DOF; }
	void SetDOF(const int dof_n)
	{
		DOF = dof_n;
	}
	long Size() const { return rows; }

	IndexType nnz() const
	{
		return DOF * DOF * size_entry_column;
	}

	IndexType* ptr;
	IndexType* col_idx;
	IndexType* entry_index;

	int GetCRSValue(double* value);
	// Print
	void Write(std::ostream& os = std::cout);
	void Write_BIN(std::ostream& os);
private:
	// Data
	double* entry;
	mutable double zero_e;
	//
	bool symmetry;
	// Topology mapping from data array to matrix. All are only pointers to the
	// correpinding members in SparseTable, and no memory are allocated for them
	long* entry_column;
	long* num_column_entries;     // number of entries of each columns in sparse
	                              // table
	long* diag_entry;
	long size_entry_column;
	long rows;
	//
	int DOF;
};
// Since the pointer to member funtions gives lower performance

} // Math_Group

#endif // sparse_table_INC
