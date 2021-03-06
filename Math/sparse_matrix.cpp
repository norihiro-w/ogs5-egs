/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "sparse_matrix.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iomanip>

#include "binarySearch.h"

namespace Math_Group
{

/*\!
 ********************************************************************
   Constructor of sparse matrix
   Arguments:
      sparse_table: Sparse graph
      dof:  Degree of freedom given by PDE
   08/2007 WW
   10/2007 WW
   02/2008 PCH Compressed Row Storage
 ********************************************************************/
CSparseMatrix::CSparseMatrix(const SparseTable& sparse_table, const int dof)
    : DOF(dof)
{
	symmetry = sparse_table.symmetry;
	size_entry_column = sparse_table.size_entry_column;
	rows = sparse_table.rows;
	// Topology mapping from data array to matrix
	// Only refer address
	entry_column = sparse_table.entry_column;
	num_column_entries = sparse_table.num_column_entries;
	diag_entry = sparse_table.diag_entry;
	// Values of all sparse entries
	entry = new double[dof * dof * size_entry_column + 1];
	entry[dof * dof * size_entry_column] = 0.;
	zero_e = 0.;

	IndexType counter_ptr = 0, counter_col_idx = 0;
	IndexType i=0, k=0, ii=0, jj=0;

	ptr = new IndexType[rows * dof + 1];
	col_idx = new IndexType[dof * dof * size_entry_column];
	entry_index = new IndexType[dof * dof * size_entry_column];

	for (ii = 0; ii < DOF; ii++)
	{
		const long row_offset_dof = ii * rows;
		for (i = 0; i < rows; i++)
		{
			ptr[row_offset_dof + i] = counter_ptr;
			const long n_columns =
				num_column_entries[i + 1] - num_column_entries[i];
			const long col_offset = num_column_entries[i];
			for (jj = 0; jj < DOF; jj++)
			{
				const long col_offset_dof = jj * rows;
				for (k = 0; k < n_columns; k++)
				{
					// column in global matrix
					col_idx[counter_col_idx] =
						col_offset_dof + entry_column[col_offset + k];

					++counter_ptr;
					++counter_col_idx;
				}
			}
		}
	}
	ptr[i + rows * (ii - 1)] = counter_ptr;
	// entry index: actual memory layout is not CRS
	long cnt = 0;
	for (int ii = 0; ii < DOF; ii++)
	{
		for (i = 0; i < rows; i++)
		{
			const long ptr0 = num_column_entries[i];
			const long ptr1 = num_column_entries[i + 1];
			for (int jj = 0; jj < DOF; jj++)
			{
				const long offset_dof = (ii * DOF + jj) * size_entry_column;
				for (long k = ptr0; k < ptr1; k++)
				{
					entry_index[cnt++] = offset_dof + k;
				}
			}
		}
	}

}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
CSparseMatrix::~CSparseMatrix()
{
	delete[] entry;
	entry = NULL;

	delete[] ptr;
	ptr = NULL;
	delete[] col_idx;
	col_idx = NULL;
	delete[] entry_index;
	entry_index = NULL;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
double& CSparseMatrix::operator()(const long i, const long j) const
{
#ifdef gDEBUG
	if (i >= rows * DOF || j >= rows * DOF)
	{
		std::cout << "\n Index exceeds the dimension of the matrix"
		          << "\n";
		abort();
	}
#endif
	long ii, jj, ir, jr, k;
	ii = i;
	jj = j;
	if (symmetry)
		if (i > j)
		{
			ii = j;
			jj = i;
		}
	ir = ii % rows;
	jr = jj % rows;
	ii /= rows;
	jj /= rows;
	//
	k = -1;

	/// Left boundary of this row: num_column_entries[ir]
	/// Right boundary of this row: num_column_entries[ir+1]
	/// Search target is jr
	k = binarySearch(entry_column,
					 jr,
					 num_column_entries[ir],
					 num_column_entries[ir + 1]);
	if (k == -1) return zero_e;

	k = (ii * DOF + jj) * size_entry_column + k;

	return entry[k];  //
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator=(const double a)
{
	long size = DOF * DOF * size_entry_column;
	for (long i = 0; i < size; i++)
		entry[i] = a;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator*=(const double a)
{
	long size = DOF * DOF * size_entry_column;
	for (long i = 0; i < size; i++)
		entry[i] *= a;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator+=(const double a)
{
	long size = DOF * DOF * size_entry_column;
	for (long i = 0; i < size; i++)
		entry[i] += a;
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator=(const CSparseMatrix& m)
{
	long size = DOF * DOF * size_entry_column;
#ifdef gDEBUG
	if (size != m.DOF * m.DOF * m.size_entry_column)
	{
		std::cout << "\n Dimensions of two matrices do not match"
		          << "\n";
		abort();
	}
#endif
	for (long i = 0; i < size; i++)
		entry[i] = m.entry[i];
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator+=(const CSparseMatrix& m)
{
	long size = DOF * DOF * size_entry_column;
#ifdef gDEBUG
	if (size != m.DOF * m.DOF * m.size_entry_column)
	{
		std::cout << "\n Dimensions of two matrices do not match"
		          << "\n";
		abort();
	}
#endif
	for (long i = 0; i < size; i++)
		entry[i] += m.entry[i];
}
/*\!
 ********************************************************************
   Desstructor of sparse matrix
   08/2007 WW
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::operator-=(const CSparseMatrix& m)
{
	long size = DOF * DOF * size_entry_column;
#ifdef gDEBUG
	if (size != m.DOF * m.DOF * m.size_entry_column)
	{
		std::cout << "\n Dimensions of two matrices do not match"
		          << "\n";
		abort();
	}
#endif
	for (long i = 0; i < size; i++)
		entry[i] -= m.entry[i];
}
/*\!
 ********************************************************************
   Output sparse matrix
   08/2007 WW
   10/2007 WW
   03/2011 WW  CRS
 ********************************************************************/
void CSparseMatrix::Write(std::ostream& os)
{
	//
	long i, k, ii, jj;
	os << "*** Non-zero entries of matrix:  "
	   << "\n";
	os.width(14);
	os.precision(8);
	//
	os << "Storage type: CRS\n";
	for (ii = 0; ii < DOF; ii++)
	{
		for (i = 0; i < rows; i++)
		{
			const long ptr0 = num_column_entries[i];
			const long ptr1 = num_column_entries[i + 1];
			for (jj = 0; jj < DOF; jj++)
			{
				const long offset_dof = (ii * DOF + jj) * size_entry_column;
				for (k = ptr0; k < ptr1; k++)
				{
					// TEST
					// if(fabs(entry[(ii*DOF+jj)*size_entry_column+counter])>DBL_MIN)
					// //DBL_EPSILON)
					os << std::setw(10) << ii * rows + i << " "
					   << std::setw(10) << jj * rows + entry_column[k]
					   << " " << std::setw(15) << entry[offset_dof + k]
					   << "\n";
				}
			}
		}
	}

}
//--------------------------------------------------------------
/*!
   \brief Write matrix to a binary file

   03.2011. WW
 */
void CSparseMatrix::Write_BIN(std::ostream& os)
{
	//
	if (DOF == 1)
	{
		os.write((char*)&rows, sizeof(long));
		os.write((char*)num_column_entries, (rows + 1) * sizeof(long));
		os.write((char*)entry_column, num_column_entries[rows] * sizeof(long));
		os.write((char*)entry, num_column_entries[rows] * sizeof(double));
	}
	else
	{
		long i, k, ii, jj, size;
		long* ptr;
		long* A_index;
		double* A_value;

		ptr = new long[DOF * rows + 1];
		size = DOF * DOF * num_column_entries[rows];
		A_index = new long[size];
		A_value = new double[size];

		long counter = 0;

		for (ii = 0; ii < DOF; ii++)
			for (i = 0; i < rows; i++)
			{
				ptr[ii * rows + i] = counter;
				for (jj = 0; jj < DOF; jj++)
					for (k = num_column_entries[i];
					     k < num_column_entries[i + 1];
					     k++)
					{
						A_index[counter] = jj * rows + entry_column[k];
						A_value[counter] =
						    entry[(ii * DOF + jj) * size_entry_column + k];
						counter++;
					}
			}
		ptr[DOF * rows] = counter;

		ii = DOF * rows;
		os.write((char*)&ii, sizeof(long));
		os.write((char*)ptr, (ii + 1) * sizeof(long));
		os.write((char*)A_index, size * sizeof(long));
		os.write((char*)A_value, size * sizeof(double));

		delete[] ptr;
		delete[] A_index;
		delete[] A_value;
	}
}

/*\!
 ********************************************************************
   Set
        A(ii,ii) = x_i,
        A(ii, j) = 0., j!=ii
        A(i, ii) = 0., i!=ii
        b_i -= A(i,k)b_k  // b_k is given
   Programm:
   10/2007 WW
   03/2011 WW  CRS storage
 ********************************************************************/
void CSparseMatrix::Diagonize(const long idiag, const double b_given, double* b)
{
	//
	double vdiag = 0.;
	long j, k, ii, jj, j0;
	long id = idiag % rows;

	ii = idiag / rows;

	const long row_end = num_column_entries[id + 1];
	/// Diagonal entry and the row where the diagonal entry exists
	j = diag_entry[id];
	vdiag = entry[(ii * DOF + ii) * size_entry_column + j];
	/// Row where the diagonal entry exists
	for (jj = 0; jj < DOF; jj++)
	{
		const long ij = (ii * DOF + jj) * size_entry_column;
		for (k = num_column_entries[id]; k < row_end; k++)
		{
			j0 = entry_column[k];
			if (id == j0 && jj == ii)  // Diagonal entry
				continue;
			entry[ij + k] = 0.;
		}
	}
#ifdef colDEBUG
	/// Clean column id
	for (i = 0; i < rows; i++)
	{
		j = binarySearch(entry_column,
						 id,
						 num_column_entries[i],
						 num_column_entries[i + 1]);
		if (j == -1) continue;
		j0 = entry_column[j];

		for (jj = 0; jj < DOF; jj++)
		{
			if (i == j0 && ii == jj) continue;
			k = (jj * DOF + ii) * size_entry_column + j;
			b[jj * rows + i] -= entry[k] * b_given;
			entry[k] = 0.;
			// Room for symmetry case
		}
	}
#endif
	b[idiag] = vdiag * b_given;
}

/********************************************************************
   Get sparse matrix values in compressed row storage
   Programm:
   02/2008 PCH
********************************************************************/
int CSparseMatrix::GetCRSValue(double* value)
{
	int success = 1;
	int i;

#pragma omp parallel for
	for (i = 0; i < size_entry_column * DOF * DOF; ++i)
		value[i] = entry[entry_index[i]];

	return success;
}

}  // Namespace
