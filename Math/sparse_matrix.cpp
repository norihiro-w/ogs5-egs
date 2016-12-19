/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*========================================================================
   GeoSys - class Matrix (Definition)
		  class vec
   Task:       Matrix object for full matrices.  If the size of matrix is
			 small, this class can do efficient matrix operation.
   Function:   See the definition below
   programming:
   22/08/2004  WW
   ==========================================================================*/

#include "sparse_matrix.h"

#include <cfloat>
#include <cmath>
#include <iomanip>
#include <algorithm>

#include "mathlib.h"

namespace Math_Group
{

////////////////////////////////////////////////////////////
#ifdef NEW_EQS

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
	max_columns = sparse_table.max_columns;
	rows = sparse_table.rows;
	storage_type = sparse_table.storage_type;  // WW
	// Topology mapping from data array to matrix
	// Only refer address
	entry_column = sparse_table.entry_column;
	num_column_entries = sparse_table.num_column_entries;
	row_index_mapping_n2o = sparse_table.row_index_mapping_n2o;
	row_index_mapping_o2n = sparse_table.row_index_mapping_o2n;
	diag_entry = sparse_table.diag_entry;
	// Values of all sparse entries
	entry = new double[dof * dof * size_entry_column + 1];
	entry[dof * dof * size_entry_column] = 0.;
	zero_e = 0.;
//
#ifdef LIS  // PCH
	IndexType counter, counter_ptr = 0, counter_col_idx = 0;
	IndexType i, k, ii, jj, J, K;
	IndexType row_in_sparse_table;

	ptr = new IndexType[rows * dof + 1];
	col_idx = new IndexType[dof * dof * size_entry_column];
	entry_index = new IndexType[dof * dof * size_entry_column];

	if (storage_type == JDS)
	{
		for (ii = 0; ii < DOF; ii++)
			for (i = 0; i < rows; i++)
			{
				// Store ptr arrary for CRS
				ptr[i + rows * ii] = counter_ptr;
				row_in_sparse_table = row_index_mapping_o2n[i];
				for (jj = 0; jj < DOF; jj++)
				{
					counter = row_in_sparse_table;
					for (k = 0; k < max_columns; k++)
					{
						if (row_in_sparse_table < num_column_entries[k])
						{
							// I = ii * rows + i; // row in global matrix
							// column in global matrix
							J = jj * rows + entry_column[counter];
							K = (ii * DOF + jj) * size_entry_column + counter;

							// Store column index for CRS
							col_idx[counter_col_idx] = J;
							entry_index[counter_col_idx] = K;

							++counter_ptr;
							++counter_col_idx;
							counter += num_column_entries[k];
						}
						else
							break;
					}
				}
			}
		ptr[i + rows * (ii - 1)] = counter_ptr;
	}
	else if (storage_type == CRS)
	{
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

#endif
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

#ifdef LIS  // PCH
	delete[] ptr;
	ptr = NULL;
	delete[] col_idx;
	col_idx = NULL;
	delete[] entry_index;
	entry_index = NULL;
#endif
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

	if (storage_type == JDS)
	{
		long row_in_parse_table, counter;
		row_in_parse_table = row_index_mapping_o2n[ir];
		counter = row_in_parse_table;
		for (k = 0; k < max_columns; k++)
		{
			if (row_in_parse_table >= num_column_entries[k]) return zero_e;
			if (entry_column[counter] == jr) break;  // Found the entry
			counter += num_column_entries[k];
		}
		if (counter >= size_entry_column) return zero_e;
		//  Zero entry;
		k = (ii * DOF + jj) * size_entry_column + counter;
	}
	else if (storage_type == CRS)
	{
		/// Left boundary of this row: num_column_entries[ir]
		/// Right boundary of this row: num_column_entries[ir+1]
		/// Search target is jr
		k = binarySearch(entry_column,
						 jr,
						 num_column_entries[ir],
						 num_column_entries[ir + 1]);
		if (k == -1) return zero_e;

		k = (ii * DOF + jj) * size_entry_column + k;
	}

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
	long i, k, ii, jj, row_in_parse_table, counter;
	os << "*** Non-zero entries of matrix:  "
	   << "\n";
	os.width(14);
	os.precision(8);
	//
	if (storage_type == CRS)
	{
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
	else if (storage_type == JDS)
	{
		os << "Storage type: JDS\n";
		for (ii = 0; ii < DOF; ii++)
			for (i = 0; i < rows; i++)
			{
				row_in_parse_table = row_index_mapping_o2n[i];
				for (jj = 0; jj < DOF; jj++)
				{
					counter = row_in_parse_table;
					for (k = 0; k < max_columns; k++)
					{
						if (row_in_parse_table < num_column_entries[k])
						{
							// TEST
							// if(fabs(entry[(ii*DOF+jj)*size_entry_column+counter])>DBL_MIN)
							// //DBL_EPSILON)
							os << std::setw(10) << ii * rows + i << " "
							   << std::setw(10)
							   << jj * rows + entry_column[counter] << " "
							   << std::setw(15)
							   << entry[(ii * DOF + jj) * size_entry_column +
										counter] << "\n";
							counter += num_column_entries[k];
						}
						else
							break;
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
	if (storage_type == JDS) return;
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

/********************************************************************
   Perform A*x
   Arguments:
	 vec_sr: M*vec_s-->vec_r
   01/2006 WW
   08/2007 WW
   10/2007 WW
   03/2011 WW      CRS storage
********************************************************************/
void CSparseMatrix::multiVec(double* vec_s, double* vec_r)
{
	long i, j, k, ii, jj, kk, ll, idof, jdof, counter;
	for (i = 0; i < rows * DOF; i++)
		vec_r[i] = 0.0;
	//
	counter = 0;
	if (DOF > 1)
	{
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1];
					 j++)
				{
					jj = entry_column[j];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							k = (idof * DOF + jdof) * size_entry_column + j;
							vec_r[kk] += entry[k] * vec_s[ll];
							if (symmetry & (kk != ll))
								vec_r[ll] += entry[k] * vec_s[kk];
						}
					}
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							j = (idof * DOF + jdof) * size_entry_column +
								counter;
							vec_r[kk] += entry[j] * vec_s[ll];
							if (symmetry & (kk != ll))
								vec_r[ll] += entry[j] * vec_s[kk];
						}
					}
					counter++;
				}
		}
	}
	else  // DOF = 1
	{
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1];
					 j++)
				{
					jj = entry_column[j];
					vec_r[ii] += entry[j] * vec_s[jj];
					if (symmetry & (ii != jj))
						vec_r[jj] += entry[j] * vec_s[ii];
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					vec_r[ii] += entry[counter] * vec_s[jj];
					if (symmetry & (ii != jj))
						vec_r[jj] += entry[counter] * vec_s[ii];
					counter++;
				}
		}
	}
}

/*\!
 ********************************************************************
   Perform A^T*x
   Arguments:
	  vec_sr: M^T*vec_s-->vec_r
   10/2010 WW
   03/2011 WW      CRS storage
 ********************************************************************/
void CSparseMatrix::Trans_MultiVec(double* vec_s, double* vec_r)
{
	long i, j, k, ii, jj, kk, ll, idof, jdof, counter;
	for (i = 0; i < rows * DOF; i++)
		vec_r[i] = 0.0;
	//
	counter = 0;
	if (DOF > 1)
	{
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1];
					 j++)
				{
					jj = entry_column[j];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							k = (idof * DOF + jdof) * size_entry_column + j;
							vec_r[ll] += entry[k] * vec_s[kk];
							if (symmetry & (kk != ll))
								vec_r[kk] += entry[k] * vec_s[ll];
						}
					}
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					for (idof = 0; idof < DOF; idof++)
					{
						kk = idof * rows + ii;
						for (jdof = 0; jdof < DOF; jdof++)
						{
							ll = jdof * rows + jj;
							j = (idof * DOF + jdof) * size_entry_column +
								counter;
							vec_r[ll] += entry[j] * vec_s[kk];
							if (symmetry & (kk != ll))
								vec_r[kk] += entry[j] * vec_s[ll];
						}
					}
					counter++;
				}
		}
	}
	else  // DOF = 1
	{
		if (storage_type == CRS)
		{
			/// ptr is num_column_entries
			for (ii = 0; ii < rows; ii++)
				for (j = num_column_entries[ii]; j < num_column_entries[ii + 1];
					 j++)
				{
					jj = entry_column[j];
					vec_r[jj] += entry[j] * vec_s[ii];
					if (symmetry & (ii != jj))
						vec_r[ii] += entry[j] * vec_s[jj];
				}
		}
		else if (storage_type == JDS)
		{
			for (k = 0; k < max_columns; k++)
				for (i = 0; i < num_column_entries[k]; i++)
				{
					ii = row_index_mapping_n2o[i];
					jj = entry_column[counter];
					vec_r[jj] += entry[counter] * vec_s[ii];
					if (symmetry & (ii != jj))
						vec_r[ii] += entry[counter] * vec_s[jj];
					counter++;
				}
		}
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

	if (storage_type == CRS)
	{
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
	}
	else if (storage_type == JDS)
	{
		const long kk = ii * DOF;
		long row_in_parse_table, counter;

		// Row is zero
		row_in_parse_table = row_index_mapping_o2n[id];
		counter = row_in_parse_table;
		for (k = 0; k < max_columns; k++)
		{
			if (row_in_parse_table < num_column_entries[k])
			{
				j0 = entry_column[counter];
				for (jj = 0; jj < DOF; jj++)
				{
					if (id == j0 && jj == ii)
					{
						vdiag = entry[(kk + jj) * size_entry_column + counter];
					}
					else
					{
						entry[(kk + jj) * size_entry_column + counter] = 0.;
					}
				}
				counter += num_column_entries[k];
			}
			else
				break;
		}
#ifdef colDEBUG
		//
		counter = 0;
		for (k = 0; k < max_columns; k++)
			for (i = 0; i < num_column_entries[k]; i++)
			{
				i0 = row_index_mapping_n2o[i];
				/*
				   if(i0 == id)
				   {
				   counter++;
				   continue;
				   }
				 */
				j0 = entry_column[counter];
				if (j0 == id)
					for (jj = 0; jj < DOF; jj++)
					{
						if (i0 == j0 && ii == jj) continue;
						j = (jj * DOF + ii) * size_entry_column + counter;
						b[jj * rows + i0] -= entry[j] * b_given;
						entry[j] = 0.;
						// Room for symmetry case
					}
				//
				counter++;
			}
#endif
	}
	b[idiag] = vdiag * b_given;
}

/*\!
 ********************************************************************
   M^{-1}*A

		  a_ij  i=j
   M = {
		  0     i!=j
   Programm:
   10/2007 WW
 ********************************************************************/
void CSparseMatrix::Precond_Jacobi(double* vec_s, double* vec_r)
{
	long i, idof;
	double diag = 0.;
	//
	if (DOF > 1)
	{
		// Although this piece of code can deal with the case
		// of DOF = 1, we also prepare a special piece of code for
		// the case of DOF = 1 just for efficiency
		for (i = 0; i < rows; i++)
			for (idof = 0; idof < DOF; idof++)
			{
				diag = entry[(idof * DOF + idof) * size_entry_column +
							 diag_entry[i]];
				if (fabs(diag) < DBL_MIN)
					//        if(fabs(diag)<DBL_EPSILON)
					diag = 1.0;
				//  std::cout<<"Diagonal entry is zero. Abort simulation!!  "
				//  <<"\n";
				vec_r[idof * rows + i] = vec_s[idof * rows + i] / diag;
			}
		//
	}
	else  // DOF = 1

		for (i = 0; i < rows; i++)
		{
			diag = entry[diag_entry[i]];
			// if(fabs(diag)<DBL_EPSILON)
			if (fabs(diag) < DBL_MIN) diag = 1.0;
			//   std::cout<<"Diagonal entry is zero. Abort simulation!!  "
			//   <<"\n";
			//
			vec_r[i] = vec_s[i] / diag;
		}
}


#ifdef LIS
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
#endif  // LIS
#endif  // NEW_EQS
///////////////////////////////////////////////////////////

}  // Namespace
