/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef sparse_table_INC
#define sparse_table_INC

#include <iostream>

namespace Math_Group
{

#ifdef NEW_EQS

/// Sparse matrix storage type
enum StorageType
{
	CRS,
	JDS
};

struct SparseTable
{
	~SparseTable();
	void Write(std::ostream& os = std::cout);

	bool symmetry;
	// Topology mapping from data array to matrix
	long* entry_column;
	long* num_column_entries;     // number of entries of each columns in sparse
								  // table
	long* row_index_mapping_n2o;  // Row index of sparse table to row index of
								  // matrix
	long* row_index_mapping_o2n;  // Inverse of last
	long* diag_entry;             // Global index to the index of  entry_column
	long size_entry_column;
	long max_columns;
	long rows;
	StorageType storage_type;  // 04.2011. WW
};

#endif

}

#endif
