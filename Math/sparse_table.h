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

	bool symmetry = false;
	// Topology mapping from data array to matrix
	long* entry_column = nullptr;
	long* num_column_entries = nullptr;     // number of entries of each columns in sparse
								  // table
	long* row_index_mapping_n2o = nullptr;  // Row index of sparse table to row index of
								  // matrix
	long* row_index_mapping_o2n = nullptr;  // Inverse of last
	long* diag_entry = nullptr;             // Global index to the index of  entry_column
	long size_entry_column = 0;
	long max_columns = 0;
	long rows = 0;
	StorageType storage_type = JDS;
};

#endif

}

#endif
