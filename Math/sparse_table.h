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

struct SparseTable
{
	~SparseTable();

	bool symmetry = false;
	// Topology mapping from data array to matrix
	long* entry_column = nullptr;
	long* num_column_entries = nullptr;     // number of entries of each columns in sparse
	                              // table
	long* diag_entry = nullptr;             // Global index to the index of  entry_column
	long size_entry_column = 0;
	long rows = 0;

	void Write(std::ostream& os = std::cout);
};

} // Math_Group

#endif // sparse_table_INC
