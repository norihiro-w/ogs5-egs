/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "sparse_table.h"

#include <iomanip>

namespace Math_Group
{

/*\!
 ********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
void SparseTable::Write(std::ostream& os)
{
	long i, k, counter = 0;

	os.width(10);
	os << "Symmetry: " << symmetry << "\n";
	os << "\n*** Row index  "
	   << "\n";

	if (storage_type == CRS)
	{
		os << "\n*** Sparse entry  "
		   << "\n";
		for (i = 0; i < rows; i++)
		{
			for (k = num_column_entries[i]; k < num_column_entries[i + 1]; k++)
				os << entry_column[k] + 1 << " ";
			os << "\n";
		}
	}
	else if (storage_type == JDS)
	{
		for (i = 0; i < rows; i++)
			os << row_index_mapping_n2o[i] + 1 << "\n";
		//
		os << "\n*** Sparse entry  "
		   << "\n";
		for (k = 0; k < max_columns; k++)
		{
			os << "--Column: " << k + 1 << "\n";
			for (i = 0; i < num_column_entries[k]; i++)
			{
				os << entry_column[counter] + 1 << "\n";
				counter++;
			}
			os << "\n";
		}
	}
}

/*\!
 ********************************************************************
   Create sparse matrix table
   08/2007 WW
   10/2007 WW
   5/2011 WW  CRS storage
 ********************************************************************/
SparseTable::~SparseTable()
{
	if (entry_column) delete[] entry_column;
	if (num_column_entries) delete[] num_column_entries;
	if (row_index_mapping_n2o) delete[] row_index_mapping_n2o;
	if (row_index_mapping_o2n) delete[] row_index_mapping_o2n;
	if (diag_entry) delete[] diag_entry;
	entry_column = NULL;
	num_column_entries = NULL;
	row_index_mapping_n2o = NULL;
	row_index_mapping_o2n = NULL;
	diag_entry = NULL;
}


}  // Namespace

