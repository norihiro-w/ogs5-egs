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
	os.width(10);
	os << "Symmetry: " << symmetry << "\n";
	os << "\n*** Row index  "
	   << "\n";

	os << "\n*** Sparse entry  "
	   << "\n";
	for (long i = 0; i < rows; i++)
	{
		for (long k = num_column_entries[i]; k < num_column_entries[i + 1]; k++)
			os << entry_column[k] + 1 << " ";
		os << "\n";
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
	entry_column = NULL;
	if (num_column_entries) delete[] num_column_entries;
	num_column_entries = NULL;
	if (diag_entry) delete[] diag_entry;
	diag_entry = NULL;
}


}  // Namespace

