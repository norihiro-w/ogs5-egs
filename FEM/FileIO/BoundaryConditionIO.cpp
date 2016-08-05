/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoundaryConditionIO.h"

#include "readNonBlankLineFromInputStream.h"

#include "GeoIO.h"
#include "ProcessIO.h"

namespace FileIO
{

void BoundaryConditionIO::write(std::ostream& out, CBoundaryCondition const& bc)
{
	// keyword
	out << "#BOUNDARY_CONDITION"
	    << "\n";

	// process and primary variable
	out << "\t$PCS_TYPE"
	    << "\n";
	out << "\t\t" << convertProcessTypeToString(bc.getProcessType()) << "\n";

	out << "\t$PRIMARY_VARIABLE"
	    << "\n";
	out << "\t\t"
	    << convertPrimaryVariableToString(bc.getProcessPrimaryVariable())
	    << "\n";

	// geometry
	out << "\t$GEO_TYPE"
	    << "\n";
	out << "\t" << bc.getGeoTypeAsString() << " " << bc.geo_name << "\n";

	// distribution type
	out << "\t$DIS_TYPE"
	    << "\n";
	out << "\t\t" << convertDisTypeToString(bc.getProcessDistributionType());
	if (bc.getProcessDistributionType() == FiniteElement::CONSTANT)
		out << "\t\t" << bc.geo_node_value << "\n";
	else if (bc.getProcessDistributionType() == FiniteElement::LINEAR)
	{
		out << "\t\t" << bc._PointsHaveDistribedBC.size() << "\n";
		for (size_t i = 0; i < bc._PointsHaveDistribedBC.size(); i++)
			out << "\t\t" << bc._PointsHaveDistribedBC[i] << "  "
			    << bc._DistribedBC[i] << "\n";
	}
	else
		out << "\n";

	// function name
	if (!bc.fct_name.empty())
	{
		out << "\t$FCT_TYPE"
		    << "\n";
		out << "\t\t" << bc.fct_name << "\n";
	}
}
}  // end namespace FileIO
