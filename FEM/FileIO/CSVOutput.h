/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef CSVOUTPUT_H
#define CSVOUTPUT_H

#include <string>

class COutput;

class CSVOutput
{
public:
	static void writeDomain(COutput* output,
	                  int timestep_number,
	                  double current_time,
	                  std::string const& baseFilename);

	static void writePoyline(COutput* output,
								   int timestep_number,
								   double /*current_time*/,
								   std::string const& baseFilename);

};

#endif  // CSVOUTPUT_H
