/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef files0_INC
#define files0_INC

#include <string>

#include "GEOObjects.h"

char* ReadString(void);

extern int ReadData(char* dateiname,
					GEOLIB::GEOObjects& geo_obj,
					std::string& unique_name);

#endif
