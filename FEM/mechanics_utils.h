/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MECHANICS_UTILS_INC
#define MECHANICS_UTILS_INC

#include "matrix_class.h"

inline double MeanStress(Math_Group::Matrix const& stress, const int gp)
{
	return 1./3.*(stress(0, gp) + stress(1, gp) + stress(2, gp));
}


#endif // MECHANICS_UTILS_INC
