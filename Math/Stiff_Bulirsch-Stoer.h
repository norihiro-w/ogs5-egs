/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef STIFF_BULIRSCH_STOER
#define STIFF_BULIRSCH_STOER

extern void derivs(double x, double y[], double dydx[], int n, long node);
extern void jacobn(double x, double y[], double dfdx[], double** dfdy, int n,
                   long node);

#endif
