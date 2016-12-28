/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object:GEO Mathlib
   Task:
   Programing:
   08/2005 CC Implementation
**************************************************************************/
#ifndef geomathlib_INC
#define geomathlib_INC

#include <cmath>
#include <limits>

#include "makros.h"
#include "Point.h"

extern double EuklVek3dDist(double* x, double* y);
extern double EuklVek3dDistCoor(
    double x1, double y1, double z1, double x2, double y2, double z2);
extern double Vek3dDistCoor(
    double x1, double y1, double z1, double x2, double y2, double z2, int norm);
extern int M3KreuzProdukt(double* vec1, double* vec2, double* vec);
extern double MBtrgVec(double* vec, long n);
extern double MSkalarprodukt(double* vec1, double* vec2, long g);
extern double M3Determinante(double* m);
extern double M4Determinante(double* m);
extern double CalcTetraederVolume(double* x, double* y, double* z);  // CC
extern double CalcPyramidVolume(double* x, double* y, double* z);    // CC
extern double CalcPrismVolume(double* x, double* y, double* z);      // CC
extern double MCalcProjectionOfPointOnPlane(
    double* pt, double* e1, double* e2, double* e3, double* proj);
extern double MCalcProjectionOfPointOnPlane(
    double* pt, double* e1, double* e2, double* e3, double* proj);
extern double MCalcDistancePointToPoint(double* pt1, double* pt2);
extern long* TOLSortNodes1(long*, double*, int);
extern int MPhi2D(double* vf, double r, double s);

/**
 * Checks if two points are within a given distance of each other
 * @param p0 The first point
 * @param p1 the second point
 * @param squaredDistance The square of the distance within which the two points
 * should be
 * @return true if p1 and p2 are within the given distance of each other, false
 * otherwise
 */
bool checkDistance(GEOLIB::Point const& p0, GEOLIB::Point const& p1,
				   double squaredDistance);

/** squared euklid norm of the vector p0 */
double sqrNrm2(const GEOLIB::Point* const p0);

/** squared dist between GEOLIB::Points p0 and p1 */
double sqrDist(const GEOLIB::Point* p0, const GEOLIB::Point* p1);

extern double MCalcDistancePointToLine(double* pt, double* l1, double* l2);
extern double MCalcProjectionOfPointOnLine(double* pt1, double* pt2,
										   double* pt3, double* pt4);
extern double MCalcDistancePointToPlane(double const* const pt, double* e1,
										double* e2, double* e3);

#endif
