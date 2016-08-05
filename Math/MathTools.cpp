/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
 * MathTools.cpp
 *
 *  Created on: Jan 13, 2010
 *      Author: TF
 */

#include "MathTools.h"

namespace MathLib
{
void crossProd(const double u[3], const double v[3], double r[3])
{
	r[0] = u[1] * v[2] - u[2] * v[1];
	r[1] = u[2] * v[0] - u[0] * v[2];
	r[2] = u[0] * v[1] - u[1] * v[0];
}

double calcProjPntToLineAndDists(const double p[3], const double a[3],
                                 const double b[3], double& lambda, double& d0)
{
	// g (lambda) = a + lambda v, v = b-a
	double v[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
	// orthogonal projection: (g(lambda)-p) * v = 0 => in order to compute
	// lambda we define a help vector u
	double u[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
	lambda = scpr(u, v, 3) / scpr(v, v, 3);

	// compute projected point
	double proj_pnt[3];
	for (size_t k(0); k < 3; k++)
		proj_pnt[k] = a[k] + lambda * v[k];

	d0 = sqrt(sqrDist(proj_pnt, a));

	return sqrt(sqrDist(p, proj_pnt));
}

float normalize(float min, float max, float val)
{
	return (val - min) / static_cast<float>(max - min);
}

double getAngle(const double p0[3], const double p1[3], const double p2[3])
{
	const double v0[3] = {p0[0] - p1[0], p0[1] - p1[1], p0[2] - p1[2]};
	const double v1[3] = {p2[0] - p1[0], p2[1] - p1[1], p2[2] - p1[2]};

	// apply Cauchy Schwarz inequality
	return acos(scpr(v0, v1, 3) /
	            (sqrt(scpr(v0, v0, 3)) * sqrt(scpr(v1, v1, 3))));
}
}  // namespace

/***************************************************************************
   GeoSys - Funktion:
           CElement::ComputeDetTri(const *double x1, const *double x2,
                                 const *double x3)
   Aufgabe:
         Compute the vulume of a triangle
   Formalparameter:
           E:
             const *double x1    : Vertex 1
             const *double x2    : Vertex 2
             const *double x3    : Vertex 3
   Programming:
   09/2004     WW        Erste Version
 **************************************************************************/
double ComputeDetTri(const double* x1, const double* x2, const double* x3)
{
	static double u[3], v[3], z[3];

	u[0] = x3[0] - x1[0];
	u[1] = x3[1] - x1[1];
	u[2] = x3[2] - x1[2];

	v[0] = x2[0] - x1[0];
	v[1] = x2[1] - x1[1];
	v[2] = x2[2] - x1[2];

	z[0] = u[1] * v[2] - u[2] * v[1];
	z[1] = u[2] * v[0] - u[0] * v[2];
	z[2] = u[0] * v[1] - u[1] * v[0];

	return 0.5 * sqrt(z[0] * z[0] + z[1] * z[1] + z[2] * z[2]);
}

/***************************************************************************
   GeoSys - Funktion:
           CElem::ComputeDetTet(const *double x1, const *double x2,
                                 const *double x3, const *double x4)
   Aufgabe:
         Compute the vulume of a tedrahedra
   Formalparameter:
           E:
             const *double x1    : Vertex 1
             const *double x2    : Vertex 2
             const *double x3    : Vertex 3
   const *double x4    : Vertex 4
   Programming:
   09/2004     WW        Erste Version
 **************************************************************************/
double ComputeDetTex(const double* x1, const double* x2, const double* x3,
                     const double* x4)
{
	return fabs((x1[0] - x4[0]) * ((x2[1] - x4[1]) * (x3[2] - x4[2]) -
	                               (x2[2] - x4[2]) * (x3[1] - x4[1])) -
	            (x1[1] - x4[1]) * ((x2[0] - x4[0]) * (x3[2] - x4[2]) -
	                               (x2[2] - x4[2]) * (x3[0] - x4[0])) +
	            (x1[2] - x4[2]) * ((x2[0] - x4[0]) * (x3[1] - x4[1]) -
	                               (x2[1] - x4[1]) * (x3[0] - x4[0]))) /
	       6.0;
}

/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   09/2005 WW Implementation
   12/2010 TF
**************************************************************************/
double NormalizeVector(double* x, size_t n)
{
	double nrm(x[0] * x[0]);
	for (size_t i = 1; i < n; i++)
		nrm += x[i] * x[i];
	double sqrt_nrm(sqrt(nrm));
	for (size_t i = 0; i < n; i++)
		x[i] /= sqrt_nrm;
	return sqrt_nrm;
}

/**************************************************************************
   MSHLib-Method:
   Task:         dim == 3
   Programing:
   09/2005 WW Implementation
**************************************************************************/
void CrossProduction(const double* x, const double* y, double* z)
{
	z[0] = x[1] * y[2] - x[2] * y[1];
	z[1] = x[2] * y[0] - x[0] * y[2];
	z[2] = x[0] * y[1] - x[1] * y[0];
}
