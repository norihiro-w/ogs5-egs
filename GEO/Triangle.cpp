/*
 * \file Triangle.cpp
 *
 *  Created on: Jun 6, 2011
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "Triangle.h"

// MathLib
#include "MathTools.h"
#include "Matrix.h"
#include "Vector3.h"
#include "GaussAlgorithm.h"

#include "geo_mathlib.h"

namespace GEOLIB
{
Triangle::Triangle(std::vector<Point*> const& pnt_vec)
    : _pnts(pnt_vec), _initialized(false), _longest_edge(0.0)
{
	_pnt_ids[0] = std::numeric_limits<size_t>::max();
	_pnt_ids[1] = std::numeric_limits<size_t>::max();
	_pnt_ids[2] = std::numeric_limits<size_t>::max();
}

Triangle::Triangle(std::vector<Point*> const& pnt_vec,
                   size_t pnt_a,
                   size_t pnt_b,
                   size_t pnt_c)
    : _pnts(pnt_vec), _initialized(true), _longest_edge(0.0)
{
	_pnt_ids[0] = pnt_a;
	_pnt_ids[1] = pnt_b;
	_pnt_ids[2] = pnt_c;
	_longest_edge = sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]]);
	double tmp(sqrDist(_pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]]));
	if (tmp > _longest_edge) _longest_edge = tmp;
	tmp = sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[2]]);
	if (tmp > _longest_edge) _longest_edge = tmp;
	_longest_edge = sqrt(_longest_edge);
}

void Triangle::setTriangle(size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
	assert(pnt_a < _pnts.size() && pnt_b < _pnts.size() &&
	       pnt_c < _pnts.size());
	_pnt_ids[0] = pnt_a;
	_pnt_ids[1] = pnt_b;
	_pnt_ids[2] = pnt_c;

	_longest_edge = sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[1]]);
	double tmp(sqrDist(_pnts[_pnt_ids[1]], _pnts[_pnt_ids[2]]));
	if (tmp > _longest_edge) _longest_edge = tmp;
	tmp = sqrDist(_pnts[_pnt_ids[0]], _pnts[_pnt_ids[2]]);
	if (tmp > _longest_edge) _longest_edge = tmp;
	_longest_edge = sqrt(_longest_edge);
}

bool Triangle::containsPoint(const double* pnt, double eps) const
{
	GEOLIB::Point const& a(*(_pnts[_pnt_ids[0]]));
	if (sqrt(MathLib::sqrDist(a.getData(), pnt)) < eps) return true;
	GEOLIB::Point const& b(*(_pnts[_pnt_ids[1]]));
	if (sqrt(MathLib::sqrDist(b.getData(), pnt)) < eps) return true;
	GEOLIB::Point const& c(*(_pnts[_pnt_ids[2]]));
	if (sqrt(MathLib::sqrDist(c.getData(), pnt)) < eps) return true;

	// check if pnt is near the edge (a,b) of the triangle
	double lambda(0.0);
	double tmp_dist(0.0);
	double proj_dist(MathLib::calcProjPntToLineAndDists(
	    pnt, a.getData(), b.getData(), lambda, tmp_dist));
	if (proj_dist / sqrt(MathLib::sqrDist(a.getData(), b.getData())) < 5e-3)
		return true;

	// check if pnt is near the edge (b,c) of the triangle
	lambda = 0.0;
	tmp_dist = 0.0;
	proj_dist = MathLib::calcProjPntToLineAndDists(
	    pnt, b.getData(), c.getData(), lambda, tmp_dist);
	if (proj_dist / sqrt(MathLib::sqrDist(b.getData(), c.getData())) < 5e-3)
		return true;

	// check if pnt is near the edge (c,a) of the triangle
	lambda = 0.0;
	tmp_dist = 0.0;
	proj_dist = MathLib::calcProjPntToLineAndDists(
	    pnt, c.getData(), a.getData(), lambda, tmp_dist);
	if (proj_dist / sqrt(MathLib::sqrDist(c.getData(), a.getData())) < 5e-3)
		return true;

	const double delta(std::numeric_limits<double>::epsilon());
	const double upper(1 + delta);

	// check special case where points of triangle have the same x-coordinate
	if (fabs(b[0] - a[0]) <= std::numeric_limits<double>::epsilon() &&
	    fabs(c[0] - a[0]) <= std::numeric_limits<double>::epsilon())
	{
		// all points of triangle have same x-coordinate
		if (fabs(pnt[0] - a[0]) / _longest_edge <= 1e-3)
		{
			// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1
			// <= 1
			MathLib::Matrix<double> mat(2, 2);
			mat(0, 0) = b[1] - a[1];
			mat(0, 1) = c[1] - a[1];
			mat(1, 0) = b[2] - a[2];
			mat(1, 1) = c[2] - a[2];
			double y[2] = {pnt[1] - a[1], pnt[2] - a[2]};

			MathLib::GaussAlgorithm<double> gauss(mat);
			gauss.execute(y);

			if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] &&
			    y[1] <= upper && y[0] + y[1] <= upper)
				return true;
			else
				return false;
		}
		else
			return false;
	}

	// check special case where points of triangle have the same y-coordinate
	if (fabs(b[1] - a[1]) <= std::numeric_limits<double>::epsilon() &&
	    fabs(c[1] - a[1]) <= std::numeric_limits<double>::epsilon())
	{
		// all points of triangle have same y-coordinate
		if (fabs(pnt[1] - a[1]) / _longest_edge <= 1e-3)
		{
			// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1
			// <= 1
			MathLib::Matrix<double> mat(2, 2);
			mat(0, 0) = b[0] - a[0];
			mat(0, 1) = c[0] - a[0];
			mat(1, 0) = b[2] - a[2];
			mat(1, 1) = c[2] - a[2];
			double y[2] = {pnt[0] - a[0], pnt[2] - a[2]};

			MathLib::GaussAlgorithm<double> gauss(mat);
			gauss.execute(y);

			if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] &&
			    y[1] <= upper && y[0] + y[1] <= upper)
				return true;
			else
				return false;
		}
		else
			return false;
	}

	// check special case where points a and b of triangle have the same x- and
	// y-coordinate
	if (fabs(b[0] - a[0]) <= std::numeric_limits<double>::epsilon() &&
	    fabs(b[1] - a[1]) <= std::numeric_limits<double>::epsilon())
	{
		// criterion: p-c = u0 * (a-c) + u1 * (b-c); 0 <= u0, u1 <= 1, u0+u1 <=
		// 1
		MathLib::Matrix<double> mat(2, 2);
		mat(0, 0) = a[0] - c[0];
		mat(0, 1) = b[0] - c[0];
		mat(1, 0) = a[2] - c[2];
		mat(1, 1) = b[2] - c[2];
		double y[2] = {pnt[0] - c[0], pnt[2] - c[2]};

		MathLib::GaussAlgorithm<double> gauss(mat);
		gauss.execute(y);

		if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] &&
		    y[1] <= upper && y[0] + y[1] <= upper)
			return true;
		else
			return false;
	}

	// check special case where points b and c of triangle have the same x- and
	// y-coordinate
	if (fabs(c[0] - b[0]) <= std::numeric_limits<double>::epsilon() &&
	    fabs(c[1] - b[1]) <= std::numeric_limits<double>::epsilon())
	{
		// criterion: p-c = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <=
		// 1
		MathLib::Matrix<double> mat(2, 2);
		mat(0, 0) = b[0] - a[0];
		mat(0, 1) = c[0] - a[0];
		mat(1, 0) = b[2] - a[2];
		mat(1, 1) = c[2] - a[2];
		double y[2] = {pnt[0] - a[0], pnt[2] - a[2]};

		MathLib::GaussAlgorithm<double> gauss(mat);
		gauss.execute(y);

		if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] &&
		    y[1] <= upper && y[0] + y[1] <= upper)
			return true;
		else
			return false;
	}

	// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
	MathLib::Matrix<double> mat(2, 2);
	mat(0, 0) = b[0] - a[0];
	mat(0, 1) = c[0] - a[0];
	mat(1, 0) = b[1] - a[1];
	mat(1, 1) = c[1] - a[1];
	double y[2] = {pnt[0] - a[0], pnt[1] - a[1]};

	MathLib::GaussAlgorithm<double> gauss(mat);
	gauss.execute(y);

	// check if the solution fulfills the third equation
	if (fabs((b[2] - a[2]) * y[0] + (c[2] - a[2]) * y[1] - (pnt[2] - a[2])) <
	    1e-3 * _longest_edge)
	{
		if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] &&
		    y[1] <= upper && y[0] + y[1] <= upper)
			return true;
		return false;
	}
	else
		return false;
}

bool Triangle::containsPoint2D(const double* pnt) const
{
	GEOLIB::Point const& a(*(_pnts[_pnt_ids[0]]));
	GEOLIB::Point const& b(*(_pnts[_pnt_ids[1]]));
	GEOLIB::Point const& c(*(_pnts[_pnt_ids[2]]));

	// criterion: p-a = u0 * (b-a) + u1 * (c-a); 0 <= u0, u1 <= 1, u0+u1 <= 1
	MathLib::Matrix<double> mat(2, 2);
	mat(0, 0) = b[0] - a[0];
	mat(0, 1) = c[0] - a[0];
	mat(1, 0) = b[1] - a[1];
	mat(1, 1) = c[1] - a[1];
	double y[2] = {pnt[0] - a[0], pnt[1] - a[1]};

	MathLib::GaussAlgorithm<double> gauss(mat);
	gauss.execute(y);

	const double delta(std::numeric_limits<double>::epsilon());
	const double upper(1 + delta);

	// check if u0 and u1 fulfills the condition (with some delta)
	if (-delta <= y[0] && y[0] <= upper && -delta <= y[1] && y[1] <= upper &&
	    y[0] + y[1] <= upper)
		return true;
	return false;
}

void getPlaneCoefficients(Triangle const& tri, double c[3])
{
	GEOLIB::Point const& p0(*(tri.getPoint(0)));
	GEOLIB::Point const& p1(*(tri.getPoint(1)));
	GEOLIB::Point const& p2(*(tri.getPoint(2)));
	MathLib::Matrix<double> mat(3, 3);
	mat(0, 0) = p0[0];
	mat(0, 1) = p0[1];
	mat(0, 2) = 1.0;
	mat(1, 0) = p1[0];
	mat(1, 1) = p1[1];
	mat(1, 2) = 1.0;
	mat(2, 0) = p2[0];
	mat(2, 1) = p2[1];
	mat(2, 2) = 1.0;
	c[0] = p0[2];
	c[1] = p1[2];
	c[2] = p2[2];

	MathLib::GaussAlgorithm<double> gauss(mat);
	gauss.execute(c);
}
}  // end namespace GEOLIB
