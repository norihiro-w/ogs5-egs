/*
 * CubicSpline.h
 *
 *  Created on: Jul 27, 2010
 *      Author: TF (moved class CubicSpline from geo_mathlib.{h.cpp})
 */

#ifndef CUBICSPLINE_H_
#define CUBICSPLINE_H_

#include <cstddef>
#include <vector>

namespace MathLib
{
class CubicSpline
{
public:
	CubicSpline(const std::vector<double>& s, const std::vector<double>& val);
	~CubicSpline();
	double interpolation(double x) const;

private:
	const size_t n;
	double* bb;
	double* cc;
	double* dd;
	std::vector<double> xx;
	std::vector<double> yy;
	void computeCoefficents();
};
}  // end namespace MathLib

#endif /* CUBICSPLINE_H_ */
