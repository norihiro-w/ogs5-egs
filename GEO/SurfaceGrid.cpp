/*
 * SurfaceGrid.cpp
 *
 *  Created on: Feb 9, 2012
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "SurfaceGrid.h"
#include "Surface.h"
#include "Triangle.h"

namespace GEOLIB
{
SurfaceGrid::SurfaceGrid(Surface const* const sfc)
    : AABB(sfc->getAABB()), _triangles_in_grid_box(NULL)
{
	double delta[3] = {0.0, 0.0, 0.0};
	for (size_t k(0); k < 3; k++)
	{
		// make the bounding box a little bit bigger,
		// such that the node with maximal coordinates fits into the grid
		_max_pnt[k] += std::abs(_max_pnt[k]) * 1e-6;
		if (fabs(_max_pnt[k]) < std::numeric_limits<double>::epsilon())
		{
			_max_pnt[k] = (_max_pnt[k] - _min_pnt[k]) * (1.0 + 1e-6);
		}
		delta[k] = _max_pnt[k] - _min_pnt[k];
	}

	if (delta[0] < std::numeric_limits<double>::epsilon())
	{
		const double max_delta(std::max(delta[1], delta[2]));
		_min_pnt[0] -= max_delta * 0.5e-3;
		_max_pnt[0] += max_delta * 0.5e-3;
		delta[0] = _max_pnt[0] - _min_pnt[0];
	}

	if (delta[1] < std::numeric_limits<double>::epsilon())
	{
		const double max_delta(std::max(delta[0], delta[2]));
		_min_pnt[1] -= max_delta * 0.5e-3;
		_max_pnt[1] += max_delta * 0.5e-3;
		delta[1] = _max_pnt[1] - _min_pnt[1];
	}

	if (delta[2] < std::numeric_limits<double>::epsilon())
	{
		const double max_delta(std::max(delta[0], delta[1]));
		_min_pnt[2] -= max_delta * 0.5e-3;
		_max_pnt[2] += max_delta * 0.5e-3;
		delta[2] = _max_pnt[2] - _min_pnt[2];
	}

	const size_t n_triangles(sfc->getNTriangles());
	const size_t n_tris_per_box(5);
	// *** condition: n_triangles / (_n_steps[0] * _n_steps[1] * _n_steps[2]) <
	// n_tris_per_box
	// *** with _n_steps[1] = _n_steps[0] * delta[1]/delta[0], _n_steps[2] =
	// _n_steps[0] * delta[2]/delta[0]
	if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() ||
	    fabs(delta[1]) < std::numeric_limits<double>::epsilon() ||
	    fabs(delta[2]) < std::numeric_limits<double>::epsilon())
	{
		// 1d case y = z = 0
		if (fabs(delta[1]) < std::numeric_limits<double>::epsilon() &&
		    fabs(delta[2]) < std::numeric_limits<double>::epsilon())
		{
			_n_steps[0] =
			    static_cast<size_t>(ceil(n_triangles / (double)n_tris_per_box));
			_n_steps[1] = 1;
			_n_steps[2] = 1;
		}
		else
		{
			// 1d case x = z = 0
			if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() &&
			    fabs(delta[2]) < std::numeric_limits<double>::epsilon())
			{
				_n_steps[0] = 1;
				_n_steps[1] = static_cast<size_t>(
				    ceil(n_triangles / (double)n_tris_per_box));
				_n_steps[2] = 1;
			}
			else
			{
				// 1d case x = y = 0
				if (fabs(delta[0]) < std::numeric_limits<double>::epsilon() &&
				    fabs(delta[1]) < std::numeric_limits<double>::epsilon())
				{
					_n_steps[0] = 1;
					_n_steps[1] = 1;
					_n_steps[2] = static_cast<size_t>(
					    ceil(n_triangles / (double)n_tris_per_box));
				}
				else
				{
					// 2d cases
					// y = 0
					if (fabs(delta[1]) < std::numeric_limits<double>::epsilon())
					{
						_n_steps[0] = static_cast<size_t>(
						    ceil(sqrt(n_triangles * delta[0] /
						              (n_tris_per_box * delta[2]))));
						_n_steps[1] = 1;
						_n_steps[2] = static_cast<size_t>(
						    ceil(_n_steps[0] * delta[2] / delta[0]));
					}
					else
					{
						// z = 0
						if (fabs(delta[2]) <
						    std::numeric_limits<double>::epsilon())
						{
							_n_steps[0] = static_cast<size_t>(
							    ceil(sqrt(n_triangles * delta[0] /
							              (n_tris_per_box * delta[1]))));
							_n_steps[1] = static_cast<size_t>(
							    ceil(_n_steps[0] * delta[1] / delta[0]));
							_n_steps[2] = 1;
						}
						else
						{
							// x = 0
							_n_steps[0] = 1;
							_n_steps[1] = static_cast<size_t>(
							    ceil(sqrt((double)n_triangles / n_tris_per_box *
							              delta[1] / delta[2])));
							_n_steps[2] = static_cast<size_t>(
							    ceil(n_triangles * delta[2] /
							         (n_tris_per_box * delta[1])));
						}
					}
				}
			}
		}
	}
	else
	{
		// 3d case
		_n_steps[0] = static_cast<size_t>(
		    ceil(pow(n_triangles * delta[0] * delta[0] /
		                 (n_tris_per_box * delta[1] * delta[2]),
		             1. / 3.)));
		_n_steps[1] = static_cast<size_t>(
		    ceil(_n_steps[0] * std::min(delta[1] / delta[0], 100.0)));
		_n_steps[2] = static_cast<size_t>(
		    ceil(_n_steps[0] * std::min(delta[2] / delta[0], 100.0)));
	}

	const size_t n_plane(_n_steps[0] * _n_steps[1]);
	_triangles_in_grid_box =
	    new std::vector<Triangle const*>[n_plane * _n_steps[2]];

	// some frequently used expressions to fill the grid vectors
	for (size_t k(0); k < 3; k++)
	{
		_step_sizes[k] = delta[k] / _n_steps[k];
		//		_inverse_step_sizes[k] = 1.0 / _step_sizes[k];
		if (_step_sizes[k] == 0)
			_inverse_step_sizes[k] = 1.0;
		else
			_inverse_step_sizes[k] = 1.0 / _step_sizes[k];
	}

	// fill the grid vectors
	size_t i_min, i_max, j_min, j_max, k_min, k_max;
	for (size_t l(0); l < n_triangles; l++)
	{
		Triangle const* const tri((*sfc)[l]);
		Point const& pnt(*(tri->getPoint(0)));
		i_min = i_max = static_cast<size_t>((pnt[0] - _min_pnt[0]) *
		                                    _inverse_step_sizes[0]);
		j_min = j_max = static_cast<size_t>((pnt[1] - _min_pnt[1]) *
		                                    _inverse_step_sizes[1]);
		k_min = k_max = static_cast<size_t>((pnt[2] - _min_pnt[2]) *
		                                    _inverse_step_sizes[2]);

		if (i_min >= _n_steps[0] || j_min >= _n_steps[1] ||
		    k_min >= _n_steps[2])
		{
			std::cout << "error computing indices "
			          << "\n";
		}

		for (size_t m(1); m < 3; m++)
		{
			Point const& pnt2(*(tri->getPoint(m)));
			const size_t i(static_cast<size_t>((pnt2[0] - _min_pnt[0]) *
			                                   _inverse_step_sizes[0]));
			const size_t j(static_cast<size_t>((pnt2[1] - _min_pnt[1]) *
			                                   _inverse_step_sizes[1]));
			const size_t k(static_cast<size_t>((pnt2[2] - _min_pnt[2]) *
			                                   _inverse_step_sizes[2]));

			if (i >= _n_steps[0] || j >= _n_steps[1] || k >= _n_steps[2])
			{
				std::cout << "error computing indices "
				          << "\n";
			}

			if (i < i_min) i_min = i;
			if (i_max < i) i_max = i;
			if (j < j_min) j_min = j;
			if (j_max < j) j_max = j;
			if (k < k_min) k_min = k;
			if (k_max < k) k_max = k;
		}

		for (size_t i(i_min); i <= i_max; i++)
		{
			for (size_t j(j_min); j <= j_max; j++)
			{
				for (size_t k(k_min); k <= k_max; k++)
				{
					_triangles_in_grid_box[i + j * _n_steps[0] + k * n_plane]
					    .push_back(tri);
				}
			}
		}
	}
}

bool SurfaceGrid::isPntInSurface(const double* pnt, double eps) const
{
	const size_t i(
	    static_cast<size_t>((pnt[0] - _min_pnt[0]) * _inverse_step_sizes[0]));
	const size_t j(
	    static_cast<size_t>((pnt[1] - _min_pnt[1]) * _inverse_step_sizes[1]));
	const size_t k(
	    static_cast<size_t>((pnt[2] - _min_pnt[2]) * _inverse_step_sizes[2]));

	if (i >= _n_steps[0] || j >= _n_steps[1] || k >= _n_steps[2])
	{
		return false;
	}

	std::vector<Triangle const*> const& triangles(
	    _triangles_in_grid_box[i + j * _n_steps[0] +
	                           k * _n_steps[0] * _n_steps[1]]);
	bool nfound(true);
	const size_t n_triangles(triangles.size());
	for (size_t kk(0); kk < n_triangles && nfound; kk++)
	{
		if (triangles[kk]->containsPoint(pnt, eps))
		{
			nfound = false;
		}
	}

	return !nfound;
}

SurfaceGrid::~SurfaceGrid()
{
	delete[] _triangles_in_grid_box;
}

}  // end namespace GEOLIB
