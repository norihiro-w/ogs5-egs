/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef element_value_INC
#define element_value_INC

#include <vector>

#include "matrix_class.h"

namespace MeshLib
{
class CElem;
}

class CRFProcess;
class COutput;

namespace FiniteElement
{
class CFiniteElementStd;

class ElementValue
{
public:
	ElementValue(CRFProcess* m_pcs, MeshLib::CElem* ele);
	~ElementValue();
	void getIPvalue_vec(const int IP, double* vec);
	void getIPvalue_vec_phase(const int IP, int phase, double* vec);
	void GetEleVelocity(double* vec);
	Math_Group::Matrix Velocity;
	Math_Group::Matrix Velocity0;

private:
	// Friend class
	friend class ::CRFProcess;
	friend class FiniteElement::CFiniteElementStd;
	friend class ::COutput;

	CRFProcess* pcs;
	Math_Group::Matrix Velocity_g;
};

}  // end namespace

extern std::vector<FiniteElement::ElementValue*> ele_gp_value;

#endif
