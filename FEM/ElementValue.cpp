/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementValue.h"

#include "MSHEnums.h"
#include "msh_elem.h"

#include "rf_pcs.h"

using namespace MeshLib;
using namespace Math_Group;

std::vector<FiniteElement::ElementValue*> ele_gp_value;

namespace FiniteElement
{

// WW 08/2007
ElementValue::ElementValue(CRFProcess* m_pcs, CElem* ele) : pcs(m_pcs)
{
	int NGPoints = 0, NGP = 0;
	int ele_dim;

	MshElemType::type ele_type = ele->GetElementType();
	ele_dim = ele->GetDimension();

	NGP = GetNumericsGaussPoints(ele_type);
	if (ele_type == MshElemType::LINE)
		// OKWW
		NGPoints = m_pcs->m_num->ele_gauss_points;
	else if (ele_type == MshElemType::TRIANGLE)
		NGPoints = 3;
	else if (ele_type == MshElemType::TETRAHEDRON)
		NGPoints = 5;  // 15;
	else
		NGPoints = (int)MathLib::fastpow(NGP, ele_dim);

	// WW Velocity.resize(m_pcs->m_msh->GetCoordinateFlag()/10, NGPoints);
	Velocity.resize(3, NGPoints);
	Velocity = 0.0;
	Velocity0.resize(3, NGPoints);
	Velocity0 = 0.0;
	// 15.3.2007 Multi-phase flow WW
	if (pcs->type == 1212 || pcs->type == 1313 || m_pcs->type == 42)
	{
		Velocity_g.resize(3, NGPoints);
		Velocity_g = 0.0;
	}
}
// WW 08/2007
void ElementValue::getIPvalue_vec(const int IP, double* vec)
{
	// SB, BG
	for (int i = 0; i < int(Velocity.Rows()); i++)
		vec[i] = Velocity(i, IP);
}
// SB, BG 09/2010
void ElementValue::getIPvalue_vec_phase(const int IP, int phase, double* vec)
{
	if (phase == 0)
		for (int i = 0; (size_t)i < Velocity.Rows(); i++)
			vec[i] = Velocity(i, IP);
	else if (phase == 10)
		for (int i = 0; (size_t)i < Velocity_g.Rows(); i++)
			vec[i] = Velocity_g(i, IP);
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2006 YD Implementation
   last modification:
**************************************************************************/
void ElementValue::GetEleVelocity(double* vec)
{
	for (int i = 0; (size_t)i < Velocity.Rows(); i++)
	{
		vec[i] = 0.0;
		for (int j = 0; (size_t)j < Velocity.Cols(); j++)
			vec[i] += Velocity(i, j);
		vec[i] /= Velocity.Cols();
	}
}
// WW
ElementValue::~ElementValue()
{
	Velocity.resize(0, 0);
	Velocity0.resize(0, 0);
	Velocity_g.resize(0, 0);
}


}  // end namespace

