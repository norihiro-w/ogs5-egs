/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ElementValueDM_INC
#define ElementValueDM_INC

#include "fem_ele.h"
#include "matrix_class.h"


namespace MeshLib
{
class CElem;
}
using namespace Math_Group;
using namespace MeshLib;

namespace FiniteElement
{

class ElementValue_DM
{
public:
	ElementValue_DM(MeshLib::CElem* ele, const int NGP, bool HM_Staggered);
	~ElementValue_DM();
	void ResetStress(bool cpl_loop);
	void Write_BIN(std::fstream& os);
	void Read_BIN(std::fstream& is);
	void ReadElementStressASCI(std::fstream& is);
	double MeanStress(const int gp)
	{
		return (*Stress)(0, gp) + (*Stress)(1, gp) + (*Stress)(2, gp);
	}

	Matrix* Stress0;  // Initial stress
	Matrix* Stress;
	Matrix* Stress_i;
	Matrix* Stress_j;
	Matrix* Strain;
	Matrix* pStrain;
	Matrix* y_surface;
	// Preconsolidation pressure
	Matrix* prep0;
	Matrix* e_i;  // Void ratio
	// Variables of single yield surface model
	Matrix* xi;    // Rotational hardening variables
	Matrix* MatP;  // Material parameters

	// Discontinuity
	double disp_j;
	double tract_j;
	bool Localized;
	Matrix* NodesOnPath;
	double* orientation;
};

}  // end namespace

extern std::vector<FiniteElement::ElementValue_DM*> ele_value_dm;
#endif
