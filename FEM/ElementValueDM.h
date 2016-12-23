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

namespace SolidProp
{
class CSolidProperties;
}
class CMediumProperties;

using namespace Math_Group;
using namespace MeshLib;

namespace FiniteElement
{

class ElementValue_DM
{
public:
	ElementValue_DM(CElem* ele, const int NGP, bool has_coupling_loop);
	~ElementValue_DM();
	void ResetStress(bool cpl_loop);
	void Write_BIN(std::fstream& os);
	void Read_BIN(std::fstream& is);
	void ReadElementStressASCI(std::fstream& is);

	Matrix* Stress0 = nullptr;  // Initial stress
	Matrix* Stress = nullptr;
	Matrix* Stress_last_ts = nullptr;
	Matrix* Stress_current_ts = nullptr;
	Matrix* dTotalStress = nullptr;

	Matrix* Strain = nullptr;
	Matrix* Strain_last_ts = nullptr;
	Matrix* pStrain = nullptr;
	Matrix* y_surface = nullptr;
	// Preconsolidation pressure
	Matrix* prep0 = nullptr;
	Matrix* e_i = nullptr;  // Void ratio
	// Variables of single yield surface model
	Matrix* xi = nullptr;    // Rotational hardening variables
	Matrix* MatP = nullptr;  // Material parameters
};

}  // end namespace

extern std::vector<FiniteElement::ElementValue_DM*> ele_value_dm;
#endif
