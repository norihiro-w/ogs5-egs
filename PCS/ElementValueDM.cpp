/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementValueDM.h"

#include <cfloat>

#include "Configure.h"

#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#if defined(NEW_EQS)
#include "equation_class.h"
#endif
#include "mathlib.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_pcs_dm.h"
#include "tools.h"


std::vector<FiniteElement::ElementValue_DM*> ele_value_dm;

namespace FiniteElement
{


/*----------------------------------------------------------------
   Class ElementValue_DM

   Allocate memory for element value
     Matrix *Mat:
 | Index  |  Paramete |
                          ----------------------
 |    0   |  alpha    |
 |    1   |  beta     |
 |    2   |  delta    |
 |    3   |  epsilon  |
 |    4   |  kappa    |
 |    5   |  gamma    |
 |    6   |  m        |
   ----------------------
   -----------------------------------------------------------------*/
ElementValue_DM::ElementValue_DM(CElem* ele, const int NGP, bool has_coupling_loop)
{
	MshElemType::type ele_type = ele->GetElementType();
	int ele_dim = ele->GetDimension();
	CSolidProperties* msp = msp_vector[ele->GetPatchIndex()];
	int Plastic = msp->Plastictity();

	int LengthBS = 4;         // Number of stress/strain components
	if (ele_dim == 2)
		LengthBS = 4;
	else if (ele_dim == 3)
		LengthBS = 6;

	int NGPoints = 0;
	if (ele_type == MshElemType::TRIANGLE)
		NGPoints = 3;
	else if (ele_type == MshElemType::TETRAHEDRON)
		NGPoints = 5;  // 15
	else if (ele_type == MshElemType::PRISM)
		NGPoints = 6;  // 9
	else
		NGPoints = MathLib::fastpow(NGP, ele_dim);

	Stress0 = new Matrix(LengthBS, NGPoints);
	dTotalStress = new Matrix(LengthBS, NGPoints);
	Stress_last_ts = new Matrix(LengthBS, NGPoints);
	if (has_coupling_loop)
		Stress_current_ts = new Matrix(LengthBS, NGPoints);
	Stress = Stress_last_ts;

	Strain = new Matrix(LengthBS, NGPoints);
	if (has_coupling_loop)
		Strain_last_ts = new Matrix(LengthBS, NGPoints);

	if (Plastic > 0)
	{
		pStrain = new Matrix(NGPoints);
		y_surface = new Matrix(NGPoints);
		*y_surface = 0.0;
		*pStrain = 0.0;
	}
	*Stress = 0.0;

	const int LengthMat = 7;  // Number of material parameter of SYS model.
	if (Plastic == 2)  // Rotational hardening model
	{
		xi = new Matrix(LengthBS - 1, NGPoints);
		MatP = new Matrix(LengthMat, NGPoints);
		*xi = 0.0;
		*MatP = 0.0;
	}
	if (Plastic == 3)  // Cam-Clay
	{
		prep0 = new Matrix(NGPoints);
		e_i = new Matrix(NGPoints);
		*prep0 = 0.0;
		*e_i = 0.0;
	}
	if (msp->CreepModel() == 1000)
	{
		xi = new Matrix(LengthBS);
		*xi = 0.0;
	}
	disp_j = 0.0;
	tract_j = 0.0;
	Localized = false;
}

ElementValue_DM::~ElementValue_DM()
{
	delete Stress0;
	delete Stress_last_ts;
	delete Stress_current_ts;
	delete dTotalStress;

	delete pStrain;
	delete y_surface;

	delete prep0;
	delete e_i;
	delete xi;
	delete MatP;

	delete NodesOnPath;
	delete orientation;

	delete Strain;
	delete Strain_last_ts;
}

void ElementValue_DM::Write_BIN(std::fstream& os)
{
	Stress0->Write_BIN(os);
	Stress_last_ts->Write_BIN(os);
	if (pStrain) pStrain->Write_BIN(os);
	if (y_surface) y_surface->Write_BIN(os);
	if (xi) xi->Write_BIN(os);
	if (MatP) MatP->Write_BIN(os);
	if (prep0) prep0->Write_BIN(os);
	if (e_i) e_i->Write_BIN(os);
	if (NodesOnPath) NodesOnPath->Write_BIN(os);
	if (orientation) os.write((char*)(orientation), sizeof(*orientation));
	os.write((char*)(&disp_j), sizeof(disp_j));
	os.write((char*)(&tract_j), sizeof(tract_j));
	os.write((char*)(&Localized), sizeof(Localized));
}

void ElementValue_DM::Read_BIN(std::fstream& is)
{
	Stress0->Read_BIN(is);
	Stress_last_ts->Read_BIN(is);
	if (pStrain) pStrain->Read_BIN(is);
	if (y_surface) y_surface->Read_BIN(is);
	if (xi) xi->Read_BIN(is);
	if (MatP) MatP->Read_BIN(is);
	if (prep0) prep0->Read_BIN(is);
	if (e_i) e_i->Read_BIN(is);
	if (NodesOnPath) NodesOnPath->Read_BIN(is);
	if (orientation) is.read((char*)(orientation), sizeof(*orientation));
	is.read((char*)(&disp_j), sizeof(disp_j));
	is.read((char*)(&tract_j), sizeof(tract_j));
	is.read((char*)(&Localized), sizeof(Localized));
}

void ElementValue_DM::ReadElementStressASCI(std::fstream& is)
{
	size_t i, j;
	size_t ns = Stress0->Rows();
	size_t nGS = Stress0->Cols();

	for (i = 0; i < ns; i++)
	{
		is >> (*Stress0)(i, 0);
		for (j = 1; j < nGS; j++)
			(*Stress0)(i, j) = (*Stress0)(i, 0);
	}

	*Stress_last_ts = *Stress0;
}

void ElementValue_DM::ResetStress(bool cpl_loop)
{
	if (Stress_current_ts == nullptr)
		return;

	if (cpl_loop)  // For coupling loop
	{
		(*Stress_current_ts) = (*Stress_last_ts);
		Stress = Stress_current_ts;
		(*Strain) = (*Strain_last_ts);
	}
	else  // Time loop
	{
		(*Stress_last_ts) = (*Stress_current_ts);
		Stress = Stress_last_ts;
		(*Strain_last_ts) = (*Strain);
	}
}


} // namespace
