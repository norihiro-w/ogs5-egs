/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_vec.h"

#include <cfloat>

#include "Configure.h"

#if defined(NEW_EQS)
#include "equation_class.h"
#endif
#include "mathlib.h"
#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#include "ElementValueDM.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_pcs_dm.h"
#include "tools.h"


#define COMP_MOL_MASS_AIR 28.96  // kg/kmol WW  28.96
#define GAS_CONSTANT 8314.41     // J/(kmol*K) WW

namespace FiniteElement
{

CFiniteElementVec::CFiniteElementVec(CRFProcessDeformation* dm_pcs,
                                     const int C_Sys_Flad,
                                     const int order)
    : CElement(C_Sys_Flad, order),
      pcs(dm_pcs),
      PressureC(NULL),
      PressureC_S(NULL),
      PressureC_S_dp(NULL),
      b_rhs(NULL)
{
	int i;

	S_Water = 1.0;
	Tem = 273.15 + 23.0;
	h_pcs = NULL;
	t_pcs = NULL;

	AuxNodal2 = NULL;
	//
	dim = pcs->m_msh->GetMaxElementDim();  // overwrite dim in CElement
	ns = 4;
	if (dim == 3) ns = 6;
	//  10.11.2010. WW
	AuxNodal0 = new double[8];
	AuxNodal = new double[8];
	AuxNodal_S0 = new double[8];
	AuxNodal_S = new double[8];
	AuxNodal1 = new double[60];

	Idx_Stress = new int[ns];
	Idx_Strain = new int[ns];
	strain_ne = new double[ns];
	stress_ne = new double[ns];
	stress0 = new double[ns];
	for (i = 0; i < 4; i++)
		NodeShift[i] = pcs->Shift[i];
	{
		// Indecex in nodal value table
		Idx_dm0[0] = pcs->GetNodeValueIndex("DISPLACEMENT_X1");
		Idx_dm0[1] = pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
		Idx_dm1[0] = Idx_dm0[0] + 1;
		Idx_dm1[1] = Idx_dm0[1] + 1;

		//     if(problem_dimension_dm==3)
		if (dim == 3)
		{
			Idx_dm0[2] = pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
			Idx_dm1[2] = Idx_dm0[2] + 1;
		}
	}

	// idx_pls =  pcs->GetNodeValueIndex("STRAIN_PLS");
	// Strain
	Idx_Strain[0] = pcs->GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = pcs->GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = pcs->GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = pcs->GetNodeValueIndex("STRAIN_XY");
	// Stress
	Idx_Stress[0] = pcs->GetNodeValueIndex("STRESS_XX");
	Idx_Stress[1] = pcs->GetNodeValueIndex("STRESS_YY");
	Idx_Stress[2] = pcs->GetNodeValueIndex("STRESS_ZZ");
	Idx_Stress[3] = pcs->GetNodeValueIndex("STRESS_XY");

	idx_S = -1;
	// Saturation;
	S_Water = 1.0;

	// For cache NW
	vec_B_matrix.resize(20);
	vec_B_matrix_T.resize(20);
	for (i = 0; i < (int)vec_B_matrix.size(); i++)
	{
		switch (dim)
		{
			case 2:
				vec_B_matrix[i] = new Matrix(4, 2);
				vec_B_matrix_T[i] = new Matrix(2, 4);
				break;
			case 3:
				vec_B_matrix[i] = new Matrix(6, 3);
				vec_B_matrix_T[i] = new Matrix(3, 6);
				break;
		}
	}

	//
	switch (dim)
	{
		case 2:
			ns = 4;
			B_matrix = new Matrix(4, 2);
			B_matrix_T = new Matrix(2, 4);
			dstress = new double[4];
			dstrain = new double[4];
			De = new Matrix(4, 4);
			ConsistDep = new Matrix(4, 4);
			AuxMatrix = new Matrix(2, 2);
			AuxMatrix2 = new Matrix(2, 4);  // NW
			Disp = new double[18];
			Temp = new double[9];
			T1 = new double[9];

			Sxx = new double[9];
			Syy = new double[9];
			Szz = new double[9];
			Sxy = new double[9];
			pstr = new double[9];

			Sxz = NULL;
			Syz = NULL;

			break;
		case 3:
			ns = 6;
			B_matrix = new Matrix(6, 3);
			B_matrix_T = new Matrix(3, 6);
			dstress = new double[6];
			dstrain = new double[6];
			De = new Matrix(6, 6);
			ConsistDep = new Matrix(6, 6);
			AuxMatrix = new Matrix(3, 3);
			AuxMatrix2 = new Matrix(3, 6);  // NW
			Disp = new double[60];
			Temp = new double[20];
			T1 = new double[20];

			Sxx = new double[20];
			Syy = new double[20];
			Szz = new double[20];
			Sxy = new double[20];
			Sxz = new double[20];
			Syz = new double[20];
			pstr = new double[20];

			// Indecex in nodal value table
			Idx_Strain[4] = pcs->GetNodeValueIndex("STRAIN_XZ");
			Idx_Strain[5] = pcs->GetNodeValueIndex("STRAIN_YZ");
			//
			Idx_Stress[4] = pcs->GetNodeValueIndex("STRESS_XZ");
			Idx_Stress[5] = pcs->GetNodeValueIndex("STRESS_YZ");
			break;
	}
	*B_matrix = 0.0;
	*B_matrix_T = 0.0;

	if (pcs->Memory_Type == 0)  // Do not store local matrices
	{
		Stiffness = new Matrix(60, 60);
		RHS = new Vector(60);
		if (H_Process) PressureC = new Matrix(60, 20);
	}
	else  // Local matrices stored and allocated the pcs.cpp;
	{
		Stiffness = NULL;
		RHS = NULL;
	}

	// Material properties
	smat = NULL;
	// Fluid coupling
	m_mfp = NULL;
	// Medium property
	m_mmp = NULL;  //

	// Coupling
	Flow_Type = -1;
	idx_P = -1;
	eleV_DM = NULL;
	idx_P0 = idx_pls = 0;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		//      if (pcs_vector[i]->pcs_type_name.find("FLOW") != string::npos) {
		// TF
		if (isFlowProcess(pcs_vector[i]->getProcessType()))
		{
			h_pcs = pcs_vector[i];
			// 25.04.2008, 04.09.2008  WW
			if (h_pcs->type == 1 || h_pcs->type == 41)
			{
				//				if (h_pcs->pcs_type_name.find("GROUND") !=
				// string::npos)
				// TF
				if (h_pcs->getProcessType() == GROUNDWATER_FLOW)
					Flow_Type = 10;
				else
					Flow_Type = 0;
				// 25.08.2005.  WW
			}
			else if (h_pcs->type == 14 || h_pcs->type == 22)
				Flow_Type = 1;
			else if (h_pcs->type == 1212 || h_pcs->type == 42)
			{
				Flow_Type = 2;  // 25.04.2008.  WW

				// 07.2011. WW
				PressureC_S = new Matrix(60, 20);
				if (pcs->m_num->nls_method == FiniteElement::NL_NEWTON &&
				    h_pcs->type == 42)  // Newton-raphson. WW
					PressureC_S_dp = new Matrix(60, 20);
			}
			// WW idx_P0 = pcs->GetNodeValueIndex("POROPRESSURE0");
			break;
			//		} else if (pcs_vector[i]->pcs_type_name.find("PS_GLOBAL") !=
			// string::npos) {
		}  // TF
		else if (pcs_vector[i]->getProcessType() == PS_GLOBAL)
		{
			h_pcs = pcs_vector[i];
			if (h_pcs->type == 1313) Flow_Type = 3;  // 05.05.2009.  PCH
			break;
		}
	}
	if (Flow_Type == 0)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
	}
	if (Flow_Type == 10)
		idx_P1 = h_pcs->GetNodeValueIndex("HEAD") + 1;
	else if (Flow_Type == 1)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		idx_P1_0 = h_pcs->GetNodeValueIndex("PRESSURE1");
		idx_S0 = h_pcs->GetNodeValueIndex("SATURATION1");
		idx_S = h_pcs->GetNodeValueIndex("SATURATION1") + 1;
	}
	else if (Flow_Type == 2)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		idx_P2 = h_pcs->GetNodeValueIndex("PRESSURE2") + 1;
		idx_S0 = h_pcs->GetNodeValueIndex("SATURATION1");
		idx_S = h_pcs->GetNodeValueIndex("SATURATION1") + 1;
		AuxNodal2 = new double[8];
	}
	else if (Flow_Type == 3)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
		idx_P2 = h_pcs->GetNodeValueIndex("PRESSURE2");
		idx_S0 = h_pcs->GetNodeValueIndex("SATURATION1");
		idx_S = idx_S0;
		idx_Snw = h_pcs->GetNodeValueIndex("SATURATION2") + 1;
		AuxNodal2 = new double[8];
	}

	for (size_t i = 0; i < pcs_vector.size(); i++)
		//      if (pcs_vector[i]->pcs_type_name.find("HEAT") != string::npos) {
		// TF
		if (pcs_vector[i]->getProcessType() == HEAT_TRANSPORT)
		{
			t_pcs = pcs_vector[i];
			break;
		}
	if (T_Flag)
	{
		idx_T0 = t_pcs->GetNodeValueIndex("TEMPERATURE1");
		idx_T1 = idx_T0 + 1;
	}
	//
	// Time unit factor
	time_unit_factor = pcs->time_unit_factor;
#if defined(USE_PETSC)     // || defined(other parallel libs)//03~04.3012. WW
	size_t size = 60;      // dim * nnodesHQ;
	row_ids = new int[size];  //> global indices of local matrix rows
	col_ids = new int[size];  //> global indices of local matrix columns
	local_idx = new int[size];  //> local index for local assemble
// local_matrix = new double[size_m * size_m]; //> local matrix
// local_vec = new double[size_m]; //> local vector
#endif
}

//  Constructor of class Element_DM
CFiniteElementVec::~CFiniteElementVec()
{
	delete B_matrix;
	delete B_matrix_T;
	delete[] dstress;
	delete[] dstrain;
	delete De;
	delete ConsistDep;
	delete AuxMatrix;
	delete AuxMatrix2;  // NW
	delete[] Disp;
	delete[] Temp;
	delete[] T1;
	delete[] Sxx;
	delete[] Syy;
	delete[] Szz;
	delete[] Sxy;
	delete[] pstr;
	delete[] Idx_Strain;
	delete[] Idx_Stress;
	delete[] strain_ne;
	delete[] stress_ne;
	delete[] stress0;
	if (Sxz) delete[] Sxz;
	if (Syz) delete[] Syz;
	if (AuxNodal2) delete[] AuxNodal2;

	if (pcs->Memory_Type == 0)  // Do not store local matrices
	{
		delete Stiffness;
		delete RHS;
		Stiffness = NULL;
		RHS = NULL;
	}

	// 11.07.2011. WW
	if (PressureC) delete PressureC;
	if (PressureC_S) delete PressureC_S;
	if (PressureC_S_dp) delete PressureC_S_dp;

	B_matrix = NULL;
	B_matrix_T = NULL;
	dstress = NULL;
	dstrain = NULL;
	De = NULL;
	ConsistDep = NULL;
	AuxMatrix = NULL;
	AuxMatrix2 = NULL;  // NW
	Disp = NULL;
	Temp = NULL;
	T1 = NULL;
	Sxx = NULL;
	Syy = NULL;
	Szz = NULL;
	Sxy = NULL;
	Sxz = NULL;
	Syz = NULL;
	pstr = NULL;
	delete[] AuxNodal;
	delete[] AuxNodal0;
	delete[] AuxNodal_S0;
	delete[] AuxNodal_S;
	delete[] AuxNodal1;
	AuxNodal = AuxNodal_S0 = AuxNodal_S = AuxNodal1 = NULL;
	//
	Idx_Strain = NULL;
	Idx_Stress = NULL;
	strain_ne = NULL;
	stress_ne = NULL;
	stress0 = NULL;
	//

	// NW
	for (int i = 0; i < (int)vec_B_matrix.size(); i++)
	{
		delete vec_B_matrix[i];
		delete vec_B_matrix_T[i];
		vec_B_matrix[i] = NULL;
		vec_B_matrix_T[i] = NULL;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::SetMaterial(const int EleIndex)

   Aufgabe:
         Set material data for local assembly
   Formalparameter:
           E:

   Programming:
   11/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::SetMaterial()
{
	//......................................................................
	// MAT group
	int MatGroup = MeshElement->GetPatchIndex();
	//......................................................................
	// MSP
	smat = msp_vector[MatGroup];
	smat->axisymmetry = pcs->m_msh->isAxisymmetry();
	// Single yield surface model
	if (smat->Plasticity_type == 2) smat->ResizeMatricesSYS(ele_dim);
	//......................................................................
	// MFP
	if (F_Flag)
	{
		m_mfp = MFPGet("LIQUID");           // YD
		if (!m_mfp) m_mfp = mfp_vector[0];  // OK
	}
	//......................................................................
	// MMP
	m_mmp = mmp_vector[MatGroup];
	//......................................................................
}

/**************************************************************************
   GeoSys - Function: SetMemory

   Aufgabe:
         Set memory for local matrices
   Programmaenderungen:
   01/2005   WW    Erste Version
   05/2005   WW    Dynamic analysis

**************************************************************************/
void CFiniteElementVec::SetMemory()
{
	int size = 0;
	ElementMatrix* EleMat = NULL;

	// Prepare local matrices
	if (pcs->Memory_Type == 0)
	{
		// If local matrices are not stored, resize the matrix
		size = dim * nnodesHQ;
		Stiffness->LimitSize(size, size);
		if (PressureC) PressureC->LimitSize(size, nnodes);
		if (PressureC_S) PressureC_S->LimitSize(size, nnodes);
		if (PressureC_S_dp) PressureC_S_dp->LimitSize(size, nnodes);
		RHS->LimitSize(size);
	}
	else
	{
		EleMat = pcs->Ele_Matrices[Index];
		Stiffness = EleMat->GetStiffness();
		RHS = EleMat->GetRHS();
		if (PressureC) PressureC = EleMat->GetCouplingMatrixA();
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: setB_Matrix(const int LocalIndex)

   Aufgabe:
          Form B matric
   Formalparameter:
           E:
             const int LocalIndex   : Local node index

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::setB_Matrix(const int LocalIndex)
{
	switch (dim)
	{
		case 2:
			// B_11, dN/dx
			(*B_matrix)(0, 0) = dshapefctHQ[LocalIndex];
			// B_12, 0.0
			(*B_matrix)(0, 1) = 0.0;

			if (axisymmetry)  // Axisymmtry
			{
				// B_21, N/r
				(*B_matrix)(1, 0) = shapefctHQ[LocalIndex] / Radius;
				// B_22, 0.0
				(*B_matrix)(1, 1) = 0.0;
				// B_31, 0.0
				(*B_matrix)(2, 0) = 0.0;
				// B_32, dN/dz
				(*B_matrix)(2, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
			}
			else
			{
				// B_21, 0.0
				(*B_matrix)(1, 0) = 0.0;
				// B_22, dN/dy
				(*B_matrix)(1, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
				// B_31, 0.0
				(*B_matrix)(2, 0) = 0.0;
				// B_32, 0.0
				(*B_matrix)(2, 1) = 0.0;
			}
			// B_41, dN/dy
			(*B_matrix)(3, 0) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_42, dN/dx
			(*B_matrix)(3, 1) = dshapefctHQ[LocalIndex];

			break;
		case 3:
			// B_11, dN/dx
			(*B_matrix)(0, 0) = dshapefctHQ[LocalIndex];
			// B_22, dN/dy
			(*B_matrix)(1, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_33, dN/dz
			(*B_matrix)(2, 2) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			//
			// B_41, dN/dy
			(*B_matrix)(3, 0) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_42, dN/dx
			(*B_matrix)(3, 1) = dshapefctHQ[LocalIndex];
			//
			// B_51, dN/dz
			(*B_matrix)(4, 0) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			// B_53, dN/dx
			(*B_matrix)(4, 2) = dshapefctHQ[LocalIndex];
			//
			// B_62, dN/dz
			(*B_matrix)(5, 1) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			// B_63, dN/dy
			(*B_matrix)(5, 2) = dshapefctHQ[nnodesHQ + LocalIndex];

			break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: setTransB_Matrix(const int LocalIndex)

   Aufgabe:
          Form the tanspose of B matric
   Formalparameter:
           E:
             const int LocalIndex   : Local node index

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::setTransB_Matrix(const int LocalIndex)
{
	setB_Matrix(LocalIndex);
	B_matrix->GetTranspose(*B_matrix_T);
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ComputeStrain()

   Aufgabe:
          Compute strains

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::ComputeStrain()
{
	int i, j = 0, k = 0;
	switch (dim)
	{
		case 2:
			for (i = 0; i < ns; i++)
				dstrain[i] = 0.0;
			if (axisymmetry)
			{
				for (i = 0; i < nnodesHQ; i++)
				{
					j = i + nnodesHQ;
					dstrain[0] += Disp[i] * dshapefctHQ[i];
					dstrain[1] += Disp[i] * shapefctHQ[i];
					dstrain[2] += Disp[j] * dshapefctHQ[j];
					dstrain[3] +=
					    Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				}
				dstrain[1] /= Radius;
			}
			else
				for (i = 0; i < nnodesHQ; i++)
				{
					j = i + nnodesHQ;
					dstrain[0] += Disp[i] * dshapefctHQ[i];
					dstrain[1] += Disp[j] * dshapefctHQ[j];
					dstrain[3] +=
					    Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				}
			break;
		case 3:
			for (i = 0; i < ns; i++)
				dstrain[i] = 0.0;
			for (i = 0; i < nnodesHQ; i++)
			{
				j = i + nnodesHQ;
				k = i + 2 * nnodesHQ;
				dstrain[0] += Disp[i] * dshapefctHQ[i];
				dstrain[1] += Disp[j] * dshapefctHQ[j];
				dstrain[2] += Disp[k] * dshapefctHQ[k];
				dstrain[3] +=
				    Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				dstrain[4] +=
				    Disp[i] * dshapefctHQ[k] + Disp[k] * dshapefctHQ[i];
				dstrain[5] +=
				    Disp[j] * dshapefctHQ[k] + Disp[k] * dshapefctHQ[j];
			}
			break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: CalDensity()

   Aufgabe:
          Caculate density of porous medium
   Programming:
   05/2005     WW        Erste Version
 **************************************************************************/
double CFiniteElementVec::CalDensity()
{
	double rho;
	// OK_MFP
	//--------------------------------------------------------------------
	// MFP fluid properties
	double density_fluid = 0.0;
	double porosity = 0.0;
	// double p_g = 0.0;
	double Sw = 0.0;
	int no_phases = (int)mfp_vector.size();
	int i = 0, phase = 0;

	rho = 0.0;
	if (F_Flag)
	{
		if ((no_phases > 0) && (no_phases > phase))
			density_fluid = m_mfp->Density();

		// OK_MMP
		//--------------------------------------------------------------------
		// MMP medium properties
		porosity = m_mmp->Porosity(this);
		// Assume solid density is constant. (*smat->data_Density)(0)
		if (smat->Density() > 0.0)
		{
			Sw = 1.0;  // JT, should be 1.0, unless multiphase (calculate below)
			           // (if unsaturated, fluid density would be negligible...
			           // so still works)
			if (Flow_Type > 0 && Flow_Type != 10)
			{
				Sw = 0.;  // WW
				for (i = 0; i < nnodes; i++)
					Sw += shapefct[i] * AuxNodal_S[i];
			}
			rho = (1. - porosity) * fabs(smat->Density()) +
			      porosity * Sw * density_fluid;

			if (Flow_Type == 2 || Flow_Type == 3)
			{
				/*
				   p_g=0.0;
				   for(i = 0; i< nnodes; i++)
				   p_g += shapefct[i]*AuxNodal2[i];
				   rho += porosity *
				   (1.0-Sw)*COMP_MOL_MASS_AIR*p_g/(GAS_CONSTANT*(Tem+273.15));
				 */
				CFluidProperties* GasProp;
				GasProp = MFPGet("GAS");
				rho += porosity * (1.0 - Sw) * GasProp->Density();
			}
		}
		else
			rho = 0.0;
	}
	else
	    // If negative value is given in the .msp file, gravity by solid is
	    // skipped
	    if (smat->Density() > 0.0)
		rho = smat->Density();
	return rho;
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: ComputeMatrix_RHS(const double fkt)
   Aufgabe:
           as the name

   Programming:
   07/2004   WW
   08/2004   OK   MFP implementation
   01/2010   NW   use chache of B,B^T matrices
 **************************************************************************/
void CFiniteElementVec::ComputeMatrix_RHS(const double fkt, const Matrix* p_D)
{
	int i, j, k, l;
	double rho, fac, fac1, fac2, dN_dx, f_buff;
	fac = fac1 = fac2 = f_buff = 0.0;
	dN_dx = 0.0;
	rho = CalDensity();
	const int nnodesHQ = this->nnodesHQ;
	const int nnodes = this->nnodes;
	const int ele_dim = this->ele_dim;
	const int ns = this->ns;

	// NW cache B, B^T
	for (i = 0; i < nnodesHQ; i++)
	{
		setTransB_Matrix(i);
		(*this->vec_B_matrix[i]) = *B_matrix;
		(*this->vec_B_matrix_T[i]) = *B_matrix_T;
	}
	Matrix* old_B_matrix = B_matrix;
	Matrix* old_B_matrix_T = B_matrix_T;

	Matrix* tmp_B_matrix = NULL;
	Matrix* tmp_B_matrix_T = NULL;
	Matrix* tmp_AuxMatrix = AuxMatrix;
	Matrix* tmp_AuxMatrix2 = AuxMatrix2;
	Matrix* tmp_Stiffness = Stiffness;

	for (i = 0; i < nnodesHQ; i++)
	{
		// NW      setTransB_Matrix(i);
		tmp_B_matrix_T = this->vec_B_matrix_T[i];
		// Local assembly of A*u=int(B^t*sigma) for Newton-Raphson method
		for (j = 0; j < ele_dim; j++)
			for (k = 0; k < ns; k++)
				(*RHS)(j* nnodesHQ + i) +=
				    (*tmp_B_matrix_T)(j, k) * (dstress[k] - stress0[k]) * fkt;
		// TEST             (*B_matrix_T)(j,k)*dstress[k]*fkt;
		if (PreLoad == 11) continue;
		(*tmp_AuxMatrix2) = 0.0;
		// NW
		tmp_B_matrix_T->multi(*p_D, *tmp_AuxMatrix2);
		for (j = 0; j < nnodesHQ; j++)
		{
			// NW          setB_Matrix(j);
			tmp_B_matrix = this->vec_B_matrix[j];
			// Compute stiffness matrix
			(*tmp_AuxMatrix) = 0.0;
			tmp_AuxMatrix2->multi(*tmp_B_matrix, *tmp_AuxMatrix);
			// NW          B_matrix_T->multi(*p_D, *B_matrix, *AuxMatrix);

			// Local assembly of stiffness matrix
			for (k = 0; k < ele_dim; k++)
				for (l = 0; l < ele_dim; l++)
					(*tmp_Stiffness)(i + k * nnodesHQ, j + l * nnodesHQ) +=
					    (*tmp_AuxMatrix)(k, l) * fkt;
		}  // loop j
	}      // loop i

	// should restore pointer NW
	B_matrix = old_B_matrix;
	B_matrix_T = old_B_matrix_T;

	//---------------------------------------------------------
	// Assemble coupling matrix
	//---------------------------------------------------------
	// LoadFactor: factor of incremental loading, prescibed in rf_pcs.cpp

	// 07.2011 WW
	if ((PressureC || PressureC_S || PressureC_S_dp) && !PreLoad)
	{
		fac = LoadFactor * fkt;

		// 07.2011. WW
		if (PressureC_S || PressureC_S_dp)
		{
			// Pressure 1
			fac2 = interpolate(AuxNodal0);
			// Saturation of phase 1
			fac1 = m_mmp->SaturationCapillaryPressureFunction(fac2);
			if (PressureC_S_dp)
				fac2 = fac1 -
				       fac2 * m_mmp->PressureSaturationDependency(fac1, true);
			// JT: dSdP now returns actual sign (<0)
		}

		if (axisymmetry)
			for (k = 0; k < nnodesHQ; k++)
			{
				for (l = 0; l < nnodes; l++)
					for (j = 0; j < ele_dim; j++)
					{
						dN_dx = dshapefctHQ[nnodesHQ * j + k];
						if (j == 0) dN_dx += shapefctHQ[k] / Radius;

						f_buff = fac * dN_dx * shapefct[l];
						(*PressureC)(nnodesHQ* j + k, l) += f_buff;
						if (PressureC_S)
							(*PressureC_S)(nnodesHQ* j + k, l) += f_buff * fac1;
						if (PressureC_S_dp)
							(*PressureC_S_dp)(nnodesHQ* j + k, l) +=
							    f_buff * fac2;
					}
			}
		else
			for (k = 0; k < nnodesHQ; k++)
			{
				for (l = 0; l < nnodes; l++)
					for (j = 0; j < ele_dim; j++)
					{
						f_buff =
						    fac * dshapefctHQ[nnodesHQ * j + k] * shapefct[l];
						(*PressureC)(nnodesHQ* j + k, l) += f_buff;
						if (PressureC_S)
							(*PressureC_S)(nnodesHQ* j + k, l) += f_buff * fac1;
						if (PressureC_S_dp)
							(*PressureC_S_dp)(nnodesHQ* j + k, l) +=
							    f_buff * fac2;
					}
			}
	}
	//---------------------------------------------------------
	// Assemble gravity force vector
	//---------------------------------------------------------
	if (rho > 0.0 && GravityForce)
	{
		// 2D, in y-direction
		// 3D, in z-direction
		i = (ele_dim - 1) * nnodesHQ;
		const double coeff = LoadFactor * rho * smat->grav_const * fkt;
		for (k = 0; k < nnodesHQ; k++)
			(*RHS)(i + k) += coeff * shapefctHQ[k];
		//        (*RHS)(i+k) += LoadFactor * rho * smat->grav_const *
		//        shapefctHQ[k] * fkt;
	}
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: LocalAssembly
   Aufgabe:
           Compute the local finite element matrices
   Formalparameter:
           E:
         const int update  : indicator to update stress and strain only

   Programming:
   06/2004   WW   Generalize for different element types as a member of class
   05/2005   WW   ...
   08/2010   WW   JFNK method
 **************************************************************************/
void CFiniteElementVec::LocalAssembly(const int update)
{
	int j;

	Index = MeshElement->GetIndex();
	SetMemory();
	SetMaterial();
	// 12.2009. WW
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];
	ns = 4;
	if (MeshElement->GetDimension() == 3) ns = 6;

#ifdef USE_PETSC
	if (MeshElement->getGhostNodeIndices())  // ghost nodes pcs->pcs_number_of_primary_nvals
	{
		act_nodes = MeshElement->getGhostNodeIndices()[0];
		act_nodes_h = MeshElement->getGhostNodeIndices()[1];
		for (int i = 0; i < act_nodes_h; i++)
			local_idx[i] = MeshElement->getGhostNodeIndices()[i + 2];
	}
	else
	{
		act_nodes = nnodes;
		act_nodes_h = nnodesHQ;
		for (int i = 0; i < nnodesHQ; i++)
			local_idx[i] = i;
	}
#endif

	// WW
	{
#ifdef NEW_EQS
		b_rhs = pcs->eqs_new->b;
#endif
	}
	(*RHS) = 0.0;
	(*Stiffness) = 0.0;
	// 07.2011. WW
	if (PressureC) (*PressureC) = 0.0;
	if (PressureC_S) (*PressureC_S) = 0.0;
	if (PressureC_S_dp) (*PressureC_S_dp) = 0.0;
	// WW
	{
		for (int i = 0; i < nnodesHQ; i++)
			eqs_number[i] = MeshElement->GetNode(i)->GetEquationIndex();
	}

	// For strain and stress extropolation all element types
	// Number of elements associated to nodes
	for (int i = 0; i < nnodes; i++)
		dbuff[i] = (double)MeshElement->GetNode(i)->getConnectedElementIDs().size();

	// Get displacement_n
	for (size_t i = 0; i < dim; i++)
		for (j = 0; j < nnodesHQ; j++)
			Disp[j + i * nnodesHQ] = pcs->GetNodeValue(nodes[j], Idx_dm0[i]);

	// Get saturation of element nodes
	if (Flow_Type > 0 && Flow_Type != 10)
	{
		for (int i = 0; i < nnodes; i++)
		{
			AuxNodal_S[i] = h_pcs->GetNodeValue(nodes[i], idx_S);
			AuxNodal_S0[i] = h_pcs->GetNodeValue(nodes[i], idx_S0);
		}
		if (Flow_Type == 2 || Flow_Type == 3)  // 09.10.2009 PCH

			for (int i = 0; i < nnodes; i++)
			{
				AuxNodal0[i] = h_pcs->GetNodeValue(nodes[i], idx_P1);
				AuxNodal2[i] = h_pcs->GetNodeValue(nodes[i], idx_P2);
			}
		// 12.03.2008 WW
		if ((Flow_Type == 1 || Flow_Type == 2) &&
		    (smat->SwellingPressureType == 3 ||
		     smat->SwellingPressureType == 4))
		{
			double fac = 1.0;
			if (Flow_Type == 1) fac = -1.0;
			for (int i = 0; i < nnodes; i++)
			{
				// Pc
				AuxNodal1[i] = fac * h_pcs->GetNodeValue(nodes[i], idx_P1 - 1);
				// dPc
				AuxNodal[i] = fac * (h_pcs->GetNodeValue(nodes[i], idx_P1) -
				                     h_pcs->GetNodeValue(nodes[i], idx_P1 - 1));
			}
		}
	}
	//

	//----------------------------------------------------
	LocalAssembly_continuum(update);

	if (update == 0)
	{
		GlobalAssembly();
		// Output matrices
		if (pcs->Write_Matrix)
		{
			(*pcs->matrix_file) << "### Element: " << MeshElement->GetGlobalIndex() << "\n";
			(*pcs->matrix_file) << "---Stiffness matrix: "
			                    << "\n";
			Stiffness->Write(*pcs->matrix_file);
			(*pcs->matrix_file) << "---RHS: "
			                    << "\n";
			RHS->Write(*pcs->matrix_file);
			if (PressureC)
			{
				(*pcs->matrix_file) << "Pressue coupling matrix: "
				                    << "\n";
				PressureC->Write(*pcs->matrix_file);
			}
			// 07.2011. WW
			if (PressureC_S)
			{
				(*pcs->matrix_file)
				    << "Saturation depedent pressue coupling matrix: "
				    << "\n";
				PressureC_S->Write(*pcs->matrix_file);
			}
			if (PressureC_S_dp)
			{
				(*pcs->matrix_file) << "Jacobi pressue coupling matrix: "
				                    << "\n";
				PressureC_S_dp->Write(*pcs->matrix_file);
			}
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
   12/2009   WW New excavtion approach
 **************************************************************************/
bool CFiniteElementVec::GlobalAssembly()
{

	GlobalAssembly_RHS();
	if (PreLoad == 11) return true;

	GlobalAssembly_Stiffness();

#ifdef USE_PETSC
	// ScreenMessage2("-> add2GlobalMatrixII()\n");
	add2GlobalMatrixII();
#endif

	return true;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_Stiffness()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
 **************************************************************************/
void CFiniteElementVec::GlobalAssembly_Stiffness()
{
	double biot = smat->biot_const;
#if defined(NEW_EQS)
	CSparseMatrix* A = pcs->eqs_new->A;
#endif

#ifndef USE_PETSC
	double f1 = 1.0;
	// Assemble stiffness matrix
	for (int i = 0; i < nnodesHQ; i++)
	{
		const long eqs_number_i = eqs_number[i];
		for (int j = 0; j < nnodesHQ; j++)
		{
			const long eqs_number_j = eqs_number[j];
			// Local assembly of stiffness matrix
			for (size_t k = 0; k < ele_dim; k++)
			{
				const long globalRowId = eqs_number_i + NodeShift[k];
				const long localRowId = i + k * nnodesHQ;
				for (size_t l = 0; l < ele_dim; l++)
				{
					double globalColId = eqs_number_j + NodeShift[l];
					double val =
					    f1 * (*Stiffness)(localRowId, j + l * nnodesHQ);
#ifdef NEW_EQS
					double& a = (*A)(globalRowId, globalColId);
#pragma omp atomic
					a += val;
//(*A)(globalRowId,globalColId) += val;
#endif
				}
			}
		}  // loop j
	}      // loop i
#endif

	// TEST OUT
	// Stiffness->Write();
	if (pcs->type / 40 != 1)  // Not monolithic scheme
		return;

	double f2 = -1.0;
	if (PressureC)
	{
		int i = 0;               // phase
		if (Flow_Type == 2)  // Multi-phase-flow
			i = 1;
		GlobalAssembly_PressureCoupling(PressureC, f2 * biot, i);
	}
	// H2: p_g- S_w*p_c
	if (PressureC_S)
		GlobalAssembly_PressureCoupling(PressureC_S, -f2 * biot, 0);
	if (PressureC_S_dp)
		GlobalAssembly_PressureCoupling(PressureC_S_dp, -f2 * biot, 0);

	// TEST OUT
	// PressureC->Write();
}
//--------------------------------------------------------------------------
/*!
   \brief Assembe the pressure coupling matrix

    to the global stiffness matrix in the monolithic scheme

   \param pCMatrix: the matrix
   \param fct: factor
   \param phase: phasse index

    07.2011. WW
 */
#if defined(USE_PETSC)  // || defined(other parallel libs)//10.3012. WW
void CFiniteElementVec::GlobalAssembly_PressureCoupling(Matrix*,
                                                        double,
                                                        const int)
{
}
#else
void CFiniteElementVec::GlobalAssembly_PressureCoupling(Matrix* pCMatrix,
                                                        double fct,
                                                        const int phase)
{
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;
	A = pcs->eqs_new->A;
#endif

	int dim_shift = dim + phase;
	// Add pressure coupling matrix to the stifness matrix
	for (int i = 0; i < nnodesHQ; i++)
	{
		for (int j = 0; j < nnodes; j++)
		{
			for (size_t k = 0; k < ele_dim; k++)
			{
#ifdef NEW_EQS
				(*A)(NodeShift[k] + eqs_number[i],
				     NodeShift[dim_shift] + eqs_number[j]) +=
				    fct * (*pCMatrix)(nnodesHQ* k + i, j);
#endif
			}
		}
	}
}
#endif

#ifdef USE_PETSC
//------------------------------------------------------
// Add the local stiff matrix to the global one
void CFiniteElementVec::add2GlobalMatrixII()
{
	int const ndof = pcs->pcs_number_of_primary_nvals;

	int m_dim, n_dim;
	double* local_matrix = NULL;
	double* local_vec = NULL;

//#define assmb_petsc_test
#ifdef assmb_petsc_test
	char rank_char[10];
	sprintf(rank_char, "%d", eqs->getMPI_Rank());
	std::string fname = FileName + rank_char + "_e_matrix.txt";
	std::ofstream os_t(fname.c_str(), std::ios::app);
	os_t << "\n==================================================\n";
#endif

	m_dim = nnodesHQ * ndof;
	n_dim = m_dim;
	local_matrix = Stiffness->getEntryArray();
	local_vec = RHS->getEntryArray();

	for (int i = 0; i < nnodesHQ; i++)
	{
		MeshLib::CNode* node = MeshElement->GetNode(i);
		const int offset = node->GetEquationIndex(true) * ndof;
		const bool isGhost = !pcs->m_msh->isNodeLocal(node->GetIndex());
		for (int k = 0; k < ndof; k++)
		{
			const int ki = k * nnodesHQ + i;
			row_ids[ki] = isGhost ? -1 : (offset + k);
			col_ids[ki] = offset + k;
		}
	}

	pcs->eqs_new->addMatrixEntries(m_dim, row_ids, n_dim, col_ids, &local_matrix[0]);

	static double temp_vec[100];
	for (int i = 0; i < m_dim; i++)
		temp_vec[i] = -local_vec[i];  // r -= RHS
	pcs->eqs_new->setArrayValues(1, m_dim, row_ids, &temp_vec[0]);

}
#endif

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_RHS()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
   05/2005   WW dyn
 **************************************************************************/
void CFiniteElementVec::GlobalAssembly_RHS()
{

	bool Residual = false;
	if (Flow_Type >= 0)
	{
		if (pcs->type / 10 == 4)  // Monolithic scheme
		{
			// If nonlinear deformation
			if (pcs_deformation > 100) Residual = true;
		}
		else  // Partitioned scheme
			Residual = true;
	}

	// Assemble coupling matrix
	// If dynamic GetNodeValue(nodes[i],idx_P0) = 0;
	if (Residual)
	{
		const double biot = smat->biot_const;
		double* p0 = pcs->GetInitialFluidPressure();  // NW should calculate (p
		                                              // - p0) because OGS
		                                              // calculates (stress -
		                                              // stress0)
		switch (Flow_Type)
		{
			case 0:  // Liquid flow
				// For monolithic scheme and liquid flow, the limit of positive
				// pressure must be removed
				for (int i = 0; i < nnodes; i++)
				{
					double val_n = h_pcs->GetNodeValue(nodes[i], idx_P1);
					//                AuxNodal[i] = LoadFactor*( val_n
					//                -Max(pcs->GetNodeValue(nodes[i],idx_P0),0.0));

					// Initial pressure should be subtracted, i.e. (p-p0)
					// because DEFORMATION
					// calculates stress balance of changes from the initial
					// stress.
					// NW 28.08.2012
					if (pcs->calcDiffFromStress0 && p0 != NULL)
					{
						val_n -= p0[nodes[i]];
					}
					AuxNodal[i] = LoadFactor * val_n;
				}
				break;
			case 10:  // Ground_flow. Will be merged to case 0
				for (int i = 0; i < nnodes; i++)
					AuxNodal[i] =
					    LoadFactor * h_pcs->GetNodeValue(nodes[i], idx_P1);
				break;
			case 1:  // Richards flow
				//
				for (int i = 0; i < nnodes; i++)
				{
					double val_n = h_pcs->GetNodeValue(nodes[i], idx_P1);
					if (biot < 0.0 && val_n < 0.0)
						AuxNodal[i] = 0.0;
					else
						AuxNodal[i] = LoadFactor * S_Water * val_n;
				}
				break;
			case 2:
			{  // Multi-phase-flow: p_g-Sw*p_c
				// 07.2011. WW
				const int dim_times_nnodesHQ(dim * nnodesHQ);
				for (int i = 0; i < dim_times_nnodesHQ; i++)
					AuxNodal1[i] = 0.0;

				if (smat->bishop_model > 0)
				{
					for (int i = 0; i < nnodes; i++)
					{
						double bishop_coef = 1.;  // bishop
						double S_e = 1.;
						switch (smat->bishop_model)
						{
							case 1:
								bishop_coef = smat->bishop_model_value;
								break;
							case 2:
								S_e = m_mmp->GetEffectiveSaturationForPerm(
								    AuxNodal_S[i], 0);
								bishop_coef =
								    pow(S_e, smat->bishop_model_value);
								break;
							default:
								break;
						}
						double val_n = 0;
						if (smat->bishop_model == 1 ||
						    smat->bishop_model == 2)  // pg-bishop*pc 05.2011 WX
							val_n = h_pcs->GetNodeValue(nodes[i], idx_P2) -
							        bishop_coef *
							            h_pcs->GetNodeValue(nodes[i], idx_P1);
						else
							val_n = h_pcs->GetNodeValue(nodes[i],
							                            idx_P2)  // pg - Sw*pc
							        -
							        AuxNodal_S[i] *
							            h_pcs->GetNodeValue(nodes[i], idx_P1);
						val_n = h_pcs->GetNodeValue(nodes[i], idx_P2) -
						        AuxNodal_S[i] *
						            h_pcs->GetNodeValue(nodes[i], idx_P1);
						if (biot < 0.0 && val_n < 0.0)
							AuxNodal[i] = 0.0;
						else
							AuxNodal[i] = val_n * LoadFactor;
					}

					PressureC->multi(AuxNodal, AuxNodal1);
				}
				else
				{
					PressureC->multi(AuxNodal2, AuxNodal1, LoadFactor);
					PressureC_S->multi(AuxNodal0, AuxNodal1, -1.0 * LoadFactor);
				}

				break;
			}
			case 3:  // Multi-phase-flow: SwPw+SgPg	// PCH 05.05.2009
				for (int i = 0; i < nnodes; i++)
				{
					double Snw = h_pcs->GetNodeValue(nodes[i], idx_Snw);
					double Sw = 1.0 - Snw;
					double Pw = h_pcs->GetNodeValue(nodes[i], idx_P1);
					double Pnw = h_pcs->GetNodeValue(nodes[i], idx_P2);
					double val_n = Sw * Pw + Snw * Pnw;
					if (biot < 0.0 && val_n < 0.0)
						AuxNodal[i] = 0.0;
					else
						AuxNodal[i] = val_n * LoadFactor;
				}
				break;
		}

		const int dim_times_nnodesHQ(dim * nnodesHQ);
		// Coupling effect to RHS
		if (Flow_Type != 2)
		{
			for (int i = 0; i < dim_times_nnodesHQ; i++)
				AuxNodal1[i] = 0.0;
			PressureC->multi(AuxNodal, AuxNodal1);
		}
		for (int i = 0; i < dim_times_nnodesHQ; i++)
			(*RHS)(i) -= fabs(biot) * AuxNodal1[i];
	}  // End if partioned


// RHS->Write();

#ifndef USE_PETSC
	//	std::cerr << "e: " << Index << "\n";
	for (size_t i = 0; i < dim; i++)
	{
		for (int j = 0; j < nnodesHQ; j++)
		{
			//		    std::cerr << eqs_number[j] << "," << NodeShift[i] << ","
			//<< i * nnodesHQ + j << "\n";
			long global_id = eqs_number[j] + NodeShift[i];
			double val = (*RHS)(i* nnodesHQ + j);
#pragma omp atomic
			b_rhs[global_id] -= val;
		}
	}
#endif
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: LocalAssembly_continuum()
   Aufgabe:
           Compute the local finite element matrices with the framework
        of continuum assumption
   Formalparameter:
           E:
         const int update  : indicator to update stress and strain only

   Programming:
   02/2000   OK   Elastische Elemente
   09/2002   WW   Local assemby of stiffness matrix of elasto-plastic
   tangential model
   Local assemby of residual
   07/2003   WW   Quadratic triangle element
   06/2004   WW   Generalize for different element types as a member of class
   12/2005   WW   Creep
 **************************************************************************/
void CFiniteElementVec::LocalAssembly_continuum(const int update)
{
	long i;

	Matrix* p_D = NULL;
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];

	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp = 0;
	gp_t = 0;
	double fkt = 0.0;

	const int PModel = smat->Plasticity_type;
	double dPhi = 0.0;  // Sclar factor for the plastic strain
	double dS = 0.0;

	double ThermalExpansion = 0.0;
	double t1 = 0.0;
	bool Strain_TCS = false;

	//
	ThermalExpansion = 0.0;
	// Thermal effect
	if (smat->Thermal_Expansion() > 0.0)
		ThermalExpansion = smat->Thermal_Expansion();

	// For swelling pressure;
	double deporo = 0.0;
	int PoroModel = m_mmp->porosity_model;
	if (PoroModel == 4)
		deporo = h_pcs->GetElementValue(
		             Index, h_pcs->GetElementValueIndex("n_sw_rate")) /
		         (double)ele_dim;

	if (T_Flag)
	{
		for (i = 0; i < nnodes; i++)
		{
			T1[i] = t_pcs->GetNodeValue(nodes[i], idx_T1);
			Temp[i] = t_pcs->GetNodeValue(nodes[i], idx_T1) -
			          t_pcs->GetNodeValue(nodes[i], idx_T0);
		}
	}
	//

	if (PModel == 1 || PModel == 10 || PModel == 11)
		smat->CalulateCoefficent_DP();
	//
	if (PModel != 3 && smat->Youngs_mode != 2)  // modified due to transverse
	                                            // isotropic elasticity: UJG
	                                            // 24.11.2009
	{
		if (smat->Youngs_mode < 10 || smat->Youngs_mode > 13)
		{
			smat->Calculate_Lame_Constant();
			smat->ElasticConsitutive(ele_dim, De);
		}
		else
			*De = *(smat->getD_tran());  // UJG/WW
	}

	if (PModel == 5) smat->CalculateCoefficent_HOEKBROWN();
	                                                         /*
	                                                            string fname=FileName+"_D.txt";
	                                                            ofstream out_f(fname.c_str());
	                                                            De->Write(out_f);
	                                                          */

	/*
	   //TEST
	   fstream oss;
	   if(update)
	   {
	   char tf_name[10];
	   #ifdef USE_MPI
	   sprintf(tf_name,"%d",myrank);
	    string fname = FileName+tf_name+".stress";
	   #else
	    string fname = FileName+".stress";
	   #endif
	   oss.open(fname.c_str(), ios::app|ios::out);
	   //    oss.open(fname.c_str(), ios::trunc|ios::out);
	   oss<<"\nElement  "<<Index<<endl;
	   oss<<endl;

	   oss<<"Diaplacement "<<endl;
	   for(i=0;i<nnodesHQ;i++)
	   {
	   oss<<nodes[i]<<"  ";
	   for(int ii=0; ii<dim; ii++)
	   oss<<Disp[ii*nnodesHQ+i]<<"  ";
	   oss<<endl;
	   }
	   oss<<"Temperature "<<endl;
	   for(i=0; i<nnodes;i++)
	   oss<<Temp[i]<<"  ";
	   oss<<endl;
	   oss.close();
	   }
	 */
	//
	if (PoroModel == 4 || T_Flag || smat->Creep_mode > 0) Strain_TCS = true;
	//
	if (smat->CreepModel() == 1000)  // HL_ODS
		smat->CleanTrBuffer_HL_ODS();
	for (i = 0; i < ns; i++)
		stress0[i] = .0;

	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeGradShapefct(2);
		ComputeShapefct(2);
		if (smat->Youngs_mode == 2)  // WW/UJG. 22.01.2009
		{
			smat->CalcYoungs_SVV(CalcStrain_v());
			smat->ElasticConsitutive(ele_dim, De);
		}
		if (smat->Youngs_mode == 3)  // AJ. 16.06.2014
		{
			smat->CalcYoungs_Drained(CalcStress_eff());
			smat->ElasticConsitutive(ele_dim, De);
		}

		ComputeStrain();
		if (update) RecordGuassStrain(gp, gp_r, gp_s, gp_t);
		if (F_Flag || T_Flag)
			ComputeShapefct(1);  // Linear order interpolation function
		//---------------------------------------------------------
		// Material properties (Integration of the stress)
		//---------------------------------------------------------
		// Initial the stress vector
		if (PModel != 3)
		{
			for (i = 0; i < ns; i++)
				dstress[i] = 0.0;
			De->multi(dstrain, dstress);
		}

		//---------------------------------------------------------
		// Integrate the stress by return mapping:
		//---------------------------------------------------------
		switch (PModel)
		{
			case -1:  // Pure elasticity
				// Non-linear elasticity: TE model. 10.03.2008. WW
				for (i = 0; i < ns; i++)
					dstress[i] += (*eleV_DM->Stress)(i, gp);
				break;
			case 1:  // Drucker-Prager model
				if (smat->StressIntegrationDP(gp, eleV_DM, dstress, dPhi,
				                              update))

					// WW DevStress = smat->devS;
					smat->ConsistentTangentialDP(ConsistDep, dPhi, ele_dim);
				break;
			case 10:  // Drucker-Prager model, direct integration. 02/06 WW
				if (smat->DirectStressIntegrationDP(gp, eleV_DM, dstress,
				                                    update))
				{
					*ConsistDep = *De;
					smat->TangentialDP(ConsistDep);
					dPhi = 1.0;
				}
				break;
			case 11:  // WX: 08.2010
			{
				double mm = 0.;  // WX:09.2010. for DP with Tension.
				switch (smat->DirectStressIntegrationDPwithTension(
				    gp, De, eleV_DM, dstress, update, mm))
				{
					case 1:
					{
						*ConsistDep = *De;
						smat->TangentialDP2(ConsistDep);
						dPhi = 1.0;
					}
					break;
					case 2:
					{
						*ConsistDep = *De;
						smat->TangentialDPwithTension(ConsistDep, mm);
						dPhi = 1.0;
					}
					break;
					case 3:
					{
						*ConsistDep = *De;
						smat->TangentialDPwithTensionCorner(ConsistDep, mm);
						dPhi = 1.0;
					}
					break;
					default:
						break;
				}
				break;
			}
			case 2:  // Rotational hardening model
				// Compute stesses and plastic multi-plier
				dPhi = 0.0;
				if (smat->CalStress_and_TangentialMatrix_SYS(
				        gp, eleV_DM, De, ConsistDep, dstress, update) > 0)
					dPhi = 1.0;
				break;
			case 3:  // Generalized Cam-Clay model
				for (i = 0; i < ns; i++)
					dstress[i] = dstrain[i];
				//
				if (smat->SwellingPressureType == 3)
				{
					double suc = interpolate(AuxNodal1);
					double dsuc = interpolate(AuxNodal);
					(*smat->data_Youngs)(7) = suc;
					(*smat->data_Youngs)(8) = dsuc;
					smat->CalStress_and_TangentialMatrix_CC_SubStep(
					    gp, eleV_DM, dstress, ConsistDep, update);
					//
					// double pn = -((*eleV_DM->Stress)(0,
					// gp)+(*eleV_DM->Stress)(1, gp)+
					//             (*eleV_DM->Stress)(2, gp))/3.0;
					// de_vsw =
					// -smat->TEPSwellingParameter(pn,suc)*dsuc/(suc+1.0e5);
				}
				/*
				   else if (smat->SwellingPressureType==4)  //Power. 07.05.2008
				   WW
				   {
				   double suc = interpolate(AuxNodal1);
				   double dsuc = interpolate(AuxNodal);
				   smat->TEPSwellingParameter_kis(suc);
				   S_Water = m_mmp->SaturationCapillaryPressureFunction(suc);
				   dS = S_Water -
				   m_mmp->SaturationCapillaryPressureFunction(suc-dsuc);
				   de_vsw = pow(S_Water, (*smat->data_Youngs)(2) )*dS;
				   }
				 */
				else
					smat->CalStress_and_TangentialMatrix_CC(
					    gp, eleV_DM, dstress, ConsistDep, update);

				dPhi = 1.0;
				break;
			case 4:  // Mohr-Coloumb	//WX:10.2010
				if (smat->DirectStressIntegrationMOHR(gp, eleV_DM, dstress,
				                                      update, De))
				{
					*ConsistDep = *De;
					// also for tension
					smat->TangentialMohrShear(ConsistDep);
					// ConsistDep->Write();
					dPhi = 1.0;
				}
				break;
				/*case 5:
				   if(smat->StressIntegrationHoekBrown(gp, eleV_DM, dstress,
				   update, De))
				   {
				   *ConsistDep = *De;
				         smat->TangentialHoekBrown(ConsistDep);		//also for
				   tension
				         //ConsistDep->Write();
				         dPhi = 1.0;
				   }
				   break;*/
		}
		// --------------------------------------------------------------------
		// Stress increment by heat, swelling, or heat
		//
		if (Strain_TCS)
		{
			if (PModel == 3) smat->ElasticConsitutive(ele_dim, De);
			for (i = 0; i < ns; i++)
				strain_ne[i] = 0.0;
			if (PoroModel == 4)  // For swelling pressure

				for (i = 0; i < 3; i++)
					strain_ne[i] -= deporo;
			//
			if (T_Flag)  // Contribution by thermal expansion
			{
				Tem = 0.0;
				t1 = 0.0;
				for (i = 0; i < nnodes; i++)
				{
					Tem += shapefct[i] * Temp[i];
					t1 += shapefct[i] * T1[i];
				}
				for (i = 0; i < 3;
				     i++)  // JT: This was commented. SHOULDN'T BE!
					strain_ne[i] -= ThermalExpansion * Tem;
			}
			// Strain increment by creep
			if (smat->Creep_mode == 1 || smat->Creep_mode == 2)
			{
				for (i = 0; i < ns; i++)
					stress_ne[i] = (*eleV_DM->Stress)(i, gp);
				smat->AddStain_by_Creep(ns, stress_ne, strain_ne, t1);
			}
			if (smat->Creep_mode == 1000)  // HL_ODS. Strain increment by creep
			{
				for (i = 0; i < ns; i++)
					stress_ne[i] = (*eleV_DM->Stress)(i, gp);
				smat->AddStain_by_HL_ODS(eleV_DM, stress_ne, strain_ne, t1);
			}
			// Stress deduced by thermal or swelling strain incremental:
			De->multi(strain_ne, dstress);
			for (i = 0; i < ns;
			     i++)  // JT: This was commented. It shouldn't be.
				dstrain[i] += strain_ne[i];
		}
		// Fluid coupling;
		S_Water = 1.0;
		if (Flow_Type > 0 && Flow_Type != 10)
			S_Water = interpolate(AuxNodal_S, 1);
		// Decovalex. Swelling pressure
		if (smat->SwellingPressureType == 1)
		{
			dS = -interpolate(AuxNodal_S0, 1);
			dS += S_Water;
			for (i = 0; i < 3; i++)
				dstress[i] -= 2.0 * S_Water * dS * smat->Max_SwellingPressure;
		}
		else if (smat->SwellingPressureType == 2)  // LBNL's model
		{
			dS = -interpolate(AuxNodal_S0, 1);
			dS += S_Water;
			for (i = 0; i < 3; i++)
				dstress[i] -= dS * smat->Max_SwellingPressure;
		}
		/*
		   else if(smat->SwellingPressureType==3||smat->SwellingPressureType==4)
		   // TEP model
		   {
		   for (i = 0; i < 3; i++)
		      strain_ne[i] = -de_vsw;
		   for (i = 3; i < ns; i++)
		      strain_ne[i] = 0.;
		   smat->ElasticConsitutive(ele_dim, De);
		   De->multi(strain_ne, dstress);
		   }
		 */
		// Assemble matrices and RHS
		if (update < 1)
		{
			//---------------------------------------------------------
			// Assemble matrices and RHS
			//---------------------------------------------------------
			if (dPhi <= 0.0)
				p_D = De;
			else
				p_D = ConsistDep;
			if (pcs->calcDiffFromStress0)
			{
				for (i = 0; i < ns; i++)
					stress0[i] = (*eleV_DM->Stress0)(i, gp);
			}
			ComputeMatrix_RHS(fkt, p_D);
		}
		else  // Update stress

			for (i = 0; i < ns; i++)
			{
				(*eleV_DM->Stress)(i, gp) = dstress[i];
				(*eleV_DM->Strain)(i, gp) += dstrain[i];
			}
	}
	// The mapping of Gauss point strain to element nodes
	if (update)
		ExtropolateGuassStrain();
	else if (smat->Creep_mode == 1000)  // HL_ODS. Strain increment by creep
		smat->AccumulateEtr_HL_ODS(eleV_DM, nGaussPoints);

	/*
	   //TEST

	      //TEST
	        if(update&&Index==0)
	        {
	           oss<<" Stress "<<endl;
	           eleV_DM->Stress->Write(oss);
	           oss.close();
	        }
	 */
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::RecordGuassValues()
   Aufgabe:
           Accumulate stress at each nodes
   Formalparameter:
           E:

   Programming:
   06/2004   WW
 **************************************************************************/
bool CFiniteElementVec::RecordGuassStrain(const int gp, const int gp_r,
                                          const int gp_s, int gp_t)
{
	int LoIndex = 0;

	//---------------------------------------------------------
	// Accumulate strains
	//---------------------------------------------------------
	switch (MeshElement->GetElementType())
	{
		case MshElemType::QUAD:  // Quadralateral
			LoIndex = GetLocalIndex(gp_r, gp_s, gp_t);
			Sxx[LoIndex] = dstrain[0];
			Syy[LoIndex] = dstrain[1];
			Sxy[LoIndex] = dstrain[3];
			Szz[LoIndex] = dstrain[2];
			break;
		case MshElemType::TRIANGLE:  // Triangle
			Sxx[gp] = dstrain[0];
			Syy[gp] = dstrain[1];
			Szz[gp] = dstrain[2];
			Sxy[gp] = dstrain[3];
			break;
		case MshElemType::HEXAHEDRON:  // Hexahedra
			LoIndex = GetLocalIndex(gp_r, gp_s, gp_t);
			if (LoIndex < 0) return false;
			Sxx[LoIndex] = dstrain[0];
			Syy[LoIndex] = dstrain[1];
			Szz[LoIndex] = dstrain[2];
			Sxy[LoIndex] = dstrain[3];
			Sxz[LoIndex] = dstrain[4];
			Syz[LoIndex] = dstrain[5];
			break;
		case MshElemType::TETRAHEDRON:  // Tedrahedra
			Sxx[gp] = dstrain[0];
			Syy[gp] = dstrain[1];
			Szz[gp] = dstrain[2];
			Sxy[gp] = dstrain[3];
			Sxz[gp] = dstrain[4];
			Syz[gp] = dstrain[5];
			break;
		case MshElemType::PYRAMID:
			Sxx[gp] = dstrain[0];
			Syy[gp] = dstrain[1];
			Szz[gp] = dstrain[2];
			Sxy[gp] = dstrain[3];
			Sxz[gp] = dstrain[4];
			Syz[gp] = dstrain[5];
			break;
		default:
			break;
			// 3D
	}
	return false;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ExtropolateGuassStrain()
   Aufgabe:
           Extropolate the Gauss point strains to nodes
   Formalparameter:
           E:

   Programming:
   06/2004   WW
   02/2007   Make it work for all 2nd variables
 **************************************************************************/
void CFiniteElementVec::ExtropolateGuassStrain()
{
	int i, j;
	//  int l1,l2,l3,l4; //, counter;
	double ESxx, ESyy, ESzz, ESxy, ESxz, ESyz;
	double avgESxx, avgESyy, avgESzz, avgESxy, avgESxz, avgESyz;
	int i_s, i_e, ish, k = 0;
	gp = 0;
	// double Area1, Area2, Tol=10e-9;

	// l1=l2=l3=l4=0;
	MshElemType::type ElementType = MeshElement->GetElementType();
	if (ElementType == MshElemType::QUAD ||
	    ElementType == MshElemType::HEXAHEDRON)
		Xi_p = CalcXi_p();

	//
	i_s = 0;
	i_e = nnodes;
	ish = 0;
	if (ElementType == MshElemType::TETRAHEDRON)  // tet
	{
		i_s = 1;
		i_e = nnodes + 1;
		ish = 1;
	}
	//---------------------------------------------------------
	// Mapping Gauss point strains to nodes and update nodes
	// strains:
	//---------------------------------------------------------
	avgESxx = avgESyy = avgESzz = avgESxy = avgESxz = avgESyz = 0.0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
	{
		// average
		avgESxx = CalcAverageGaussPointValues(Sxx);
		avgESyy = CalcAverageGaussPointValues(Syy);
		avgESzz = CalcAverageGaussPointValues(Szz);
		avgESxy = CalcAverageGaussPointValues(Sxy);
		avgESxz = CalcAverageGaussPointValues(Sxz);
		avgESyz = CalcAverageGaussPointValues(Syz);
	}

	for (i = 0; i < nnodes; i++)
	{
		ESxx = ESyy = ESzz = ESxy = ESxz = ESyz = 0.0;

		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			SetExtropoGaussPoints(i);
			ComputeShapefct(1);  // Linear interpolation function
			//
			for (j = i_s; j < i_e; j++)
			{
				k = j - ish;
				ESxx += Sxx[j] * shapefct[k];
				ESyy += Syy[j] * shapefct[k];
				ESxy += Sxy[j] * shapefct[k];
				ESzz += Szz[j] * shapefct[k];
				if (ele_dim == 3)
				{
					ESxz += Sxz[j] * shapefct[k];
					ESyz += Syz[j] * shapefct[k];
				}
			}
		}
		else if (this->GetExtrapoMethod() ==
		         ExtrapolationMethod::EXTRAPO_AVERAGE)
		{
			// average
			ESxx = avgESxx;
			ESyy = avgESyy;
			ESxy = avgESxy;
			ESzz = avgESzz;
			if (ele_dim == 3)
			{
				ESxz = avgESxz;
				ESyz = avgESyz;
			}
		}

		// Average value of the contribution of all neighbor elements
		ESxx /= dbuff[i];
		ESyy /= dbuff[i];
		ESxy /= dbuff[i];
		ESzz /= dbuff[i];

		ESxx += pcs->GetNodeValue(nodes[i], Idx_Strain[0]);
		ESyy += pcs->GetNodeValue(nodes[i], Idx_Strain[1]);
		ESzz += pcs->GetNodeValue(nodes[i], Idx_Strain[2]);
		ESxy += pcs->GetNodeValue(nodes[i], Idx_Strain[3]);

		pcs->SetNodeValue(nodes[i], Idx_Strain[0], ESxx);
		pcs->SetNodeValue(nodes[i], Idx_Strain[1], ESyy);
		pcs->SetNodeValue(nodes[i], Idx_Strain[2], ESzz);
		pcs->SetNodeValue(nodes[i], Idx_Strain[3], ESxy);

		if (ele_dim == 3)
		{
			ESxz /= dbuff[i];
			ESyz /= dbuff[i];
			ESxz += pcs->GetNodeValue(nodes[i], Idx_Strain[4]);
			ESyz += pcs->GetNodeValue(nodes[i], Idx_Strain[5]);
			//
			pcs->SetNodeValue(nodes[i], Idx_Strain[4], ESxz);
			pcs->SetNodeValue(nodes[i], Idx_Strain[5], ESyz);
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ExtropolateGuassStress()
   Aufgabe:
           Extropolate the Gauss point strains to nodes
   Formalparameter:
           E:

   Programming:
   06/2004   WW
   03/2007   WW  Generize for all 2nd variables
 **************************************************************************/
void CFiniteElementVec::ExtropolateGuassStress()
{
	int i, j, gp_r, gp_s, gp_t;
	// int l1,l2,l3,l4; //, counter;
	double ESxx, ESyy, ESzz, ESxy, ESxz, ESyz, Pls;
	double avgESxx, avgESyy, avgESzz, avgESxy, avgESxz, avgESyz, avgPls;
	int i_s, i_e, ish, k = 0;
	MshElemType::type ElementType = MeshElement->GetElementType();
	long node_i = 0;
	// For strain and stress extropolation all element types
	// Number of elements associated to nodes
	for (i = 0; i < nnodes; i++)
		dbuff[i] =
		    (double)MeshElement->GetNode(i)->getConnectedElementIDs().size();
	//
	gp = gp_r = gp_s = gp_t = 0;
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];
	if (eleV_DM->pStrain)  // 08.02.2008 WW
		idx_pls = pcs->GetNodeValueIndex("STRAIN_PLS");
	//
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		if (ElementType == MshElemType::QUAD ||
		    ElementType == MshElemType::HEXAHEDRON)
		{
			SetGaussPoint(gp, gp_r, gp_s, gp_t);
			i = GetLocalIndex(gp_r, gp_s, gp_t);
			if (i == -1) continue;
		}
		else
			i = gp;

		Sxx[i] = (*eleV_DM->Stress)(0, gp);
		Syy[i] = (*eleV_DM->Stress)(1, gp);
		Szz[i] = (*eleV_DM->Stress)(2, gp);
		Sxy[i] = (*eleV_DM->Stress)(3, gp);
		if (eleV_DM->pStrain)
			pstr[i] = (*eleV_DM->pStrain)(gp);
		else
			pstr[i] = 0.0;  // 08.02.2008 WW
		if (ele_dim == 3)
		{
			Sxz[i] = (*eleV_DM->Stress)(4, gp);
			Syz[i] = (*eleV_DM->Stress)(5, gp);
		}
	}
	//
	if (ElementType == MshElemType::QUAD ||
	    ElementType == MshElemType::HEXAHEDRON)
		Xi_p = CalcXi_p();

	//
	i_s = 0;
	i_e = nnodes;
	ish = 0;
	if (ElementType == MshElemType::TETRAHEDRON)  // tet
	{
		i_s = 1;
		i_e = nnodes + 1;
		ish = 1;
	}
	//---------------------------------------------------------
	// Mapping Gauss point strains to nodes and update nodes
	// strains:
	//---------------------------------------------------------
	avgESxx = avgESyy = avgESzz = avgESxy = avgESxz = avgESyz = avgPls = 0.0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
	{
		// average
		avgESxx = CalcAverageGaussPointValues(Sxx);
		avgESyy = CalcAverageGaussPointValues(Syy);
		avgESzz = CalcAverageGaussPointValues(Szz);
		avgESxy = CalcAverageGaussPointValues(Sxy);
		avgESxz = CalcAverageGaussPointValues(Sxz);
		avgESyz = CalcAverageGaussPointValues(Syz);
		avgPls = CalcAverageGaussPointValues(pstr);
	}

	for (i = 0; i < nnodes; i++)
	{
		ESxx = ESyy = ESzz = ESxy = ESxz = ESyz = Pls = 0.0;

		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			//
			SetExtropoGaussPoints(i);
			//
			ComputeShapefct(1);  // Linear interpolation function
			//
			for (j = i_s; j < i_e; j++)
			{
				k = j - ish;
				ESxx += Sxx[j] * shapefct[k];
				ESyy += Syy[j] * shapefct[k];
				ESxy += Sxy[j] * shapefct[k];
				ESzz += Szz[j] * shapefct[k];
				Pls += pstr[j] * shapefct[k];
				if (ele_dim == 3)
				{
					ESxz += Sxz[j] * shapefct[k];
					ESyz += Syz[j] * shapefct[k];
				}
			}
		}
		else if (this->GetExtrapoMethod() ==
		         ExtrapolationMethod::EXTRAPO_AVERAGE)
		{
			// average
			ESxx = avgESxx;
			ESyy = avgESyy;
			ESxy = avgESxy;
			ESzz = avgESzz;
			Pls = avgPls;
			if (ele_dim == 3)
			{
				ESxz = avgESxz;
				ESyz = avgESyz;
			}
		}

		// Average value of the contribution of ell neighbor elements
		ESxx /= dbuff[i];
		ESyy /= dbuff[i];
		ESxy /= dbuff[i];
		ESzz /= dbuff[i];
		Pls /= dbuff[i];
		//
		node_i = nodes[i];
		ESxx += pcs->GetNodeValue(node_i, Idx_Stress[0]);
		ESyy += pcs->GetNodeValue(node_i, Idx_Stress[1]);
		ESzz += pcs->GetNodeValue(node_i, Idx_Stress[2]);
		ESxy += pcs->GetNodeValue(node_i, Idx_Stress[3]);
		if (eleV_DM->pStrain)  // 08.02.2008 WW
			Pls += pcs->GetNodeValue(node_i, idx_pls);

		pcs->SetNodeValue(node_i, Idx_Stress[0], ESxx);
		pcs->SetNodeValue(node_i, Idx_Stress[1], ESyy);
		pcs->SetNodeValue(node_i, Idx_Stress[2], ESzz);
		pcs->SetNodeValue(node_i, Idx_Stress[3], ESxy);
		if (eleV_DM->pStrain)  // 08.02.2008 WW
			pcs->SetNodeValue(node_i, idx_pls, fabs(Pls));

		if (ele_dim == 3)
		{
			ESxz /= dbuff[i];
			ESyz /= dbuff[i];

			ESxz += pcs->GetNodeValue(node_i, Idx_Stress[4]);
			ESyz += pcs->GetNodeValue(node_i, Idx_Stress[5]);

			pcs->SetNodeValue(node_i, Idx_Stress[4], ESxz);
			pcs->SetNodeValue(node_i, Idx_Stress[5], ESyz);
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
      CFiniteElementVec::ComputePrincipleStresses(double *dEStress)
      (2D only)
   Aufgabe:

   Formalparameter:
         E:
   const double *Stresses: Stresses

   Return: Angle of maxium principle stress component to x direction
   sig2<sig1

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
double CFiniteElementVec::ComputePrincipleStresses(const double* Stresses, double* pr_stress)
{
	double prin_ang, var;
	// Angle of the principle plane
	if (fabs(Stresses[0] - Stresses[1]) < MKleinsteZahl) return 0.0;
	prin_ang = atan(2.0 * Stresses[3] / (Stresses[0] - Stresses[1]));
	// Principle stress 1
	pr_stress[0] = 0.5 * (Stresses[0] + Stresses[1]) +
	               0.5 * (Stresses[0] - Stresses[1]) * cos(prin_ang) +
	               Stresses[3] * sin(prin_ang);
	// Principle stress 2
	pr_stress[1] = 0.5 * (Stresses[0] + Stresses[1]) -
	               0.5 * (Stresses[0] - Stresses[1]) * cos(prin_ang) -
	               Stresses[3] * sin(prin_ang);
	pr_stress[2] = Stresses[2];

	prin_ang *= 0.5;

	// Jump direction
	if (pr_stress[1] >= pr_stress[0])
	{
		prin_ang += 0.5 * pai;
		var = pr_stress[0];
		pr_stress[0] = pr_stress[1];
		pr_stress[1] = var;
	}
	// if(pr_stress[0]<0.0) prin_ang += pai;
	return prin_ang;
}

/***************************************************************************
   GeoSys - Funktion:
            CFiniteElementVec:: CalcStrain_v()
   Aufgabe:
           Calculate effictive strain at Gauss points
   Formalparameter:
           E:

   Programming:
   01/2009   WW/UWG
 **************************************************************************/
double CFiniteElementVec::CalcStrain_v()
{
	int i;
	for (i = 0; i < ns; i++)
		dstrain[i] = 0.0;
	//
	for (i = 0; i < nnodesHQ; i++)
		dstrain[i] +=
		    pcs->GetNodeValue(nodes[i], Idx_Strain[0]) * shapefctHQ[i];
	double val = 0;
	for (i = 0; i < 3; i++)
		val += dstrain[i] * dstrain[i];
	for (i = 3; i < ns; i++)
		val += 0.5 * dstrain[i] * dstrain[i];

	return sqrt(2.0 * val / 3.);
}
/***************************************************************************
   GeoSys - Funktion:
            CFiniteElementVec:: CalcStress_eff()
   Aufgabe:
           Calculate effective stress at Gauss points
   Formalparameter:
           E:

   Programming:
   06/2014   AJ
 **************************************************************************/
double CFiniteElementVec::CalcStress_eff()
{
	size_t j;
	double val = 0.;
	ElementValue_DM* ele_value = ele_value_dm[MeshElement->GetIndex()];

	for (j = 0; j < dim; j++)
		val -= (*ele_value->Stress)(j, gp);

	val /= dim;

	return val;
}
}
// end namespace FiniteElement
