/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_std.h"

#include <cfloat>

#include "Configure.h"
#include "memory.h"

#ifdef NEW_EQS
#include "equation_class.h"
#endif
#include "mathlib.h"
#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#include "ElementMatrix.h"
#include "ElementValue.h"
#include "eos.h"
#include "fem_ele_vec.h"
#include "mechanics_utils.h"
#include "rfmat_cp.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_pcs_dm.h"
#include "SparseMatrixDOK.h"


extern double gravity_constant;     // TEST, must be put in input file
#define COMP_MOL_MASS_AIR 28.96     // kg/kmol WW  28.96
#define COMP_MOL_MASS_WATER 18.016  // WW 18.016
#define GAS_CONSTANT 8314.41        // J/(kmol*K) WW
#define GAS_CONSTANT_V 461.5        // WW
#define T_KILVIN_ZERO 273.15        // WW

using namespace std;

namespace FiniteElement
{
//========================================================================
// Element calculation
//========================================================================

EnumProcessType getEnumProcessType(CRFProcess* pcs)
{
	EnumProcessType PcsType = EnumProcessType_INVALID;
	char pcsT = (convertProcessTypeToString(pcs->getProcessType()))[0];
	if (pcs->getProcessType() == AIR_FLOW)  // OK
		pcsT = 'A';
	else if (pcs->getProcessType() == MULTI_PHASE_FLOW)
		pcsT = 'V';  // Non-isothermal multiphase flow
	else if (pcs->getProcessType() == DEFORMATION_H2)  // 09.08.2010 WW
		pcsT = 'V';  // Non-isothermal multiphase flow
	else if (pcs->getProcessType() == DEFORMATION_FLOW)  // NW
		pcsT = 'L';

	switch (pcsT)
	{
		default:
			break;
		case 'L':  // Liquid flow
			PcsType = L;
			break;
		case 'U':  // Unconfined flow
			PcsType = U;
			break;
		case 'G':  // Groundwater flow
			PcsType = G;
			break;
		case 'T':  // Two-phase flow
			PcsType = T;
			break;
		case 'C':  // Componental flow
			PcsType = C;
			break;
		case 'H':  // heat transport
			PcsType = H;
			break;
		case 'M':  // Mass transport
			PcsType = M;
			break;
		case 'R':  // OK4104 Richards flow
			PcsType = R;
			break;
		case 'A':  // Air (gas) flow
			PcsType = A;
			break;
		case 'F':  // Fluid Momentum Process
			PcsType =
			    R;  // R should include L if the eqn of R is written right.
			break;
		case 'V':  // 24.02.2007 WW
			PcsType = V;
			break;
		case 'P':  // 04.03.2009 PCH
			PcsType = P;
			break;
	}

	if (pcs->getProcessType() == FiniteElement::TH_MONOLITHIC) PcsType = TH;
	return PcsType;
}

/**************************************************************************
   GeoSys - Function: Constructor
   Programmaenderungen:
   01/2005   WW    Erste Version
**************************************************************************/
CFiniteElementStd::CFiniteElementStd(CRFProcess* Pcs, const int C_Sys_Flad,
                                     const int order,
                                     const EnumProcessType pcsType)
    : CElement(C_Sys_Flad, order),
      phase(0),
      comp(0),
      SolidProp(NULL),
      FluidProp(NULL),
      MediaProp(NULL),
      pcs(Pcs),
      dm_pcs(NULL),
      HEAD_Flag(false)
{
	int i;
	int size_m = 64;  // 25.2.2007
	string name2;
	char name1[MAX_ZEILE];
	cpl_pcs = NULL;
	// 27.2.2007 WW
	newton_raphson = false;
	// WW
	if (pcs->m_num->nls_method_name.compare("NEWTON_RAPHSON") == 0)
		newton_raphson = true;
	Mass = NULL;
	Mass2 = NULL;
	Laplace = NULL;
	Advection = NULL;
	Storage = NULL;
	Content = NULL;
	StrainCoupling = NULL;
	RHS = NULL;
	FCT_MassL = NULL;  // NW

	//
	edlluse = edttuse = NULL;
	idx_vel_disp = NULL;  // WW
	weight_func = NULL;   // WW
	idx_vel = new int[3];
	NodalVal1 = new double[size_m];
	NodalVal2 = new double[size_m];
	NodalVal3 = new double[size_m];
	NodalVal4 = new double[size_m];
	NodalValC = new double[size_m];
	NodalValC1 = new double[size_m];
	NodalVal_Sat = new double[size_m];
	NodalVal_SatNW = new double[size_m];
	NodalVal_p2 = new double[size_m];
	NodalVal_p20 = new double[size_m];  // AKS
	NodalVal_T0 = new double[size_m];   // AKS/NB
	NodalVal_T1 = new double[size_m];   // AKS/NB
	NodalVal_p0 = new double[size_m];   // AKS/NB
	NodalVal_p1 = new double[size_m];   // AKS/NB
	mat = new double[9];                // WW
#if 0
	fem1 = NULL;
	fem2 = NULL;
#endif
	// NW
	switch (C_Sys_Flad / 10)
	{
		case 1:
			weight_func = new double[2];
			break;
		case 2:
			weight_func = new double[4];
			break;
		case 3:
			weight_func = new double[8];
			break;
	}
	//
	// 27.2.2007. GravityMatrix = NULL;
	eqs_rhs = NULL;  // 08.2006 WW
	//
	// 12.12.2007 WW
	for (i = 0; i < 4; i++)
		NodeShift[i] = 0;
	//
	idx_vel_disp = new int[3];

	dm_pcs = NULL;
	heat_phase_change = false;

	idx_vel_disp[0] = idx_vel_disp[1] = idx_vel_disp[2] = -1;
	idx_pres = -1;

	idxS = idx3 = -1;
	if (pcs->primary_variable_name.compare("HEAD") == 0) HEAD_Flag = true;
	// SB4218 added
	string pcs_primary = pcs->pcs_primary_function_name[0];
	if (pcs_primary.compare("HEAD") == 0) HEAD_Flag = true;
	for (i = 0; i < 9; i++)
		mat[i] = 0.0;

	idx0 = idx1 = 0;  // column index in the node value data
	LocalShift = 0;

	if (pcsType == EnumProcessType_INVALID)
		PcsType = getEnumProcessType(pcs);
	else
		PcsType = pcsType;

	switch (PcsType)
	{
		case L:  // Liquid flow
			idx0 = pcs->GetNodeValueIndex("PRESSURE1");
			idx1 = idx0 + 1;
			idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
			idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
			idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
			break;
		case G:  // Groundwater flow
			idx0 = pcs->GetNodeValueIndex("HEAD");
			idx1 = idx0 + 1;
			// WW
			idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
			// WW
			idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
			// WW
			idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
			break;
		case H:  // heat transport
			idx0 = pcs->GetNodeValueIndex("TEMPERATURE1");
			idx1 = idx0 + 1;
			break;
		case M:  // Mass transport
			sprintf(name1, "%s", pcs->pcs_primary_function_name[0]);
			name2 = name1;
			idx0 = pcs->GetNodeValueIndex(name2);
			idx1 = idx0 + 1;
			break;
		case R:  // OK4104 Richards flow
			// 02.2.2007 GravityMatrix = new  SymMatrix(size_m);
			idx0 = pcs->GetNodeValueIndex("PRESSURE1");
			idx1 = idx0 + 1;
			idxS = pcs->GetNodeValueIndex("SATURATION1") + 1;
			idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
			idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
			idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
			if ((int)pcs->dof > 1)  // Dual porosity model. WW
			{
				idxp20 = pcs->GetNodeValueIndex("PRESSURE2");
				idxp21 = idxp20 + 1;
				// WW
				Advection = new Matrix(size_m, size_m);
				// 12.12.2007 WW
				for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
					NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
			}
			break;
		case A:  // Air (gas) flow
			idx0 = pcs->GetNodeValueIndex("PRESSURE1");
			idx1 = idx0 + 1;  // OK
			break;
		case V:  // 24.02.2007 WW
			     // // 02.2.2007 GravityMatrix = new  SymMatrix(size_m);
			// 12.12.2007 WW
			for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
				NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
			//
			idx0 = pcs->GetNodeValueIndex("PRESSURE1");
			idx1 = idx0 + 1;
			idxp20 = pcs->GetNodeValueIndex("PRESSURE2");
			idxp21 = idxp20 + 1;
			idxS = pcs->GetNodeValueIndex("SATURATION1") + 1;
			idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
			idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
			idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
			size_m = 40;
			break;
		case P:  // 04.03.2009 PCH
			for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
				NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
			//
			idx0 = pcs->GetNodeValueIndex("PRESSURE1");
			idx1 = idx0 + 1;
			idxSn0 = pcs->GetNodeValueIndex("SATURATION2");
			idxSn1 = idxSn0 + 1;
			idxS = pcs->GetNodeValueIndex("SATURATION1") + 1;
			idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
			idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
			idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
			size_m = 40;
			break;
		case TH:
			for (i = 0; i < pcs->pcs_number_of_primary_nvals; i++)
				NodeShift[i] = i * pcs->m_msh->GetNodesNumber(false);
#if 0
		fem1 = new CFiniteElementStd(pcs,C_Sys_Flad,order,L);
		fem2 = new CFiniteElementStd(pcs,C_Sys_Flad,order,H);
#else
			idxp0 = pcs->GetNodeValueIndex("PRESSURE1");
			idxp1 = idxp0 + 1;
			idxT0 = pcs->GetNodeValueIndex("TEMPERATURE1");
			idxT1 = idxT0 + 1;
			idx_vel[0] = pcs->GetNodeValueIndex("VELOCITY_X1");
			idx_vel[1] = pcs->GetNodeValueIndex("VELOCITY_Y1");
			idx_vel[2] = pcs->GetNodeValueIndex("VELOCITY_Z1");
#endif
			size_m = 40;
			break;
		case U:  // Unconfined flow
		case T:  // Two-phase flow
		case C:  // Componental flow
		case F:  // Fluid Momentum Process
			break;
		default:
#if 0
		//WW GravityMatrix = new  SymMatrix(size_m);
		{
			idx0 = pcs->GetNodeValueIndex("PRESSURE1");
			idx1 = idx0 + 1;
		}
#endif
			break;
	}

	if (pcs->Memory_Type == 0)  // Do not store local matrices
	{
		// 04.03.2009 PCH
		if (PcsType == V || PcsType == P)
			Mass2 = new Matrix(size_m, size_m);
		else
			Mass = new Matrix(size_m, size_m);
		Laplace = new Matrix(size_m, size_m);
		if (PcsType == H || PcsType == M || PcsType == A)
		{
			Advection = new Matrix(size_m, size_m);
			Storage = new Matrix(size_m, size_m);
			Content = new Matrix(size_m, size_m);
		}
		else if (PcsType == TH)
		{
			Advection = new Matrix(size_m, size_m);
		}
		if (D_Flag) StrainCoupling = new Matrix(size_m, 60);
		RHS = new Vector(size_m);
	}
	//
	StiffMatrix = new Matrix(size_m, size_m);
	AuxMatrix = new Matrix(size_m, size_m);
	AuxMatrix1 = new Matrix(size_m, size_m);

	if (this->pcs->m_num->fct_method > 0)  // NW
		FCT_MassL = new DiagonalMatrix(size_m);
	time_unit_factor = pcs->time_unit_factor;

	check_matrices = true;
	//
	SolidProp1 = NULL;
	MediaProp1 = NULL;
	flag_cpl_pcs = false;  // OK
	// size_m changed
	NodalVal = new double[size_m];
	NodalVal0 = new double[size_m];

#if defined(USE_PETSC)       // || defined(other parallel libs)//03~04.3012. WW
	row_ids = new int[size_m];  //> global indices of local matrix rows
	col_ids = new int[size_m];  //> global indices of local matrix columns
#endif
	add2global = false;
	dof_index = 0;
	drho_gw_dT = .0;
	dSdp = .0;
	GasProp = NULL;
	index = 0;
	M_g = 0;
	mfp_pcs = NULL;
	p_gw = .0;
	PG = PG0 = PG2 = PG20 = 0;
	poro = 0;
	rho_g = rho_ga = rho_gw = rhow = 0;
	Sw = 0;
	TG = TG0 = 0;
	tort = 0;
}

/**************************************************************************
   GeoSys - Function: Destructor
   Programmaenderungen:
   01/2005   WW    Erste Version
**************************************************************************/
// Destructor
CFiniteElementStd::~CFiniteElementStd()
{
	//  02.2.2007 if(GravityMatrix) delete GravityMatrix;
	// 02.2.2007  GravityMatrix = NULL;

	if (pcs->Memory_Type == 0)  // Do not store local matrices
	{
		if (Mass) delete Mass;
		if (Mass2) delete Mass2;
		if (Laplace) delete Laplace;
		if (Advection) delete Advection;
		if (Storage) delete Storage;
		if (Content) delete Content;
		if (StrainCoupling) delete StrainCoupling;
		if (RHS) delete RHS;
		if (FCT_MassL) delete FCT_MassL;

		Mass = NULL;
		Laplace = NULL;
		Advection = NULL;
		Storage = NULL;
		Content = NULL;
		StrainCoupling = NULL;
		RHS = NULL;
		FCT_MassL = NULL;
	}

	delete StiffMatrix;
	delete AuxMatrix;
	delete AuxMatrix1;
	if (edlluse) delete[] edlluse;
	if (edttuse) delete[] edttuse;

	StiffMatrix = NULL;
	AuxMatrix = NULL;
	AuxMatrix1 = NULL;
	// 27.2.2007 WW
	delete[] NodalVal;
	delete[] NodalVal0;
	delete[] NodalVal1;
	delete[] NodalVal2;
	delete[] NodalVal3;
	delete[] NodalVal4;
	delete[] NodalValC;
	delete[] NodalValC1;
	delete[] NodalVal_Sat;
	delete[] NodalVal_SatNW;
	delete[] NodalVal_p2;
	delete[] mat;
	if (idx_vel_disp) delete[] idx_vel_disp;
	delete[] idx_vel;
	delete[] NodalVal_p20;  // AKS
	delete[] NodalVal_T0;
	delete[] NodalVal_T1;
	delete[] NodalVal_p0;
	delete[] NodalVal_p1;
	// NW
	if (weight_func) delete[] weight_func;  // Remove bug. WW
	weight_func = NULL;
}
/**************************************************************************
   GeoSys - Function: SetMemory

   Aufgabe:
         Set memory for local matrices
   Programmaenderungen:
   01/2005   WW    Erste Version

**************************************************************************/
void CFiniteElementStd::SetMemory()
{
	int Size = nnodes;
	if (PcsType == V || PcsType == P || PcsType == TH)  // 4.3.2009 PCH
		Size *= 2;
	ElementMatrix* EleMat = NULL;
	// Prepare local matrices
	// If local matrices are not stored, resize the matrix
	if (pcs->Memory_Type == 0)
	{
		if (PcsType == V || PcsType == P)  // 04.3.2009 PCH
			Mass2->LimitSize(Size, Size);
		else
			Mass->LimitSize(Size, Size);  // Mass->LimitSize(nnodes); //
		                                  // unsymmetric in case of Upwinding
		Laplace->LimitSize(Size, Size);
		if (PcsType == H || PcsType == M || PcsType == A)
		{
			Advection->LimitSize(Size, Size);  // SB4200
			Storage->LimitSize(Size, Size);    // SB4200
			Content->LimitSize(Size, Size);    // SB4209
		}
		else if (PcsType == TH)
		{
			Advection->LimitSize(Size, Size);
		}

		if (PcsType == R && pcs->type == 22)  // dual-porosity. WW
			Advection->LimitSize(Size, Size);
		if (D_Flag > 0) StrainCoupling->LimitSize(Size, dim * nnodesHQ);
		RHS->LimitSize(Size);
	}
	else
	{
		EleMat = pcs->Ele_Matrices[Index];
		// if(PcsType==V) //24.2.2007 WW
		// Mass2 = EleMat->GetMass2();
		Mass = EleMat->GetMass();
		Laplace = EleMat->GetLaplace();
		// Advection, Storage, Content SB4200
		if (PcsType == M || PcsType == H)
		{
			Advection = EleMat->GetAdvection();
			Storage = EleMat->GetStorage();
			Content = EleMat->GetContent();
		}
		RHS = EleMat->GetRHS();
		if (D_Flag > 0) StrainCoupling = EleMat->GetCouplingMatrixB();
		if (D_Flag == 41) LocalShift = dim * nnodesHQ;
	}

	// 25.2.2007.WW if(GravityMatrix) GravityMatrix->LimitSize(nnodes);

	StiffMatrix->LimitSize(Size, Size);
	AuxMatrix->LimitSize(Size, Size);
	AuxMatrix1->LimitSize(Size, Size);
	if (this->pcs->m_num->fct_method > 0)  // NW
		FCT_MassL->LimitSize(Size);
}

/**************************************************************************
   GeoSys - Function: ConfigureCoupling

   Aufgabe:
         Set coupling information for local fem calculation
   Programmaenderungen:
   01/2005   WW    Erste Version
   02/2007   WW    Multi phase flow
    03/2009   PCH	 PS_GLOBAL

**************************************************************************/
void CFiniteElementStd::ConfigureCoupling(CRFProcess* pcs, const int* Shift)
{
	// ProcessType pcs_type (pcs->getProcessType());
	EnumProcessType CplPcsType = getEnumProcessType(pcs);

	if (D_Flag > 0)
	{
		for (size_t i = 0; i < pcs_vector.size(); i++)
			if (isDeformationProcess(pcs_vector[i]->getProcessType()))
			{
				dm_pcs = (CRFProcessDeformation*)pcs_vector[i];
				break;
			}
		Idx_dm0[0] = dm_pcs->GetNodeValueIndex("DISPLACEMENT_X1");
		Idx_dm0[1] = dm_pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
		Idx_dm1[0] = Idx_dm0[0] + 1;
		Idx_dm1[1] = Idx_dm0[1] + 1;
		//     if(problem_dimension_dm==3)
		if (dim == 3)
		{
			Idx_dm0[2] = dm_pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
			Idx_dm1[2] = Idx_dm0[2] + 1;
		}
		if (dm_pcs->type / 10 == 4)
			for (size_t i = 0; i < pcs->GetPrimaryVNumber(); i++)
				NodeShift[i] = Shift[i];
	}

	switch (CplPcsType)
	{
		default:
			if (T_Flag)
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case L:  // Liquid flow
			if (T_Flag)
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case U:  // Unconfined flow
			break;
		case G:  // Groundwater flow
			if (T_Flag)
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case T:  // Two-phase flow
			if (pcs->pcs_type_number == 0)
			{
				cpl_pcs = pcs_vector[pcs->pcs_number + 1];
				idx_c0 = cpl_pcs->GetNodeValueIndex("SATURATION2");
				idx_c1 = idx_c0 + 1;
			}
			else if (pcs->pcs_type_number == 1)
			{
				cpl_pcs = pcs_vector[pcs->pcs_number - 1];
				idx_c0 = cpl_pcs->GetNodeValueIndex("PRESSURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
#if 0
	case S:                             // Multi-phase flow. 24.2.2007 WW
		if(C_Flag)                //if(PCSGet("HEAT_TRANSPORT"))
		{
			cpl_pcs = PCSGet("MASS_TRANSPORT");
			idx_c0 = cpl_pcs->GetProcessComponentNumber();
			idx_c1 = idx_c0 + 1;
		}
		break;
#endif
		case C:  // Componental flow
			break;
		case H:  // heat transport
			// SB CMCD this needs to be fixed
			cpl_pcs = PCSGet("GROUNDWATER_FLOW");
			if (cpl_pcs)  // WW
			{
				idx_c0 = cpl_pcs->GetNodeValueIndex("HEAD");
				idx_c1 = idx_c0 + 1;
			}
			else
			{
				cpl_pcs = PCSGet("LIQUID_FLOW");
				if (cpl_pcs == NULL)
				{
					// OK
					cpl_pcs = PCSGet("RICHARDS_FLOW");
					if (cpl_pcs)
						// WW
						idxS = cpl_pcs->GetNodeValueIndex("SATURATION1") + 1;
				}
				if (cpl_pcs == NULL)
				{
					// 24.042.2004 WW
					cpl_pcs = PCSGet("MULTI_PHASE_FLOW");
					if (cpl_pcs)
						// WW
						idxS = cpl_pcs->GetNodeValueIndex("SATURATION1") + 1;
				}
				if (cpl_pcs == NULL)  // 23.02.2009 NB 4.9.05
				{
					cpl_pcs = PCSGet("TWO_PHASE_FLOW");
					if (cpl_pcs)
						idxS = cpl_pcs->GetNodeValueIndex("SATURATION1") + 1;
				}
				if (cpl_pcs == NULL)  // 23.02.2009 NB 4.9.05

					cpl_pcs = PCSGet("AIR_FLOW");  // 23.01.2009 NB

				if (cpl_pcs)  // MX
				{
					idx_c0 = cpl_pcs->GetNodeValueIndex("PRESSURE1");
					idx_c1 = idx_c0 + 1;
				}
			}
			break;
		case M:  // Mass transport
			if (C_Flag && T_Flag)
			{
				if (cpl_pcs == NULL)
				{
					cpl_pcs = PCSGetFlow();
					if (cpl_pcs)
						idx_c0 = cpl_pcs->GetNodeValueIndex("PRESSURE1");
					idx_c1 = idx_c0 + 1;
				}
			}  // Mass transport
			else if (T_Flag)
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case R:          // Richards flow
			if (T_Flag)  // if(PCSGet("HEAT_TRANSPORT"))
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case V:          // Multi-phase flow. 24.2.2007 WW
			if (T_Flag)  // if(PCSGet("HEAT_TRANSPORT"))
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case A:          // Gas flow
			if (T_Flag)  // NB 23.01.2009 4.9.05
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
		case P:  // Multi-phase flow. 03.03.2009 PCH
			if (T_Flag)
			{
				cpl_pcs = PCSGet("HEAT_TRANSPORT");
				idx_c0 = cpl_pcs->GetNodeValueIndex("TEMPERATURE1");
				idx_c1 = idx_c0 + 1;
			}
			break;
	}
}

/*************************************************************************
   FEMLib-Function:
   Task: Set material pointers to the current element
   01/2005 WW Implementation
   03/2005 OK MultiMSH
   11/2005 YD Set cursor of gas
   06/2009 OK MMP test not here (time consuming)
   01/2010 NW Set geo_area here
 **************************************************************************/
void CFiniteElementStd::SetMaterial(int /*phase*/)
{
	//	phase = 0;
	//----------------------------------------------------------------------
	// MMP
	int mmp_index = 0;
	long group = MeshElement->GetPatchIndex();
	mmp_index = group;
	// Single continua thermal:
	if (msp_vector.size() > 0)
	{
		// TODO: For multidomain meshes accessing msp_vector can result in SEGV.
		// Either use checked \c at() access or ensure
		// mmp_index < msp_vector.size().
		SolidProp = msp_vector[mmp_index];
		// KR not used: SolidProp->m_pcs = pcs;   //NW
		SolidProp->Fem_Ele_Std = this;  // CMCD for Decovalex
	}

	if (pcs->type == 22)  // WW/YD
	{
		if (pcs->GetContinnumType() == 0)  // Matrix //WW
			mmp_index = 2 * group;
		else  // fracture //WW
			mmp_index = 2 * group + 1;
	}
	// TODO: For multidomain meshes accessing mmp_vector can result in SEGV,
	// like above msp_vector[].
	MediaProp = mmp_vector[mmp_index];
	MediaProp->m_pcs = pcs;
	MediaProp->Fem_Ele_Std = this;
	MeshElement->SetFluxArea(MediaProp->geo_area);  // NW

	if (MediaProp->storage_model == 7)  // 29.11.2011. WW
		SolidProp->Calculate_Lame_Constant();

	//----------------------------------------------------------------------
	// MSP
	// If dual thermal:
	/*
	   if(msp_vector.size()>0)
	   {
	   SolidProp = msp_vector[mmp_index];
	   SolidProp->Fem_Ele_Std = this;//CMCD for Decovalex
	   }
	 */
	if (pcs->type == 22)  // WW
	{
		if (pcs->GetContinnumType() == 0)  // Matrix //WW
			mmp_index = 2 * group + 1;
		else  // fracture //WW
			mmp_index = 2 * group;
		MediaProp1 = mmp_vector[mmp_index];
		MediaProp1->m_pcs = pcs;
		MediaProp1->Fem_Ele_Std = this;
		//----------------------------------------------------------------------
		// MSP
		// If dual thermal:
		/*
		   if(msp_vector.size()>0)
		   {
		   SolidProp1 = msp_vector[mmp_index];
		   SolidProp1->Fem_Ele_Std = this;//CMCD for Decovalex
		   }
		 */
	}
	//----------------------------------------------------------------------
	// MFP
	/* Comment out - NW
	   if(PCSGet("LIQUID_FLOW")){
	    FluidProp = MFPGet("LIQUID");
	    if(!FluidProp)
	      cout << "Warning: LIQUID was not found in fluid properties." << endl;
	   }
	 */
	if (mfp_vector.size() > 0)
	{
		FluidProp = mfp_vector[0];
		FluidProp->Fem_Ele_Std = this;
	}
#ifndef OGS_ONLY_TH
	// 03.2009 PCH
	if ((PCSGet("RICHARDS_FLOW") && PCSGet("HEAT_TRANSPORT")) ||
	    pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || pcs->getProcessType() == FiniteElement::PS_GLOBAL || pcs->type == 42)
	{
		FluidProp = MFPGet("LIQUID");
		FluidProp->Fem_Ele_Std = this;
		// FluidProp = mfp_vector[0];
		GasProp = MFPGet("GAS");
		if (GasProp) GasProp->Fem_Ele_Std = this;
	}
#endif
	//	if (SolidProp && FluidProp) //NW
	//		if(D_Flag > 0  && FluidProp->Density() > MKleinsteZahl)
	//			SolidProp->Calculate_Lame_Constant();
	//----------------------------------------------------------------------
	// MCP
	//----------------------------------------------------------------------
}


/**************************************************************************
   FEMLib-Method:
   Task: Calculate nodal enthalpy
   Programming: WW 09/2005
**************************************************************************/
void CFiniteElementStd::CalNodalEnthalpy()
{
	int i;
	double temp, dT;
	for (i = 0; i < nnodes; i++)
	{
		heat_phase_change = SolidProp->CheckTemperature_in_PhaseChange(
		    NodalVal0[i], NodalVal1[i]);
		if (heat_phase_change) break;
	}
	if (!heat_phase_change) return;
	// Calculate node enthalpy
	for (i = 0; i < nnodes; i++)
	{
		NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS);
		SetCenterGP();
		temp = FluidProp->Density() *
		       MediaProp->Porosity(Index, pcs->m_num->ls_theta) *
		       NodalVal_Sat[i];
		// Enthalpy
		dT = 0.0;
		NodalVal2[i] = SolidProp->Enthalpy(NodalVal0[i], temp);
		if (fabs(NodalVal1[i] - NodalVal0[i]) < 1.0e-8) dT = 1.0e-8;
		NodalVal3[i] = SolidProp->Enthalpy(NodalVal1[i] + dT, temp);
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   08/2005 OK Gas flow
   10/2005 YD/OK: general concept for heat capacity
   11/2005 CMCD Heat capacity function included in mmp
   01/2007 OK Two-phase flow
**************************************************************************/
double CFiniteElementStd::CalCoefMass()
{
	int Index = MeshElement->GetIndex();
	double val = 0.0;
	double humi = 1.0;
	double rhov = 0.0;
	double arg[2];
	double biot_val, poro_val = 0.0, rho_val = 0.0, Se;
	int tr_phase = 0;         // SB, BG
	double saturation = 0.0;  // SB, BG
	CompProperties* m_cp = NULL;
	double drho_dp_rho = 0.0;  // drho/dp * 1/rho

	if (pcs->m_num->ele_mass_lumping) ComputeShapefct(1);
	switch (PcsType)
	{
		default:
			std::cout << "Fatal error in CalCoefMass: No valid PCS type"
			          << "\n";
			break;
		case L:  // Liquid flow
			// Is this really needed?
			val = MediaProp->StorageFunction(Index, unit, pcs->m_num->ls_theta);

			// get drho/dp/rho from material model or direct input
			if (FluidProp->compressibility_model_pressure > 0)
			{
				rho_val = FluidProp->Density();
				arg[0] = interpolate(NodalVal1);   //   p
				arg[1] = interpolate(NodalValC1);  //   T
				drho_dp_rho = FluidProp->drhodP(arg) / rho_val;
			}
			else
				drho_dp_rho = FluidProp->drho_dp;

			// JT 2010, needed storage term and fluid compressibility...
			// We derive here the storage at constant strain, or the inverse of
			// Biot's "M" coefficient
			// Assumptions are the most general possible::  Invarience under
			// "pi" (Detournay & Cheng) loading.
			// Se = 1/M = poro/Kf + (alpha-poro)/Ks    ::    Cf = 1/Kf = 1/rho *
			// drho/dp    ::    alpha = 1 - K/Ks
			// Second term (of Se) below vanishes for incompressible grains
			// WW if(D_Flag > 0  && rho_val > MKleinsteZahl)
			if (dm_pcs && MediaProp->storage_model == 7)
			{
				biot_val = SolidProp->biot_const;
				poro_val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
				val = 0.;

				if (fabs(SolidProp->K) < DBL_MIN)  // WW 29.09.2011
				{
					if (SolidProp->Youngs_mode < 10 ||
					    SolidProp->Youngs_mode > 13)  // JM,WX: 2013
						SolidProp->K = SolidProp->E / 3 /
						               (1 - 2 * SolidProp->PoissonRatio);
					else
					{
						double E_av;   // average Youngs modulus
						double nu_av;  // average Poisson ratio
						double nu_ai;  // Poisson ratio perpendicular to the
						               // plane of isotropie, due to strain in
						               // the plane of isotropie
						double nu_ia;  // Poisson ratio in the plane of
						               // isotropie, due to strain perpendicular
						               // to the plane of isotropie
						double nu_i;   // Poisson ratio in the plane of isotropy

						E_av = 2. / 3. * (*SolidProp->data_Youngs)(0) +
						       1. / 3. * (*SolidProp->data_Youngs)(1);

						nu_ia = (*SolidProp->data_Youngs)(2);
						nu_ai =
						    nu_ia * (*SolidProp->data_Youngs)(1) /
						    (*SolidProp->data_Youngs)(0);  //  nu_ai=nu_ia*Ea/Ei

						nu_i = SolidProp->Poisson_Ratio();
						//           12     13    21   23   31    32
						//           ai     ai    ia   ii   ia    ii
						nu_av = 1. / 3. * (nu_ai + nu_ia + nu_i);

						SolidProp->K = E_av / 3 / (1 - 2 * nu_av);
					}
				}
				val += poro_val * drho_dp_rho +
				       (biot_val - poro_val) * (1.0 - biot_val) / SolidProp->K;
				// Will handle the dual porosity version later...
			}
			else if (MediaProp->storage_model == 9)
			{
				val = 0.0;
				poro_val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
				rho_val = FluidProp->Density();
				FluidProp->compressibility_model_pressure = 0;
				drho_dp_rho = FluidProp->drhodP(NULL) / rho_val;
				val += poro_val * drho_dp_rho;
				// if (fabs(SolidProp->K) < DBL_MIN)
				// 	SolidProp->Calculate_Lame_Constant();
				// val += (SolidProp->biot_const - poro_val)*(1.0 - biot_val) / SolidProp->K;
			}
			else if (MediaProp->storage_model != 1)
			{
				poro_val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
				val += poro_val * drho_dp_rho;
			}

			// AS,WX: 08.2012 storage function eff stress
			if (MediaProp->storage_effstress_model > 0)
			{
				double storage_effstress = 1.;
				CFiniteElementStd* h_fem;
				h_fem = this;
				storage_effstress =
				    MediaProp->StorageFunctionEffStress(Index, nnodes, h_fem);
				val *= storage_effstress;
			}

			if (M_Process && pcs->use_total_stress_coupling && !pcs->disable_deformation_coupling_in_flow)
			{
				if (fabs(SolidProp->K) < DBL_MIN)
					SolidProp->Calculate_Lame_Constant();
				val += std::pow(SolidProp->biot_const, 2) / SolidProp->K;
			}

			val /= time_unit_factor;
			break;
		case U:  // Unconfined flow
			break;
		case G:  // MB now Groundwater flow
			if (MediaProp->unconfined_flow_group > 0)  // OK
				val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
			else
				val = MediaProp->StorageFunction(Index, unit,
				                                 pcs->m_num->ls_theta);
			break;
		case T:  // Two-phase flow
			// val = (1/rho*n*d_rho/d_p*S + Se*S )
			if (pcs->pcs_type_number == 0)
			{
				// PCH cpl_pcs gives a funny process number.
				// It is just the opposite of the phase. So, I get the value the
				// other way around.
				idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");
				for (int i = 0; i < nnodes; i++)
					NodalVal_Sat[i] = cpl_pcs->GetNodeValue(nodes[i], idxS + 1);
				Sw = 1.0 - interpolate(NodalVal_Sat);
				// Is this really needed?
				val = MediaProp->StorageFunction(Index, unit,
				                                 pcs->m_num->ls_theta) *
				      MMax(0., Sw);

				// JT 2010, generalized poroelastic storage. See single phase
				// version in case "L".
				// Se = 1/M = poro/Kf + (alpha-poro)/Ks
				rho_val = FluidProp->Density();
				if (D_Flag > 0 && rho_val > MKleinsteZahl)
				{
					biot_val = SolidProp->biot_const;
					poro_val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
					Se =
					    poro_val * FluidProp->drho_dp +
					    (biot_val - poro_val) * (1.0 - biot_val) / SolidProp->K;
					// The poroelastic portion
					val += Se * MMax(0., Sw);
				}

				// If Partial-Pressure-Based model
				if (pcs->PartialPS == 1)
				{
					// Let's get dPcdSw and apply to the mat_fac
					CMediumProperties* m_mmp = NULL;
					m_mmp = mmp_vector[0];

					// dSedPc always return positive numbers for default case
					// However, the value should be negative analytically.
					double dSwdPc =
					    m_mmp->PressureSaturationDependency(Sw, true);
					val += poro_val * dSwdPc;
				}
			}
			if (pcs->pcs_type_number == 1)
				val = MediaProp->Porosity(Index, pcs->m_num->ls_theta) *
				      MediaProp->geo_area;
			// PCH: geo_area is only used for 1 and 2 dimensions.
			// This is originated from Old RockFlow that handles 1, 2, and 3
			// dimensional
			// elements in seperate functions that show inconsistency in
			// handling geo_area.
			// For example, geo_area is never used for 3D even in Old RockFlow.
			break;
		case C:  // Componental flow
			// OK comp = m_pcs->pcs_type_number;
			// OK coefficient = MPCCalcStorativityNumber(ele,phase,comp,gp);
			break;
		//....................................................................
		case H:  // Heat transport
			TG = interpolate(NodalVal1);
			val = MediaProp->HeatCapacity(Index, pcs->m_num->ls_theta, this);
			val /= time_unit_factor;
			break;
		//....................................................................
		case M:  // Mass transport //SB4200
			// Porosity
			val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
			// SB Transport in both phases
			tr_phase = cp_vec[this->pcs->pcs_component_number]->transport_phase;
			// Multi phase transport of components
			saturation = PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type,
			                                          "SATURATION1", 1);
			if (tr_phase == 0)  // Water phase
				val *= saturation;
			else if (tr_phase == 10)  // non wetting phase
				val *= (1.0 - saturation);
			m_cp = cp_vec[pcs->pcs_component_number];
			// Retardation Factor
			val *= m_cp->CalcElementRetardationFactorNew(Index, unit, pcs);
			break;
		case R:  // Richards
			Sw = 1.0;
			dSdp = 0.;
			PG = interpolate(NodalVal1);  // 12.02.2007.  Important! WW

			if (PG <
			    0.0)  // JM skip to call these two functions in saturated case
			{
				Sw = MediaProp->SaturationCapillaryPressureFunction(-PG);
				dSdp = -MediaProp->PressureSaturationDependency(
				    Sw, true);  // JT: dSdp now returns actual sign (i.e. <0)
			}

			poro = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
			rhow = FluidProp->Density();

			if (FluidProp->compressibility_model_pressure > 0)
			{  // drho_dp from rho-p-T relation
				arg[0] = PG;
				arg[1] = interpolate(NodalValC1);  // T
				drho_dp_rho = FluidProp->drhodP(arg) / rhow;
			}
			else
				drho_dp_rho = FluidProp->drho_dp;

			// Storativity
			val =
			    MediaProp->StorageFunction(Index, unit, pcs->m_num->ls_theta) *
			    Sw;

			// Fluid compressibility
			if (rhow > 0.0) val += poro * Sw * drho_dp_rho;
			// Capillarity
			if (PG < 0.0)  // dSdp gives always a value>0, even if p>0!
				val += poro * dSdp;
			// WW
			if (MediaProp->heat_diffusion_model == 1)
			{
				//           PG = fabs(interpolate(NodalVal1));
				TG = interpolate(NodalValC) + T_KILVIN_ZERO;
				// Rv = GAS_CONSTANT;
				humi = exp(PG / (GAS_CONSTANT_V * TG * rhow));
				rhov = humi * FluidProp->vaporDensity(TG);
				//
				val -= poro * rhov * dSdp / rhow;
				val += (1.0 - Sw) * poro * rhov /
				       (rhow * rhow * GAS_CONSTANT_V * TG);
			}
			break;
		case F:  // Fluid Momentum
			val = 1.0;
			break;
		case A:  // Air (gas) flow
			val = MediaProp->Porosity(Index, pcs->m_num->ls_theta) /
			      interpolate(NodalVal1);
			break;
	}
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoefStorage()
{
	int Index = MeshElement->GetIndex();
	double val = 0.0;
	double saturation = 0.0;      // SB, BG
	int tr_phase = 0;             // SB, BG
	CompProperties* m_cp = NULL;  // CMCD
	// CompProperties *m_cp = cp_vec[pcs->pcs_component_number]; //SB4200
	switch (PcsType)
	{
		default:
			std::cout << "Fatal error in CalCoefStorage: No valid PCS type"
			          << "\n";
			break;
		case L:  // Liquid flow
			break;
		case U:  // Unconfined flow
			break;
		case G:  // MB now Groundwater flow
			break;
		case T:  // Two-phase flow
			break;
		case C:  // Componental flow
			break;
		case H:  // heat transport
			val = 0.0;
			break;
		case M:  // Mass transport //SB4200
			// CMCD
			m_cp = cp_vec[pcs->pcs_component_number];
			// Porosity
			val = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
			// SB, BG
			tr_phase = cp_vec[this->pcs->pcs_component_number]->transport_phase;
			// Multi phase transport of components
			saturation = PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type,
			                                          "SATURATION1", 1);
			if (tr_phase == 0)  // Water phase
				val *= saturation;
			else if (tr_phase == 10)        // non wetting phase
				val *= (1.0 - saturation);  // SB, BG
			val *= m_cp->CalcElementDecayRateNew(Index, pcs);
			// Retardation Factor
			val *= m_cp->CalcElementRetardationFactorNew(Index, unit, pcs);
			break;
		case R:  // Richards
			break;
		case F:  // Fluid Momentum
			break;
		case A:  // Air (gas) flow
			break;
	}
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Content matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoefContent()
{
	int Index = MeshElement->GetIndex();
	double val = 0.0;
	double dS = 0.0;
	double nodeval0, nodeval1, porval0, porval1;
	int tr_phase = 0;  // SB, BG
	// CompProperties *m_cp = NULL; //SB4200
	string name;

	switch (PcsType)
	{
		default:
			std::cout << "Fatal error in CalCoefContent: No valid PCS type"
			          << "\n";
			break;
		case L:  // Liquid flow
			break;
		case U:  // Unconfined flow
			break;
		case G:  // MB now Groundwater flow
			break;
		case T:  // Two-phase flow
			break;
		case C:  // Componental flow
			break;
		case H:  // heat transport
			break;
		case M:  // Mass transport //SB4200
		{
// kg44 added changing Porosity for GEMS coupling

#ifdef GEM_REACT
			// kg44 for GEMS coupling this should be updated to arbitrary flow
			// processes
			porval0 = PCSGetFlow()->GetElementValue(
			    Index,
			    PCSGetFlow()->GetElementValueIndex(
			        "POROSITY"));  // for GEMS we need old and new porosity!
			porval1 = PCSGetFlow()->GetElementValue(
			    Index, PCSGetFlow()->GetElementValueIndex("POROSITY") + 1);
#else
			porval0 = MediaProp->Porosity(
			    Index, pcs->m_num->ls_theta);  // get porosity..this is the
			                                   // "old" behaviour without GEMS
			                                   // coupling
			porval1 = porval0;                 // mimic no porosity change
#endif
			// Get saturation change:
			// Get saturation change, depending on phase // SB, BG
			tr_phase = cp_vec[this->pcs->pcs_component_number]->transport_phase;
			nodeval0 = PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type,
			                                        "SATURATION1", 0);
			nodeval1 = PCSGetEleMeanNodeSecondary_2(Index, pcs->flow_pcs_type,
			                                        "SATURATION1", 1);
			if (tr_phase == 10)
			{
				nodeval0 = 1.0 - nodeval0;
				nodeval1 = 1.0 - nodeval1;
			}  // SB, BG
			dS = porval1 * nodeval1 -
			     porval0 * nodeval0;  // 1/dt accounted for in assemble function
			//		if(Index == 195) cout << val << "Sat_old = " << nodeval0 <<
			//",
			// Sa_new: "<< nodeval1<< ", dS: " << dS << endl;
			val = dS;
			break;
		}
		case R:  // Richards
			break;
		case F:  // Fluid Momentum
			break;
		case A:  // Air (gas) flow
			break;
	}
	return val;
}
/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix
   Programing:
   01/2005 WW/OK Implementation
   02/2005 OK Richards flow
   03/2005 WW Heat transport
   06/2005 OK Overland flow based on CalcEle2DQuad_OF by MB
   07/2005 WW Change for geometry element object
   08/2005 OK Air (gas) flow
   01/2007 OK Two-phase flow
   10/2008 PCH Two-phase flow modified
**************************************************************************/
void CFiniteElementStd::CalCoefLaplace(bool Gravity, int ip)
{
	double dens_arg[3];  // AKS
	double mat_fac = 1.0;
	double Dpv = 0.0;
	double poro = 0.0;
	double tort = 0.0;
	double humi = 1.0;
	double rhow = 0.0;
	double* tensor = NULL;
	double w[3];
	int Index = MeshElement->GetIndex();
	double k_rel;
	ComputeShapefct(1);  //  12.3.2007 WW
	// double variables[3];                  //OK4709
	int tr_phase = 0;            // SB, BG
	double perm_effstress = 1.;  // AS:08.2012
	// WX:12.2012 perm depends on p or strain, same as CalCoefLaplace2
	CFiniteElementStd* h_fem;
	h_fem = this;
	double fac_perm = 1.0;

	// For nodal value interpolation
	//======================================================================
	switch (PcsType)
	{
		default:
			break;
		case TH:
		case L:  // Liquid flow
			tensor = MediaProp->PermeabilityTensor(Index);
			// AS:08.2012 permeability function eff stress
			if (MediaProp->permeability_effstress_model > 0)
			{
				CFiniteElementStd* h_fem;
				h_fem = this;
				perm_effstress = MediaProp->PermeabilityFunctionEffStress(
				    Index, nnodes, h_fem);
			}

			// if (ele_dim != dim)
			if (dim > MediaProp->geo_dimension)
			{
				Matrix local_tensor(dim, dim);
				Matrix temp_tensor(dim, dim);
				if (MeshElement->getTransformTensor() == NULL)
				{
					std::cout << "***Error: Geometric dimension in MMP is not "
					             "consistent with element."
					          << "\n";
					exit(0);
				}
				Matrix const& transform_tensor(*MeshElement->getTransformTensor());
				Matrix t_transform_tensor(transform_tensor);
				transform_tensor.GetTranspose(t_transform_tensor);
				Matrix global_tensor(dim, dim);
				for (size_t i = 0; i < ele_dim; i++)
					for (size_t j = 0; j < ele_dim; j++)
						local_tensor(i, j) = tensor[j + i * ele_dim];
				// cout << "K':" << endl; local_tensor.Write();
				local_tensor.multi(t_transform_tensor, temp_tensor);
				for (size_t i = 0; i < dim; i++)
					for (size_t j = 0; j < dim; j++)
						for (size_t k = 0; k < dim; k++)
							global_tensor(i, j) +=
							    transform_tensor(i, k) *
							    temp_tensor(k, j);
				// cout << "K:" << endl; global_tensor.Write();
				for (size_t i = 0; i < dim; i++)
					for (size_t j = 0; j < dim; j++)
						tensor[dim * i + j] = global_tensor(i, j);
			}
#if 0
		variables[0] = interpolate(NodalVal1); //OK4709 pressure
		if(T_Flag)
			variables[1] = interpolate(NodalValC); //OK4709 temperature
		else
			variables[1] = 15;//WX
		//OK4709
		mat_fac = FluidProp->Viscosity(variables);
#else
			mat_fac = FluidProp->Viscosity();
#endif
			// OK4709 mat_fac = FluidProp->Viscosity();
			if (gravity_constant < MKleinsteZahl)  // HEAD version
				mat_fac = 1.0;
			if (HEAD_Flag) mat_fac = 1.0;
			// Modified LBNL model WW
			if (MediaProp->permeability_stress_mode > 1)
			{
				if (cpl_pcs)
					TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
				else
					TG = 296.0;
				MediaProp->CalStressPermeabilityFactor(w, TG);
				for (size_t i = 0; i < dim; i++)
					tensor[i * dim + i] *= w[i];
			}
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] / mat_fac *
				         perm_effstress;  // AS:perm. dependent eff stress.

			break;
		case G:  // Groundwater flow
			     /* SB4218 - moved to ->PermeabilityTensor(Index);
			             if(MediaProp->permeability_model==2){ //?efficiency
			               for(i=0;i<(int)pcs->m_msh->mat_names_vector.size();i++){
			                 if(pcs->m_msh->mat_names_vector[i].compare("PERMEABILITY")==0)
			                   break;
			               }

			               mat_fac = MeshElement->mat_vector(i);
			               mat_fac /= FluidProp->Viscosity();
			              for(i=0; i<dim; i++) //WW
			                 mat[i*dim+i] = mat_fac;
			        }
			        else{
			      */
			tensor = MediaProp->PermeabilityTensor(Index);
			// TK/NW 10.10.2011
			if (dim > MediaProp->geo_dimension)
			{
				Matrix local_tensor(dim, dim);
				Matrix temp_tensor(dim, dim);
				if (MeshElement->getTransformTensor() == NULL)
				{
					std::cout << "***Error: Geometric dimension in MMP is not "
					             "consistent with element."
					          << "\n";
					exit(0);
				}
				Matrix const& transform_tensor(*MeshElement->getTransformTensor());
				Matrix t_transform_tensor(transform_tensor);
				transform_tensor.GetTranspose(t_transform_tensor);
				Matrix global_tensor(dim, dim);
				for (size_t i = 0; i < ele_dim; i++)
					for (size_t j = 0; j < ele_dim; j++)
						local_tensor(i, j) = tensor[j + i * ele_dim];
				// cout << "K':" << endl; local_tensor.Write();
				local_tensor.multi(t_transform_tensor, temp_tensor);
				for (size_t i = 0; i < dim; i++)
				{
					for (size_t j = 0; j < dim; j++)
						for (size_t k = 0; k < dim; k++)
							global_tensor(i, j) +=
							    transform_tensor(i, k) *
							    temp_tensor(k, j);
				}
				// cout << "K:" << endl; global_tensor.Write();
				for (size_t i = 0; i < dim; i++)
				{
					for (size_t j = 0; j < dim; j++)
					{
						tensor[dim * i + j] = global_tensor(i, j);
					}
				}
			}
			// TK/NW 10.10.2011
			for (size_t i = 0; i < dim * dim; i++)
				// 16.10.2009 .WW
				mat[i] = tensor[i] * time_unit_factor;
			break;
		//..................................................................
		case T:  // Two-phase flow
			// PCH Rewriting...
			// PCH Laplace mat_fac is accounted for two phases here.
			// thought to be related to the reference pressure.
			tensor = MediaProp->PermeabilityTensor(Index);
			if (pcs->pcs_type_number == 0)
			{
				if (!Gravity)
				{
					// PCH Laplace mat_fac is accounted for two phases here.
					// thought to be related to the reference pressure.
					int numOfPhases = 2;
					double mat_fac = 0.0;
					for (int p = 0; p < numOfPhases; ++p)
					{
						// PCH Check if Capillary term is on
						if (pcs->ML_Cap == 1) p = 1;

						idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");

						for (int i = 0; i < nnodes; i++)
							NodalVal_Sat[i] =
							    cpl_pcs->GetNodeValue(nodes[i], idxS + 1);

						// Whatever phase, we use Sw for getting relative
						// permeabilities
						Sw = 1.0 - interpolate(NodalVal_Sat);
						k_rel =
						    MediaProp->PermeabilitySaturationFunction(Sw, p);

						// Note here mat_fac is += meaning adding two phases
						mat_fac += time_unit_factor * k_rel /
						           mfp_vector[p]->Viscosity();
						// If Partial-Pressure-Based model
						if (pcs->PartialPS == 1) p = 1;
					}
					for (size_t i = 0; i < dim * dim; i++)
						mat[i] = tensor[i] * mat_fac;
				}
				else  // Here is only active for Cal_Velocity
				{
					// This is to calculate velocity
					// WW					int numOfPhases = 2;
					double mat_fac = 0.0;

					idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");

					for (int i = 0; i < nnodes; i++)
						NodalVal_Sat[i] =
						    cpl_pcs->GetNodeValue(nodes[i], idxS + 1);

					Sw = 1.0 - interpolate(NodalVal_Sat);
					k_rel = MediaProp->PermeabilitySaturationFunction(
					    Sw, pcs->pcs_type_number);

					// Note here mat_fac is += meaning adding two phases
					mat_fac = time_unit_factor * k_rel /
					          mfp_vector[pcs->pcs_type_number]->Viscosity();

					for (size_t i = 0; i < dim * dim; i++)
						mat[i] = tensor[i] * mat_fac;
				}
			}
			else if (pcs->pcs_type_number == 1)
			{
				// PCH Check if Capillary term is on
				if (pcs->ML_Cap == 1)
				{
					int phase = pcs->pcs_type_number;

					idxS = pcs->GetNodeValueIndex("SATURATION2");
					for (int i = 0; i < nnodes; i++)
						NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS + 1);
					Sw = 1.0 - interpolate(NodalVal_Sat);
					k_rel =
					    MediaProp->PermeabilitySaturationFunction(Sw, phase);

					// Here only the second phase accounted.
					// Let's get dPcdSw and apply to the mat_fac
					CMediumProperties* m_mmp = NULL;
					m_mmp = mmp_vector[0];
					double dPcdSw = 0.0;
					if (m_mmp->capillary_pressure_values[0] < MKleinsteZahl)
						dPcdSw = 0.0;
					else
						dPcdSw = -m_mmp->PressureSaturationDependency(
						    Sw, false);  // JT: now returns correct sign.
					mat_fac = time_unit_factor * k_rel /
					          mfp_vector[phase]->Viscosity() * (-dPcdSw);
					for (size_t i = 0; i < dim * dim; i++)
						mat[i] = tensor[i] * mat_fac;
				}
				else
				{
					int phase = pcs->pcs_type_number;

					idxS = pcs->GetNodeValueIndex("SATURATION2");
					for (int i = 0; i < nnodes; i++)
						NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS + 1);
					Sw = 1.0 - interpolate(NodalVal_Sat);
					k_rel =
					    MediaProp->PermeabilitySaturationFunction(Sw, phase);

					// Here only the second phase accounted.
					mat_fac = time_unit_factor * k_rel /
					          mfp_vector[phase]->Viscosity();
					for (size_t i = 0; i < dim * dim; i++)
						mat[i] = tensor[i] * mat_fac;
				}
			}

			break;
		//..................................................................
		case C:  // Componental flow
			break;
		case H:  // heat transport
			if (SolidProp->GetConductModel() ==
			    2)  // Boiling model. DECOVALEX THM2
			{
				TG = interpolate(NodalVal1);
				for (size_t i = 0; i < dim * dim; i++)
					mat[i] = 0.0;
				for (size_t i = 0; i < dim; i++)
					mat[i * dim + i] = SolidProp->Heat_Conductivity(TG);
			}
			// DECOVALEX THM1 or Curce 12.09. WW
			else if (SolidProp->GetConductModel() == 3 ||
			         SolidProp->GetConductModel() == 4)
			{
				// WW
				PG = interpolate(NodalValC1);
				if (cpl_pcs->type != 1212) PG *= -1.0;
				Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
				for (size_t i = 0; i < dim * dim; i++)
					mat[i] = 0.0;
				mat_fac = SolidProp->Heat_Conductivity(Sw);
				for (size_t i = 0; i < dim; i++)
					mat[i * dim + i] = mat_fac;
			}
			// WW        else if(SolidProp->GetCapacityModel()==1 &&
			// MediaProp->heat_diffusion_model == 273){
			else if (SolidProp->GetConductModel() == 1)
			{
				TG = interpolate(NodalVal1);
				tensor = MediaProp->HeatDispersionTensorNew(ip);
				for (size_t i = 0; i < dim * dim; i++)
					mat[i] = tensor[i];
			}
			else
			{
				tensor = MediaProp->HeatConductivityTensor(Index);
				for (size_t i = 0; i < dim * dim; i++)
					mat[i] = tensor[i];  // mat[i*dim+i] = tensor[i];
			}
			break;
		case M:             // Mass transport
			mat_fac = 1.0;  // MediaProp->Porosity(Index,pcs->m_num->ls_theta);
			                // // porosity now included in
			                // MassDispersionTensorNew()
			// Get transport phase of component, to obtain correct velocities in
			// dispersion tensor
			tr_phase = cp_vec[this->pcs->pcs_component_number]->transport_phase;
			// SB, BG
			tensor = MediaProp->MassDispersionTensorNew(ip, tr_phase);
			// CB
			// SB->CB I think this does not belong here
			// mat_fac *= PCSGetEleMeanNodeSecondary_2(Index,
			// pcs->flow_pcs_type, "SATURATION1", 1);
			// if(PCSGet("RICHARDS_FLOW"))
			//	    mat_fac *= PCSGetEleMeanNodeSecondary(Index,
			//"RICHARDS_FLOW", "SATURATION1", 1);
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
		//------------------------------------------------------------------
		case R:  // Richards flow
			// The following line only applies when Fluid Momentum is on
			PG = interpolate(NodalVal1);  // 05.01.07 WW
			// 05.01.07 WW
			Sw = MediaProp->SaturationCapillaryPressureFunction(-PG);

			if (MediaProp->permeability_pressure_model > 0)  // 12.2012. WX
				fac_perm = MediaProp->PermeabilityFunctionPressure(Index, PG);
			if (MediaProp->permeability_strain_model > 0)  // 12.2012 WX
				fac_perm *=
				    MediaProp->PermeabilityFunctionStrain(Index, nnodes, h_fem);

			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = time_unit_factor *
			          MediaProp->PermeabilitySaturationFunction(Sw, 0) /
			          FluidProp->Viscosity();
			// Modified LBNL model WW
			if (MediaProp->permeability_stress_mode > 1)
			{
				if (cpl_pcs)
					TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
				else
					TG = 296.0;
				MediaProp->CalStressPermeabilityFactor(w, TG);
				for (size_t i = 0; i < dim; i++)
					tensor[i * dim + i] *= w[i];
			}
			//
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * fac_perm;  // WX:12.2012

			if (MediaProp->heat_diffusion_model == 1 && !Gravity)
			{
				rhow = FluidProp->Density();
				// PG = fabs(interpolate(NodalVal1));
				TG = interpolate(NodalValC) + T_KILVIN_ZERO;
				poro = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
				tort = MediaProp->TortuosityFunction(Index, unit,
				                                     pcs->m_num->ls_theta);
				// Rv = GAS_CONSTANT;
				humi = exp(PG / (GAS_CONSTANT_V * TG * rhow));
				//
				Dpv = 2.16e-5 * tort * (1 - Sw) * poro *
				      pow(TG / T_KILVIN_ZERO, 1.8);
				Dpv *= time_unit_factor * FluidProp->vaporDensity(TG) * humi /
				       (GAS_CONSTANT_V * rhow * TG);
				for (size_t i = 0; i < dim; i++)
					mat[i * dim + i] += Dpv / rhow;
			}
			break;
		//------------------------------------------------------------------
		case A:  // Air flow
			dens_arg[0] = interpolate(NodalVal1);
			dens_arg[1] = interpolate(NodalValC1) + T_KILVIN_ZERO;
			dens_arg[2] = Index;
			double vis = FluidProp->Viscosity(dens_arg);
			mat_fac = vis;
			tensor = MediaProp->PermeabilityTensor(Index);
			// WX:09.2011
			fac_perm = 1.;
			if (MediaProp->permeability_pressure_model > 0)
			{
				fac_perm =
				    MediaProp->PermeabilityFunctionPressure(Index, dens_arg[0]);
				mat_fac /= fac_perm;
			}
			if (MediaProp->permeability_strain_model > 0)
			{
				fac_perm =
				    MediaProp->PermeabilityFunctionStrain(Index, nnodes, h_fem);
				mat_fac /= fac_perm;
			}

			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] / mat_fac;
			break;
			//------------------------------------------------------------------
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for advection matrix
   Programing:
   01/2005 WW/OK Implementation
   03/2005 WW Heat transport
   07/2005 WW Change for geometry element object
   09/2005 SB
   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoefAdvection()
{
	double val = 0.0;
	double dens_arg[3];  // AKS
	// OK long Index = MeshElement->GetIndex();
	//----------------------------------------------------------------------
	switch (PcsType)
	{
		default:
			cout << "Fatal error in CalCoefAdvection: No valid PCS type"
			     << "\n";
			break;
		case L:  // Liquid flow
			break;
		case U:  // Unconfined flow
			break;
		case G:  // MB now Groundwater flow
			break;
		case T:  // Two-phase flow
			break;
		case C:  // Componental flow
			break;
		case H:  // heat transport
			if (FluidProp->density_model == 14 &&
			    MediaProp->heat_diffusion_model == 1 && cpl_pcs)
			{
				dens_arg[0] = interpolate(NodalValC1);
				dens_arg[1] = interpolate(NodalVal1) + T_KILVIN_ZERO;
				dens_arg[2] = Index;
				val = FluidProp->SpecificHeatCapacity(dens_arg) *
				      FluidProp->Density(dens_arg);
			}
			else
				val = FluidProp->SpecificHeatCapacity() * FluidProp->Density();
			break;
		case M:  // Mass transport //SB4200
			val =
			    1.0 *
			    time_unit_factor;  //*MediaProp->Porosity(Index,pcs->m_num->ls_theta);
			                       //// Porosity;
			break;
		case R:  // Richards
			break;
		case F:  // Fluid Momentum
			break;
		case A:                                  // Air (gas) flow
			val = 1.0 / interpolate(NodalVal1);  // 1/p
			break;
	}
	return val;
}

/***************************************************************************
   GeoSys - Function: CalCoefStrainCouping

   Aufgabe:
      Calculate coefficient for StrainCouping matrix
   Programmaenderungen:
   01/2005   WW/OK    Erste Version
   07/2005 WW Change for geometry element object
 **************************************************************************/
double CFiniteElementStd::CalCoefStrainCouping(const int phase)
{
	double val = 0.0;
	/*
	   double r = unit[0];
	   double s = unit[1];
	   double t = unit[2];
	 */

	switch (PcsType)
	{
		default:
			break;
		case L:  // Liquid flow
			//
			val = 1.0;
			break;
		case U:  // Unconfined flow
		case G:  // Groundwater
		case T:  // Two-phase flow
		case C:  // Componental flow
			break;
		case R:                                // Richard flow
			return interpolate(NodalVal_Sat);  // Water saturation
			break;
		case V:
			if (phase == 0)
			{
				PG = interpolate(NodalVal1);
				Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
				val = Sw;
			}
			else
				val = 1. - Sw;
			return val;
			break;
	}
	return val;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcMass
   Aufgabe:
           Compute mass matrix, i.e. int (N.mat.N). Linear interpolation

   Programming:
   01/2005   WW
   02/2005 OK GEO factor
   01/2010 NW SUPG
 **************************************************************************/
void CFiniteElementStd::CalcMass()
{
	int i, j;  // OK411 k;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, mat_fac;
	// Material
	mat_fac = 1.0;
	double vel[3];  // NW
//  int indice = MeshElement->GetIndex();
//  int phase = pcs->pcs_type_number;

#ifndef OGS_ONLY_TH
	int upwind_method = pcs->m_num->ele_upwind_method;
	double alpha[3] = {}, summand[8] = {};
	if (PcsType == T)
		if (upwind_method > 0)
		{
			// CB 11/07 this is to provide the velocity at the element center of
			// gravity
			// call to his function here is also required for upwinding in
			// CalcCoefLaplace
			Cal_Velocity_2();
			UpwindAlphaMass(alpha);  // CB 160507
		}
#endif

	ElementValue* gp_ele = ele_gp_value[Index];  // NW

	//----------------------------------------------------------------------
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		// if(PcsType==T)
		//{
		//  if((upwind_method == 1) || (upwind_method == 2))
		//      UpwindUnitCoord(phase, gp, indice); // phase 0
		//}
		ComputeShapefct(1);                   // Linear interpolation function
		if (pcs->m_num->ele_supg_method > 0)  // NW
			ComputeGradShapefct(1);           // Linear interpolation function

		// Material
		mat_fac = CalCoefMass();
		// if(Index < 0) cout << "mat_fac in CalCoeffMass: " << mat_fac << endl;
		// GEO factor
		mat_fac *= fkt;
// Calculate mass matrix
#ifndef OGS_ONLY_TH
		if (PcsType == T)
		{
			// upwinding: addiere SUPG-Summanden auf shapefct entsprechend Fkt.
			// Mphi2D_SPG
			if (pcs->m_num->ele_upwind_method > 0)
				UpwindSummandMass(gp, gp_r, gp_s, gp_t, alpha, summand);
			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					(*Mass)(i, j) +=
					    mat_fac * (shapefct[i] + summand[i]) * shapefct[j];
		}
		else
#endif
		{
			for (i = 0; i < nnodes; i++)
			{
				for (j = 0; j < nnodes; j++)
				{
					if (pcs->m_num->ele_supg_method == 0)
						if (j > i) continue;
					(*Mass)(i, j) += mat_fac * shapefct[i] * shapefct[j];
				}
			}

			if (pcs->m_num->ele_supg_method > 0)
			{
				vel[0] = gp_ele->Velocity(0, gp);
				vel[1] = gp_ele->Velocity(1, gp);
				vel[2] = gp_ele->Velocity(2, gp);

				double tau = 0;
				CalcSUPGWeightingFunction(vel, gp, tau, weight_func);

				// tau*({v}[dN])^T*[N]
				for (i = 0; i < nnodes; i++)
					for (j = 0; j < nnodes; j++)
						(*Mass)(i, j) +=
						    mat_fac * tau * weight_func[i] * shapefct[j];
			}
		}  // end else
	}      // loop gauss points

	if (PcsType != T && pcs->m_num->ele_supg_method == 0)
	{
		for (i = 0; i < nnodes; i++)
			for (j = 0; j < nnodes; j++)
				if (j > i) (*Mass)(i, j) = (*Mass)(j, i);
	}
	// Test Output
	// Mass->Write();
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcLumpedMass
   Aufgabe:
           Compute lumped mass matrix, i.e. int (N.mat.N). Linear interpolation

   Programming:
   01/2005   WW
   02/2005 OK GEO factor
   02/2007   WW Multi-phase flow
   05/2007   WW Axismmetry volume
   01/2010   NW geometrical area
 **************************************************************************/
void CFiniteElementStd::CalcLumpedMass()
{
	int i, gp_r, gp_s, gp_t;
	double factor, vol = 0.0;
	gp = 0;
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	// Initialize
	(*Mass) = 0.0;
	// Volume
	if (axisymmetry)
	{  // This calculation should be done in CompleteMesh.
		// However, in order not to destroy the concise of the code,
		// it is put here. Anyway it is computational cheap. WW
		vol = 0.0;
		for (gp = 0; gp < nGaussPoints; gp++)
			//---------------------------------------------------------
			//  Get local coordinates and weights
			//  Compute Jacobian matrix and its determinate
			//---------------------------------------------------------
			vol += GetGaussData(gp, gp_r, gp_s, gp_t);
	}
	else
		// NW multiply geo_area
		vol = MeshElement->GetVolume() * MeshElement->GetFluxArea();
	// Center of the reference element
	SetCenterGP();
	factor = CalCoefMass();
	pcs->timebuffer = factor;  // Tim Control "Neumann"
	factor *= vol / (double)nnodes;
	for (i = 0; i < nnodes; i++)
		(*Mass)(i, i) = factor;
//
#ifdef otherLumpedMass
	int i, j;
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt;
	//----------------------------------------------------------------------
	for (i = 0; i < nnodes; i++)
	{
		for (j = 0; j < ele_dim; j++)
			x2buff[j] = nodes_xyz[j * nnodes + i];
		UnitCoordinates(x2buff);
		fkt = GetGaussData(i, gp_r, gp_s, gp_t) * CalCoefMass();
		(*Mass)(i + in * nnodes, i + jn * nnodes) += fkt;
	}
//----------------------------------------------------------------------
#endif
	// TEST OUT
	// Mass->Write();
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcStorage
   Aufgabe:
           Compute mass matrix, i.e. int (N.mat.N). Linear interpolation

   Programming:
   01/2005   WW
   02/2005 OK GEO factor
 **************************************************************************/
void CFiniteElementStd::CalcStorage()
{
	int i, j;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, mat_fac;
	// Material
	mat_fac = 1.0;
	//----------------------------------------------------------------------
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeShapefct(1);  // Linear interpolation function
		// Material
		mat_fac = CalCoefStorage();
		// GEO factor
		fkt *= mat_fac;
		for (i = 0; i < nnodes; i++)
			for (j = 0; j < nnodes; j++)
				(*Storage)(i, j) += fkt * shapefct[i] * shapefct[j];
	}
	// TEST OUTPUT
	//  if(Index == 195){cout << "Storage Matrix: " << "\n"; Storage->Write(); }
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcContent
   Aufgabe:
           Compute Content matrix, i.e. int (N.mat.N). Linear interpolation

   Programming:
   01/2005   WW
   02/2005 OK GEO factor
 **************************************************************************/
void CFiniteElementStd::CalcContent()
{
	int i, j;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, mat_fac;
	// Material
	mat_fac = 1.0;
	//----------------------------------------------------------------------
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeShapefct(1);  // Linear interpolation function
		// Material
		mat_fac = CalCoefContent();
		// GEO factor
		fkt *= mat_fac;
		for (i = 0; i < nnodes; i++)
			for (j = 0; j < nnodes; j++)
				(*Content)(i, j) += fkt * shapefct[i] * shapefct[j];
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcLaplace
   Aufgabe:
           Compute mass matrix, i.e. int (gradN.mat.gradN). Linear interpolation

   Programming:
   01/2005   WW
   02/2005 OK GEO factor
   02/2007 WW Multi-phase
    03/2009 PCH PS_GLOBAL for Multiphase flow
 **************************************************************************/
void CFiniteElementStd::CalcLaplace()
{
	int i, j, l, in, jn;
	size_t k;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp_t = 0;
	double fkt;
	int dof_n = 1;

	// 03.03 2009 PCH
	if (PcsType == V || PcsType == P) dof_n = 2;

	//----------------------------------------------------------------------
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
		ComputeGradShapefct(1);  // Linear interpolation function
                                 // Calculate mass matrix
#ifndef OGS_ONLY_TH
		double water_depth = 1.0;
		// The following "if" is done by WW
		if (PcsType == G && MediaProp->unconfined_flow_group == 1 &&
			MeshElement->GetDimension() == 2 && !pcs->m_msh->hasCrossSection())
		{
			water_depth = 0.0;
			for (i = 0; i < nnodes; i++)
				water_depth +=
				    (pcs->GetNodeValue(nodes[i], idx1) - Z[i]) * shapefct[i];
		}
		fkt *= water_depth;
#endif
		//---------------------------------------------------------

		for (in = 0; in < dof_n; in++)
		{
			//			const int ishd = in * dof_n ;
			const int ish = in * nnodes;
			for (jn = 0; jn < dof_n; jn++)
			{
				// Material
				if (dof_n == 1) CalCoefLaplace(false, gp);
#ifndef OGS_ONLY_TH
				else if (dof_n == 2)
				{
					if (PcsType == V)
						CalCoefLaplace2(false, in * dof_n + jn);
					else if (PcsType == P)
						CalCoefLaplacePSGLOBAL(false, in * dof_n + jn);
				}
#endif
				const int jsh = jn * nnodes;
				for (i = 0; i < nnodes; i++)
				{
					const int iish = i + ish;
					for (j = 0; j < nnodes; j++)
					{
						const int jjsh = j + jsh;
						//  if(j>i) continue;
						for (k = 0; k < dim; k++)
						{
							const int ksh = k * nnodes + i;
							const int km = dim * k;
							for (l = 0; l < (int)dim; l++)
							{
								(*Laplace)(iish, jjsh) +=
								    fkt * dshapefct[ksh] * mat[km + l] *
								    dshapefct[l * nnodes + j];
							}
						}
					}  // j: nodes
				}      // i: nodes
			}
		}  //	//TEST OUTPUT
	}
	// Laplace->Write();
}
/**************************************************************************
   FEMLib-Method:
   10/2006 YD Implementation
   01/2007 WW Fundamental changes
**************************************************************************/
void CFiniteElementStd::Assemble_DualTransfer()
{
	int i, j;
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double W, fkt, mat_fac = 0.;
#if defined(NEW_EQS)
	CSparseMatrix* A = pcs->eqs_new->getA();
#endif

	// Inintialize
	//-------------------------- WW
	W = pcs->continuum_vector[pcs->GetContinnumType()];
	//
	for (i = 0; i < nnodes; i++)
	{
		// Pressure 1
		NodalVal3[i] = pcs->GetNodeValue(nodes[i], idx1);
		// Pressure 2
		NodalVal4[i] = pcs->GetNodeValue(nodes[i], idxp21);
	}
	(*Advection) = 0.0;
	//---------------------------------------------------------
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determination
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Material
		ComputeShapefct(1);  // Moved here by NW 25.10.2011
		mat_fac = CalcCoefDualTransfer();
		mat_fac *= fkt;
		// Calculate mass matrix
		for (i = 0; i < nnodes; i++)
			for (j = 0; j < nnodes; j++)
				(*Advection)(i, j) += mat_fac * shapefct[i] * shapefct[j];
	}
// Add local matrix to global matrix
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	// 15.02.2007 WW
	long cshift = pcs->m_msh->GetNodesNumber(false);
#endif
	double fm = 1.0 / W;
	//
	if (pcs->continuum == 0)
	{
#ifndef USE_PETSC
		double ff = 1.0 / (1.0 - W);
		if (MediaProp->transfer_coefficient < 0.0)  // for LBNL
			ff = 1.0;
#endif
		for (int i = 0; i < nnodes; i++)
		{
			for (int j = 0; j < nnodes; j++)
			{
#ifdef NEW_EQS
				(*A)(eqs_number[i], eqs_number[j] + cshift) +=
				    -fm * (*Advection)(i, j);
				(*A)(eqs_number[i] + cshift, eqs_number[j]) +=
				    -ff * (*Advection)(i, j);
#endif
			}
		}
	}
	else if (MediaProp->transfer_coefficient < 0.0)  // for LBNL
		fm = 1.0;
	//
	(*Advection) *= fm;
	(*Laplace) += (*Advection);
	//
	//-------------------------- WW
}
/**************************************************************************
   FEMLib-Method:
   10/2006 YD Implementation
   01/2007 WW Fundamental changes
**************************************************************************/
double CFiniteElementStd::CalcCoefDualTransfer()
{
	double Sm = 0.0, Sf = 0.0, ExFac = 0.0;
	double pm = 0.0, pf = 0.0, matrix_conductivity, val = 0;
	// double* permeability;
	double* permeability = NULL;
	//-------------------------------------------WW
	CMediumProperties* m_matrix = NULL;
	CMediumProperties* f_matrix = NULL;
	if (pcs->GetContinnumType() == 0)
	{
		m_matrix = MediaProp;
		f_matrix = MediaProp1;
	}
	else  // fracture //WW
	{
		m_matrix = MediaProp1;
		f_matrix = MediaProp;
	}
	//-------------------------------------------WW
	switch (PcsType)
	{
		default:
			break;
		case R:
			pm = interpolate(NodalVal3);
			pf = interpolate(NodalVal4);
			// Matrix
			Sm = m_matrix->SaturationCapillaryPressureFunction(-pm);
			// Fracture
			Sf = f_matrix->SaturationCapillaryPressureFunction(-pf);
			permeability = m_matrix->PermeabilityTensor(Index);
			ExFac = m_matrix->transfer_coefficient;
			// Dual by van Genuchten
			if (ExFac > 0.0)
				matrix_conductivity =
				    0.5 * (m_matrix->PermeabilitySaturationFunction(Sm, 0) +
				           m_matrix->PermeabilitySaturationFunction(Sf, 0)) /
				    FluidProp->Viscosity();

			else  // by LBNL. WW
			{
				double Sf_e =
				    f_matrix->GetEffectiveSaturationForPerm(Sf, phase);
				matrix_conductivity =
				    Sf_e * m_matrix->PermeabilitySaturationFunction(Sm, 0) /
				    FluidProp->Viscosity();
				ExFac *= -1.0;
			}
			//
			val = time_unit_factor * permeability[0] * matrix_conductivity *
			      ExFac;
			break;
		//---------------------------------------------------------
		case H:

			break;
	}
	return val;
}

// SB4200
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcAdvection
   Aufgabe:  Calculate the advection matrix

   Programming:
   01/2005   WW
   02/2005   OK GEO factor
   09/2005   SB - adapted to advection
   03/2007   WW - Fluid advection with multiphase flow
   05/2008   WW - General densty for multiphase flow
   01/2010   NW - SUPG
 **************************************************************************/
void CFiniteElementStd::CalcAdvection()
{
	int i, j;
	int gp_r = 0, gp_s = 0, gp_t;
	double fkt, mat_factor = 0.0;
	double vel[3];
#ifndef OGS_ONLY_TH
	double dens_aug[3];
	CFluidProperties* m_mfp_g = NULL;
	bool multiphase = false;
	// 18.02.2008, 04.09.2008 WW
	if (!cpl_pcs && (pcs->type != 2) && (pcs->type != 5)) return;
	if (cpl_pcs && cpl_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
	{
		multiphase = true;
		m_mfp_g = mfp_vector[1];
		GasProp = MFPGet("GAS");
	}
#endif
	ElementValue* gp_ele = ele_gp_value[Index];
#ifndef OGS_ONLY_TH
	CRFProcess* pcs_fluid_momentum = PCSGet("FLUID_MOMENTUM");
#endif

	// Initial values
	gp_t = 0;
#ifndef OGS_ONLY_TH
	(*Advection) = 0.0;
#endif

	//----------------------------------------------------------------------
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		//---------------------------------------------------------
		// Compute geometry
		ComputeGradShapefct(
		    1);  // Linear interpolation function....dNJ-1....var dshapefct
		ComputeShapefct(1);  // Linear interpolation N....var shapefct
		//---------------------------------------------------------
		mat_factor = CalCoefAdvection();  // this should be called after
		                                  // calculating shape functions. NW
		// Velocity
		vel[0] = mat_factor * gp_ele->Velocity(0, gp);
		vel[1] = mat_factor * gp_ele->Velocity(1, gp);
		vel[2] = mat_factor * gp_ele->Velocity(2, gp);
#ifndef OGS_ONLY_TH
		// If component is in non - wetting phase, as designated by
		// transport_phase == 10 // SB, BG
		if (cp_vec.size() > 0 && this->pcs->pcs_component_number >= 0)
			// // SB, BG
			if (cp_vec[this->pcs->pcs_component_number]->transport_phase == 10)
			{
				vel[0] = mat_factor * gp_ele->Velocity_g(0, gp);
				vel[1] = mat_factor * gp_ele->Velocity_g(1, gp);
				vel[2] = mat_factor * gp_ele->Velocity_g(2, gp);
			}            // SB, BG
		if (multiphase)  // 02/2007 WW
		{
			PG2 = interpolate(NodalVal_p2);
			PG = interpolate(NodalValC1);
			TG = interpolate(NodalVal1) + T_KILVIN_ZERO;
			rhow = FluidProp->Density();
			rho_gw =
			    FluidProp->vaporDensity(TG) *
			    exp(-PG * COMP_MOL_MASS_WATER / (rhow * GAS_CONSTANT * TG));
			p_gw = rho_gw * GAS_CONSTANT * TG / COMP_MOL_MASS_WATER;
			dens_aug[0] = PG2 - p_gw;
			dens_aug[1] = TG;
			// 29.05.2008. WW/ 2 Dec 2010 AKS
			rho_g = rho_gw + GasProp->Density(dens_aug);
			mat_factor = rho_g * m_mfp_g->SpecificHeatCapacity();
			vel[0] += mat_factor * gp_ele->Velocity_g(0, gp);
			vel[1] += mat_factor * gp_ele->Velocity_g(1, gp);
			vel[2] += mat_factor * gp_ele->Velocity_g(2, gp);
		}

		// Velocity by Fluid_Momentum - 13.11.2009  PCH
		if (pcs_fluid_momentum)
		{
			CRFProcess* m_pcs = pcs_fluid_momentum;

			vel[0] = mat_factor *
			         m_pcs->GetElementValue(
			             index, m_pcs->GetElementValueIndex("VELOCITY1_X") + 1);
			vel[1] = mat_factor *
			         m_pcs->GetElementValue(
			             index, m_pcs->GetElementValueIndex("VELOCITY1_Y") + 1);
			vel[2] = mat_factor *
			         m_pcs->GetElementValue(
			             index, m_pcs->GetElementValueIndex("VELOCITY1_Z") + 1);
		}
#endif
		for (i = 0; i < nnodes; i++)
			for (j = 0; j < nnodes; j++)
				for (size_t k = 0; k < dim; k++)
					(*Advection)(i, j) +=
					    fkt * shapefct[i] * vel[k] * dshapefct[k * nnodes + j];

		if (pcs->m_num->ele_supg_method > 0)
		{
			vel[0] = gp_ele->Velocity(0, gp);
			vel[1] = gp_ele->Velocity(1, gp);
			vel[2] = gp_ele->Velocity(2, gp);
#ifndef OGS_ONLY_TH
			if (pcs_fluid_momentum)
			{
				CRFProcess* m_pcs = pcs_fluid_momentum;

				vel[0] = m_pcs->GetElementValue(
				    index, m_pcs->GetElementValueIndex("VELOCITY1_X") + 1);
				vel[1] = m_pcs->GetElementValue(
				    index, m_pcs->GetElementValueIndex("VELOCITY1_Y") + 1);
				vel[2] = m_pcs->GetElementValue(
				    index, m_pcs->GetElementValueIndex("VELOCITY1_Z") + 1);
			}
#endif

			double tau = 0;
			CalcSUPGWeightingFunction(vel, gp, tau, weight_func);

			// Calculate mat_factor*tau*({v}[dN])^T*({v}[dN])
			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					(*Advection)(i, j) += fkt * mat_factor * tau *
					                      weight_func[i] * weight_func[j];
		}
	}
	// TEST OUTPUT
	// cout << "Advection Matrix: " << endl; Advection->Write();
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcAdvection
   Aufgabe:  Calculate the advection matrix

   Programming:
   12/2005   WW
***************************************************************************/
void CFiniteElementStd::CalcRHS_by_ThermalDiffusion()
{
	int i, j;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp = 0;
	double fkt;
	double Dv = 0.0;
	double Dtv = 0.0;
	double poro = 0.0;
	double tort = 0.0;
	double humi = 1.0;
	double rhov = 0.0;
	double drdT = 0.0;
	double beta = 0.0;
	// 12.12.2007 WW
	long cshift = 0;
	if (pcs->dof > 1) cshift = NodeShift[pcs->continuum];

	(*Laplace) = 0.0;
	(*Mass) = 0.0;
	//----------------------------------------------------------------------
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
		ComputeGradShapefct(1);  // Linear interpolation function
		//---------------------------------------------------------
		// Material
		//	  if(FluidProp->diffusion_model=273)
		//    {
		ComputeShapefct(1);
		double rhow = FluidProp->Density();
		PG = interpolate(NodalVal1);
		TG = interpolate(NodalValC) + T_KILVIN_ZERO;
		// WW
		Sw = MediaProp->SaturationCapillaryPressureFunction(-PG);
		poro = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
		tort = MediaProp->TortuosityFunction(Index, unit, pcs->m_num->ls_theta);
		beta = poro *
		       MediaProp->StorageFunction(Index, unit, pcs->m_num->ls_theta) *
		       Sw;
		// Rv = GAS_CONSTANT;
		humi = exp(PG / (GAS_CONSTANT_V * TG * rhow));
		Dv = MediaProp->base_heat_diffusion_coefficient * tort * (1 - Sw) *
		     poro * pow(TG / T_KILVIN_ZERO, 1.8);
		rhov = humi * FluidProp->vaporDensity(TG);
		drdT = (FluidProp->vaporDensity_derivative(TG) * humi -
		        rhov * PG / (GAS_CONSTANT_V * rhow * TG * TG)) /
		       rhow;
		Dtv = time_unit_factor * Dv * drdT;

		//    }
		//---------------------------------------------------------
		// Calculate a Laplace
		for (i = 0; i < nnodes; i++)
			for (j = 0; j < nnodes; j++)
			{
				if (j > i) continue;
				for (size_t k = 0; k < dim; k++)
					(*Laplace)(i, j) += fkt * Dtv * dshapefct[k * nnodes + i] *
					                    dshapefct[k * nnodes + j];
				(*Mass)(i, j) += fkt * poro * (beta + (1.0 - Sw) * drdT) *
				                 shapefct[i] * shapefct[j];
			}
	}
	// Symmetry
	for (i = 0; i < nnodes; i++)
		for (j = 0; j < nnodes; j++)
		{
			if (j <= i) continue;
			(*Laplace)(i, j) = (*Laplace)(j, i);
		}
	cshift += NodeShift[problem_dimension_dm];

	for (i = 0; i < nnodes; i++)
	{
		for (j = 0; j < nnodes; j++)
		{
			(*RHS)(i) -= (*Laplace)(i, j) * (NodalValC[j] + T_KILVIN_ZERO);
			(*RHS)(i) += (*Mass)(i, j) * (NodalValC1[j] - NodalValC[j]) / dt;
		}
		eqs_rhs[cshift + eqs_number[i]] += (*RHS)(i);
	}

	// TEST OUTPUT
	// Laplace->Write();
	// Mass->Write();
	// RHS->Write();
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: Coordinates for high order nodes
   Aufgabe:
           Compute the strain couping matrix

   Programming:
   02/2007   WW
 **************************************************************************/
void CFiniteElementStd::SetHighOrderNodes()
{
	int i = 0;
	setOrder(2);
	// Swap cordinates in case of (x, 0.0, z) only for 2D problem
	if (coordinate_system % 10 == 2)  // Z has number
	{
		switch (dim)
		{
			case 1:
				for (i = 0; i < nNodes; i++)
				{
					//				X[i] = MeshElement->GetNode(i)->Z();
					//				Y[i] = MeshElement->GetNode(i)->Y();
					//				Z[i] = MeshElement->GetNode(i)->X();
					double const* const coords(
					    MeshElement->GetNode(i)->getData());
					X[i] = coords[2];
					Y[i] = coords[1];
					Z[i] = coords[0];
				}
				break;
			case 2:
				for (i = 0; i < nNodes; i++)
				{
					//				X[i] = MeshElement->GetNode(i)->X();
					//				Y[i] = MeshElement->GetNode(i)->Z();
					//				Z[i] = MeshElement->GetNode(i)->Y();
					double const* const coords(
					    MeshElement->GetNode(i)->getData());
					X[i] = coords[0];
					Y[i] = coords[2];
					Z[i] = coords[1];
				}
				break;
			case 3:
				for (i = nnodes; i < nnodesHQ; i++)
				{
					//				X[i] = MeshElement->GetNode(i)->X();
					//				Y[i] = MeshElement->GetNode(i)->Y();
					//				Z[i] = MeshElement->GetNode(i)->Z();
					double const* const coords(
					    MeshElement->GetNode(i)->getData());
					X[i] = coords[0];
					Y[i] = coords[1];
					Z[i] = coords[2];
				}
				break;
		}
	}
	else
	{
		if (dim == 1 || dim == 2)
			for (i = nnodes; i < nnodesHQ; i++)
			{
				//				X[i] = MeshElement->GetNode(i)->X();
				//				Y[i] = MeshElement->GetNode(i)->Y();
				//				Z[i] = MeshElement->GetNode(i)->Z();
				double const* const coords(MeshElement->GetNode(i)->getData());
				X[i] = coords[0];
				Y[i] = coords[1];
				Z[i] = coords[2];
			}
	}
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcStrainCoupling
   Aufgabe:
           Compute the strain couping matrix

   Programming:
   01/2005   WW
 **************************************************************************/
void CFiniteElementStd::CalcStrainCoupling(int phase)
{
	(*StrainCoupling) = 0.0;

	SetHighOrderNodes();
	for (int gp = 0; gp < nGaussPoints; gp++)
	{
		int gp_r, gp_s, gp_t;
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		ComputeGradShapefct(2);
		ComputeShapefct(2);
		ComputeShapefct(1);
		//
		fkt *= CalCoefStrainCouping(phase);
		for (size_t i = 0; i < dim; i++)
		{
			for (int k = 0; k < nnodes; k++)
			{
				for (int l = 0; l < nnodesHQ; l++)
				{
					int kl = nnodesHQ * i + l;
					double du = dshapefctHQ[kl];
					if (i == 0 && axisymmetry)
						du += shapefctHQ[l] / Radius;
					(*StrainCoupling)(k, kl) += shapefct[k] * du * fkt;
				}
			}
		}
	}

	setOrder(1);
	// StrainCoupling->Write();

	if (pcs->matrix_file)
	{
		(*pcs->matrix_file) << "---Strain couping matrix: " << endl;
		StrainCoupling->Write(*pcs->matrix_file);
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: Assemby_Gravity
   Aufgabe:
           Assemble the contribution of gravity to RHS in Darcy flow
           to the global system

   Programming:
   01/2005   WW/OK
   08/2006   WW Re-implement
   02/2007   WW Multi-phase flow
 **************************************************************************/
// Local assembly
void CFiniteElementStd::Assemble_Gravity()
{
	// int Index = MeshElement->GetIndex();
	if ((coordinate_system) % 10 != 2)  // NW: exclude (!axisymmetry)

		// 27.2.2007 WW (*GravityMatrix) = 0.0;
		return;
	int i, ii;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp_t = 0;
	double fkt, rho;     //, dens_arg[3];         //, rich_f;
#if !defined(USE_PETSC)  // && defined(other parallel libs)//03~04.3012. WW
	double k_rel_iteration = 1.0;
#endif
	// GEO
	// NW  double geo_fac = MediaProp->geo_area;
	if (!FluidProp->CheckGravityCalculation()) return;
	long cshift = 0;  // WW
	//
	//
	int dof_n = 1;  // 27.2.2007 WW
	if (PcsType == V || PcsType == P) dof_n = 2;

	// WW 05.01.07
	cshift = 0;
	if (pcs->dof > 1) cshift = NodeShift[pcs->continuum];

	// rich_f = 1.0;
	// if(PcsType==R) rich_f = -1.0; //WW

	//	k_rel_iteration = 1.0;

	for (i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;

	// (*GravityMatrix) = 0.0;
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determination
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Moved from CalCoefLaplace(). 12.3.2007 WW
		// Material
		rho = FluidProp->Density();
		if (gravity_constant < MKleinsteZahl)  // HEAD version
			rho = 1.0;
		else if (HEAD_Flag)
			rho = 1.0;
		else
			rho *= gravity_constant;
		fkt *= rho;  //*rich_f;
		//
		for (ii = 0; ii < dof_n; ii++)
		{
			if (dof_n == 1)
			{
				if (PcsType == T)
					CalCoefLaplace(false);
				else
					CalCoefLaplace(true);
			}

			if (dof_n == 2)
			{
				if (PcsType == V)
					CalCoefLaplace2(true, ii * dof_n + 1);
				else if (PcsType == P)
					CalCoefLaplacePSGLOBAL(true, ii * dof_n);
			}
			// Calculate mass matrix
			const int iinn = ii * nnodes;  // 19.06.2012. WW
			for (i = 0; i < nnodes; i++)
			{
				const int ipiinn = iinn + i;  // 19.06.2012. WW
				for (size_t k = 0; k < dim; k++)
				{
					NodalVal[ipiinn] -= fkt * dshapefct[k * nnodes + i] *
					                    mat[dim * k + dim - 1];
				}
			}
		}
	}
	//
	/// 02.2011. WW
	int dm_shift = 0;
	if (pcs->type / 10 == 4) dm_shift = problem_dimension_dm;
	int ii_sh = 0;
	for (ii = 0; ii < dof_n; ii++)  // 07.02.07 WW
	{
		if (pcs->type == 22)  // Dual porosity model. 06.2011. WW
			cshift += NodeShift[ii + dm_shift];
		else
			cshift = NodeShift[ii + dm_shift];
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && defined(other parallel libs)//03~04.3012. WW
			if (add2global)
				eqs_rhs[cshift + eqs_number[i]] +=
				    k_rel_iteration * NodalVal[i + ii_sh];
// NW not necessary to multiply geo_area(geo_fac) here. It's already multiplied
// in ComputeJacobian() through fkt.
//          eqs_rhs[cshift + eqs_number[i]]
//                  += k_rel_iteration* geo_fac*NodalVal[i+ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) += NodalVal[i + ii_sh];
		}
	}
	// TEST OUTPUT
	// RHS->Write();
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: Velocity calulation

   Programming:  WW
   08/2005
   03/2007   WW  Multi-phase flow
   01/2010   NW  Fix multi-dimensional case
   02/2010   WW  Fix a bug in velocity of the first phase
 **************************************************************************/
// Local assembly
void CFiniteElementStd::Cal_Velocity()
{
	int k;
	static double vel[3], vel_g[3];
	double dens_arg[3];  // AKS
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	double coef = 0.0;
	int dof_n = 1;
	if (PcsType == V || PcsType == P) dof_n = 2;
	//
	gp_t = 0;

	// Get room in the memory for local matrices
	SetMemory();
	// Set material
	SetMaterial();

	ElementValue* gp_ele = ele_gp_value[Index];

	// gp_ele->Velocity = 0.0; // CB commented and inserted below due to
	// conflict with transport calculation, needs velocities
	// Loop over Gauss points
	k = (coordinate_system) % 10;

	if (PcsType == TH) idx1 = idxp1;

	if (PcsType == T)  // WW/CB
	{
		if (pcs->pcs_type_number == 0)
		{
			// gas pressure
			idx1 = pcs->GetNodeValueIndex("PRESSURE1") + 1;
			for (int i = 0; i < nnodes; i++)
				NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
		}
		else if (pcs->pcs_type_number == 1)
		{
			idxp21 = pcs->GetNodeValueIndex("PRESSURE_CAP");
			// gas pressure
			idx1 = cpl_pcs->GetNodeValueIndex("PRESSURE1") + 1;
			gp_ele = ele_gp_value[Index + (long)pcs->m_msh->ele_vector.size()];
			for (int i = 0; i < nnodes; i++)
				// P_l = P_g - P_cap
				NodalVal[i] = cpl_pcs->GetNodeValue(nodes[i], idx1) -
				              pcs->GetNodeValue(nodes[i], idxp21);
		}
	}
	else
		// This should be enough for Vw in PS_GLOBAL as well,
		// since the first primary variable is Pw.
		for (int i = 0; i < nnodes; i++)
		{
			NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
			NodalVal1[i] = NodalVal[i];
		}
	//
	if (PcsType == V)
	{
		gp_ele->Velocity_g = 0.0;  // WW
		for (int i = 0; i < nnodes; i++)
		{
			// 02.2010. WW
			NodalVal2[i] = pcs->GetNodeValue(nodes[i], idxp21);
			NodalVal[i] = NodalVal2[i] - NodalVal[i];
		}
	}
	if (PcsType == P)
	{
		gp_ele->Velocity_g = 0.0;  // PCH
		// Just get Pnw, which is the secondary variable in PS_GLOBAL
		int idx_pn = pcs->GetNodeValueIndex("PRESSURE2");
		for (int i = 0; i < nnodes; i++)
			NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx_pn);
	}
	Matrix tmp_gp_velocity(gp_ele->Velocity);
	tmp_gp_velocity = 0.0;
	// gp_ele->Velocity = 0.0;                     // CB inserted here and
	// commented above due to conflict with transport calculation, needs
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determination
		//---------------------------------------------------------
		GetGaussData(gp, gp_r, gp_s, gp_t);

		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Moved from CalCoefLaplace(). 12.3.2007 WW
		// WW/CB
		if ((PcsType == T) && (pcs->pcs_type_number == 1)) flag_cpl_pcs = true;
		// Material
		if (dof_n == 1)
			CalCoefLaplace(true);
		else if (dof_n == 2 && PcsType == V)  // PCH 05.2009
			CalCoefLaplace2(true, 0);
		else if (dof_n == 2 && PcsType == P)  // PCH 05.2009
			CalCoefLaplacePSGLOBAL(true, 0);
		else if (dof_n == 2 && PcsType == TH)
			CalCoefLaplace(true);
		// WW/CB
		if ((PcsType == T) && (pcs->pcs_type_number == 1)) flag_cpl_pcs = false;
		// Velocity
		for (size_t i = 0; i < dim; i++)
		{
			vel[i] = 0.0;
			for (int j = 0; j < nnodes; j++)
				vel[i] += NodalVal[j] * dshapefct[i * nnodes + j];
			//			 vel[i] += fabs(NodalVal[j])*dshapefct[i*nnodes+j];
		}
		if (PcsType == V || PcsType == P)  // PCH 05.2009

			for (size_t i = 0; i < dim; i++)
			{
				vel_g[i] = 0.0;
				for (int j = 0; j < nnodes; j++)
					// Change   NodalVal2 to NodalVal1. 02.2010. WW
					vel_g[i] += NodalVal2[j] * dshapefct[i * nnodes + j];
			}
		// Gravity term
		// NW
		if (k == 2 && (!HEAD_Flag) && FluidProp->CheckGravityCalculation())
		{
			if (FluidProp->density_model == 14)
			{
				dens_arg[0] = interpolate(NodalVal1);
				dens_arg[1] = interpolate(NodalValC) + T_KILVIN_ZERO;
				dens_arg[2] = Index;
				coef = gravity_constant * FluidProp->Density(dens_arg);
			}
			coef = gravity_constant * FluidProp->Density();
			if (dim == 3 && ele_dim == 2)
			{
				vel[dim - 1] += coef;  // NW local permeability tensor is
				                       // already transformed to global one in
				                       // CalCoefLaplace()
				if (PcsType == V || PcsType == P)
				{
					auto& tensor(*MeshElement->getTransformTensor());
					for (size_t i = 0; i < dim; i++)
						for (size_t j = 0; j < ele_dim; j++)
						{
							if (PcsType == V)
								vel_g[i] +=
								    rho_g * gravity_constant *
								    tensor(i, k) *
								    tensor(2, k);
							if (PcsType == P)  // PCH 05.2009
								vel_g[i] +=
								    coef * GasProp->Density() /
								    FluidProp->Density() *
								    tensor(i, k) *
								    tensor(2, k);
						}
				}
			}  // To be correctted
			else
			{
				if (PcsType == V)
				{
					vel[dim - 1] += coef;
					vel_g[dim - 1] += gravity_constant * rho_ga;
				}
				else if (PcsType == P)  // PCH 05.2009
				{
					// vel[dim-1] -= coef;
					// CB_merge_0513 ?? gravity term
					vel[dim - 1] += coef;  // CB I think this should be added
					vel_g[dim - 1] += gravity_constant * GasProp->Density();
				}
				else
					vel[dim - 1] += coef;
			}
		}
		//
		if (PcsType == V)
		{
			for (size_t i = 0; i < dim; i++)  // 02.2010. WW
			{
				for (size_t j = 0; j < dim; j++)
					tmp_gp_velocity(i, gp) +=
					    mat[dim * i + j] * vel[j] / time_unit_factor;
				// gp_ele->Velocity(i, gp) +=
				// mat[dim*i+j]*vel[j]/time_unit_factor;
			}
			CalCoefLaplace2(true, 3);
			double coef_tmp;  // WX:08.2010.
			coef_tmp = rhow / rho_ga;
			for (size_t i = 0; i < dim; i++)
				for (size_t j = 0; j < dim; j++)
					gp_ele->Velocity_g(i, gp) -= coef_tmp * mat[dim * i + j] *
					                             vel_g[j] / time_unit_factor;
			// WX:modified.08.2010
		}
		else  // 02.2010. WW
		{
			for (size_t i = 0; i < dim; i++)
			{
				for (size_t j = 0; j < dim; j++)
					//              gp_ele->Velocity(i, gp) -=
					//              mat[dim*i+j]*vel[j];  // unit as that given
					//              in input file
					// SI Unit
					tmp_gp_velocity(i, gp) -=
					    mat[dim * i + j] * vel[j] / time_unit_factor;
				// gp_ele->Velocity(i, gp) -=
				// mat[dim*i+j]*vel[j]/time_unit_factor;
			}
		}
		if (PcsType == P)  // PCH 05.2009
		{
			// Juse use the coefficient of PSGLOBAL Pressure-based velocity (4)
			CalCoefLaplacePSGLOBAL(true, 4);
			for (size_t i = 0; i < dim; i++)
				for (size_t j = 0; j < dim; j++)
					gp_ele->Velocity_g(i, gp) -=
					    mat[dim * i + j] * vel_g[j] / time_unit_factor;
		}
		//
	}
	gp_ele->Velocity = tmp_gp_velocity;
//
#if 0
	if(pcs->Write_Matrix)
	{
		(*pcs->matrix_file) << "### Element: " << Index << "\n";
		(*pcs->matrix_file) << "---Velocity of water " << "\n";
		gp_ele->Velocity.Write(*pcs->matrix_file);
		if(gp_ele->Velocity_g.Size() > 0)
		{
			(*pcs->matrix_file) << "---Velocity of gas " << "\n";
			gp_ele->Velocity_g.Write(*pcs->matrix_file);
		}
	}
	// gp_ele->Velocity.Write();
#endif
}

/***************************************************************************
   GeoSys - Funktion: Cal_GP_Velocity_FM
   CFiniteElementStd:: Velocity calulation in gauss points from
   node velocities obtained by fluid momentum for one element

   Programming:  SB
   09/2009	first version
 **************************************************************************/
void CFiniteElementStd::Cal_GP_Velocity_FM(int* i_ind)
{
	int i;
	/* //WW
	   static double vel_g_old[3]=
	   {
	   0.0,0.0,0.0
	   }
	 */
	double vel_g[3] = {0.0, 0.0, 0.0};
	// ---- Gauss integral
	// WW int gp_r=0, gp_s=0, gp_t=0;
	// WW double fkt=0.0;                             //OK411 coef = 0.0
	int i_idx;
	// Get fluid_momentum process
	CRFProcess* m_pcs_fm = PCSGet("FLUID_MOMENTUM");

	ElementValue* gp_ele = ele_gp_value[Index];

	// Gauss point loop
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		// Get gauss point data
		// GetGaussData(gp, gp_r, gp_s, gp_t);
		// WW fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute the shape function for interpolation within element
		ComputeShapefct(1);

		// Save former gp velocity
		// WW  for(i_dim=0;i_dim<dim;i_dim++) vel_g_old[i_dim] =
		// gp_ele->Velocity(i_dim,gp);

		// Interpolate velocity from nodes to gauss point for all three velocity
		// components
		for (size_t i_dim = 0; i_dim < dim; i_dim++)
		{
			// Get  velocities from FLUID_MOMENTUM process in element nodes:
			i_idx = i_ind[i_dim];
			for (i = 0; i < nnodes; i++)
			{
				NodalVal[i] = m_pcs_fm->GetNodeValue(nodes[i], i_idx);
				// dirty fix for permebility to conductivity
				NodalVal[i] = NodalVal[i] / gravity_constant / 1000.0 * 0.001;
			}
			vel_g[i_dim] = interpolate(NodalVal);
		}  // end for dim

		// Set gauss point velocity
		for (size_t i_dim = 0; i_dim < dim; i_dim++)
			gp_ele->Velocity(i_dim, gp) = vel_g[i_dim];

		/*        // Write out differences:
		     if((Index < 100)&&(Index > 0)&&(gp < 3)){
		     cout << " Element: " << Index << ", GP: " << gp << ": ";
		     cout << "vel_fem: " ;
		     for(i_dim=0;i_dim<dim;i_dim++) cout << vel_g_old[i_dim] << "  ";
		   //     cout << "vel_FM: " ;
		   //	  for(i_dim=0;i_dim<dim;i_dim++) cout << vel_g[i_dim] << "  ";
		   //	  cout << "vel_diff: " ;
		   //	  for(i_dim=0;i_dim<dim;i_dim++) cout <<
		   vel_g_old[i_dim]-vel_g[i_dim] << "  ";
		     cout << endl;
		     }
		 */
	}  // end gauss point loop

	// Output
	// gp_ele->Velocity.Write();
}

/***************************************************************************
   GeoSys - Funktion: InterpolatePropertyToGausspoint
   CFiniteElementStd:: necessary for using precalculated density and viscosity
 BG, 11/2010

   Programming:  BG
   11/2010	first version
 **************************************************************************/
double CFiniteElementStd::InterpolatePropertyToGausspoint(int GPIndex,
                                                          CRFProcess* m_pcs,
                                                          int Variableindex)
{
	(void)GPIndex;  // unused
	// double fkt = 0.0;
	// int gp_r=0, gp_s=0, gp_t=0;
	// ElementValue* gp_ele = ele_gp_value[Index];
	int i;
	double variable;
	int size_m = 20;  // assigned to the value in CFiniteElementStd(CRFProcess
	                  // *Pcs, const int C_Sys_Flad, const int order=1);
	double* NodalVal_BG;

	NodalVal_BG = new double[size_m];  // BG
	// Get gauss point data
	// GetGaussData(gp, gp_r, gp_s, gp_t);
	// fkt = GetGaussData(GPIndex, gp_r, gp_s, gp_t);
	// Compute the shape function for interpolation within element
	// ComputeShapefct(1);
	// read density from nodes
	for (i = 0; i < nnodes; i++)
		NodalVal_BG[i] = m_pcs->GetNodeValue(nodes[i], Variableindex);
	// Interpolate density from nodes to gauss point
	variable = interpolate(NodalVal_BG);

	return variable;
}

/***************************************************************************
   GeoSys - Funktion: Get_Element_Velocity
   CFiniteElementStd:: Velocity calulation in gauss points from
   node velocities obtained by DUMUX or ECLIPSE

   Programming:  BG
   08/2010	first version
 **************************************************************************/
double CFiniteElementStd::Get_Element_Velocity(int Index, CRFProcess* /*m_pcs*/,
                                               int phase_index, int dimension)
{
	ostringstream temp;
	string tempstring;
	double velocity[3];

	ElementValue* gp_ele = ele_gp_value[Index];

	// Gauss point loop
	velocity[dimension] = 0;
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		// for(i_dim = 0; i_dim < dim; i_dim++)
		//{
		// velocity[i_dim] = 0;
		//}

		// Compute the shape function for interpolation within element
		ComputeShapefct(1);

		// Get gp velocity
		if (phase_index == 0)
		{
			// cout << " Water " << Index;
			velocity[dimension] += gp_ele->Velocity(dimension, gp);
			// cout << " " << gp_ele->Velocity(dimension, gp);
		}
		else
		{
			// cout << " Gas " << Index;
			velocity[dimension] += gp_ele->Velocity_g(dimension, gp);
			// cout << " " << gp_ele->Velocity(dimension, gp);
		}
		// for(i_dim = 0; i_dim < dim; i_dim++)
		//{
		//	if (phase_index == 0)
		//		velocity[i_dim] += gp_ele->Velocity(i_dim, gp);
		//	else
		//		velocity[i_dim] += gp_ele->Velocity_g(i_dim, gp);
		//}
	}
	// cout << endl;
	// cout << gp_ele->Velocity(dimension, 0) << " " <<
	// gp_ele->Velocity(dimension, 1) << " " << gp_ele->Velocity(dimension, 2)
	// << " " << gp_ele->Velocity(dimension, 3) << " " <<
	// gp_ele->Velocity(dimension, 4) << " " << gp_ele->Velocity(dimension, 5)
	// << " " << gp_ele->Velocity(dimension, 6) << " " <<
	// gp_ele->Velocity(dimension, 7) << " " << gp_ele->Velocity(dimension, 8)
	// << " " << gp_ele->Velocity(dimension, 9) << " " <<
	// gp_ele->Velocity(dimension, 10) << " " << gp_ele->Velocity(dimension,11)
	// << " " << gp_ele->Velocity(dimension, 12) << " " <<
	// gp_ele->Velocity(dimension, 13) << " " << gp_ele->Velocity(dimension, 14)
	// << endl;

	// Calculate average element velocity
	velocity[dimension] /= nGaussPoints;
	// for(i_dim = 0; i_dim < dim; i_dim++)
	//{
	//	velocity[i_dim] /= nGaussPoints;
	//}

	return velocity[dimension];
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: Assemby_Gravity
   Aufgabe:
           Assemble the contribution of known gradient of hydraulic head or
         pressure and gravity to RHS in Darcy flow
           to the global system

   Programming:
   05/2005   PCH
   09/2005   PCH
 **************************************************************************/
// Local assembly
void CFiniteElementStd::AssembleRHS(int dimension)
{
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp_t = 0;
	double fkt, fktG, rho;

	// Declare two known properties on node
	// Since I declare these variables locally, the object of Vec should handle
	// destruction nicely
	// when this local function is done so that I don't bother with memory leak.

	// Initialize Pressure from the value already computed previously.
	CRFProcess* m_pcs = NULL;
	for (size_t i = 0; i < pcs_vector.size(); ++i)
	{
		m_pcs = pcs_vector[i];
		//		if(m_pcs->pcs_type_name.find("LIQUID_FLOW")!=string::npos) // TF
		if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
		{
			PcsType = L;
			break;
			//		} else if (m_pcs->pcs_type_name.find("RICHARDS_FLOW") !=
			// string::npos) { // TF
		}
		else if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
		{
			PcsType = R;
			break;
			//		} else if (m_pcs->pcs_type_name.find("GROUNDWATER_FLOW") //
			// TF
		}
		else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		{
			PcsType = G;
			break;
		}
	}
	// Update the process for proper coefficient calculation.
	pcs = m_pcs;
	int nidx1;
	//	if (!(m_pcs->pcs_type_name.find("GROUNDWATER_FLOW") != string::npos)) //
	// TF
	if (!(m_pcs->getProcessType() == GROUNDWATER_FLOW))
		nidx1 = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
	else  // then, this is GROUNDWATER_FLOW
	{
		nidx1 = m_pcs->GetNodeValueIndex("HEAD") + 1;
		HEAD_Flag = 1;
		PcsType = G;
	}

	for (int i = 0; i < nnodes; ++i)
	{
		NodalVal[i] = 0.0;
		NodalVal1[i] = m_pcs->GetNodeValue(nodes[i], nidx1);
		NodalVal2[i] = 0.0;
	}

	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determination
		//---------------------------------------------------------
		fktG = fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function

		// Material
		CalCoefLaplace(true);

		// Calculate vector that computes dNj/dx*Ni*Pressure(j)
		// These index are very important.
		rho = FluidProp->Density();
		// Let's get the viscosity too.
		// 09/2010 TF compiler complains about unused value
		//		CFluidProperties *FluidProp = mfp_vector[0];
		if (gravity_constant < MKleinsteZahl)  // HEAD version
			rho = 1.0;
		else if (HEAD_Flag)
			// FS/WW 21.05.2010
			// fkt = fkt*rho * gravity_constant/FluidProp->Viscosity();
			rho = 1.0;
		else
			rho *= gravity_constant;
		//			rho *= gravity_constant/FluidProp->Viscosity();		// This
		//seems
		// to divide viscosity two times. Thus, wrong.

		fktG *= rho;
		for (int i = 0; i < nnodes; i++)
			for (int j = 0; j < nnodes; j++)
				for (size_t k = 0; k < dim; k++)
				{
					NodalVal[i] -= fkt * dshapefct[k * nnodes + j]
					               // NW  dshapefct[dimension*nnodes+j] ->
					               // dshapefct[k*nnodes+j]
					               *
					               mat[dim * dimension + k] * shapefct[i] *
					               NodalVal1[j];
					//	*************************************
					// FS/WW 21.05.2010
					if (HEAD_Flag) continue;
					//***************************************
					NodalVal2[i] += fktG * dshapefct[k * nnodes + j]
					                // NW  dshapefct[dimension*nnodes+j] ->
					                // dshapefct[k*nnodes+j]
					                *
					                mat[dim * dimension + k] * shapefct[i] *
									MeshElement->GetNode(j)->getData()[2];
				}
	}

	// Just influence when it's the gravitational direction in the case of
	// Liquid_Flow
	// Thus, it needs one more switch to tell Liquid_Flow and Groundwater_Flow.
	int IsGroundwaterIntheProcesses = 0;
	for (size_t i = 0; i < pcs_vector.size(); ++i)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->pcs_type_name.find("GROUNDWATER_FLOW") !=
		//string::npos)
		//// TF
		if (m_pcs->getProcessType() == GROUNDWATER_FLOW)
			IsGroundwaterIntheProcesses = 1;
	}

	// Checking the coordinateflag for proper solution.
	int checkZaxis = 0;
	int coordinateflag = pcs->m_msh->GetCoordinateFlag();
	if ((coordinateflag == 12) || (coordinateflag == 22 && dimension == 1) ||
	    (coordinateflag == 32 && dimension == 2))
		checkZaxis = 1;  // Then, this gotta be z axis.

	// Compansate the gravity term along Z direction
	if (checkZaxis && IsGroundwaterIntheProcesses == 0)
		for (int i = 0; i < nnodes; i++)
			NodalVal[i] -= NodalVal2[i];

	// Store the influence into the global vectors.
#if defined(NEW_EQS)
	m_pcs = PCSGet("FLUID_MOMENTUM");
	for (int i = 0; i < nnodes; i++)
	{
		m_pcs->eqs_new->getRHS()[eqs_number[i]] += NodalVal[i];
	}
#endif
	// OK. Let's add gravity term that incorporates the density coupling term.
	// This is convenient. The function is already written in RF.
	// Assemble_Gravity();
}

void CFiniteElementStd::AssembleParabolicEquation()
{
	int i, ii;
	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	double relax0, relax1, pcs_time_step, dt_inverse;
	long dm_shift = 0, cshift = 0;  // WW 05.01.07
	bool H2_mono = false;           // 15. 07.2011. WW
	if (PcsType == V || PcsType == P) H2_mono = true;
	//
	// JT2012: Get the time step of this process! Now dt can be independently
	// controlled.
	pcs_time_step = pcs->Tim->time_step_length;
	dt_inverse = 1.0 / pcs_time_step;  // (also, no need to check minimum. It is
	                                   // handeled in Tim.

	// WW 05.01.07
	relax0 = pcs->m_num->nls_relaxation;  // WW

	relax1 = 1.0;
	if (relax0 < DBL_MIN) relax0 = 1.0;
	relax1 = 1.0 - relax0;
	//
	if (pcs->dof > 1) cshift = NodeShift[pcs->continuum];
	if (pcs->type / 10 == 4) dm_shift = problem_dimension_dm;
	//----------------------------------------------------------------------
	double fac1, fac2;
	//----------------------------------------------------------------------
	// Initialize.
	// if (pcs->Memory_Type==2) skip the these initialization
	if (H2_mono)
		(*Mass2) = 0.0;
	else
		(*Mass) = 0.0;
	(*Laplace) = 0.0;
	//----------------------------------------------------------------------
	// GEO
	// double geo_fac = MediaProp->geo_area;
	//----------------------------------------------------------------------
	// Calculate matrices
	// Mass matrix..........................................................
	if (this->pcs->tim_type == FiniteElement::TIM_TRANSIENT)
	{
		if (PcsType == V)  // WW
		{
			if (pcs->m_num->ele_mass_lumping)
				CalcLumpedMass2();
			else
				CalcMass2();
		}
		else if (PcsType == P)  // PCH
		{
			if (pcs->m_num->ele_mass_lumping)
				CalcLumpedMassPSGLOBAL();
			else
				CalcMassPSGLOBAL();
		}
		else
		{
			if (pcs->m_num->ele_mass_lumping)
				CalcLumpedMass();
			else
				CalcMass();
		}
	}

	// Laplace matrix.......................................................
	CalcLaplace();
	if (RD_Flag)  // YD /WW
		Assemble_DualTransfer();
	if (pcs->Tim->time_control_type == TimeControlType::NEUMANN)
		pcs->timebuffer /= mat[0];  // YD
	//======================================================================
	// Assemble global matrix
	//----------------------------------------------------------------------
	// Assemble local left matrix:
	// [C]/dt + theta [K] non_linear_function for static problems
	{
		fac1 = dt_inverse;
		fac2 = relax0;  // unterrelaxation WW theta* non_linear_function_iter;
		                // //*geo_fac;
	}

	// Mass matrix
	if (pcs->PartialPS != 1)  // PCH if not partial-pressure-based
	{
		if (H2_mono)
			*StiffMatrix = *Mass2;
		else
			*StiffMatrix = *Mass;
		(*StiffMatrix) *= fac1;
	}
	// Laplace matrix
	// PCH to reduce PDE to ODE in Saturation model
	// PCH: If equation 2 in Two-phase flow.
	if (pcs->pcs_type_number == 1 && pcs->ML_Cap != 0)
	{  // then, Laplace equation is no need. Only solve for ODE
	}
	else
	{
		*AuxMatrix = *Laplace;
	}

	(*AuxMatrix) *= fac2;
	*StiffMatrix += *AuxMatrix;

	// if(index >= 199){
	// Mass->Write();
	// Laplace->Write();
	// StiffMatrix->Write();
	//}
	//----------------------------------------------------------------------
	// Add local matrix to global matrix
	/// Initialize temporal vector
	for (i = 0; i < nnodes; i++)
		NodalVal[i] = 0.0;
	if (PcsType == V)  // For DOF>1: 27.2.2007 WW
		for (i = 0; i < nnodes; i++)
			NodalVal[i + nnodes] = 0.0;
	if (pcs->m_num->nls_method == FiniteElement::NL_NEWTON)
		StiffMatrix->multi(NodalVal1, NodalVal, -1.0);

	//{
	//----------------------------------------------------------------------
	// Add local matrix to global matrix
	// add2GlobalMatrixII(); //TN - added again 07/2013
	//}
	//======================================================================
	// Assemble local RHS vector:
	// ( [C]/dt - (1.0-theta) [K] non_linear_function ) u0  for static problems
	{
		fac1 = dt_inverse;
		fac2 = relax1;  // Unerrelaxation. WW  (1.0-theta) *
		                // non_linear_function_t0; //*geo_fac;
	}

	// Mass - Storage
	if (pcs->PartialPS != 1)  // PCH if not partial-pressure-based
	{
		if (PcsType == V || PcsType == P)  // PCH
			*AuxMatrix1 = *Mass2;
		else
			*AuxMatrix1 = *Mass;
		(*AuxMatrix1) *= fac1;
	}
	// Laplace - Diffusion
	// Laplace matrix
	// PCH: If equation 2 in Two-phase flow.
	if (pcs->pcs_type_number == 1 && pcs->ML_Cap != 0)
	{  // then, Laplace equation is no need. Only solve for ODE
	}
	else
	{
		*AuxMatrix = *Laplace;
	}
	(*AuxMatrix) *= fac2;
	*AuxMatrix1 -= *AuxMatrix;
	// 07.01.07 WW
	int idx = idx0;
	if (pcs->continuum == 1) idx = idxp20;
	for (i = 0; i < nnodes; i++)
		NodalVal0[i] = pcs->GetNodeValue(nodes[i], idx);

	if (PcsType == P)  // For DOF>1:

		for (i = 0; i < nnodes; i++)
		{
			NodalVal[i] = 0.0;
			NodalVal0[i + nnodes] = pcs->GetNodeValue(nodes[i], idxSn0);
			NodalVal[i + nnodes] = 0.0;
		}
	AuxMatrix1->multi(NodalVal0, NodalVal);
	// PCH: Type III (Cauchy boundary conditions) in need, it should be added
	// here.

	//
	if (H2_mono)
	{
		int nDF = 2;
		for (ii = 0; ii < nDF; ii++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			i_sh = NodeShift[ii + dm_shift];
#endif
			ii_sh = ii * nnodes;
			for (i = 0; i < nnodes; i++)
			{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
				eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i + ii_sh];
#endif
				(*RHS)(i + LocalShift + ii_sh) += NodalVal[i + ii_sh];
			}
		}
	}
	else
	{
		cshift += NodeShift[dm_shift];
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[cshift + eqs_number[i]] += NodalVal[i];
#endif
			(*RHS)(i + LocalShift) += NodalVal[i];
		}
	}
	//
}
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
                        /*!
                           \brief Add the local stiff matrix to the global one

                           04.2012. WW
                         */
//------------------------------------------------------
void CFiniteElementStd::add2GlobalMatrixII(bool updateA, bool updateRHS)
{
	if (pcs->m_num->petsc_split_fields)
	{
		add2GlobalMatrixII_Split(updateA, updateRHS);
		return;
	}

	const int dof = (pcs->GetContinnumType() == 1) ? 1 : pcs->pcs_number_of_primary_nvals;

	const int m_dim = nnodes * dof;
	const int n_dim = m_dim;
	double* local_matrix = StiffMatrix->getEntryArray();
	double* local_vec = RHS->getEntryArray();

	for (int i = 0; i < nnodes; i++)
	{
		MeshLib::CNode* node = MeshElement->GetNode(i);
		const int offset = node->GetEquationIndex() * dof;
		const bool isGhost = !pcs->m_msh->isNodeLocal(node->GetIndex());
		for (int k = 0; k < dof; k++)
		{
			const int ki = k * nnodes + i;
			row_ids[ki] = isGhost ? -1 : (offset + k);
			col_ids[ki] = offset + k;
		}
	}

//	std::stringstream ss;
//	ss  << "\nElement: local=" << MeshElement->GetIndex() << ", global=" << MeshElement->GetGlobalIndex();
//	ss  << "\nrow_ids = ";
//	for (int i=0; i<nnodes*dof; i++)
//		ss << row_ids[i] << " ";
//	ss  << "\ncol_ids = ";
//	for (int i=0; i<nnodes*dof; i++)
//		ss << col_ids[i] << " ";
//	ScreenMessage2("%s\n", ss.str().data());

	petsc_group::PETScLinearSolver* eqs = pcs->eqs_new;
	if (updateA)
		eqs->addMatrixEntries(m_dim, row_ids, n_dim, col_ids, local_matrix);
	if (updateRHS)
		eqs->setArrayValues(1, m_dim, row_ids, local_vec);
}

void CFiniteElementStd::add2GlobalMatrixII_Split(bool updateA, bool updateRHS)
{
	assert(pcs->m_num->petsc_split_fields);

	const int dof = pcs->pcs_number_of_primary_nvals;
	const int c_nnodes = nnodes;

	petsc_group::PETScLinearSolver* eqs = pcs->eqs_new;
	const int m_dim = nnodes;
	const int n_dim = nnodes;
	for (int i = 0; i < n_dim; i++)
		col_ids[i] = MeshElement->GetNode(i)->GetEquationIndex();
	for (int i = 0; i < m_dim; i++)
		row_ids[i] = MeshElement->GetNode(i)->GetEquationIndex();
	double const* const loc_cpl_mat = StiffMatrix->getEntryArray();
	double const* const loc_cpl_rhs = RHS->getEntryArray();
	const unsigned n_cpl_mat_columns = nnodes * dof;

#if 0
	if (myrank==1 && hasGhostNodes) {
		bool found = false;
		for (int i=0; i<nnodes; i++)
			if (MeshElement->GetNode(i)->GetEquationIndex()==292) found = true;
		if (found) {
			std::cout << "-> Index: " << Index << "\n";
			std::cout << "-> Nodes: \n";
			for (int i=0; i<nnodes; i++)
				std::cout << MeshElement->GetNode(i)->GetEquationIndex() << " ";
			std::cout << "\n";
			std::cout << "-> u0,u1: \n";
			for (int i=0; i<nnodes; i++) {
				std::cout << NodalVal_p0[i] << " "<< NodalVal_p1[i] << " "<< NodalVal_T0[i] << " "<< NodalVal_T1[i] << "\n";
			}
			std::cout << "\n";
			std::cout << "-> Velocity: \n";
			ele_gp_value[Index]->Velocity.Write();
			std::cout << "-> StiffMatrix: \n";
			StiffMatrix->Write(std::cout);
			std::cout << "-> RHS: \n";
			RHS->Write(std::cout);
			std::cout << "-> act_nodes =" << act_nodes << "\n";
			std::cout << "-> local_idx: \n";
			for (int i = 0; i < act_nodes; i++)
				std::cout << local_idx[i] << " ";
			std::cout << "\n";
			std::cout << "-> idxm: \n";
			for (int i = 0; i < m_dim; i++)
				std::cout << idxm[i] << " ";
			std::cout << "\n";
			for (int ii=0; ii<dof; ii++)
			{
//				ii= 1;
				const int ii_offset = ii*n_dim*n_cpl_mat_columns;
#if 0
				for (int jj=0; jj<dof; jj++)
				{
//					jj = 0;
					const int jj_offset = jj*nnodes;
					std::cout << "sub mat " << ii << jj << ":\n";
					for (int i = 0; i < m_dim; i++)
					{
						const int i_offest = ii_offset + local_idx[i]*nnodes*dof + jj_offset;
						for (int j = 0; j < nnodes; j++) {
							assert(i_offest + j < StiffMatrix->Rows()*StiffMatrix->Cols());
							std::cout << loc_cpl_mat[i_offest + j] << " ";
						}
						std::cout << "\n";
					}
					std::cout << "\n";
//					break;
				}
#endif
				std::cout << "sub rhs " << ii << ":\n";
				for (int i = 0; i < m_dim; i++)
					std::cout << loc_cpl_rhs[local_idx[i]+ ii*nnodes] << " ";
				std::cout << "\n";
//				break;
			}
			std::cout << "\n";
		}
	}
#endif

	//	if (myrank==0) { // updateRHS
	//		std::cout << "# Element " << Index << "\n";
	//		StiffMatrix->Write();
	//		RHS->Write();
	//	}

	PetscInt ierr;
	for (int ii = 0; ii < dof; ii++)
	{
		if (updateA)
		{
			const int ii_offset = ii * n_dim * n_cpl_mat_columns;
			for (int jj = 0; jj < dof; jj++)
			{
				(*Laplace) = .0;
				double* local_matrix =
				    Laplace->getEntryArray();  // Temporary use
				const int jj_offset = jj * c_nnodes;
				for (int i = 0; i < m_dim; i++)
				{
					const int i_offest =
						ii_offset + i * c_nnodes * dof + jj_offset;
					for (int j = 0; j < c_nnodes; j++)
					{
						// assert(i_offest + j <
						// StiffMatrix->Rows()*StiffMatrix->Cols());
						local_matrix[i * c_nnodes + j] =
						    loc_cpl_mat[i_offest + j];
					}
				}
				ierr = MatSetValues(eqs->vec_subA[ii * dof + jj], m_dim, row_ids,
									n_dim, col_ids, local_matrix, ADD_VALUES);
				CHKERRABORT(PETSC_COMM_WORLD, ierr);
			}
		}

		if (updateRHS)
		{
			double* local_vec = NodalVal;
			for (int i = 0; i < c_nnodes; i++)
				local_vec[i] = .0;
			for (int i = 0; i < m_dim; i++)
				local_vec[i] = loc_cpl_rhs[i + ii * c_nnodes];

			ierr = VecSetValues(eqs->vec_subRHS[ii], m_dim, row_ids, local_vec,
			                    ADD_VALUES);
			CHKERRABORT(PETSC_COMM_WORLD, ierr);
		}  // update RHS
	}
}

#else
//------------------------------------------------------
/*!
   \brief Add the local stiff matrix to the global one

   22.06.2011. WW
 */
//------------------------------------------------------
void CFiniteElementStd::add2GlobalMatrixII(const int block_cols)
{
	if (!add2global) return;

	const int cshift = (pcs->dof > 1) ? NodeShift[pcs->continuum] : 0;
	const int dm_shift = (pcs->type / 10 == 4) ? problem_dimension_dm : 0;

	CSparseMatrix* A = pcs->eqs_new->getA();

	// For DOF>1:
	if (PcsType == V || PcsType == P || PcsType == TH)
	{
		const int nDF = 2;
		for (int ii = 0; ii < nDF; ii++)
		{
			const int i_sh = NodeShift[ii + dm_shift];
			const int ii_sh = ii * nnodes;
			for (int jj = 0; jj < block_cols; jj++)
			{
				const int j_sh = NodeShift[jj + dm_shift];
				const int jj_sh = jj * nnodes;
				for (int i = 0; i < nnodes; i++)
				{
					const int kk = i_sh + eqs_number[i];  // 02.2011. WW
					for (int j = 0; j < nnodes; j++)
					{
						(*A)(kk, j_sh + eqs_number[j]) +=
						    (*StiffMatrix)(i + ii_sh, j + jj_sh);
					}
				}
			}
		}
	}
	else
	{
		const int cshift_dm = cshift + NodeShift[dm_shift];  // WW 05.01.07
		for (int i = 0; i < nnodes; i++)
		{
			const int kk = cshift_dm + eqs_number[i];  // 02.2011. WW
			for (int j = 0; j < nnodes; j++)
			{
				(*A)(kk, cshift_dm + eqs_number[j]) += (*StiffMatrix)(i, j);
			}
		}
	}

	if (pcs->matrix_file)
	{
		(*pcs->matrix_file) << "Stiffness: " << endl;
		StiffMatrix->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << endl;
	}
}
#endif
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2010 NW Implementation
   05/2013 NW Support MPI with PETSc
**************************************************************************/
void CFiniteElementStd::CalcFEM_FCT()
{
	const double dt_inverse = 1.0 / dt;
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;  // WW
	A = pcs->eqs_new->getA();
#endif

	//----------------------------------------------------------------------
	// Construct lumped mass matrix
	//----------------------------------------------------------------------
	// assemble local matrix
	(*FCT_MassL) = 0.0;
	for (int i = 0; i < nnodes; i++)
		for (int j = 0; j < nnodes; j++)
			(*FCT_MassL)(i) += (*Mass)(i, j);
	// add into a global diagonal vector
	Math_Group::Vector* ML = this->pcs->Gl_ML;
	for (int i = 0; i < nnodes; i++)
	{
#ifdef USE_PETSC
		long node_i_id = MeshElement->GetNode(i)->GetEquationIndex();
#else
		long node_i_id = this->MeshElement->GetNodeIndex(i);
#endif
		(*ML)(node_i_id) += (*FCT_MassL)(i);
	}
	//----------------------------------------------------------------------
	// Initialize FCT flux with consistent mass matrix: f_ij = m_ij
	//----------------------------------------------------------------------
	Math_Group::SparseMatrixDOK* FCT_Flux = this->pcs->FCT_AFlux;
	for (int i = 0; i < nnodes; i++)
	{
		long node_i_id = this->MeshElement->GetNodeIndex(i);
		//    for (j=i; j<nnodes; j++) {
		for (int j = i + 1; j < nnodes; j++)  // symmetric
		{
			long node_j_id = this->MeshElement->GetNodeIndex(j);
			double v = (*this->Mass)(i, j);
#ifdef USE_PETSC
			if (v == .0) v = (*this->Mass)(j, i);  // look for inner nodes
#endif
			(*FCT_Flux)(node_i_id, node_j_id) += v;
			(*FCT_Flux)(node_j_id, node_i_id) += v;
		}
	}

	//----------------------------------------------------------------------
	// calculate transport operator K
	//----------------------------------------------------------------------
	// local K
	*AuxMatrix = *Laplace;
	*AuxMatrix += *Advection;
	*AuxMatrix += *Storage;

#ifdef USE_PETSC
	// store K (global)
	for (int i = 0; i < nnodes; i++)
	{
		long glob_i = MeshElement->GetNode(i)->GetEquationIndex();
		for (int j = 0; j < nnodes; j++)
		{
			long glob_j = MeshElement->GetNode(j)->GetEquationIndex();
			(*this->pcs->FCT_K)(glob_i, glob_j) += (*AuxMatrix)(i, j);
		}
	}
#else
	// Add K matrix to a global coefficient matrix
	for (int i = 0; i < nnodes; i++)
	{
		for (int j = 0; j < nnodes; j++)
		{
#ifdef NEW_EQS
			(*A)(NodeShift[problem_dimension_dm] + eqs_number[i],
			     NodeShift[problem_dimension_dm] + eqs_number[j]) +=
			    (*AuxMatrix)(i, j);
#endif
		}
	}
#endif

//----------------------------------------------------------------------
// Setup local coefficient matrix and RHS vector
//----------------------------------------------------------------------
#ifdef USE_PETSC
	const double theta = pcs->m_num->ls_theta;
	// A=1/dt*ML + theta*K
	*AuxMatrix *= theta;
	*StiffMatrix = *FCT_MassL;
	*StiffMatrix *= dt_inverse;
	*StiffMatrix += *AuxMatrix;  // StiffMatrix is later added into a global
	                             // matrix in add2GlobalMatrixII()

	// rhs=(1/dt*ML-(1-theta)*K)u^n
	*AuxMatrix1 = *FCT_MassL;
	*AuxMatrix1 *= dt_inverse;
	*AuxMatrix = *Laplace;
	*AuxMatrix += *Advection;
	*AuxMatrix += *Storage;
	*AuxMatrix *= -(1.0 - theta);
	*AuxMatrix1 += *AuxMatrix;
	for (int i = 0; i < nnodes; i++)
	{
		NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx0);
		NodalVal[i] = 0.0;
	}
	AuxMatrix1->multi(NodalVal1, NodalVal);
	for (int i = 0; i < nnodes; i++)
		(*RHS)(i + LocalShift) += NodalVal[i];  // RHS is later added into a
                                                // global RHS in
                                                // add2GlobalMatrixII()

#else
	// assemble part of RHS: b_i += 1/dt * ml_i * u_i^n
	double fac_mass = dt_inverse;  //*geo_fac;
	for (int i = 0; i < nnodes; i++)
		NodalVal[i] =
		    fac_mass * (*FCT_MassL)(i)*pcs->GetNodeValue(nodes[i], idx0);
	for (int i = 0; i < nnodes; i++)
	{
		eqs_rhs[NodeShift[problem_dimension_dm] + eqs_number[i]] += NodalVal[i];
		(*RHS)(i + LocalShift) += NodalVal[i];
	}
#endif
}

void CFiniteElementStd::AssembleMixedHyperbolicParabolicEquation()
{
	int i;
#ifndef USE_PETSC
	int j;
#endif
	double pcs_time_step, dt_inverse;
	ElementMatrix* EleMat = NULL;  // SB-3
	// NUM
	double theta = pcs->m_num->ls_theta;  // OK
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;  // WW
	A = pcs->eqs_new->getA();
#endif

	// JT2012: Get the time step of this process! Now dt can be independently
	// controlled
	pcs_time_step = pcs->Tim->time_step_length;
	dt_inverse = 1.0 / pcs_time_step;  // (also, no need to check minimum. It is
	                                   // handeled in Tim.
	//
	//----------------------------------------------------------------------
	unit[0] = unit[1] = unit[2] = 0.0;
	// Non-linearities
	//  double non_linear_function_iter = 1.0; //OK
	//  MediaProp->NonlinearFlowFunction(Index,unit,theta);
	//  double non_linear_function_t0   = 1.0; //OK
	//  MediaProp->NonlinearFlowFunction(Index,unit,0.0);
	double fac_mass, fac_laplace, fac_advection;
#ifndef OGS_ONLY_TH
	double fac_storage, fac_content;
#endif
	// if(((aktueller_zeitschritt==1)||(pcs->tim_type_name.compare("TRANSIENT")==0))){
	// //SB-3
	// SB-3
	if (((aktueller_zeitschritt == 1) || (pcs->Memory_Type == 0)))
	{
		// Initialize.
		(*Mass) = 0.0;
		(*Laplace) = 0.0;
		(*Advection) = 0.0;
#ifndef OGS_ONLY_TH
		(*Storage) = 0.0;
		(*Content) = 0.0;
#endif
		//----------------------------------------------------------------------
		// GEO
		// double geo_fac = MediaProp->geo_area;
		//----------------------------------------------------------------------
		// Calculate matrices
		// Mass matrix..........................................................
		// NW
		if (this->pcs->tim_type == FiniteElement::TIM_TRANSIENT)
		{
			if (pcs->m_num->ele_mass_lumping)
				CalcLumpedMass();
			else
				CalcMass();
		}
		// Laplace matrix.......................................................
		CalcLaplace();
		// Advection matrix.....................................................
		CalcAdvection();
#ifndef OGS_ONLY_TH
		// Calc Storage Matrix for decay
		CalcStorage();
		// Calc Content Matrix for  saturation changes
		CalcContent();
#endif

		// Store matrices to memory for steady state element matrices     //SB-3
		if (pcs->Memory_Type > 0)
		{
			EleMat = pcs->Ele_Matrices[Index];
			EleMat->SetMass_notsym(Mass);
			EleMat->SetLaplace(Laplace);
			EleMat->SetAdvection(Advection);
#ifndef OGS_ONLY_TH
			EleMat->SetStorage(Storage);
			EleMat->SetContent(Content);
#endif
		}
	}  // SB-3
	else
	{
		if (Index < 1)
			cout << "        Skipping calculation of element matrices " << endl;
		// Get Element Matrices
		EleMat = pcs->Ele_Matrices[Index];
		Mass = EleMat->GetMass_notsym();
		Laplace = EleMat->GetLaplace();
		Advection = EleMat->GetAdvection();
#ifndef OGS_ONLY_TH
		Storage = EleMat->GetStorage();
		Content = EleMat->GetContent();
#endif
	}  // pcs->tim_type    //SB-3
	//======================================================================
	// Assemble global matrix
	//----------------------------------------------------------------------
	// Assemble local left matrix:
	// [C]/dt + theta [K] non_linear_function for static problems

	fac_mass = dt_inverse;  //*geo_fac;
	fac_laplace = theta;    //* non_linear_function_iter; //*geo_fac;
	fac_advection = theta;
#ifndef OGS_ONLY_TH
	fac_storage = theta;
	fac_content = theta * dt_inverse;
#endif

	if (this->pcs->femFCTmode)  // NW

		this->CalcFEM_FCT();
	else
	{
		// Mass matrix
		*StiffMatrix = *Mass;
		(*StiffMatrix) *= fac_mass;
		// Laplace matrix
		*AuxMatrix = *Laplace;
		(*AuxMatrix) *= fac_laplace;
		*StiffMatrix += *AuxMatrix;
		// Advection matrix
		*AuxMatrix = *Advection;
		(*AuxMatrix) *= fac_advection;
		*StiffMatrix += *AuxMatrix;
#ifndef OGS_ONLY_TH
		// Storage matrix
		*AuxMatrix = *Storage;
		(*AuxMatrix) *= fac_storage;
		*StiffMatrix += *AuxMatrix;
// Content matrix
//*AuxMatrix      = *Content;		//SB, BG; Korrektur Stofftransport bei
// Mehrphasenströmung
//(*AuxMatrix)   *= fac_content;
//*StiffMatrix   += *AuxMatrix; // SB, BG
#endif
		if (pcs->matrix_file)
		{
			(*pcs->matrix_file) << "Stiffness: " << endl;
			StiffMatrix->Write(*pcs->matrix_file);
			(*pcs->matrix_file) << endl;
		}

#ifdef NEW_EQS
		//----------------------------------------------------------------------
		// Add local matrix to global matrix
		if (add2global)
			for (i = 0; i < nnodes; i++)
			{
				for (j = 0; j < nnodes; j++)
				{
					(*A)(NodeShift[problem_dimension_dm] + eqs_number[i],
					     NodeShift[problem_dimension_dm] + eqs_number[j]) +=
					    (*StiffMatrix)(i, j);
				}
			}
#endif
		//======================================================================
		// Assemble local RHS vector:
		// ( [C]/dt - (1.0-theta) [K] non_linear_function ) u0  for static
		// problems

		fac_mass = dt_inverse;         //*geo_fac;
		fac_laplace = -(1.0 - theta);  // * non_linear_function_t0; //*geo_fac;
		fac_advection = -(1.0 - theta);
#ifndef OGS_ONLY_TH
		fac_storage = -(1.0 - theta);  //*lambda
		fac_content = -(1.0 - theta) * dt_inverse;
#endif

		// Mass - Storage
		*AuxMatrix1 = *Mass;
		(*AuxMatrix1) *= fac_mass;
		// Laplace - Diffusion
		*AuxMatrix = *Laplace;
		(*AuxMatrix) *= fac_laplace;
		*AuxMatrix1 += *AuxMatrix;
		// Advection
		*AuxMatrix = *Advection;
		(*AuxMatrix) *= fac_advection;
		*AuxMatrix1 += *AuxMatrix;
#ifndef OGS_ONLY_TH
		// Storage
		*AuxMatrix = *Storage;
		(*AuxMatrix) *= fac_storage;
		*AuxMatrix1 += *AuxMatrix;
		// Content
		*AuxMatrix = *Content;
		(*AuxMatrix) *= fac_content;
		*AuxMatrix1 += *AuxMatrix;
#endif

		for (i = 0; i < nnodes; i++)
		{
			NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx0);
			NodalVal[i] = 0.0;
		}
		AuxMatrix1->multi(
		    NodalVal1,
		    NodalVal);  // AuxMatrix1 times vector NodalVal1 = NodalVal
		//----------------------------------------------------------------------
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			if (add2global)
				eqs_rhs[NodeShift[problem_dimension_dm] + eqs_number[i]] +=
				    NodalVal[i];
#endif
			(*RHS)(i + LocalShift) += NodalVal[i];
		}
	}  // end: femFCTmode
	   //----------------------------------------------------------------------
	   // Debug output
	   /*
	      if(Index < 10){
	      cout << " Element Number " << Index << endl;
	      cout << " Mass matrix" << endl;
	      Mass->Write();
	      cout << " Advection matrix" << endl;
	      Advection->Write();
	      cout << " Dispersion matrix" << endl;
	      Laplace->Write();
	      cout << " Storage matrix" << endl;
	      Storage->Write();
	      cout << " Content matrix" << endl;
	      Content->Write();
	      cout << " Left matrix" << endl;
	      StiffMatrix->Write();
	      cout << " Right matrix" << endl;
	      AuxMatrix1->Write();
	      cout << "RHS: " << endl ;
	      for (i=0;i<nnodes; i++) cout << "| " << NodalVal[i] << " |" << endl;
	      cout << " initial concentrations" << endl;
	      for (i=0;i<nnodes; i++) cout << "| " << NodalVal1[i] << " |" << endl;
	      //	cout << " RHS vector: " << endl;
	      //	for (i=0;i<nnodes; i++) cout << "| " <<    (double)(*RHS)(i+LocalShift)
	      << " |" << endl;
	      }
	    */
}

void CFiniteElementStd::Assemble_totalStressCPL(const int phase)
{
	if (this->pcs->tim_type == FiniteElement::TIM_STEADY)
		return;

	bool updateRHS = false;
	if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION)
	{
		updateRHS = true;
	}
	else if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW
			 && FiniteElement::isNewtonKind(dm_pcs->m_num->nls_method))
	{
		updateRHS = true;
	}

	ElementValue_DM const* ev_dm = ele_value_dm[Index];
	for (int k = 0; k < nnodes; k++)
		NodalVal[k] = 0.0;

	// r += -Np^T biot/K*d(TotalStress_v)/dt
	for (int gp = 0; gp < nGaussPoints; gp++)
	{
		int gp_r, gp_s, gp_t;
		double fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		ComputeShapefct(1);

		//
		//double gp_p = interpolate(NodalVal1);

		// increment of volumetric total stress
		double dSv = MeanStress(*ev_dm->dTotalStress, gp);
		//dSv -= SolidProp->biot_const * gp_p;

		// alpha/K*dSv/dt
		fkt *= - SolidProp->biot_const / SolidProp->K * dSv / dt;

		//
		for (int k = 0; k < nnodes; k++)
		{
			(*RHS)(k) += shapefct[k] * fkt;
			NodalVal[k] += shapefct[k] * fkt;
		}
	}

	if (updateRHS)
	{
#ifdef NEW_EQS
		// add to global RHS
		int shift_index = 0;
		if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
			shift_index = problem_dimension_dm + phase;

		for (int i = 0; i < nnodes; i++)
		{
			eqs_rhs[NodeShift[shift_index] + eqs_number[i]] += NodalVal[i];
		}
#endif
	}

	//TODO MONO
	if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
	{
		assert("Not supported yet.");
		ogsAbort(1);
		//Assemble_strainCPL_Matrix(fac, phase);
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: Assemby_strainCPL
   Aufgabe:
           Assemble local metrices of strain coupling
           to the global system

   Programming:
   01/2005   WW/OK
   05/2005   WW dyn
   07/2005   WW Change due to geometry element object
   06/2011   WW for multi-phase flow
 **************************************************************************/
void CFiniteElementStd::Assemble_strainCPL(const int phase)
{
	if (this->pcs->tim_type == FiniteElement::TIM_STEADY ||
	    this->pcs->disable_deformation_coupling_in_flow)
		return;

	if (pcs->use_total_stress_coupling)
	{
		Assemble_totalStressCPL(phase);
		return;
	}

	bool updateRHS = false;
	if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION)
	{
		updateRHS = true;
	}
	else if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW
			 && FiniteElement::isNewtonKind(dm_pcs->m_num->nls_method))
	{
		updateRHS = true;
	}

#if 0
	if (MediaProp->storage_model == 7)
		fac *= MediaProp->storage_model_values[0];
	else
		fac *= fabs(SolidProp->biot_const);  // WX:11.2012. biot coeff is
	                                         // needed, in some case biot is
	                                         // defined negative
#endif

	//
	CalcStrainCoupling(phase);

	// prepare node values of du
	for (int i = nnodes; i < nnodesHQ; i++)
		nodes[i] = MeshElement->GetNodeIndex(i);

	if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION)
	{
		for (int i = 0; i < nnodesHQ; i++)
		{
			NodalVal2[i] = (dm_pcs->GetNodeValue(nodes[i], Idx_dm1[0]) -
						dm_pcs->GetNodeValue(nodes[i], Idx_dm0[0]));
			NodalVal3[i] = (dm_pcs->GetNodeValue(nodes[i], Idx_dm1[1]) -
						dm_pcs->GetNodeValue(nodes[i], Idx_dm0[1]));
			if (dim == 3)  // 3D.
				NodalVal4[i] = (dm_pcs->GetNodeValue(nodes[i], Idx_dm1[2]) -
							dm_pcs->GetNodeValue(nodes[i], Idx_dm0[2]));
		}
	}
	else if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW
			 && FiniteElement::isNewtonKind(dm_pcs->m_num->nls_method))
	{
		// du is stored in u_0
		for (int i = 0; i < nnodesHQ; i++)
		{
			NodalVal2[i] = pcs->GetNodeValue(nodes[i], Idx_dm0[0]);
			NodalVal3[i] =  pcs->GetNodeValue(nodes[i], Idx_dm0[1]);
			if (dim == 3)
				NodalVal4[i] = pcs->GetNodeValue(nodes[i], Idx_dm0[2]);
		}
	}

#if 0
	if (Index < 5)
	{
		std::cout << "-- u0\n";
		for (int i = 0; i < nnodesHQ; i++)
			std::cout << dm_pcs->GetNodeValue(nodes[i], Idx_dm0[0]) << " ";
		std::cout << "\n";
		for (int i = 0; i < nnodesHQ; i++)
			std::cout << dm_pcs->GetNodeValue(nodes[i], Idx_dm0[1]) << " ";
		std::cout << "\n";

		std::cout << "-- u1\n";
		for (int i = 0; i < nnodesHQ; i++)
			std::cout << dm_pcs->GetNodeValue(nodes[i], Idx_dm1[0]) << " ";
		std::cout << "\n";
		for (int i = 0; i < nnodesHQ; i++)
			std::cout << dm_pcs->GetNodeValue(nodes[i], Idx_dm1[1]) << " ";
		std::cout << "\n";

		std::cout << "-- du\n";
		for (int i = 0; i < nnodesHQ; i++)
			std::cout << NodalVal2[i] << " ";
		std::cout << "\n";
		for (int i = 0; i < nnodesHQ; i++)
			std::cout << NodalVal3[i] << " ";
		std::cout << "\n";
	}
#endif

	double fac = 1.0 / dt;
	fac *= SolidProp->biot_const;

	if (updateRHS)
	{
		// -biot*du/dt
		for (int i = 0; i < nnodesHQ; i++)
		{
			NodalVal2[i] *= -fac;
			NodalVal3[i] *= -fac;
			if (dim == 3)
				NodalVal4[i] *= -fac;
		}

		// N_p^T B * (-biot*du/dt)
		for (int i = 0; i < nnodes; i++)
		{
			NodalVal[i] = 0.0;
			for (int j = 0; j < nnodesHQ; j++)
			{
				NodalVal[i] += (*StrainCoupling)(i, j) * NodalVal2[j];
				NodalVal[i] += (*StrainCoupling)(i, j + nnodesHQ) * NodalVal3[j];
				if (dim == 3)
					NodalVal[i] += (*StrainCoupling)(i, j + 2 * nnodesHQ) * NodalVal4[j];
			}
		}

		// Add to local RHS
		for (int i = 0; i < nnodes; i++)
			(*RHS)(i + LocalShift) += NodalVal[i];

#ifdef NEW_EQS
		// add to global RHS
		int shift_index = 0;
		if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
			shift_index = problem_dimension_dm + phase;

		for (int i = 0; i < nnodes; i++)
		{
			eqs_rhs[NodeShift[shift_index] + eqs_number[i]] += NodalVal[i];
		}
#endif
	}

	if (dm_pcs->getProcessType() == FiniteElement::DEFORMATION_FLOW)
		Assemble_strainCPL_Matrix(fac, phase);
}

//**************************************************************************
/*!
   \brief Assemble the local strain coupling matrix to the golbal one

   28.11.2011 WW
 */
//**************************************************************************
#if defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
void CFiniteElementStd::Assemble_strainCPL_Matrix(const double /*fac*/,
                                                  const int /*phase*/)
{
#else
void CFiniteElementStd::Assemble_strainCPL_Matrix(const double fac,
                                                  const int phase)
{
	// TODO
	int i, j;
	int shift_index;
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;
	A = pcs->eqs_new->getA();
#endif
	// if Richard, StrainCoupling should be multiplied with -1.
	shift_index = problem_dimension_dm + phase;
	for (i = 0; i < nnodes; i++)
	{
		for (j = 0; j < nnodesHQ; j++)
		{
#ifdef NEW_EQS
			(*A)(NodeShift[shift_index] + eqs_number[i],
			     eqs_number[j] + NodeShift[0]) += (*StrainCoupling)(i, j) * fac;
			(*A)(NodeShift[shift_index] + eqs_number[i],
			     eqs_number[j] + NodeShift[1]) +=
			    (*StrainCoupling)(i, j + nnodesHQ) * fac;
			if (problem_dimension_dm == 3)
				(*A)(NodeShift[shift_index] + eqs_number[i],
				     eqs_number[j] + NodeShift[2]) +=
				    (*StrainCoupling)(i, j + 2 * nnodesHQ) * fac;
#endif
		}
	}
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: Assemble local mass matrices to the global system
   Programing:
   05/2005 PCH Implementation
**************************************************************************/
void CFiniteElementStd::AssembleMassMatrix(int option)
{
	// Calculate matrices
	// Mass matrix..........................................................
	// ---- Gauss integral
	int gp;
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt;  // WW ,mat_fac;
// Material
// WW mat_fac = 1.0;

#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;  // PCH
	A = pcs->eqs_new->getA();
#endif
	//----------------------------------------------------------------------
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		// Compute geometry
		ComputeShapefct(1);  // Linear interpolation function

		if (option == 0)  // The consistent method

			// Calculate mass matrix
			for (int i = 0; i < nnodes; i++)
				for (int j = 0; j < nnodes; j++)
					//			if(j>i) continue;
					(*Mass)(i, j) += fkt * shapefct[i] * shapefct[j];

		else if (option == 1)  // The lumped method

			// Calculate mass matrix
			for (int i = 0; i < nnodes; i++)
				for (int j = 0; j < nnodes; j++)
					(*Mass)(i, i) += fkt * shapefct[i] * shapefct[j];
	}

//----------------------------------------------------------------------
// Add local matrix to global matrix
#ifndef USE_PETSC
	if (PcsType == V || PcsType == P)  // For DOF>1: 03.03.2009 PCH
	{
		int ii_sh, jj_sh;
		long i_sh, j_sh = 0;
		for (int ii = 0; ii < pcs->dof; ii++)
		{
			i_sh = NodeShift[ii];
			ii_sh = ii * nnodes;
			for (int jj = 0; jj < pcs->dof; jj++)
			{
				j_sh = NodeShift[jj];
				jj_sh = jj * nnodes;
				for (int i = 0; i < nnodes; i++)
				{
					for (int j = 0; j < nnodes; j++)
					{
#ifdef NEW_EQS
						(*A)(i_sh + eqs_number[i], j_sh + eqs_number[j]) +=
						    (*Mass)(i + ii_sh, j + jj_sh);
#endif
					}
				}
			}
		}
	}
	else
	{
		int cshift = 0;
		// WW 05.01.07
		cshift += NodeShift[problem_dimension_dm];
		for (int i = 0; i < nnodes; i++)
		{
			for (int j = 0; j < nnodes; j++)
			{
#ifdef NEW_EQS
				(*A)(cshift + eqs_number[i], cshift + eqs_number[j]) +=
				    (*Mass)(i, j);
#endif
			}
		}
	}
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: Config material and knowns for local assembly
   Programing:
   08/2008 WW Implementation
**************************************************************************/
void CFiniteElementStd::Config()
{
	int i, nn;
	//----------------------------------------------------------------------
	index = Index;
	//----------------------------------------------------------------------
	nn = nnodes;
	// ?2WW
	// ?2WW
	if (pcs->type / 10 == 4 || pcs->type == 4) nn = nnodesHQ;
//----------------------------------------------------------------------
// For DDC WW
#ifdef NEW_EQS
	eqs_rhs = pcs->eqs_new->getRHS();
#endif

#if defined(NEW_EQS)
	for (i = 0; i < nn; i++)
		eqs_number[i] = MeshElement->GetNode(i)->GetEquationIndex();
#endif
	//----------------------------------------------------------------------
	// Get room in the memory for local matrices
	SetMemory();
	//----------------------------------------------------------------------
	// Set material
	SetMaterial();
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	// ?2WW
	if (D_Flag)
		dm_pcs = (CRFProcessDeformation*)PCSGetDeformation();
	//----------------------------------------------------------------------
	// Initialize RHS
	if (pcs->Memory_Type > 0)
		for (i = LocalShift; (size_t)i < RHS->Size(); i++)
			(*RHS)(i) = 0.0;
	else
		(*RHS) = 0.0;
	//----------------------------------------------------------------------
	// Node value of the previous time step
	int idx00 = idx0;  //----------WW 05.01.07
	int idx11 = idx1;
	if (pcs->GetContinnumType() == 1)
	{
		idx00 = idxp20;
		idx11 = idxp21;
	}
	if (PcsType != TH)
	{
		for (i = 0; i < nnodes; i++)
		{
			NodalVal0[i] = pcs->GetNodeValue(nodes[i], idx00);
			NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx11);
		}  //----------WW 05.01.07
	}
	if (PcsType == V)  // 25.2.2007
		for (i = 0; i < nnodes; i++)
		{
			NodalVal_p2[i] = pcs->GetNodeValue(nodes[i], idxp21);
			NodalVal0[i + nnodes] = pcs->GetNodeValue(nodes[i], idxp20);
			NodalVal1[i + nnodes] = pcs->GetNodeValue(nodes[i], idxp21);
		}
	if (PcsType == P)
		for (i = 0; i < nnodes; i++)
			NodalVal_SatNW[i] = pcs->GetNodeValue(nodes[i], idxSn1);
	//----------WW 05.01.07
	if (cpl_pcs)  // ?2WW: flags are necessary

		for (i = 0; i < nnodes; i++)
		{
			NodalValC[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c0);
			NodalValC1[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c1);
			if (cpl_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || cpl_pcs->type == 42)
			{
				NodalVal_p2[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c1 + 2);
				NodalVal_p20[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c0 + 2);
			}
			if (cpl_pcs->type == 1111)
			{
				NodalVal_T0[i] =
				    cpl_pcs->GetNodeValue(nodes[i], idx_c0 + 2) + T_KILVIN_ZERO;
				NodalVal_T1[i] =
				    cpl_pcs->GetNodeValue(nodes[i], idx_c1 + 2) + T_KILVIN_ZERO;
			}
		}
}
/**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices to the global system
   Programing:
   01/2005 WW Implementation
   02/2005 OK Richards flow
   02/2005 WW Matrix output
   03/2005 WW Heat transport
   04/2005 OK MSH
   05/2005 OK regional PCS
   08/2005 OK Air (gas) flow
   10/2005 OK DDC
   06/2005 WW Adjustment in DDC
   07/2007 WW Nonisothermal multi-phase flow
   10/2007 OK Two-phase flow
   08/2008 WW Extract the configuration of material properties and knowns as
   a single function
**************************************************************************/
void CFiniteElementStd::Assembly(bool updateA, bool updateRHS,
                                 const bool add2global)
{
	this->add2global = add2global;
	int i;
	Config();  // 26.08.2008

	// If output matrices and vectors. 07.2011. WW
	if (pcs->Write_Matrix)
		(*pcs->matrix_file) << "### Element: " << Index << "\n";

	//======================================================================
	switch (PcsType)
	{
		//....................................................................
		case L:  // Liquid flow
			AssembleParabolicEquation();
			Assemble_Gravity();
			Assemble_RHS_LIQUIDFLOW();
			if (dm_pcs)
				Assemble_strainCPL();
			add2GlobalMatrixII();
			break;
		//....................................................................
		// case U: // Unconfined flow  //  part of Groundwater flow mmp keyword
		// ($UNCONFINED)
		//....................................................................
		case G:  // Groundwater flow
			AssembleParabolicEquation();
			// RHS->Write();
			if (dm_pcs) Assemble_strainCPL();
			add2GlobalMatrixII();
			break;
		//....................................................................
		case T:  // Two-phase flow
			if (pcs->pcs_type_number == 0)
			{
				// Start partial-pressure-based model
				pcs->PartialPS = 0;

				AssembleParabolicEquation();

				if (pcs->PartialPS == 1)  // If it is partial-pressure-based
					AssembleRHSVector();
				//			PrintTheSetOfElementMatrices("Pressure1");
				AssembleCapillaryEffect();
				Assemble_Gravity_Multiphase();
			}
			else if (pcs->pcs_type_number == 1)
			{
				// Turn off the partial-pressure-based model for Snw equation
				pcs->PartialPS = 0;

				pcs->ML_Cap = 0;
				AssembleParabolicEquation();
				pcs->ML_Cap = 0;

				AssembleRHSVector();
				Assemble_Gravity_Multiphase();
			}
			add2GlobalMatrixII();
			break;
		//....................................................................
		case C:  // Componental flow
			for (i = 0; i < nnodes; i++)
				NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS);
			break;
		//....................................................................
		case H:                         // Heat transport
			heat_phase_change = false;  // ?2WW
			//  if(SolidProp->GetCapacityModel()==2) // Boiling model
			//    CalNodalEnthalpy();
			// CMCD4213
			AssembleMixedHyperbolicParabolicEquation();
			if (FluidProp->density_model == 14 &&
			    MediaProp->heat_diffusion_model == 1 && cpl_pcs)
				Assemble_RHS_HEAT_TRANSPORT();  // This include when need
			                                    // pressure terms n dp/dt +
			                                    // nv.Nabla p//AKS
			if (MediaProp->evaporation == 647)
				Assemble_RHS_HEAT_TRANSPORT2();  // AKS

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
			add2GlobalMatrixII();
#endif
			break;
		//....................................................................
		case M:  // Mass transport
			// SB4200
			AssembleMixedHyperbolicParabolicEquation();
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
			add2GlobalMatrixII();
#endif
			break;
		//....................................................................
		case R:  // Richards flow
			if (MediaProp->heat_diffusion_model == 1)
				CalcRHS_by_ThermalDiffusion();
			AssembleParabolicEquation();  // OK
			Assemble_Gravity();
			Assemble_RHS_LIQUIDFLOW();  // JM  (thermal expansion fluid)
			if (dm_pcs) Assemble_strainCPL();

			if (pcs->m_num->nls_method ==
			    FiniteElement::NL_NEWTON)  // Newton-Raphson. 07.2011. WW
				ComputeAdditionalJacobi_Richards();
			add2GlobalMatrixII();
			break;
		//....................................................................
		case F:  // Fluid Momentum - Assembly handled in Assembly in
			     // Fluid_Momentum file
			break;
		//....................................................................
		case A:  // Air (gas) flow
			// To account advection like term nv.Nabla p
			AssembleMixedHyperbolicParabolicEquation();
			// AKS
			if (MediaProp->heat_diffusion_model == 1 && cpl_pcs)
				Assemble_RHS_AIR_FLOW();  // n*drho/dt + Nabla.[rho*k/mu rho
                                          // g]//AKS
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
			add2GlobalMatrixII();
#endif
			break;
		case V:
			// TEST                   dm_pcs = NULL;
			// Multi-phase flow 24.02.2007 WW
			AssembleParabolicEquation();
			Assemble_Gravity();
			if (cpl_pcs && MediaProp->heat_diffusion_model == 1)
				Assemble_RHS_T_MPhaseFlow();
			if (dm_pcs) Assemble_RHS_M();
			if (pcs->m_num->nls_method ==
			    FiniteElement::NL_NEWTON)  // Newton-Raphson. 06.2011. WW
			{
				ComputeAdditionalJacobi_H2();

				if (dm_pcs && dm_pcs->type == 42)
				{
					(*StrainCoupling) = 0.0;
					CalcStrainCoupling(0);
					Assemble_strainCPL_Matrix(1.0, 0);  // Phase 0

					(*StrainCoupling) = 0.0;
					CalcStrainCoupling(1);
					Assemble_strainCPL_Matrix(1.0, 1);  // Phase 1
				}
			}
			add2GlobalMatrixII();
			break;

		case P:  // PS_GLOBAL for Multi-phase flow 03.03 2009 PCH
			AssembleParabolicEquation();
			PrintTheSetOfElementMatrices("Laplace");
			if (pcs->num_type_name.find("DirectPc") != string::npos)
			{
				Assemble_RHS_Pc();
				PrintTheSetOfElementMatrices("RHS_Pc");
			}
			Assemble_Gravity();

			if (dm_pcs) Assemble_RHS_M();
			Assemble_RHS_T_PSGlobal();
			add2GlobalMatrixII();
			break;
		case TH:
			//AssembleTHEquation(updateA, updateRHS);
			if (updateA)
				AssembleTHJacobian();
			if (updateRHS)
				AssembleTHResidual();
#ifdef USE_PETSC
			add2GlobalMatrixII(updateA, updateRHS);
#else
			if (updateA)
				add2GlobalMatrixII();
#endif
			break;
		//....................................................................
		default:
			cout << "Fatal error: No valid PCS type" << endl;
			break;
	}

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Output matrices
#if 1
	if (pcs->Write_Matrix)
	{
		std::ostream& os = *pcs->matrix_file;
#else
	if (Index < 3)
	{
		std::ostream& os = std::cout;
		os << "### Element: " << Index << "\n";
#endif
//		os << "---Mass matrix: "
//		   << "\n";
//		if (Mass)
//			Mass->Write(os);
//		else if (Mass2)
//			Mass2->Write(os);
//		os << "---Laplacian matrix: "
//		   << "\n";
//		Laplace->Write(os);
//		if (Advection)
//		{
//			os << "---Advective matrix: "
//			   << "\n";
//			Advection->Write(os);
//		}
//		os << "---RHS: "
//		   << "\n";
//		RHS->Write(os);
//		os << "---U0: "
//		   << "\n";
//		for (int i = 0; i < nnodes; i++)
//			os << "| " << NodalVal0[i] << " | "
//			   << "\n";
//		os << "\n";
		if (updateRHS)
		{
			os << "---U: \n";
			os.setf(std::ios::scientific, std::ios::floatfield);
			os.precision(12);
			for (int i = 0; i < nnodes; i++)
				os << pcs->GetNodeValue(nodes[i], idxp1) << "\n";
			for (int i = 0; i < nnodes; i++)
				os << pcs->GetNodeValue(nodes[i], idxT1) << "\n";
			os << "\n";
			os << "---RHS: "
			   << "\n";
			RHS->Write(os);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices to the global system
   Programing:
   01/2005 WW Implementation
   02/2005 OK Richards flow
   02/2005 WW Matrix output
   03/2005 WW Heat transport
   08/2005 PCH for Fluid_Momentum
   last modification:
**************************************************************************/
void CFiniteElementStd::Assembly(int option, int dimension)
{
	int i, nn;

	nn = nnodes;
	// PCH should check the following line carefully.
	if (pcs->type / 10 == 4 || pcs->type == 4)
		nn = nnodesHQ;

#if defined(NEW_EQS)
	eqs_rhs = pcs->eqs_new->getRHS();
#endif

	for (i = 0; i < nn; i++)
		eqs_number[i] = MeshElement->GetNode(i)->GetEquationIndex();

	// Get room in the memory for local matrices
	SetMemory();

	// Set material
	SetMaterial();

	// Initialize.
	// if (pcs->Memory_Type==2) skip the these initialization
	(*Mass) = 0.0;
	(*Laplace) = 0.0;
	if (pcs->Memory_Type > 0)
		for (i = LocalShift; (size_t)i < RHS->Size(); i++)
			(*RHS)(i) = 0.0;
	else
		(*RHS) = 0.0;

	// Fluid Momentum
	AssembleMassMatrix(option);  // This is exactly same with CalcMass().
	AssembleRHS(dimension);
	// Output matrices
	if (pcs->Write_Matrix)
	{
		for (i = 0; i < nnodes; i++)
			(*RHS)(i) = NodalVal[i];
		(*pcs->matrix_file) << "### Element: " << Index << "\n";
		(*pcs->matrix_file) << "---Mass matrix: "
		                    << "\n";
		Mass->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "---Laplacian matrix: "
		                    << "\n";
		Laplace->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "---RHS: "
		                    << "\n";
		RHS->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "\n";
		(*pcs->matrix_file) << "Stiffness: "
		                    << "\n";
		StiffMatrix->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "\n";
	}
}
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   18/02/2006 WW Implementation
**************************************************************************/
void CFiniteElementStd::ExtropolateGauss(CRFProcess* m_pcs, const int idof)
{
	int i, j, gp, gp_r, gp_s, gp_t, idx_v2 = 0;
	int i_s, i_e, ish;
	double EV, EV1 = 0.0, varx = 0.0;
	//
	MshElemType::type ElementType = MeshElement->GetElementType();

	// Multi-phase flow 03.2009 PCH
	if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || m_pcs->getProcessType() == FiniteElement::PS_GLOBAL || m_pcs->type == 42)
	{
		switch (idof)
		{
			case 0:
				idx_v2 = m_pcs->GetNodeValueIndex("VELOCITY_X2");
				break;
			case 1:
				idx_v2 = m_pcs->GetNodeValueIndex("VELOCITY_Y2");
				break;
			case 2:
				idx_v2 = m_pcs->GetNodeValueIndex("VELOCITY_Z2");
				break;
		}
	}
	// For strain and stress extrapolation all element types
	// Number of elements associated to nodes
	for (i = 0; i < nnodes; i++)
		dbuff[i] =
		    (double)MeshElement->GetNode(i)->getConnectedElementIDs().size();
	//
	gp_r = gp_s = gp_t = gp = 0;
	ElementValue* gp_ele = ele_gp_value[Index];
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

		NodalVal1[i] = gp_ele->Velocity(idof, gp) * time_unit_factor;
		//
		//
		// PCH 05.2009
		if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || m_pcs->getProcessType() == FiniteElement::PS_GLOBAL || m_pcs->type == 42)
			NodalVal2[i] = gp_ele->Velocity_g(idof, gp) * time_unit_factor;
	}

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
	double avgEV = .0;
	double avgEV1 = .0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
	{
		// average
		avgEV = CalcAverageGaussPointValues(NodalVal1);
		if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || m_pcs->getProcessType() == FiniteElement::PS_GLOBAL)
			avgEV1 = CalcAverageGaussPointValues(NodalVal2);
	}

	for (i = 0; i < nnodes; i++)
	{
		EV = EV1 = varx = 0.0;

		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			SetExtropoGaussPoints(i);
			//
			ComputeShapefct(1);  // Linear interpolation function
			for (j = i_s; j < i_e; j++)
				EV += NodalVal1[j] * shapefct[j - ish];
		}
		else if (this->GetExtrapoMethod() ==
		         ExtrapolationMethod::EXTRAPO_AVERAGE)
			// average
			EV = avgEV;
		// for(j=i_s; j<i_e; j++)
		// EV += NodalVal1[j];
		// EV /=(i_e-i_s);	//WX:09.2010. Use average value for nodes.
		// Average value of the contribution of ell neighbor elements
		EV /= dbuff[i];
		EV += m_pcs->GetNodeValue(nodes[i], idx_vel[idof]);
		m_pcs->SetNodeValue(nodes[i], idx_vel[idof], EV);
		//
		// Multi-phase flow PCH 05.2009
		if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || m_pcs->getProcessType() == FiniteElement::PS_GLOBAL || m_pcs->type == 42)
		{
			// Calculate values at nodes
			if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
				for (j = i_s; j < i_e; j++)
					EV1 += NodalVal2[j] * shapefct[j - ish];
			else if (this->GetExtrapoMethod() ==
			         ExtrapolationMethod::EXTRAPO_AVERAGE)
				// average
				EV1 = avgEV1;
			// for(j=i_s; j<i_e; j++)
			// EV += NodalVal1[j];
			// EV /=(i_e-i_s);	//WX:09.2010. Use average value for nodes.
			//
			EV1 /= dbuff[i];
			EV1 += m_pcs->GetNodeValue(nodes[i], idx_v2);
			m_pcs->SetNodeValue(nodes[i], idx_v2, EV1);
		}
		//
	}
}

/***********************************************************************
   27.03.2007 WW
***********************************************************************/
void CFiniteElementStd::CalcSatution()
{
	int i, j, gp, gp_r, gp_s, gp_t, idx_cp, idx_S;
	int i_s, i_e, ish;
	//  int l1,l2,l3,l4; //, counter;
	double sign, eS = 0.0;
	// CB_merge_0513
	double* tens = NULL;
	int Index;
	Index = MeshElement->GetIndex();

	MshElemType::type ElementType = MeshElement->GetElementType();
	//----------------------------------------------------------------------
	// Media
	int mmp_index = 0;
	long group = MeshElement->GetPatchIndex();
	mmp_index = group;
	//
	if (pcs->type == 22)
	{
		if (pcs->GetContinnumType() == 0)  // Matrix //WW
			mmp_index = 2 * group;
		else  // fracture //WW
			mmp_index = 2 * group + 1;
	}
	MediaProp = mmp_vector[mmp_index];
	MediaProp->m_pcs = pcs;
	MediaProp->Fem_Ele_Std = this;
	// CB_merge_0513
	tens = MediaProp->PermeabilityTensor(Index);
	//
	sign = -1.0;
	idx_cp = pcs->GetNodeValueIndex("PRESSURE1") + 1;
	idx_S = pcs->GetNodeValueIndex("SATURATION1", true);
	// Dual Richards
	if (pcs->type == 22 && pcs->GetContinnumType() == 1)
	{
		idx_cp = pcs->GetNodeValueIndex("PRESSURE2") + 1;
		idx_S = pcs->GetNodeValueIndex("SATURATION2") + 1;
	}
	if (pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW || pcs->type == 42) sign = 1.0;
	//
	for (i = 0; i < nnodes; i++)
	{
		// Number of elements associated to nodes
		dbuff[i] =
		    (double)MeshElement->GetNode(i)->getConnectedElementIDs().size();
		// pressure
		NodalVal0[i] = sign * pcs->GetNodeValue(nodes[i], idx_cp);
	}
	//
	gp_r = gp_s = gp_t = gp = 0;
	//
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		SetGaussPoint(gp, gp_r, gp_s, gp_t);
		if (ElementType == MshElemType::QUAD ||
		    ElementType == MshElemType::HEXAHEDRON)
		{
			i = GetLocalIndex(gp_r, gp_s, gp_t);
			if (i == -1) continue;
		}
		else
			i = gp;
		//
		if (i > nnodes) continue;
		ComputeShapefct(1);
		//
		// CB_merge_0513 in case of het K, store local K
		MediaProp->local_permeability = tens[0];
		PG = interpolate(NodalVal0);
		NodalVal_Sat[i] = MediaProp->SaturationCapillaryPressureFunction(PG);
	}

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
	double avgSat = .0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
		// average
		avgSat = CalcAverageGaussPointValues(NodalVal_Sat);
	for (i = 0; i < nnodes; i++)
	{
		eS = 0.0;
		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			SetExtropoGaussPoints(i);
			//
			ComputeShapefct(1);  // Linear interpolation function
			for (j = i_s; j < i_e; j++)
				eS += NodalVal_Sat[j] * shapefct[j - ish];
		}
		else if (this->GetExtrapoMethod() ==
		         ExtrapolationMethod::EXTRAPO_AVERAGE)
			eS = avgSat;
		// Average value of the contribution of ell neighbor elements
		eS /= dbuff[i];
		eS += pcs->GetNodeValue(nodes[i], idx_S);
		// In case the node is on the material interface
		if (eS > 1.0) eS = 1.0;
		//
		pcs->SetNodeValue(nodes[i], idx_S, eS);
	}
}
/**************************************************************************
   FEMLib-Method:
   Task: Caculate material parameter at element nodes for output
   Programing:
   04/2007 WW Implementation
**************************************************************************/
void CFiniteElementStd::CalcNodeMatParatemer()
{
	int i, gp_r, gp_s, gp_t, idx_perm[3], idxp = 0;
	int i_s, i_e, ish;
	double w[3], nval = 0.0;
	//
	MshElemType::type ElementType = MeshElement->GetElementType();
	//----------------------------------------------------------------------
	gp = 0;
	index = Index;
	w[0] = w[1] = w[2] = 1.0;
	//----------------------------------------------------------------------
	setOrder(1);
	// Set material
	SetMaterial();
	//----------------------------------------------------------------------
	// Node value of the previous time step
	int idx11 = idx1;
	if (pcs->GetContinnumType() == 1) idx11 = idxp21;
	for (i = 0; i < nnodes; i++)
		NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx11);
	if (PcsType == V)
		for (i = 0; i < nnodes; i++)
			NodalVal_p2[i] = pcs->GetNodeValue(nodes[i], idxp21);
	if (PcsType == P)  // 4.3.2009 PCH

		for (i = 0; i < nnodes; i++)
			NodalVal_SatNW[i] = pcs->GetNodeValue(nodes[i], idxSn1);
	if (cpl_pcs)
		for (i = 0; i < nnodes; i++)
		{
			NodalValC[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c0);
			NodalValC1[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c1);
			if (cpl_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
				NodalVal_p2[i] = cpl_pcs->GetNodeValue(nodes[i], idx_c1 + 2);
			// AKS
			NodalVal_p20[i] = pcs->GetNodeValue(nodes[i], idx_c0 + 2);
		}
	//
	if ((pcs->additioanl2ndvar_print > 0) && (pcs->additioanl2ndvar_print < 3))
	{
		idx_perm[0] = pcs->GetNodeValueIndex("PERMEABILITY_X1");
		idx_perm[1] = pcs->GetNodeValueIndex("PERMEABILITY_Y1");
		if (dim == 3)  // 3D
			idx_perm[2] = pcs->GetNodeValueIndex("PERMEABILITY_Z1");
	}
	if (pcs->additioanl2ndvar_print > 1)
		idxp = pcs->GetNodeValueIndex("POROSITY");
	// Number of elements associated to nodes
	for (i = 0; i < nnodes; i++)
		dbuff[i] =
		    (double)MeshElement->GetNode(i)->getConnectedElementIDs().size();
	//
	gp_r = gp_s = gp_t = gp = 0;
	//
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		SetGaussPoint(gp, gp_r, gp_s, gp_t);
		if (ElementType == MshElemType::QUAD ||
		    ElementType == MshElemType::HEXAHEDRON)
		{
			i = GetLocalIndex(gp_r, gp_s, gp_t);
			if (i == -1) continue;
		}
		else
			i = gp;
		//
		if (i > nnodes) continue;
		ComputeShapefct(1);
		PG = interpolate(NodalVal1);
		//
		if ((pcs->additioanl2ndvar_print > 0) &&
		    (pcs->additioanl2ndvar_print < 3))
		{
			double* tensor = MediaProp->PermeabilityTensor(Index);
			// Modified LBNL model
			if (MediaProp->permeability_stress_mode == 2 ||
			    MediaProp->permeability_stress_mode == 3)
			{
				if (cpl_pcs)
					TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
				else
					TG = 293.15;
				MediaProp->CalStressPermeabilityFactor(w, TG);
				for (size_t j = 0; j < dim; j++)
					tensor[j * dim + j] *= w[j];
			}
			NodalVal2[i] = tensor[0];                          // w[0];
			NodalVal3[i] = tensor[dim + 1];                    // w[1]; //
			if (dim == 3) NodalVal4[i] = tensor[2 * dim + 2];  // w[2]; //
		}
		// Porosity
		if (pcs->additioanl2ndvar_print > 1)
			// MediaProp->Porosity(this);
			NodalVal0[i] = MediaProp->Porosity(MeshElement->GetIndex(), 1.0);
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
	double avgW[3] = {};
	double avgVal = 0.0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
	{
		// average
		if ((pcs->additioanl2ndvar_print > 0) &&
		    (pcs->additioanl2ndvar_print < 3))
		{
			avgW[0] = CalcAverageGaussPointValues(NodalVal2);
			avgW[1] = CalcAverageGaussPointValues(NodalVal3);
			avgW[2] = CalcAverageGaussPointValues(NodalVal4);
		}
		if (pcs->additioanl2ndvar_print > 1)
			avgVal = CalcAverageGaussPointValues(NodalVal0);
	}
	for (i = 0; i < nnodes; i++)
	{
		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			SetExtropoGaussPoints(i);
			//
			ComputeShapefct(1);  // Linear interpolation function
		}
		if ((pcs->additioanl2ndvar_print > 0) &&
		    (pcs->additioanl2ndvar_print < 3))
		{
			w[0] = w[1] = w[2] = 0.0;
			if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
				for (int j = i_s; j < i_e; j++)
				{
					w[0] += NodalVal2[j] * shapefct[j - ish];
					w[1] += NodalVal3[j] * shapefct[j - ish];
					if (dim == 3) w[2] += NodalVal4[j] * shapefct[j - ish];
				}
			else if (this->GetExtrapoMethod() ==
			         ExtrapolationMethod::EXTRAPO_AVERAGE)
				for (size_t k = 0; k < dim; k++)
					w[k] = avgW[k];
			// Average value of the contribution of ell neighbor elements
			for (size_t k = 0; k < dim; k++)
			{
				w[k] /= dbuff[i];
				w[k] += pcs->GetNodeValue(nodes[i], idx_perm[k]);
				//
				pcs->SetNodeValue(nodes[i], idx_perm[k], w[k]);
			}
		}
		if (pcs->additioanl2ndvar_print > 1)
		{
			nval = 0.0;
			if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
				for (int j = i_s; j < i_e; j++)
					nval += NodalVal0[j] * shapefct[j - ish];
			else if (this->GetExtrapoMethod() ==
			         ExtrapolationMethod::EXTRAPO_AVERAGE)
				nval = avgVal;
			nval /= dbuff[i];
			nval += pcs->GetNodeValue(nodes[i], idxp);
			//
			pcs->SetNodeValue(nodes[i], idxp, nval);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   01/2006 OK Implementation
**************************************************************************/
// void CFiniteElementStd::AssembleLHSMatrix()
void CFiniteElementStd::AssembleParabolicEquationRHSVector()
{
	int i;
	//----------------------------------------------------------------------
	// TIM
	double dt_inverse = 0.0;
	dt_inverse = 1.0 / dt;
	//----------------------------------------------------------------------
	// Initialize
	// if (pcs->Memory_Type==2) skip the these initialization
	(*Mass) = 0.0;
	(*Laplace) = 0.0;
	//----------------------------------------------------------------------
	// Calculate matrices
	// Mass matrix..........................................................
	if (pcs->m_num->ele_mass_lumping)
		CalcLumpedMass();
	else
		CalcMass();
	// Laplace matrix.......................................................
	CalcLaplace();
	//----------------------------------------------------------------------
	// Assemble local LHS matrix:
	// [C]/dt + theta [K]
	// Mass matrix
	*StiffMatrix = *Mass;
	(*StiffMatrix) *= dt_inverse;
	// Laplace matrix
	*AuxMatrix = *Laplace;
	*StiffMatrix += *AuxMatrix;
	//----------------------------------------------------------------------
	for (i = 0; i < nnodes; i++)
	{
		NodalVal1[i] = pcs->GetNodeValue(nodes[i], idx1);
		NodalVal[i] = 0.0;
	}
	//----------------------------------------------------------------------
	StiffMatrix->multi(NodalVal1, NodalVal);
//----------------------------------------------------------------------
#ifdef NEW_EQS
	eqs_rhs = pcs->eqs_new->getRHS();
	for (i = 0; i < nnodes; i++)
	{
		eqs_number[i] = MeshElement->GetNode(i)->GetEquationIndex();
		eqs_rhs[eqs_number[i]] += NodalVal[i];
	}
#endif
	//----------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate right hand terms temperature coupled term and body force
   Programing:
   05/2010 AKS Implementation
   last modification:
**************************************************************************/

double CFiniteElementStd::CalCoef_RHS_AIR_FLOW(int dof_index)
{
	double val = 0.0;
	int Index = MeshElement->GetIndex();
	PG = interpolate(NodalVal1);
	TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
	TG0 = interpolate(NodalValC) + T_KILVIN_ZERO;
	switch (dof_index)
	{
		case 0:
			val = -MediaProp->Porosity(Index, pcs->m_num->ls_theta) / TG;
			val *= (TG - TG0);
			break;

		case 1:
			val = -1.0 / TG;
			break;
	}
	return val;
}
/**************************************************************************
   FEMLib-Method:
   Task: Calculate RHS of pressure coupled term
   Programing:
   05/2010 AKS Implementation
   last modification:
**************************************************************************/

double CFiniteElementStd::CalCoef_RHS_HEAT_TRANSPORT(int dof_index)
{
	double val = 0.0, rho_g = 0.0, rho_0 = 0.0;
	int Index = MeshElement->GetIndex();
	double dens_arg[3];
	dens_arg[0] = interpolate(NodalValC1);
	dens_arg[1] = interpolate(NodalVal1) + T_KILVIN_ZERO;
	dens_arg[2] = Index;
	rho_g = FluidProp->Density(dens_arg);
	dens_arg[0] = 4.0e6;
	dens_arg[1] = 120 + T_KILVIN_ZERO;
	rho_0 = FluidProp->Density(dens_arg);

	switch (dof_index)
	{
		case 0:
			val = (interpolate(NodalValC1) - interpolate(NodalValC)) *
			      MediaProp->Porosity(Index, pcs->m_num->ls_theta) * rho_g /
			      rho_0;
			break;

		case 1:
			val = rho_g / rho_0;
			val -= 1.0;  // term coresponding to the 'Viscour dissipation'
			break;
	}
	return val;
}

/***************************************************************************
   GeoSys - Funktion:
          Assemble_RHS_LIQUIDFLOW
   Programming:
   11/2012   NW
 **************************************************************************/
void CFiniteElementStd::Assemble_RHS_LIQUIDFLOW()
{
	return;
	if (!isTemperatureCoupling()) return;
	if ((FluidProp->drho_dT == .0 &&
	     (FluidProp->density_model < 8 || FluidProp->density_model > 14)) &&
	    SolidProp->Thermal_Expansion() == .0)
		return;
	if (pcs->tim_type == FiniteElement::TIM_STEADY) return;

#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	int dm_shift = 0;
	if (pcs->type / 10 == 4) dm_shift = problem_dimension_dm;
#endif
	//----------------------------------------------------------------------
	for (int i = 0; i < nnodes; i++)
		NodalVal[i] = 0.0;
	//======================================================================
	// Loop over Gauss points
	int gp_r = 0, gp_s = 0, gp_t = 0;
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		const double gp_fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeShapefct(1);  // Linear interpolation function
		//---------------------------------------------------------
		//  Evaluate variables
		//---------------------------------------------------------
		const double T_n = interpolate(NodalValC);
		const double T_n1 = interpolate(NodalValC1);
		const double dT = T_n1 - T_n;
		//---------------------------------------------------------
		//  Evaluate material property
		//---------------------------------------------------------
		const double poro = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
		double alpha_T_s =
		    3. *
		    SolidProp
		        ->Thermal_Expansion();  // multiply 3 for volumetrix expression
		Sw = 1.0;
		double alpha_T_l;
		if (FluidProp->density_model > 7 && FluidProp->density_model < 15)
		{
			double arg[2];
			arg[0] = interpolate(NodalVal1);   // p
			arg[1] = interpolate(NodalValC1);  // T
			alpha_T_l = -FluidProp->drhodT(arg) / FluidProp->Density();
		}
		else
			alpha_T_l =
			    -FluidProp
			         ->drho_dT;  // negative sign is required due to OGS input

		if (PcsType == R)
		{
			// for Richards:
			PG = interpolate(NodalVal1);
			if (PG < 0.0)
			{
				if (FluidProp->drho_dT_unsaturated)
					Sw = MediaProp->SaturationCapillaryPressureFunction(-PG);
				else
					alpha_T_l = alpha_T_s = 0.0;
			}
		}
		const double eff_thermal_expansion =
		    (SolidProp->biot_const - poro) * alpha_T_s + poro * Sw * alpha_T_l;
		//---------------------------------------------------------
		//  Compute RHS+=int{N^T alpha_T dT/dt}
		//---------------------------------------------------------
		const double fac =
		    eff_thermal_expansion * dT / dt / time_unit_factor;  // WX:bug fixed
		for (int i = 0; i < nnodes; i++)
			NodalVal[i] += gp_fkt * fac * shapefct[i];
	}
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	int i_sh = NodeShift[dm_shift];
#endif
	for (int i = 0; i < nnodes; i++)
	{
#ifndef USE_PETSC
		eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i];
#endif
		(*RHS)(i + LocalShift) += NodalVal[i];
	}
}

/***************************************************************************
   GeoSys - Funktion:
   Assemble_RHS_AIR_FLOW
   Programming:
    05/2010 AKS
 **************************************************************************/

void CFiniteElementStd::Assemble_RHS_AIR_FLOW()
{
	int j, ii;  // KR,idxd;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;  // KR ,z_sum;
	double vel[3];                     // KR,rhoz[3];
	double fkt, fac, mat_fac, fluid_density;
	double dens_arg[3];  // 08.05.2008 WW
	double* tensor = NULL;
	// KR CFEMesh* m_msh;
	int GravityOn = 1;  // Initialized to be on
	// If no gravity, then set GravityOn to be zero.
	if ((coordinate_system) % 10 != 2 && (!axisymmetry)) GravityOn = 0;
	// Material
	int dof_n = 1;
	//----------------------------------------------------------------------
	for (int i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function
		ElementValue* gp_ele = ele_gp_value[Index];

		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			fac = CalCoef_RHS_AIR_FLOW(ii) / dt;

			for (int i = 0; i < nnodes; i++)
				NodalVal[i + ii * nnodes] += fac * fkt * shapefct[i];
		}

		// grad T
		for (ii = 0; ii < dof_n; ii++)
		{
			fac = CalCoef_RHS_AIR_FLOW(ii + 1);
			// Velocity
			vel[0] = fac * gp_ele->Velocity(0, gp);
			vel[1] = fac * gp_ele->Velocity(1, gp);
			vel[2] = fac * gp_ele->Velocity(2, gp);

			for (int i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i + ii * nnodes] +=
						    fkt * shapefct[i] * vel[k] *
						    dshapefct[k * nnodes + j] *
						    (NodalValC1[j] + T_KILVIN_ZERO);
		}

		// Body force term
		if (GravityOn)
		{
			dens_arg[0] = interpolate(NodalVal1);
			dens_arg[1] = interpolate(NodalValC1) + T_KILVIN_ZERO;
			dens_arg[2] = Index;
			fluid_density = FluidProp->Density(dens_arg);
			mat_fac = FluidProp->Viscosity(dens_arg);
			tensor = MediaProp->PermeabilityTensor(Index);
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] / mat_fac;
			for (ii = 0; ii < dof_n; ii++)
			{
				for (int i = 0; i < nnodes; i++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i + ii * nnodes] -=
						    fkt * fluid_density * gravity_constant *
						    mat[dim * k + dim - 1] * dshapefct[k * nnodes + i];
			}
		}
	}
	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	for (ii = 0; ii < pcs->dof; ii++)
	{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
		i_sh = NodeShift[ii];
#endif
		ii_sh = ii * nnodes;
		for (int i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i + ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) -= NodalVal[i + ii_sh];
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
   Assemble_RHS_HEAT_TRANSPORT: This include when need pressure terms n dp/dt +
 nv.Nabla p
   Programming:
   05/2010   AKS
 **************************************************************************/

void CFiniteElementStd::Assemble_RHS_HEAT_TRANSPORT()
{
	int i, j, ii;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double vel[3];
	double fkt = 0.0, fac = 0.0;
	// Material
	int dof_n = 1;
	//----------------------------------------------------------------------
	for (i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function
		ElementValue* gp_ele = ele_gp_value[Index];

		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			fac = CalCoef_RHS_HEAT_TRANSPORT(ii) / dt;

			for (i = 0; i < nnodes; i++)
				NodalVal[i + ii * nnodes] += fac * fkt * shapefct[i];
		}

		// grad P

		for (ii = 0; ii < dof_n; ii++)
		{
			fac = CalCoef_RHS_HEAT_TRANSPORT(ii + 1);
			// Velocity
			vel[0] = fac * gp_ele->Velocity(0, gp);
			vel[1] = fac * gp_ele->Velocity(1, gp);
			vel[2] = fac * gp_ele->Velocity(2, gp);

			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i + ii * nnodes] +=
						    fkt * vel[k] * shapefct[i] *
						    dshapefct[k * nnodes + j] * NodalValC1[j];
		}
	}
	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	for (ii = 0; ii < pcs->dof; ii++)
	{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
		i_sh = NodeShift[ii];
#endif
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i + ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) -= NodalVal[i + ii_sh];
		}
	}
}
/**************************************************************************
   FEMLib-Method:
   Task: Calculate RHS of pressure coupled term
   Programing:
   05/2010 AKS Implementation
   last modification:
**************************************************************************/

double CFiniteElementStd::CalCoef_RHS_HEAT_TRANSPORT2(int dof_index)
{
	// TF unused variable - comment fix compile warning
	//      ElementValue* gp_ele = ele_gp_value[Index];
	double* tensor = NULL;
	double val = 0.0, mat_fac;
	// TF unused variable - comment fix compile warning
	//      double Tc=647.096;
	double H_vap = 0.0, dens_arg[3];
	ComputeShapefct(1);
	PG = interpolate(NodalValC1);
	PG2 = interpolate(NodalVal_p2);
	TG = interpolate(NodalVal1) + T_KILVIN_ZERO;
	PG0 = interpolate(NodalValC);
	PG20 = interpolate(NodalVal_p20);
	dens_arg[1] = TG;
	dens_arg[0] = PG2 - PG;
	rhow = FluidProp->Density(dens_arg);
	Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
	dSdp = MediaProp->PressureSaturationDependency(
	    Sw, true);  // JT: now returns correct sign.
	poro = MediaProp->Porosity(Index, pcs->m_num->ls_theta);
	if (MediaProp->evaporation == 647)
		H_vap = -2257000;  // pow((Tc - TG),0.38)*2.5397E+5;//It is specific you
	                       // can change thi value as you chaning fluid from
	                       // water
	for (size_t i = 0; i < dim * dim; i++)
		mat[i] = 0.0;
	switch (dof_index)
	{
		case 0:
			val = H_vap * rhow * poro * dSdp;
			val *= (PG - PG0);
			return val;
			break;

		case 1:
			val = 0.0;
			return val;
			break;

		case 2:
			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = H_vap * rhow *
			          MediaProp->PermeabilitySaturationFunction(Sw, 0) /
			          FluidProp->Viscosity();
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;

		case 3:
			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = -H_vap * rhow *
			          MediaProp->PermeabilitySaturationFunction(Sw, 0) /
			          FluidProp->Viscosity();
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
	}
	return 0.;  // WW
}
/***************************************************************************
   GeoSys - Funktion:
   Assemble_RHS_HEAT_TRANSPORT2: This include when need pressure terms n dp/dt +
 nv.Nabla p
   Programming:
   05/2010   AKS
 **************************************************************************/

void CFiniteElementStd::Assemble_RHS_HEAT_TRANSPORT2()
{
	int i, j, ii;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	// TF unused variable - comment fix compile warning
	//      double *tensor = NULL,
	double dens_arg[3];
	// TF unused variable - comment fix compile warning
	//      double H_vap=0;
	// TF unused variable - comment fix compile warning
	//      double Tc=647.096;
	double fkt = 0.0, fac = 0.0;  // WW,mat_fac;
	// Material
	int dof_n = 1;
	//----------------------------------------------------------------------
	for (i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function
		// TF unused variable - comment fix compile warning
		//         ElementValue* gp_ele = ele_gp_value[Index];
		int dof_n = 2;
		int GravityOn = 1;  // Initialized to be on
		// If no gravity, then set GravityOn to be zero.
		if ((coordinate_system) % 10 != 2 && (!axisymmetry)) GravityOn = 0;
		TG = interpolate(NodalVal1) + T_KILVIN_ZERO;
		PG = interpolate(NodalValC1);
		PG2 = interpolate(NodalVal_p2);
		dens_arg[1] = TG;
		dens_arg[0] = PG2 - PG;

		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			fac = fkt * CalCoef_RHS_HEAT_TRANSPORT2(ii) / dt;
			// Calculate THS
			for (i = 0; i < nnodes; i++)
				NodalVal[i] += fac * shapefct[i];
		}

		// grad pc
		for (ii = 0; ii < dof_n - 1; ii++)
		{
			// Material
			CalCoef_RHS_HEAT_TRANSPORT2(ii + dof_n);
			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i] +=
						    fkt * mat[dim * k + k] * dshapefct[k * nnodes + i] *
						    dshapefct[k * nnodes + j] * NodalValC1[j];
		}

		// grad pc
		for (ii = 0; ii < dof_n - 1; ii++)
		{
			// Material
			CalCoef_RHS_HEAT_TRANSPORT2(ii + dof_n + 1);
			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i] +=
						    fkt * mat[dim * k + k] * dshapefct[k * nnodes + i] *
						    dshapefct[k * nnodes + j] * NodalVal_p2[j];
		}

		// gravity
		if (GravityOn)
		{
			CalCoef_RHS_HEAT_TRANSPORT2(2);
			for (i = 0; i < nnodes; i++)
				for (size_t k = 0; k < dim; k++)
					NodalVal[i] -= fkt * mat[dim * k + dim - 1] *
					               FluidProp->Density(dens_arg) *
					               gravity_constant * dshapefct[k * nnodes + i];
		}
	}
	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	for (ii = 0; ii < pcs->dof; ii++)
	{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
		i_sh = NodeShift[ii];
#endif
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i + ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) -= NodalVal[i + ii_sh];
		}
	}
}

/**************************************************************************
   PCSLib-Method:
   01/2007 OK Implementation
   02/2009 PCH modified to handle known Snw in the mass matrix
**************************************************************************/
void CFiniteElementStd::AssembleRHSVector()
{
	int i;
	int idx_fv = 0, idx_pw = 0, idx_pc = 0;
	double NodalVal_FV[20];
	// OK411 double FV;
	CRFProcess* pcs_p = NULL;
	CRFProcess* pcs_s = NULL;
	//----------------------------------------------------------------------
	// Initializations
	for (i = 0; i < nnodes; i++)
		NodalVal[i] = 0.0;

	// TF fixed warning -Wunused-but-set-variable
	//      double temp[8];

	switch (PcsType)
	{
		//....................................................................
		case T:                       // Two-phase flow
			if (pcs->PartialPS == 0)  // If not partial-pressure-based
				(*Laplace) = 0.0;
			else
				(*Mass) = 0.0;
			break;
		default:
			break;
			//....................................................................
	}
	//----------------------------------------------------------------------
	// Field variables
	switch (PcsType)
	{
		//....................................................................
		case T:  // Two-phase flow
			if (pcs->PartialPS == 0)
			{
				pcs_p = pcs_vector[0];
				pcs_s = pcs_vector[1];

				idx_pw = pcs_p->GetNodeValueIndex("PRESSURE1");
				idx_pc = pcs_p->GetNodeValueIndex("PRESSURE_CAP");
				idx_fv = pcs_s->GetNodeValueIndex("SATURATION2");
				// WW CMediumProperties *m_mmp = NULL;
				// WW m_mmp = mmp_vector[0];
				for (i = 0; i < nnodes; i++)
				{
					Sw = 1.0 - pcs_s->GetNodeValue(nodes[i], idx_fv + 1);
					double Pw = pcs_p->GetNodeValue(nodes[i], idx_pw + 1);
					double Pc = pcs_p->GetNodeValue(nodes[i], idx_pc + 1);
					if (pcs->ML_Cap == 0)  // If ODE method for Snw,
						NodalVal_FV[i] = -(Pw + Pc);
					else  // If PDE method,
						NodalVal_FV[i] = -Pw;
				}
			}
			else
			{
			}
			break;
		default:
			break;
			//....................................................................
	}

	//----------------------------------------------------------------------
	// Element matrices
	switch (PcsType)
	{
		//....................................................................
		case T:  // Two-phase flow
			if (pcs->PartialPS == 0)
				CalcLaplace();
			else
				CalcMass();
			break;
		default:
			break;
			//....................................................................
	}
	//----------------------------------------------------------------------
	// Calc RHS contribution
	switch (PcsType)
	{
		//....................................................................
		case T:  // Two-phase flow
			if (pcs->PartialPS == 0)
				Laplace->multi(NodalVal_FV, NodalVal);
			else
				Mass->multi(NodalVal_FV, NodalVal);
			break;
		default:
			break;
			//....................................................................
	}

	// TF fixed warning -Wunused-but-set-variable
	//      for(i=0;i<nnodes;i++)
	//         temp[i]=NodalVal[i];

	//----------------------------------------------------------------------
	// Store RHS contribution
	for (i = 0; i < nnodes; i++)
	{
// CB 04008
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
                        // TODO
#elif defined(NEW_EQS)
		pcs->eqs_new->getRHS()[NodeShift[problem_dimension_dm] + eqs_number[i]] +=
		    NodalVal[i];
#else
		pcs->eqs->b[NodeShift[problem_dimension_dm] + eqs_number[i]] +=
		    NodalVal[i];
#endif
		(*RHS)(i + LocalShift) += NodalVal[i];
	}
	//----------------------------------------------------------------------
	// RHS->Write();
}

/**************************************************************************
   PCSLib-Method:
   02/2009 PCH Implementation
**************************************************************************/
void CFiniteElementStd::PrintTheSetOfElementMatrices(std::string mark)
{
	// Output matrices
	if (pcs->Write_Matrix)
	{
		(*pcs->matrix_file) << "### Mark: " << mark << "\n";

		(*pcs->matrix_file) << "### Element: " << Index << "\n";
		(*pcs->matrix_file) << "---Mass matrix: "
		                    << "\n";
		if (Mass)
			Mass->Write(*pcs->matrix_file);
		else if (Mass2)
			Mass2->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "---Laplacian matrix: "
		                    << "\n";
		Laplace->Write(*pcs->matrix_file);

		(*pcs->matrix_file) << "---AuxMatrix1 matrix: "
		                    << "\n";
		AuxMatrix1->Write(*pcs->matrix_file);  // PCH for debug
		if (Advection)
		{
			// CMCD
			(*pcs->matrix_file) << "---Advective matrix: "
			                    << "\n";
			Advection->Write(*pcs->matrix_file);
		}
		if (StrainCoupling)
		{
			(*pcs->matrix_file) << "---Strain couping matrix: "
			                    << "\n";
			StrainCoupling->Write(*pcs->matrix_file);
		}
		(*pcs->matrix_file) << "---RHS: "
		                    << "\n";
		RHS->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "\n";
		(*pcs->matrix_file) << "Stiffness: "
		                    << "\n";
		StiffMatrix->Write(*pcs->matrix_file);
		(*pcs->matrix_file) << "\n";
	}
}

}  // end namespace

