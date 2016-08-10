/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_pcs_dm.h"

#include <cfloat>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>

#include "makros.h"
#include "display.h"
#include "StringTools.h"

#include "geo_sfc.h"

#include "msh_elem.h"

#include "mathlib.h"
#if defined(NEW_EQS)
#include "equation_class.h"
#endif
#include "FEMEnums.h"
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
#ifdef USE_PETSC
#include "PETSC/PETScLinearSolver.h"
#endif
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_node.h"
#include "rf_msp_new.h"
#include "rf_out_new.h"
#include "rf_pcs.h"
#include "rf_st_new.h"
#include "rf_tim_new.h"
#include "tools.h"


double LoadFactor = 1.0;
double Tolerance_global_Newton = 0.0;
double Tolerance_Local_Newton = 0.0;
int enhanced_strain_dm = 0;
int number_of_load_steps = 1;
int problem_dimension_dm = 0;
int PreLoad = 0;
bool GravityForce = true;

bool Localizing = false;  // for tracing localization
// Last discontinuity element correponding to SeedElement

using namespace std;

vector<DisElement*> LastElement(0);
vector<long> ElementOnPath(0);  // Element on the discontinuity path

using FiniteElement::CFiniteElementVec;
using FiniteElement::CFiniteElementStd;
using FiniteElement::ElementValue_DM;
using SolidProp::CSolidProperties;
using Math_Group::Matrix;

CRFProcessDeformation::CRFProcessDeformation()
    : CRFProcess(),
      fem_dm(NULL),
      ARRAY(NULL),
      p0(NULL),
      counter(0),
      InitialNormR0(0.0)

{
	norm_du0_pre_cpl_itr = 0.0;
	idata_type = none;
	_isInitialStressNonZero = false;  // NW
	InitialNormDU0 = 0.0;
}

CRFProcessDeformation::~CRFProcessDeformation()
{
	if (ARRAY) delete[] ARRAY;
	if (p0) delete[] p0;
	if (fem_dm) delete fem_dm;

	fem_dm = NULL;
	ARRAY = NULL;
	// Release memory for element variables
	// alle stationaeren Matrizen etc. berechnen
	long i;
	// Write Gauss stress // TEST for excavation analysis
	//   if(reload==1)
	// if(reload == 1 || reload == 3)
	if (idata_type == write_all_binary || idata_type == read_write)
	{
		WriteGaussPointStress();
		if (type == 41)  // mono-deformation-liquid
			WriteSolution();
	}
	//
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		delete ele_value_dm[i];
		ele_value_dm[i] = NULL;
	}
	if (enhanced_strain_dm > 0)
	{
		while (ele_value_dm.size() > 0)
			ele_value_dm.pop_back();
		for (i = 0; i < (long)LastElement.size(); i++)
		{
			DisElement* disEle = LastElement[i];
			delete disEle->InterFace;
			delete disEle;
			disEle = NULL;
		}
		while (LastElement.size() > 0)
			LastElement.pop_back();
	}
}

/*************************************************************************
   Task: Initilization for deformation process
   Programming:
   05/2003 OK/WW Implementation
   08/2003 WW   Some changes for monolithic scheme
   08/2004 WW   Changes based on PCSCreateMProcess(obsolete)
   last modified: WW
 **************************************************************************/
void CRFProcessDeformation::Initialization()
{
	//-- NW 25.10.2011
	// this section has to be executed at latest before calling InitGauss()
	// Control for reading and writing solution
	if (reload == 1) idata_type = write_all_binary;
	if (reload == 2) idata_type = read_all_binary;
	if (reload == 3) idata_type = read_write;
	//--
	if ((reload == 2 || reload == 3) && calcDiffFromStress0)
		_isInitialStressNonZero = true;  // NW

	// Local assembliers
	// An instaniate of CFiniteElementVec
	int i, Axisymm = 1;  // ani-axisymmetry
	//
	if (m_msh->isAxisymmetry()) Axisymm = -1;  // Axisymmetry is true
	fem_dm = new CFiniteElementVec(this, Axisymm * m_msh->GetCoordinateFlag());
	fem_dm->SetGaussPointNumber(m_num->ele_gauss_points);
	//
	// Monolithic scheme
	if (type / 10 == 4)
		fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());
	//
	pcs_number_deformation = pcs_number;
	//
	if (m_num)
	{
		Tolerance_Local_Newton = m_num->nls_plasticity_local_tolerance;
		Tolerance_global_Newton = m_num->nls_error_tolerance[0];
	}
	//

	// Initialize material transform tensor for tansverse isotropic elasticity
	// UJG/WW. 25.11.2009
	for (i = 0; i < (int)msp_vector.size(); i++)
		msp_vector[i]->CalculateTransformMatrixFromNormalVector(
		    problem_dimension_dm);

	if (!msp_vector.size())
	{
		std::cout << "***ERROR: MSP data not found!"
		          << "\n";
		return;
	}
	InitialMBuffer();
	InitGauss();


	// Initial pressure is stored to evaluate pressure difference from the
	// initial
	// because DEFORMATION calculates stress balance of changes from the initial
	// stress
	// NW 28.08.2012
	if (_isInitialStressNonZero)
	{
		std::cout << "->Initial stress is given."
		          << "\n";
		CRFProcess* h_pcs = PCSGet("LIQUID_FLOW");
		if (h_pcs)
		{  // NW
			assert(p0 == NULL);
			std::cout << "->Found LIQUID_FLOW. Store initial liquid pressure."
			          << "\n";
			int n_nodes = m_msh->GetNodesNumber(false);
			p0 = new double[n_nodes];
			const int id_p0 = h_pcs->GetNodeValueIndex("PRESSURE1");
			for (i = 0; i < n_nodes; i++)
			{
				p0[i] = h_pcs->GetNodeValue(i, id_p0);
			}
		}
	}

}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::InitialMBuffer
   Task:  Initialize the temporarily used variables
   Programming:
   12/2003 WW
 **************************************************************************/
void CRFProcessDeformation::InitialMBuffer()
{
	if (!msp_vector.size())
	{
		cout << "No .msp file.   " << endl;
		abort();
	}

	size_t bufferSize(0);
	bool HM_Stagered = false;
	if (GetObjType() == 4)
	{
		bufferSize = GetPrimaryVNumber() * m_msh->GetNodesNumber(true);
		if (H_Process) HM_Stagered = true;
	}
	else if (GetObjType() == 41)
		bufferSize = (GetPrimaryVNumber() - 1) * m_msh->GetNodesNumber(true) +
		             m_msh->GetNodesNumber(false);
	else if (GetObjType() == 42)
		bufferSize = (GetPrimaryVNumber() - 2) * m_msh->GetNodesNumber(true) +
		             2 * m_msh->GetNodesNumber(false);

	// Allocate memory for  temporal array
	ARRAY = new double[bufferSize];

	// Allocate memory for element variables
	MeshLib::CElem* elem = NULL;
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		//       if (elem->GetMark()) // Marked for use
		//       {
		ElementValue_DM* ele_val =
		    new ElementValue_DM(elem, m_num->ele_gauss_points, HM_Stagered);
		ele_value_dm.push_back(ele_val);
		//       }
	}
}

double CRFProcessDeformation::getNormOfDisplacements()
{
#ifdef USE_PETSC
	const int g_nnodes = m_msh->getNumNodesLocal_Q();
	const int size = g_nnodes * pcs_number_of_primary_nvals;
	vector<int> ix(size);
	vector<double> val(size);
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		const int nidx0 = GetNodeValueIndex(pcs_primary_function_name[i]);
		for (int j = 0; j < g_nnodes; j++)
		{
			int ish = pcs_number_of_primary_nvals * j + i;
			ix[ish] =
			    pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[j] +
			    i;
			val[ish] = GetNodeValue(j, nidx0 + 1);
		}
	}
	eqs_new->setArrayValues(0, size, &ix[0], &val[0], INSERT_VALUES);
	eqs_new->AssembleUnkowns_PETSc();
	double norm_u_k1 = eqs_new->GetVecNormX();
#else
	const int g_nnodes = m_msh->GetNodesNumber(true);
	//const int size = g_nnodes * pcs_number_of_primary_nvals;
	double val = .0;
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		const int nidx0 = GetNodeValueIndex(pcs_primary_function_name[i]);
		for (int j = 0; j < g_nnodes; j++)
		{
			val += std::pow(GetNodeValue(j, nidx0 + 1), 2.0);
		}
	}
	double norm_u_k1 = std::sqrt(val);
#endif
	return norm_u_k1;
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  Solve plastic deformation by generalized Newton-Raphson method
   Programming:
   02/2003 OK Implementation
   05/2003 WW Polymorphism function by OK
   last modified: 23.05.2003
 **************************************************************************/
double CRFProcessDeformation::Execute(int loop_process_number)
{
	ScreenMessage("\n================================================\n");
	ScreenMessage("->Process %d: %s\n", loop_process_number,
	              convertProcessTypeToString(getProcessType()).c_str());
	ScreenMessage("================================================\n");

	clock_t dm_time;
	dm_time = -clock();

	counter++;  // Times of this method  to be called
	// For pure elesticity
	const bool isLinearProblem = (pcs_deformation <= 100);

	// setup mesh
	m_msh->SwitchOnQuadraticNodes(true);
	if (hasAnyProcessDeactivatedSubdomains || Deactivated_SubDomain.size() > 0)
		CheckMarkedElement();

#ifdef NEW_EQS
	eqs_new->ConfigNumerics(m_num);
#endif

	if (this->first_coupling_iteration)
		StoreLastSolution();  // u_n-->temp

	//  Reset stress for each coupling step when partitioned scheme is applied
	//  to HM
	if (H_Process && (type / 10 != 4)) ResetCouplingStep();

	//
	// Compute the maximum ratio of load increment and
	//   predict the number of load steps
	// ---------------------------------------------------------------
	// Compute the ratio of the current load to initial yield load
	// ---------------------------------------------------------------
	number_of_load_steps = 1;
	if (type / 10 == 4)  // For monolithic scheme
		number_of_load_steps = 1;
	LoadFactor = 1.0;
	double damping = 1.0;
	for (int l = 1; l <= number_of_load_steps; l++)
	{
		// Initialize incremental displacement: w=0
		InitializeNewtonSteps();
		double NormDU = 0.0, NormR = 0.0;
		double ErrorR = 0.0, ErrorU = 0.0;
		const int MaxIteration =
		    isLinearProblem ? 1 : m_num->nls_max_iterations;
		if (!isLinearProblem)
		{
			ErrorR = ErrorU = NormR = NormDU = 1.0e+8;
			ScreenMessage("Starting loading step %d/%d. Load factor: %g\n", l,
			              number_of_load_steps, LoadFactor);
			ScreenMessage("------------------------------------------------\n");
		}

		// Begin Newton-Raphson steps
		ite_steps = 0;
		while (ite_steps < MaxIteration)
		{
			ite_steps++;
			ScreenMessage("-->Starting Newton-Raphson iteration: %d/%d\n",
			              ite_steps, MaxIteration);
			ScreenMessage("------------------------------------------------\n");
// Refresh solver
#if defined(USE_PETSC)
			eqs_new->Initialize();
#elif defined(NEW_EQS)
			eqs_new->Initialize();
#endif

			// Assemble and solve system equation
			ScreenMessage("Assembling equation system...\n");
			GlobalAssembly();

			// init solution vector
			if (isLinearProblem && type != 41)
#if defined(USE_PETSC)
				InitializeRHS_with_u0();
#else
				SetInitialGuess_EQS_VEC();
#endif

			ScreenMessage("Calling linear solver...\n");
			// Linear solver
#if defined(USE_PETSC)
			//			eqs_new->EQSV_Viewer("eqs" +
			// number2str(aktueller_zeitschritt) + "b");
			eqs_new->Solver();
			eqs_new->MappingSolution();
#elif defined(NEW_EQS)
			bool compress_eqs =
				(type / 10 == 4 || this->Deactivated_SubDomain.size() > 0);
			eqs_new->Solver(this->m_num, compress_eqs);
#endif

			if (!isLinearProblem)
			{
// Get norm of residual vector, solution increment
#if defined(USE_PETSC)
				NormR = eqs_new->GetVecNormRHS();
				NormDU = eqs_new->GetVecNormX();
#elif defined(NEW_EQS)
				NormR = eqs_new->NormRHS();
				NormDU = eqs_new->NormX();
#endif

				if (ite_steps == 1)
				{
					if (this->first_coupling_iteration)
					{
						InitialNormDU_coupling = NormDU;
						norm_du0_pre_cpl_itr = .0;
					}
					InitialNormR0 = NormR;
					InitialNormDU0 = NormDU;
				}

				// calculate errors
				ErrorR = NormR / (InitialNormR0 == 0 ? 1 : InitialNormR0);
				ErrorU = NormDU / (InitialNormDU0 == 0 ? 1 : InitialNormDU0);

				// Compute damping for Newton-Raphson step
				damping = 1.0;

#if 0
				if(ErrorR / Error1 > 1.0e-1 || ErrorU / ErrorU1 > 1.0e-1)
					damping = 0.5;
#endif

				//
				ScreenMessage("-->End of Newton-Raphson iteration: %d/%d\n",
				              ite_steps, MaxIteration);
				ScreenMessage(
				    "   NR-Error  RHS Norm 0  RHS Norm    Unknowns Norm  "
				    "Damping\n");
				ScreenMessage("   %8.2e  %8.2e    %8.2e    %8.2e       %8.2e\n",
				              ErrorR, InitialNormR0, NormR, NormDU, damping);
				ScreenMessage(
				    "------------------------------------------------\n");

				if (ErrorR > 100.0 && ite_steps > 1)
				{
					ScreenMessage(
					    "***Attention: Newton-Raphson step is diverged. "
					    "Programme halt!\n");
					exit(1);
				}
				//				if(InitialNormR0 < 10 * Tolerance_global_Newton)
				//					break;
				//				if(NormR < 0.001 * InitialNormR0)
				//					break;
				if (ErrorR <= Tolerance_global_Newton)
				{
					if (ite_steps == 1) UpdateIterativeStep(damping, 0);
					break;
				}
			}

			UpdateIterativeStep(damping, 0);  // w = w+dw
		}                                     // Newton-Raphson iteration

		// Update stresses
		UpdateStress();

		// Update displacements, u=u+w for the Newton-Raphson
		// u1 = u0 for linear problems
		UpdateIterativeStep(1.0, 1);
	}  // Load step

	// Determine the discontinuity surface if enhanced strain methods is on.
	if (enhanced_strain_dm > 0) Trace_Discontinuity();

	// Recovery the old solution.  Temp --> u_n	for flow proccess
	RecoverSolution();

	// get coupling error
	const double norm_u_k1 =
	    getNormOfDisplacements();  // InitialNormDU_coupling
	const double cpl_abs_error =
	    std::abs(InitialNormDU0 - norm_du0_pre_cpl_itr) /
	    (norm_u_k1 == 0 ? 1 : norm_u_k1);
	cpl_max_relative_error = cpl_abs_error / m_num->cpl_error_tolerance[0];
	cpl_num_dof_errors = 1;
	ScreenMessage("   ||u^k+1||=%g, ||du^k+1||-||du^k||=%g\n", norm_u_k1,
	              std::abs(InitialNormDU0 - norm_du0_pre_cpl_itr));

	// store current du0
	norm_du0_pre_cpl_itr = InitialNormDU0;

#ifdef NEW_EQS  // WW
	// Also allocate temporary memory for linear solver. WW
	eqs_new->Clean();
#endif

	//
	dm_time += clock();
	ScreenMessage("PCS error: %g\n", cpl_max_relative_error);
	ScreenMessage("CPU time elapsed in deformation: %g s\n",
	              (double)dm_time / CLOCKS_PER_SEC);
	ScreenMessage("------------------------------------------------\n");

	return cpl_max_relative_error;
}

/**************************************************************************
   ROCKFLOW - Funktion: InitializeStress

   Aufgabe:
   Initilize all Gausss values and others

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - const int NodesOfEelement:

   Ergebnis:
   - void -

   Programmaenderungen:
   01/2003  WW  Erste Version
   09/2007  WW  Parallelize the released load for the excavation modeling
   letzte Aenderung:
**************************************************************************/
void CRFProcessDeformation::InitGauss(void)
{
	const int LenMat = 7;
	size_t i;
	int j, k, gp, NGS, MatGroup, n_dom;
	int PModel = 1;
	int gp_r = 0, gp_s = 0, gp_t = 0;
	//  double z=0.0;
	double xyz[3];
	static double Strs[6];
	ElementValue_DM* eleV_DM = NULL;
	CSolidProperties* SMat = NULL;
	CInitialCondition* m_ic = NULL;
	std::vector<CInitialCondition*> stress_ic(6);

	// double M_cam = 0.0;
	double pc0 = 0.0;
	double OCR = 1.0;
	n_dom = k = 0;

	int Idx_Strain[9];

	int NS = 4;
	Idx_Strain[0] = GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = GetNodeValueIndex("STRAIN_XY");

	if (problem_dimension_dm == 3)
	{
		NS = 6;
		Idx_Strain[4] = GetNodeValueIndex("STRAIN_XZ");
		Idx_Strain[5] = GetNodeValueIndex("STRAIN_YZ");
	}
	Idx_Strain[NS] = GetNodeValueIndex("STRAIN_PLS");

	for (j = 0; j < NS; j++)
		stress_ic[j] = NULL;
	for (j = 0; j < (long)ic_vector.size(); j++)
	{
		m_ic = ic_vector[j];
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_XX)
			stress_ic[0] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_YY)
			stress_ic[1] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_ZZ)
			stress_ic[2] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_XY)
			stress_ic[3] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_XZ)
			stress_ic[4] = m_ic;
		if (m_ic->getProcessPrimaryVariable() == FiniteElement::STRESS_YZ)
			stress_ic[5] = m_ic;
	}
	int ccounter = 0;
	for (j = 0; j < NS; j++)
		if (stress_ic[j]) ccounter++;
	if (ccounter > 0) reload = -1000;

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
		for (j = 0; j < NS + 1; j++)
			SetNodeValue(i, Idx_Strain[j], 0.0);
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			MatGroup = elem->GetPatchIndex();
			SMat = msp_vector[MatGroup];
			elem->SetOrder(true);
			fem_dm->ConfigElement(elem);
			eleV_DM = ele_value_dm[i];
			*(eleV_DM->Stress0) = 0.0;
			*(eleV_DM->Stress) = 0.0;
			PModel = SMat->Plasticity_type;

			for (j = 3; j < fem_dm->ns; j++)
				Strs[j] = 0.0;

			if (PModel == 2) *(eleV_DM->xi) = 0.0;

			if (PModel == 3)
			{
				// WW M_cam = (*SMat->data_Plasticity)(0);
				pc0 = (*SMat->data_Plasticity)(
				    3);  // The initial preconsolidation pressure
				         // Void ratio
				*(eleV_DM->e_i) = (*SMat->data_Plasticity)(4);
				OCR = (*SMat->data_Plasticity)(5);  // Over consolidation ratio
				for (j = 0; j < 3; j++)
					Strs[j] = (*SMat->data_Plasticity)(6 + j);

				/*
				   g_s = GetSolidDensity(i);
				   if(g_s<=0.0)
				   {
				   printf("\n !!! Input error. Gravity density should not be
				   less than zero with Cam-Clay model\n  ");
				   abort();
				   }

				   if(EleType== TriEle) // Triangle
				   nh = 6;
				   // Set soil profile. Cam-Clay. Step 2
				   for (j = 0; j < nh; j++)
				   h_node[j]=GetNodeY(element_nodes[j]); //Note: for 3D, should
				   be Z
				 */
			}

			//
			// if 2D //ToDo: Set option for 3D
			// Loop over Gauss points
			NGS = fem_dm->GetNumGaussPoints();
			// WW NGSS = fem_dm->GetNumGaussSamples();

			for (gp = 0; gp < NGS; gp++)
			{
				if (ccounter > 0)
				{
					fem_dm->GetGaussData(gp, gp_r, gp_s, gp_t);
					fem_dm->ComputeShapefct(2);
					fem_dm->RealCoordinates(xyz);
					for (j = 0; j < NS; j++)
					{
						m_ic = stress_ic[j];
						if (!m_ic) continue;
						n_dom = m_ic->GetNumDom();
						for (k = 0; k < n_dom; k++)
						{
							if (MatGroup != m_ic->GetDomain(k)) continue;
							(*eleV_DM->Stress)(j, gp) =
							    m_ic->getLinearFunction()->getValue(
							        m_ic->GetDomain(k), xyz[0], xyz[1], xyz[2]);
							(*eleV_DM->Stress0)(j, gp) =
							    (*eleV_DM->Stress)(j, gp);
						}
					}
				}
				else
				{
					switch (PModel)
					{
						case 2:  // Weimar's model
							// Initial stress_xx, yy,zz
							for (j = 0; j < 3; j++)
								(*eleV_DM->Stress)(j, gp) =
								    (*SMat->data_Plasticity)(20 + j);
							break;
						case 3:  // Cam-Clay
							for (j = 0; j < 3; j++)
								(*eleV_DM->Stress0)(j, gp) = Strs[j];
							(*eleV_DM->Stress) = (*eleV_DM->Stress0);
							break;
					}
				}
				if (eleV_DM->Stress_j)
					(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
				//
				switch (PModel)
				{
					case 2:  // Weimar's model
						for (j = 0; j < LenMat; j++)
							(*eleV_DM->MatP)(j, gp) =
							    (*SMat->data_Plasticity)(j);
						break;
					case 3:          // Cam-Clay
						pc0 *= OCR;  /// TEST
						(*eleV_DM->prep0)(gp) = pc0;
						break;
				}
				//
			}
// Initial condition by LBNL
////////////////////////////////////////////////////////
//#define  EXCAVATION
#ifdef EXCAVATION
			int gp_r, gp_s, gp_t;
			double z = 0.0;
			double xyz[3];
			for (gp = 0; gp < NGS; gp++)
			{
				fem_dm->GetGaussData(gp, gp_r, gp_s, gp_t);
				fem_dm->ComputeShapefct(2);
				fem_dm->RealCoordinates(xyz);
				/*
				   //THM2
				   z = 250.0-xyz[1];
				   (*eleV_DM->Stress)(1, gp) = -2360*9.81*z;
				   (*eleV_DM->Stress)(2, gp) = 0.5*(*eleV_DM->Stress)(1, gp);
				   (*eleV_DM->Stress)(0, gp) = 0.6*(*eleV_DM->Stress)(1, gp);
				 */

				// THM1
				z = 500 - xyz[2];  // 3D xyz[1]; //2D
				(*eleV_DM->Stress)(2, gp) = -(0.02 * z + 0.6) * 1.0e6;
				(*eleV_DM->Stress)(1, gp) = -2700 * 9.81 * z;
				(*eleV_DM->Stress)(0, gp) = -(0.055 * z + 4.6) * 1.0e6;

				if (eleV_DM->Stress_j)
					(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
			}
#endif
			////////////////////////////////////////////////////////
			elem->SetOrder(false);
		}
	}
	// Reload the stress results of the previous simulation
	// if(reload >= 2)
	if (idata_type == read_all_binary || idata_type == read_write)
	{
		ReadGaussPointStress();
		if (type == 41)  // mono-deformation-liquid
			ReadSolution();
	}
	// For excavation simulation. Moved here on 05.09.2007 WW
	if (num_type_name.find("EXCAVATION") != 0) Extropolation_GaussValue();
	//
}
/*************************************************************************
   ROCKFLOW - Function: Calculations of initial stress and released load
   Programming:
   09/2007 WW
 **************************************************************************/
void CRFProcessDeformation::CreateInitialState4Excavation()
{
	size_t i;
	int j;
	int Idx_Strain[9];
	int NS = 4;
	if (num_type_name.find("EXCAVATION") != 0) return;
	//
	Idx_Strain[0] = GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = GetNodeValueIndex("STRAIN_XY");

	if (problem_dimension_dm == 3)
	{
		NS = 6;
		Idx_Strain[4] = GetNodeValueIndex("STRAIN_XZ");
		Idx_Strain[5] = GetNodeValueIndex("STRAIN_YZ");
	}
	Idx_Strain[NS] = GetNodeValueIndex("STRAIN_PLS");
	// For excavation simulation. Moved here on 05.09.2007 WW
	if ((idata_type == write_all_binary || idata_type == none) &&
	    reload != -1000)
	//	if(reload < 2 && reload != -1000)
	{
		GravityForce = true;
		cout << "\n ***Excavation simulation: 1. Establish initial stress "
		        "profile..." << endl;
		counter = 0;
		Execute(0);
	}
	else
		UpdateInitialStress(true);  // s0 = 0
	//
	Extropolation_GaussValue();
	//
	cout << "\n ***Excavation simulation: 2. Excavating..." << endl;
	counter = 0;
	InitializeNewtonSteps(true);
	GravityForce = false;
	//
	ReleaseLoadingByExcavation();
	// GravityForce = true;
	UpdateInitialStress(false);  // s-->s0
	m_msh->ConnectedElements2Node();
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
		for (j = 0; j < NS + 1; j++)
			SetNodeValue(i, Idx_Strain[j], 0.0);

	if (reload == -1000) reload = 1;
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::InitializeStress_EachCouplingStep()
   Programming:
   12/2005 WW
 **************************************************************************/
void CRFProcessDeformation::ResetCouplingStep()
{
	long i, e;
	int j;
	long number_of_nodes;
	long shift = 0;
	ElementValue_DM* eleV_DM = NULL;
	for (e = 0; e < (long)m_msh->ele_vector.size(); e++)
		if (m_msh->ele_vector[e]->GetMark())
		{
			eleV_DM = ele_value_dm[e];
			eleV_DM->ResetStress(true);
		}
	shift = 0;
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
			SetNodeValue(j, p_var_index[i], ARRAY[shift + j]);
		shift += number_of_nodes;
	}
}
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::InitializeStress_EachCouplingStep()
   Programming:
   12/2005 WW
 **************************************************************************/
void CRFProcessDeformation::ResetTimeStep()
{
	long e;
	ElementValue_DM* eleV_DM = NULL;
	for (e = 0; e < (long)m_msh->ele_vector.size(); e++)
		if (m_msh->ele_vector[e]->GetMark())
		{
			eleV_DM = ele_value_dm[e];
			eleV_DM->ResetStress(false);
		}
}

/*************************************************************************
   ROCKFLOW - Funktion: TransferNodeValuesToVectorLinearSolver

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E LINEAR_SOLVER *ls: Zeiger auf eine Instanz vom Typ LINEAR_SOLVER.

   Programmaenderungen:
   02/2000   OK   aus TransferNodeValuesToVectorLinearSolver abgeleitet
   07/2005   WW  aus  TransferNodeValuesToVectorLinearSolver(OK)
   11/2010   WW  Modification for H2M
*************************************************************************/
void CRFProcessDeformation::SetInitialGuess_EQS_VEC()
{
	int i;
	long j, v_idx = 0;
	long number_of_nodes;
	long shift = 0;
	double* eqs_x = NULL;
#if defined(USE_PETSC)  // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)
	eqs_x = eqs_new->x;
#else
	eqs_x = eqs->x;
#endif
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		v_idx = p_var_index[i];
		if (i < problem_dimension_dm)
		{
			v_idx--;
			for (j = 0; j < number_of_nodes; j++)
				eqs_x[shift + j] = GetNodeValue(j, v_idx);
		}
		else
			for (j = 0; j < number_of_nodes; j++)
				eqs_x[shift + j] = 0.;
		shift += number_of_nodes;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: UpdateIterativeStep(LINEAR_SOLVER * ls, const int Type)

   Aufgabe:
   Update solution in Newton-Raphson procedure

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver
   const double damp : damping for Newton-Raphson method
   const int type    : 0,  update w=w+dw
                       1,  update u=u+w

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
void CRFProcessDeformation::UpdateIterativeStep(const double damp,
                                                const int u_type)
{
	int i, j;
	long shift = 0;
	long number_of_nodes;
	int ColIndex = 0;
	double* eqs_x = NULL;

#if defined(USE_PETSC)  // || defined (other parallel solver lib). 04.2012 WW
	eqs_x = eqs_new->GetGlobalSolution();
#elif defined(NEW_EQS)
	eqs_x = eqs_new->x;
#else
	eqs_x = eqs->x;
#endif

	//
	for (i = 0; i < problem_dimension_dm; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		//
		ColIndex = p_var_index[i] - 1;
		///  Update Newton step: w = w+dw
		if (u_type == 0)
		{
			for (j = 0; j < number_of_nodes; j++)
			{
#ifdef USE_PETSC
				long k = m_msh->Eqs2Global_NodeIndex[j] *
				             pcs_number_of_primary_nvals +
				         i;
				SetNodeValue(j, ColIndex,
				             GetNodeValue(j, ColIndex) + eqs_x[k] * damp);
#else
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) +
				                              eqs_x[j + shift] * damp);
#endif
			}
			shift += number_of_nodes;
		}
		else
			for (j = 0; j < number_of_nodes; j++)
			{
				SetNodeValue(j, ColIndex + 1, GetNodeValue(j, ColIndex + 1) +
				                                  GetNodeValue(j, ColIndex));
			}
	}

	// if(type == 42&&m_num->nls_method>0)         //H2M, Newton-Raphson.
	// 06.09.2010. WW
	if (type / 10 == 4)  // H2M, HM. 28.09.2011. WW
	{
		/// $p_{n+1}=p_{n+1}+\Delta p$ is already performed when type = 0
		if (u_type == 1) return;

		for (i = problem_dimension_dm; i < pcs_number_of_primary_nvals; i++)
		{
			number_of_nodes = num_nodes_p_var[i];
			//
			ColIndex = p_var_index[i];

			for (j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) +
				                              eqs_x[j + shift] * damp);

			shift += number_of_nodes;
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: InitializeNewtonSteps(LINEAR_SOLVER * ls)

   Aufgabe:
   Initialize the incremental unknows in Newton-Raphson procedure

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E:
   LINEAR_SOLVER * ls: linear solver
   const int type    : 0,  update w=0 (u0=0)
                       1,  update u=0 (u1=0)

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
   06/2007   WW   Rewrite
**************************************************************************/
void CRFProcessDeformation::InitializeNewtonSteps(const bool ini_excav)
{
	long i, j;
	long number_of_nodes;
	int col0, Col = 0, start, end;
	//
	//
	start = 0;
	end = pcs_number_of_primary_nvals;
	//

	/// u_0 = 0
	if (type == 42)  // H2M
		end = problem_dimension_dm;

	/// Dynamic: plus p_0 = 0
	if (type == 41)
	{
		// p_1 = 0
		for (i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			Col = p_var_index[i];
			col0 = Col - 1;
			number_of_nodes = num_nodes_p_var[i];
			if (i < problem_dimension_dm)
				for (j = 0; j < number_of_nodes; j++)
					// SetNodeValue(j, Col, 0.0);
					SetNodeValue(j, col0, 0.0);

			else
			{
				if (FiniteElement::isNewtonKind(
				        m_num->nls_method))  // If newton. 29.09.2011. WW
					continue;

				for (j = 0; j < number_of_nodes; j++)
					SetNodeValue(j, Col, 0.0);
			}
		}
	}
	else  // non HM monolithic
	{
		for (i = start; i < end; i++)
		{
			Col = p_var_index[i] - 1;
			number_of_nodes = num_nodes_p_var[i];
			for (j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, Col, 0.0);
		}
	}
	/// Excavation: plus u_1 = 0;
	if (ini_excav)
		// p_1 = 0
		for (i = 0; i < problem_dimension_dm; i++)
		{
			Col = p_var_index[i];
			number_of_nodes = num_nodes_p_var[i];
			for (j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, Col, 0.0);
		}
}

/**************************************************************************
   ROCKFLOW - Funktion: NormOfUpdatedNewton

   Aufgabe:
   Compute the norm of Newton increment

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   12/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
double CRFProcessDeformation::NormOfUpdatedNewton()
{
	int i, j;
	long number_of_nodes;
	double NormW = 0.0;
	double val;
	int Colshift = 1;
	//
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
		{
			val = GetNodeValue(j, p_var_index[i] - Colshift);
			NormW += val * val;
		}
	}
	return sqrt(NormW);
}

/**************************************************************************
   ROCKFLOW - Funktion: StoreDisplacement

   Aufgabe:
   Copy the displacement of the previous time interval to a vector
   temporarily

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/

void CRFProcessDeformation::StoreLastSolution(const int ty)
{
	int i, j;
	long number_of_nodes;
	long shift = 0;

	// Displacement
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
			ARRAY[shift + j] = GetNodeValue(j, p_var_index[i] - ty);
		shift += number_of_nodes;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: RetrieveDisplacement(LINEAR_SOLVER * ls)

   Aufgabe:
   Retrive the displacement from the temporary array
   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
void CRFProcessDeformation::RecoverSolution(const int ty)
{
	int i, j, idx;
	long number_of_nodes;
	int Colshift = 1;
	long shift = 0;
	double tem = 0.0;

	int start, end;

	start = 0;
	end = pcs_number_of_primary_nvals;

	// If monolithic scheme for p-u coupling,  p_i-->p_0 only
	if (pcs_deformation % 11 == 0 && ty > 0)
	{
		start = problem_dimension_dm;
		for (i = 0; i < start; i++)
			shift += num_nodes_p_var[i];

		// TODO: end = problem_dimension_dm;
	}
	for (i = start; i < end; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		idx = p_var_index[i] - Colshift;
		for (j = 0; j < number_of_nodes; j++)
		{
			if (ty < 2)
			{
				if (ty == 1) tem = GetNodeValue(j, idx);
				SetNodeValue(j, idx, ARRAY[shift + j]);
				if (ty == 1) ARRAY[shift + j] = tem;
			}
			else if (ty == 2)
			{
				tem = ARRAY[shift + j];
				ARRAY[shift + j] = GetNodeValue(j, idx);
				SetNodeValue(j, idx, tem);
			}
		}
		shift += number_of_nodes;
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: NormOfDisp

   Aufgabe:
   Compute the norm of  u_{n+1}
   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: LINEAR_SOLVER * ls: linear solver

   Ergebnis:
   - double - Eucleadian Norm

   Programmaenderungen:
   10/2002   WW   Erste Version
   11/2007   WW   Change to fit the new equation class
**************************************************************************/
double CRFProcessDeformation::NormOfDisp()
{
	int i, j;
	long number_of_nodes;
	double Norm1 = 0.0;
	//
	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		number_of_nodes = num_nodes_p_var[i];
		for (j = 0; j < number_of_nodes; j++)
			Norm1 += GetNodeValue(j, p_var_index[i]) *
			         GetNodeValue(j, p_var_index[i]);
	}
	return Norm1;
}


/**************************************************************************
   ROCKFLOW - Funktion: MaxiumLoadRatio

   Aufgabe:
   Calculate the muxium effective stress, Smax, of all Gauss points.
   (For 2-D 9 nodes element only up to now). Then compute the maxium ration
   by:

   Smax/Y0(initial yield stress)

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - const int NodesOfEelement:

   Ergebnis:
   - void -

   Programmaenderungen:
   01/2003  WW  Erste Version

   letzte Aenderung:

**************************************************************************/
//#define Modified_B_matrix
double CRFProcessDeformation::CaclMaxiumLoadRatio(void)
{
	int j, gp, gp_r, gp_s;  //, gp_t;
	int PModel = 1;
	long i = 0;
	double* dstrain;

	double S0 = 0.0, p0 = 0.0;
	double MaxS = 0.000001;
	double EffS = 0.0;

	MeshLib::CElem* elem = NULL;
	ElementValue_DM* eleV_DM = NULL;
	CSolidProperties* SMat = NULL;

	// Weimar's model
	Matrix* Mat = NULL;
	double II = 0.0;
	double III = 0.0;

	double PRatio = 0.0;
	const double MaxR = 20.0;

	int NGS, NGPS;

	// gp_t = 0;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			fem_dm->ConfigElement(elem);
			fem_dm->SetMaterial();
			eleV_DM = ele_value_dm[i];
			SMat = fem_dm->smat;
			SMat->axisymmetry = m_msh->isAxisymmetry();
			PModel = SMat->Plasticity_type;
			//
			switch (PModel)
			{
				case 1:
#ifdef RFW_FRACTURE
					SMat->Calculate_Lame_Constant(elem);
#endif
#ifndef RFW_FRACTURE
					SMat->Calculate_Lame_Constant();
#endif
					SMat->ElasticConsitutive(fem_dm->Dim(), fem_dm->De);
					SMat->CalulateCoefficent_DP();
					S0 = MSqrt2Over3 * SMat->BetaN * SMat->Y0;
					break;
				case 2:
#ifdef RFW_FRACTURE
					SMat->Calculate_Lame_Constant(elem);
#endif
#ifndef RFW_FRACTURE
					SMat->Calculate_Lame_Constant();
#endif
					SMat->ElasticConsitutive(fem_dm->Dim(), fem_dm->De);
					Mat = eleV_DM->MatP;
					break;
				case 3:
					Mat = SMat->data_Plasticity;
					S0 = (*Mat)(3);
					break;
			}
			NGS = fem_dm->GetNumGaussPoints();
			NGPS = fem_dm->GetNumGaussSamples();
			//
			for (gp = 0; gp < NGS; gp++)
			{
				switch (elem->GetElementType())
				{
					case MshElemType::TRIANGLE:  // Triangle
						SamplePointTriHQ(gp, fem_dm->unit);
						break;
					case MshElemType::QUAD:  // Quadralateral
						gp_r = (int)(gp / NGPS);
						gp_s = gp % NGPS;
						fem_dm->unit[0] = MXPGaussPkt(NGPS, gp_r);
						fem_dm->unit[1] = MXPGaussPkt(NGPS, gp_s);
						break;
					default:
						std::cerr
						    << "CRFProcessDeformation::CaclMaxiumLoadRatio "
						       "MshElemType not handled"
						    << "\n";
				}
				fem_dm->computeJacobian(2);
				fem_dm->ComputeGradShapefct(2);
				fem_dm->ComputeStrain();

				dstrain = fem_dm->GetStrain();

				if (PModel == 3)  // Cam-Clay
				{
					p0 =
					    ((*eleV_DM->Stress)(0, gp) + (*eleV_DM->Stress)(1, gp) +
					     (*eleV_DM->Stress)(2, gp)) /
					    3.0;
					// Swelling index: (*SMat->data_Plasticity)(2)
					if (fabs(p0) < MKleinsteZahl)
						// The initial preconsolidation pressure
						p0 = (*SMat->data_Plasticity)(3);

					SMat->K = (1.0 + (*eleV_DM->e_i)(gp)) * fabs(p0) /
					          (*SMat->data_Plasticity)(2);
					SMat->G = 1.5 * SMat->K * (1 - 2.0 * SMat->PoissonRatio) /
					          (1 + SMat->PoissonRatio);
					SMat->Lambda = SMat->K - 2.0 * SMat->G / 3.0;
					SMat->ElasticConsitutive(fem_dm->Dim(), fem_dm->De);
				}

				// Stress of the previous time step
				for (j = 0; j < fem_dm->ns; j++)
					fem_dm->dstress[j] = (*eleV_DM->Stress)(j, gp);

				// Compute try stress, stress incremental:
				fem_dm->De->multi(dstrain, fem_dm->dstress);

				p0 = DeviatoricStress(fem_dm->dstress) / 3.0;

				switch (PModel)
				{
					case 1:  // Drucker-Prager model
						EffS = sqrt(TensorMutiplication2(fem_dm->dstress,
						                                 fem_dm->dstress,
						                                 fem_dm->Dim())) +
						       3.0 * SMat->Al * p0;

						if (EffS > S0 && EffS > MaxS &&
						    fabs(S0) > MKleinsteZahl)
						{
							MaxS = EffS;
							PRatio = MaxS / S0;
						}
						break;

					case 2:  // Single yield surface
						// Compute try stress, stress incremental:
						II = TensorMutiplication2(
						    fem_dm->dstress, fem_dm->dstress, fem_dm->Dim());
						III = TensorMutiplication3(fem_dm->dstress,
						                           fem_dm->dstress,
						                           fem_dm->dstress,
						                           fem_dm->Dim());
						p0 *= 3.0;
						EffS =
						    sqrt(II * pow(1.0 + (*Mat)(5) * III / pow(II, 1.5),
						                  (*Mat)(6)) +
						         0.5 * (*Mat)(0) * p0 * p0 +
						         (*Mat)(2) * (*Mat)(2) * p0 * p0 * p0 * p0) +
						    (*Mat)(1) * p0 + (*Mat)(3) * p0* p0;

						if (EffS > (*Mat)(4))
						{
							if ((*Mat)(4) > 0.0)
							{
								if (EffS > MaxS) MaxS = EffS;
								PRatio = MaxS / (*Mat)(4);
								if (PRatio > MaxR) PRatio = MaxR;
							}
							else
								PRatio = EffS;
						}
						break;

					case 3:  // Cam-Clay
						II = 1.5 * TensorMutiplication2(fem_dm->dstress,
						                                fem_dm->dstress,
						                                fem_dm->Dim());
						if (S0 > 0.0)
						{
							EffS = II / (p0 * (*Mat)(0) * (*Mat)(0)) + p0;
							if (EffS > S0)
								PRatio = EffS / S0;
							else
								PRatio = 1.0;
						}
						else
							PRatio = 1.0;
						break;
				}
			}
		}
	}

	return PRatio;
}

/**************************************************************************
   ROCKFLOW - Funktion: Extropolation_GaussValue

   Aufgabe:
   Calculate the stresses of element nodes using the values at Gauss points.

   Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - const int NodesOfEelement:

   Ergebnis:
   - void -

   Programmaenderungen:
   10/2002  WW  Erste Version
   07/2003  WW  Extroplolation in quadraitc triangle element is added
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::Extropolation_GaussValue()
{
	int k, NS;
	long i = 0;
	int Idx_Stress[7];
	const long LowOrderNodes = m_msh->GetNodesNumber(false);
	MeshLib::CElem* elem = NULL;

	// Clean nodal stresses
	NS = 4;
	Idx_Stress[0] = GetNodeValueIndex("STRESS_XX");
	Idx_Stress[1] = GetNodeValueIndex("STRESS_YY");
	Idx_Stress[2] = GetNodeValueIndex("STRESS_ZZ");
	Idx_Stress[3] = GetNodeValueIndex("STRESS_XY");
	if (problem_dimension_dm == 3)
	{
		NS = 6;
		Idx_Stress[4] = GetNodeValueIndex("STRESS_XZ");
		Idx_Stress[5] = GetNodeValueIndex("STRESS_YZ");
	}
	Idx_Stress[NS] = GetNodeValueIndex("STRAIN_PLS");
	NS++;
	for (i = 0; i < LowOrderNodes; i++)
		for (k = 0; k < NS; k++)
			SetNodeValue(i, Idx_Stress[k], 0.0);

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			fem_dm->ConfigElement(elem);
			fem_dm->SetMaterial();
			//         eval_DM = ele_value_dm[i];
			// TEST        (*eval_DM->Stress) += (*eval_DM->Stress0);
			fem_dm->ExtropolateGuassStress();
			// TEST        if(!update)
			//           (*eval_DM->Stress) -= (*eval_DM->Stress0);
		}
	}
}

/*--------------------------------------------------------------------------
   Trace discontinuity path. Belong to Geometry
   --------------------------------------------------------------------------*/
void CRFProcessDeformation::Trace_Discontinuity()
{
	long k, l;
	int i, nn, Size, bFaces, bFacesCounter, intP;
	int b_node_counter;
	int locEleFound, numf, nPathNodes;
	bool thisLoop;  //, neighborSeed;

	// Element value
	ElementValue_DM* eleV_DM;

	static double xn[20], yn[20], zn[20];
	static double xa[3], xb[3], n[3], ts[3];
	double v1, v2;

	int FNodes0[8];
	std::vector<long> SeedElement;

	locEleFound = 0;
	nPathNodes = 2;  // 2D element
	bFaces = 0;

	intP = 0;

	// Check all element for bifurcation
	MeshLib::CElem* elem = NULL;
	for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
	{
		elem = m_msh->ele_vector[l];
		if (elem->GetMark())  // Marked for use

			if (fem_dm->LocalAssembly_CheckLocalization(elem))
			{
				locEleFound++;
				// If this is first bifurcated element, call them as seeds
				if (!Localizing) SeedElement.push_back(l);
			}
	}

	if (locEleFound > 0 && !Localizing)  // Bifurcation inception
	{
		// TEST
		// mesh1   de =23;
		// mesh2_iregular de = 76
		// mesh coarst de = 5;
		// mesh quad de=23
		// crack tri de=0

		/*
		   SeedElement.clear();
		   int de = 39; //64; //itri  //39; //Quad //72; //rtri crack
		   SeedElement.push_back(de);
		 */

		// TEST

		// Determine the seed element
		Size = (long)SeedElement.size();
		for (l = 0; l < Size; l++)
		{
			k = SeedElement[l];
			elem = m_msh->ele_vector[k];

			numf = elem->GetFacesNumber();
			eleV_DM = ele_value_dm[k];

			// If seed element are neighbor. Choose one
			/*//TEST
			   neighborSeed = false;
			   for(m=0; m<Size; m++)
			   {
			   if(m==l) continue;
			   for(i=0; i<numf; i++)
			   {
			       if(neighbor[i]==SeedElement[m])
			       {
			          neighborSeed = true;
			          break;
			   }
			   }
			   }
			   if(neighborSeed)
			   {
			   delete eleV_DM->orientation;
			   eleV_DM->orientation = NULL;
			   continue;
			   }*/
			//

			nn = elem->GetNodesNumber(true);
			for (i = 0; i < nn; i++)
			{
				// Coordinates of all element nodes
				//               xn[i] = elem->nodes[i]->X();
				//               yn[i] = elem->nodes[i]->Y();
				//               zn[i] = elem->nodes[i]->Z();
				double const* const coords(elem->GetNode(i)->getData());
				xn[i] = coords[0];
				yn[i] = coords[1];
				zn[i] = coords[2];
			}
			// Elements which have one only boundary face are chosen as seed
			// element
			bFaces = -1;
			bFacesCounter = 0;
			for (i = 0; i < numf; i++)
				if (elem->GetNeighbor(i)->GetDimension() != elem->GetDimension())
				{
					bFaces = i;
					bFacesCounter++;
				}

			// Elements which have only one boundary face or one boundary node
			// are chosen as seed element
			if (bFacesCounter != 1)
			{
				//
				b_node_counter = 0;
				for (i = 0; i < elem->GetVertexNumber(); i++)
				{
					bFaces = i;
					b_node_counter++;
				}
				if (b_node_counter != 1)
				{
					eleV_DM->Localized = false;
					delete eleV_DM->orientation;
					eleV_DM->orientation = NULL;
					continue;
				}
			}

			fem_dm->ConfigElement(elem);
			// 2D
			elem->GetElementFaceNodes(bFaces, FNodes0);
			if (elem->GetElementType() == MshElemType::QUAD ||
			    elem->GetElementType() == MshElemType::TRIANGLE)
				nPathNodes = 2;
			// Locate memory for points on the path of this element
			eleV_DM->NodesOnPath = new Matrix(3, nPathNodes);
			*eleV_DM->NodesOnPath = 0.0;
			if (nPathNodes == 2)  // 2D
			{
				// Departure point
				(*eleV_DM->NodesOnPath)(0, 0) =
				    0.5 * (xn[FNodes0[0]] + xn[FNodes0[1]]);
				(*eleV_DM->NodesOnPath)(1, 0) =
				    0.5 * (yn[FNodes0[0]] + yn[FNodes0[1]]);
				(*eleV_DM->NodesOnPath)(2, 0) =
				    0.5 * (zn[FNodes0[0]] + zn[FNodes0[1]]);

				xa[0] = (*eleV_DM->NodesOnPath)(0, 0);
				xa[1] = (*eleV_DM->NodesOnPath)(1, 0);
				xa[2] = (*eleV_DM->NodesOnPath)(2, 0);

				// Check oreintation again.
				ts[0] = xn[FNodes0[1]] - xn[FNodes0[0]];
				ts[1] = yn[FNodes0[1]] - yn[FNodes0[0]];
				// ts[2] = zn[FNodes0[1]]-zn[FNodes0[0]];
				v1 = sqrt(ts[0] * ts[0] + ts[1] * ts[1]);
				ts[0] /= v1;
				ts[1] /= v1;

				n[0] = cos(eleV_DM->orientation[0]);
				n[1] = sin(eleV_DM->orientation[0]);
				v1 = n[0] * ts[0] + n[1] * ts[1];
				n[0] = cos(eleV_DM->orientation[1]);
				n[1] = sin(eleV_DM->orientation[1]);
				v2 = n[0] * ts[0] + n[1] * ts[1];
				if (fabs(v2) > fabs(v1))
				{
					v1 = eleV_DM->orientation[0];
					eleV_DM->orientation[0] = eleV_DM->orientation[1];
					eleV_DM->orientation[1] = v1;
				}

				intP = fem_dm->IntersectionPoint(bFaces, xa, xb);

				(*eleV_DM->NodesOnPath)(0, 1) = xb[0];
				(*eleV_DM->NodesOnPath)(1, 1) = xb[1];
				(*eleV_DM->NodesOnPath)(2, 1) = xb[2];
			}

			// Last element to this seed
			DisElement* disEle = new DisElement;
			disEle->NumInterFace = 1;
			disEle->ElementIndex = k;
			disEle->InterFace = new int[1];
			disEle->InterFace[0] = intP;
			LastElement.push_back(disEle);
			ElementOnPath.push_back(k);
		}

		Localizing = true;
	}

	// Seek path from the last bifurcated element of corrsponding seeds.
	if (Localizing)
	{
		Size = (long)LastElement.size();
		for (i = 0; i < Size; i++)
		{
			thisLoop = true;
			while (thisLoop)
			{
				nn = MarkBifurcatedNeighbor(i);
				if (nn < 0) break;
				ElementOnPath.push_back(nn);
			}
		}

		Size = (long)ElementOnPath.size();
		for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
		{
			elem = m_msh->ele_vector[l];
			if (elem->GetMark())  // Marked for use
			{
				eleV_DM = ele_value_dm[l];
				if (eleV_DM->Localized)
				{
					thisLoop = false;
					for (k = 0; k < Size; k++)
						if (l == ElementOnPath[k])
						{
							thisLoop = true;
							break;
						}
					if (!thisLoop)
					{
						eleV_DM->Localized = false;
						delete eleV_DM->orientation;
						eleV_DM->orientation = NULL;
					}
				}
			}
		}
	}
}

// WW
long CRFProcessDeformation::MarkBifurcatedNeighbor(const int PathIndex)
{
	int j;
	int f1, f2, nb, numf1;
	long index, Extended;
	bool adjacent, onPath;
	ElementValue_DM* eleV_DM, *eleV_DM1;
	DisElement* disEle;
	static double n1[2], n2[2], xA[3], xB[3];
	// WW static int Face_node[8];                    // Only 2D
	MeshLib::CElem* elem;
	MeshLib::CElem* elem1;

	double pd1, pd2;

	Extended = -1;
	// 2D only
	disEle = LastElement[PathIndex];
	index = disEle->ElementIndex;
	f1 = disEle->InterFace[0];

	elem = m_msh->ele_vector[index];
	eleV_DM = ele_value_dm[index];

	// numf = elem->GetFacesNumber();

	n1[0] = cos(eleV_DM->orientation[0]);
	n1[1] = sin(eleV_DM->orientation[0]);

	xA[0] = (*eleV_DM->NodesOnPath)(0, 1);
	xA[1] = (*eleV_DM->NodesOnPath)(1, 1);
	xA[2] = (*eleV_DM->NodesOnPath)(2, 1);

	// nfnode = elem->GetElementFaceNodes(f1, Face_node);

	// Check discintinuity path goes to which neighbor
	elem1 = elem->GetNeighbor(f1);
	// Boundary reached
	if (elem1->GetDimension() != elem->GetDimension()) return -1;
	nb = elem1->GetIndex();
	// Check if the element is already in discontinuity line/surface
	onPath = false;
	for (j = 0; j < (int)ElementOnPath.size(); j++)
		if (nb == ElementOnPath[j])
		{
			onPath = true;
			break;
		}

	// If has neighbor and it is not on the discontinuity surface.
	if (!onPath)  // Has neighbor
		if (ele_value_dm[nb]->Localized)
		{
			// TEST OUT
			// cout <<" element on track  " <<nb<<endl;

			adjacent = false;
			numf1 = elem1->GetFacesNumber();
			eleV_DM1 = ele_value_dm[nb];
			fem_dm->ConfigElement(elem1);
			// Search faces of neighbor's neighbors
			for (j = 0; j < numf1; j++)
			{
				// Neighbor is a face on surface
				if (elem1->GetNeighbor(j)->GetDimension() !=
				    elem1->GetDimension())
					continue;
				if ((size_t)index != elem1->GetNeighbor(j)->GetIndex()) continue;
				{
					adjacent = true;
					Extended = nb;
					// Choose a smooth direction
					n2[0] = cos(eleV_DM1->orientation[0]);
					n2[1] = sin(eleV_DM1->orientation[0]);
					pd1 = n1[0] * n2[0] + n1[1] * n2[1];
					n2[0] = cos(eleV_DM1->orientation[1]);
					n2[1] = sin(eleV_DM1->orientation[1]);
					pd2 = n1[0] * n2[0] + n1[1] * n2[1];
					if (pd2 > pd1)  // Always use the first entry of orientation
					{
						// Swap the values
						pd1 = eleV_DM1->orientation[1];
						eleV_DM1->orientation[1] = eleV_DM1->orientation[0];
						eleV_DM1->orientation[0] = pd1;
					}
					eleV_DM1->NodesOnPath = new Matrix(3, 2);
					*eleV_DM1->NodesOnPath = 0.0;
					// Get another intersection point
					f2 = fem_dm->IntersectionPoint(j, xA, xB);
					(*eleV_DM1->NodesOnPath)(0, 0) = xA[0];
					(*eleV_DM1->NodesOnPath)(1, 0) = xA[1];
					(*eleV_DM1->NodesOnPath)(2, 0) = xA[2];
					(*eleV_DM1->NodesOnPath)(0, 1) = xB[0];
					(*eleV_DM1->NodesOnPath)(1, 1) = xB[1];
					(*eleV_DM1->NodesOnPath)(2, 1) = xB[2];

					// Last element
					disEle->ElementIndex = nb;
					disEle->NumInterFace = 1;
					disEle->InterFace[0] = f2;
					break;
				}
			}

			// If not on the discontinuity surface
			// Release memory
			if (!adjacent)
			{
				delete eleV_DM1->orientation;
				eleV_DM1->orientation = NULL;
				eleV_DM1->Localized = false;
			}
		}

	// If true, discontinuity extended to its neighbor
	return Extended;
}

/**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices and RHS for each element
   Programing:
   02/2005 WW
**************************************************************************/
void CRFProcessDeformation::GlobalAssembly()
{
	// WW
	{
		GlobalAssembly_DM();

		if (type / 10 == 4)
		{  // p-u monolithic scheme

			// if(!fem_dm->dynamic)   ///
			//  RecoverSolution(1);  // p_i-->p_0
			// 2.
			// Assemble pressure eqs
			// Changes for OpenMP
			GlobalAssembly_std(true);
#if 0
            const size_t n_nodes_linear = m_msh->GetNodesNumber(false);
            const size_t n_nodes_quard = m_msh->GetNodesNumber(true);
            const size_t offset_H = problem_dimension_dm * n_nodes_quard;
			if (this->eqs_new->size_A > offset_H + n_nodes_linear)
			{
	            // set dummy diagonal entry of rows corresponding to unused quadratic nodes for H
	            std::cout << "set dummy diagonal entry of rows corresponding to unused quadratic nodes for H\n";
	            std::cout << "-> Linear nodes = " << n_nodes_linear << ", Quadratic nodes = " << n_nodes_quard << "\n";
	            std::cout << "-> Constrain equation index from " << offset_H +  n_nodes_linear << " to " << offset_H + n_nodes_quard << "\n";
	            for (size_t i=n_nodes_linear; i<n_nodes_quard; i++) {
	                (*this->eqs_new->A)(offset_H+i,offset_H+i)=1.0;
	            }
			}
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03~04.3012. WW
            eqs_new->EQSV_Viewer("eqs" + number2str(aktueller_zeitschritt) + "a");
#endif
#endif
		}
// if(!fem_dm->dynamic)
//   RecoverSolution(2);  // p_i-->p_0

//----------------------------------------------------------------------
//
// {			 MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); // abort();}

// DumpEqs("rf_pcs1.txt");

#if 0
            {
		   ofstream Dum(std::string("eqs_after_assembly.txt").c_str(), ios::out); // WW
		   this->eqs_new->Write(Dum);
		   Dum.close();
            }
#endif
		// Apply Neumann BC
		IncorporateSourceTerms();
// DumpEqs("rf_pcs2.txt");

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012.
		ScreenMessage2d("assemble PETSc matrix and vectors...\n");
		eqs_new->AssembleUnkowns_PETSc();
		eqs_new->AssembleRHS_PETSc();
		eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);
//		eqs_new->EQSV_Viewer("eqs_after_assembl");
#endif

		// {			MXDumpGLS("rf_pcs1.txt",1,eqs->b,eqs->x); // abort();}
		//#if defined(USE_PETSC)  // || defined(other parallel
		// libs)//03~04.3012.
		////		eqs_new->EQSV_Viewer("eqs_after_ST");
		//		eqs_new->AssembleRHS_PETSc();
		//#endif

		// Apply Dirchlete bounday condition
		IncorporateBoundaryConditions();
//  {			 MXDumpGLS("rf_pcs_dm1.txt",1,eqs->b,eqs->x);  //abort();}
//

#if 0
            {
           ofstream Dum(std::string("eqs_after_BCST.txt").c_str(), ios::out); // WW
           this->eqs_new->Write(Dum);
           Dum.close();
            }
#endif

#define atest_dump
#ifdef test_dump
		string fname = FileName + "rf_pcs_omp.txt";
		ofstream Dum1(fname.c_str(), ios::out);  // WW
		eqs_new->Write(Dum1);
		Dum1.close();  //   abort();
#endif

#define atest_bin_dump
#ifdef test_bin_dump  // WW
		string fname = FileName + ".eqiation_binary.bin";

		ofstream Dum1(fname.data(), ios::out | ios::binary | ios::trunc);
		if (Dum1.good()) eqs_new->Write_BIN(Dum1);
		Dum1.close();
#endif
		//
	}
	ScreenMessage("Global assembly is done\n");
}

/*!  \brief Assembe matrix and vectors
      for deformation process

      24.11.2010. WW
 */
void CRFProcessDeformation::GlobalAssembly_DM()
{
	long i;
	MeshLib::CElem* elem = NULL;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark())  // Marked for use
			continue;

		elem->SetOrder(true);
		fem_dm->ConfigElement(elem);
		fem_dm->LocalAssembly(0);
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Update stresses and straines at each Gauss points
   Argument:
   Programing:
   02/2005 WW
   06/2005 WW  Parallelization
**************************************************************************/
void CRFProcessDeformation::UpdateStress()
{
	long i;
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			elem->SetOrder(true);
			fem_dm->ConfigElement(elem);
			fem_dm->LocalAssembly(1);
		}
	}
	//}
}

// Coupling
/*
   void CRFProcessDeformation::ConfigureCoupling()
   {
    fem_dm->ConfigureCoupling(this, Shift);
   }
 */

std::string CRFProcessDeformation::GetGaussPointStressFileName()
{
	std::string m_file_name = FileName;
#if defined(USE_PETSC)
	m_file_name += "_rank" + number2str(myrank);
#endif
	m_file_name += ".sts";
	return m_file_name;
}

/**************************************************************************
   ROCKFLOW - Funktion: WriteGaussPointStress()

   Aufgabe:
   Write Gauss point stresses to a file

   Programmaenderungen:
   03/2005  WW  Erste Version
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::WriteGaussPointStress()
{
	long i;
	string StressFileName = GetGaussPointStressFileName();
	fstream file_stress(StressFileName.data(), ios::binary | ios::out);
	ElementValue_DM* eleV_DM = NULL;

	long ActiveElements = 0;
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
			ActiveElements++;
	}
	file_stress.write((char*)(&ActiveElements), sizeof(ActiveElements));
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			eleV_DM = ele_value_dm[i];
			file_stress.write((char*)(&i), sizeof(i));
			//          *eleV_DM->Stress_i += *eleV_DM->Stress0;
			// TEST           *eleV_DM->Stress0 = 0.0;
			*eleV_DM->Stress0 = *eleV_DM->Stress_i;
			eleV_DM->Write_BIN(file_stress);
		}
	}
	//
	file_stress.close();
}
/**************************************************************************
   ROCKFLOW - Funktion: ReadGaussPointStress()

   Aufgabe:
   Read Gauss point stresses

   Programmaenderungen:
   03/2005  WW  Erste Version
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::ReadGaussPointStress()
{
	long i, index, ActiveElements;
	string StressFileName = GetGaussPointStressFileName();
	fstream file_stress(StressFileName.data(), ios::binary | ios::in);
	ElementValue_DM* eleV_DM = NULL;
	//
	file_stress.read((char*)(&ActiveElements), sizeof(ActiveElements));
	for (i = 0; i < ActiveElements; i++)
	{
		file_stress.read((char*)(&index), sizeof(index));
		eleV_DM = ele_value_dm[index];
		eleV_DM->Read_BIN(file_stress);
		(*eleV_DM->Stress0) = (*eleV_DM->Stress);
		if (eleV_DM->Stress_j) (*eleV_DM->Stress_j) = (*eleV_DM->Stress);
	}
	//
	file_stress.close();
}

/**************************************************************************
   ROCKFLOW - Funktion: ReadGaussPointStress()

   Aufgabe:
   Read element-wise stress data

   Programmaenderungen:
   10/2011  WW  Erste Version
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::ReadElementStress()
{
	long i, index, ActiveElements;
	string StressFileName = FileName + ".ele_stress.asc";
	fstream file_stress(StressFileName.data());
	ElementValue_DM* eleV_DM = NULL;
	//
	file_stress >> ActiveElements;
	for (i = 0; i < ActiveElements; i++)
	{
		file_stress >> index;
		eleV_DM = ele_value_dm[index];
		eleV_DM->ReadElementStressASCI(file_stress);
		(*eleV_DM->Stress) = (*eleV_DM->Stress0);
		if (eleV_DM->Stress_j) (*eleV_DM->Stress_j) = (*eleV_DM->Stress);
	}
	//
	file_stress.close();
}

/**************************************************************************
   ROCKFLOW - Funktion: ReleaseLoadingByExcavation()

   Aufgabe:
   Compute the nodal forces produced by excavated body

   Programmaenderungen:
   04/2005  WW  Erste Version
   09/2007  WW  Set as a boundary condition
   letzte Aenderung:

**************************************************************************/
void CRFProcessDeformation::ReleaseLoadingByExcavation()
{
	long i, actElements;
	int j, k, l, SizeSt, SizeSubD;
	ElementValue_DM* ele_val = NULL;

	std::vector<int> ExcavDomainIndex;
	std::vector<long> NodesOnCaveSurface;

	CSourceTerm* m_st = NULL;
	SizeSt = (int)st_vector.size();
	bool exist = false;
	double* eqs_b = NULL;

#if defined(USE_PETSC)  // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)
	eqs_b = eqs_new->b;
#else
	eqs_b = eqs->b;
#endif

	for (k = 0; k < SizeSt; k++)
	{
		m_st = st_vector[k];
		if (m_st->getProcessPrimaryVariable() == FiniteElement::EXCAVATION)
		{
			// ---- 16.01.2009 WW
			exist = false;

			for (j = k + 1; j < SizeSt; j++)
				if (m_st->getSubDomainIndex() ==
				    st_vector[j]->getSubDomainIndex())
				{
					exist = true;
					break;
				}
			if (!exist) ExcavDomainIndex.push_back(m_st->getSubDomainIndex());
		}
	}
	SizeSubD = (int)ExcavDomainIndex.size();
	if (SizeSubD == 0) return;  // 05.09.2007 WW
	exist = false;              // 16.02
	// 1. De-active host domain to be exvacated
	actElements = 0;
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		elem->SetMark(false);
		for (k = 0; k < SizeSubD; k++)
			if (elem->GetPatchIndex() ==
			    static_cast<size_t>(ExcavDomainIndex[k]))
				elem->SetMark(true);
		if (elem->GetMark()) actElements++;
	}
	if (actElements == 0)
	{
		cout << "No element specified for excavation. Please check data in .st "
		        "file " << endl;
		abort();
	}
// 2. Compute the released node loading

	for (i = 0; i < 4; i++)  // In case the domain decomposition is employed
		fem_dm->NodeShift[i] = Shift[i];
	//
	PreLoad = 11;
	LoadFactor = 1.0;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			fem_dm->ConfigElement(elem);
			fem_dm->LocalAssembly(0);
			ele_val = ele_value_dm[i];
			// Clear stresses in excavated domain
			(*ele_val->Stress0) = 0.0;
			(*ele_val->Stress) = 0.0;
			if (ele_val->Stress_j) (*ele_val->Stress_j) = 0.0;
		}
	}

	// 3 --------------------------------------------------------
	// Store the released loads to source term buffer
	long number_of_nodes;
	CNodeValue* m_node_value = NULL;
	std::vector<long> nodes_vector(0);

	number_of_nodes = 0;
	RecordNodeVSize((long)st_node_value.size());

	// TEST
	st_node_value.clear();
	//

	for (k = 0; k < SizeSt; k++)
	{
		// Get nodes on cave surface
		m_st = st_vector[k];
		if (m_st->getProcessPrimaryVariable() != FiniteElement::EXCAVATION)
			continue;
		if (m_st->getGeoType() == GEOLIB::POLYLINE)
		{
			CGLPolyline* m_polyline(GEOGetPLYByName(m_st->getGeoName()));

			// reset the min edge length of mesh
			double mesh_min_edge_length(m_msh->getMinEdgeLength());
			m_msh->setMinEdgeLength(m_polyline->epsilon);

			if (m_st->getGeoObj())
			{
				m_msh->GetNODOnPLY(
				    static_cast<const GEOLIB::Polyline*>(m_st->getGeoObj()),
				    nodes_vector);
				// reset min edge length of mesh
				m_msh->setMinEdgeLength(mesh_min_edge_length);
			}
			m_msh->setMinEdgeLength(mesh_min_edge_length);
		}
		if (m_st->getGeoType() == GEOLIB::SURFACE)
		{
			// CC 10/05
			Surface* m_surface = GEOGetSFCByName(m_st->getGeoName());
			//			 07/2010 TF ToDo: to do away with the global vector
			// surface_vector
			//			                  fetch the geometry from CFEMesh
			//			Surface *m_surface
			//(surface_vector[m_st->getGeoObjIdx()]);
			if (m_surface)
			{
				if (m_surface->type == 100)
					m_msh->GetNodesOnCylindricalSurface(m_surface,
					                                    nodes_vector);
				else
					m_msh->GetNODOnSFC_PLY(m_surface, nodes_vector);
			}
		}
		// Set released node forces from eqs->b;
		number_of_nodes = (int)nodes_vector.size();
		for (j = 0; j < problem_dimension_dm; j++)
			for (i = 0; i < number_of_nodes; i++)
			{
				m_node_value = new CNodeValue();
				m_node_value->msh_node_number = nodes_vector[i] + Shift[j];
				m_node_value->geo_node_number = nodes_vector[i];
				m_node_value->node_value =
				    -eqs_b[m_node_value->geo_node_number + Shift[j]];
				m_node_value->CurveIndex = m_st->CurveIndex;
				// Each node only take once
				exist = false;
				for (l = 0; l < (int)st_node_value[k].size(); l++)
					if (st_node_value[k][l]->msh_node_number ==
					    m_node_value->msh_node_number)
					{
						exist = true;
						break;
					}
				if (!exist) st_node_value[k].push_back(m_node_value);
			}
	}
	//
	// Deactivate the subdomains to be excavated
	Deactivated_SubDomain.resize(SizeSubD);
	for (j = 0; j < SizeSubD; j++)
		Deactivated_SubDomain[j] = ExcavDomainIndex[j];

	// Activate the host domain for excavtion analysis
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark()) elem->SetMark(true);
	}
	PreLoad = 1;
	// TEST OUTPUT
	//   {MXDumpGLS("rf_pcs.txt",1,eqs->b,eqs->x);  abort();}
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::UpdateInitialStress()
   Task:  Compute number of element neighbors to a node
   Dim : Default=2
   Programming:
   12/2003 WW
 **************************************************************************/
void CRFProcessDeformation::UpdateInitialStress(bool ZeroInitialS)
{
	long i;
	ElementValue_DM* eval_DM;

	// Over all elements
	MeshLib::CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			eval_DM = ele_value_dm[i];
			if (ZeroInitialS)
				(*eval_DM->Stress0) = 0.0;
			else
				(*eval_DM->Stress0) = (*eval_DM->Stress);
		}
	}
}
/**************************************************************************
   GEOSYS - Funktion: CalcBC_or_SecondaryVariable_Dynamics(bool BC);
   Programmaenderungen:
   05/2005  WW  Erste Version
   letzte Aenderung:
   09/2011 TF substituted pow by fastpow
**************************************************************************/
bool CRFProcessDeformation::CalcBC_or_SecondaryVariable_Dynamics(bool BC)
{
	const char* function_name[7];
	size_t i;
	long j;
	double v, bc_value, time_fac = 1.0;

	std::vector<int> bc_type;
	long bc_msh_node;
	long bc_eqs_index;
	int interp_method = 0;
	int curve, valid = 0;
	int idx_disp[3], idx_vel[3], idx_acc[3], idx_acc0[3];
	int idx_pre, idx_dpre, idx_dpre0;
	int nv, k;

	size_t Size = m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false);
	CBoundaryCondition* m_bc = NULL;
	bc_type.resize(Size);

	v = 0.0;
	// 0: not given
	// 1, 2, 3: x,y, or z is given
	for (size_t i = 0; i < Size; i++)
		bc_type[i] = 0;

	idx_dpre0 = GetNodeValueIndex("PRESSURE_RATE1");
	idx_dpre = idx_dpre0 + 1;
	idx_pre = GetNodeValueIndex("PRESSURE1");

	nv = 2 * problem_dimension_dm + 1;
	if (m_msh->GetCoordinateFlag() / 10 == 2)  // 2D
	{
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
		idx_disp[0] = GetNodeValueIndex("DISPLACEMENT_X1");
		idx_disp[1] = GetNodeValueIndex("DISPLACEMENT_Y1");
		idx_vel[0] = GetNodeValueIndex("VELOCITY_DM_X");
		idx_vel[1] = GetNodeValueIndex("VELOCITY_DM_Y");
		idx_acc0[0] = GetNodeValueIndex("ACCELERATION_X1");
		idx_acc0[1] = GetNodeValueIndex("ACCELERATION_Y1");
		idx_acc[0] = idx_acc0[0] + 1;
		idx_acc[1] = idx_acc0[1] + 1;
	}
	else if (m_msh->GetCoordinateFlag() / 10 == 3)  // 3D
	{
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
		idx_disp[0] = GetNodeValueIndex("DISPLACEMENT_X1");
		idx_disp[1] = GetNodeValueIndex("DISPLACEMENT_Y1");
		idx_disp[2] = GetNodeValueIndex("DISPLACEMENT_Z1");
		idx_vel[0] = GetNodeValueIndex("VELOCITY_DM_X");
		idx_vel[1] = GetNodeValueIndex("VELOCITY_DM_Y");
		idx_vel[2] = GetNodeValueIndex("VELOCITY_DM_Z");
		idx_acc0[0] = GetNodeValueIndex("ACCELERATION_X1");
		idx_acc0[1] = GetNodeValueIndex("ACCELERATION_Y1");
		idx_acc0[2] = GetNodeValueIndex("ACCELERATION_Z1");
		for (k = 0; k < 3; k++)
			idx_acc[k] = idx_acc0[k] + 1;
	}

	//
	for (size_t i = 0; i < bc_node_value.size(); i++)
	{
		CBoundaryConditionNode* m_bc_node = bc_node_value[i];
		m_bc = bc_node[i];
		for (j = 0; j < nv; j++)
			if (convertPrimaryVariableToString(
			        m_bc->getProcessPrimaryVariable())
			        .compare(function_name[j]) == 0)
				break;
		if (j == nv)
		{
			cout << "No such primary variable found in "
			        "CalcBC_or_SecondaryVariable_Dynamics." << endl;
			abort();
		}
		bc_msh_node = m_bc_node->geo_node_number;
		if (!m_msh->nod_vector[bc_msh_node]->GetMark()) continue;
		if (bc_msh_node >= 0)
		{
			curve = m_bc_node->CurveIndex;
			if (curve > 0)
			{
				time_fac =
				    GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
				if (!valid) continue;
			}
			else
				time_fac = 1.0;
			bc_value = time_fac * m_bc_node->node_value;
			bc_eqs_index = m_msh->nod_vector[bc_msh_node]->GetEquationIndex();

			if (BC)
			{
				if (j < problem_dimension_dm)  // da
				{
					bc_eqs_index += Shift[j];
// da = v = 0.0;
#if defined(USE_PETSC)  // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)  // WW
					eqs_new->SetKnownX_i(bc_eqs_index, 0.);
#else
					MXRandbed(bc_eqs_index, 0.0, eqs->b);
#endif
				}
				else if (j == nv - 1)  // P
				{
					bc_eqs_index += Shift[problem_dimension_dm];
// da = v = 0.0;
#if defined(USE_PETSC)  // || defined (other parallel solver lib). 04.2012 WW
// TODO
#elif defined(NEW_EQS)  // WW
					eqs_new->SetKnownX_i(bc_eqs_index, 0.);
#else
					MXRandbed(bc_eqs_index, 0.0, eqs->b);
#endif
				}
			}
			else
			{
				// Bit operator
				if (!(bc_type[bc_eqs_index] & (int)MathLib::fastpow(2, j)))
					bc_type[bc_eqs_index] += (int)MathLib::fastpow(2, j);
				if (j < problem_dimension_dm)  // Disp
				{
					SetNodeValue(bc_eqs_index, idx_disp[j], bc_value);
					SetNodeValue(bc_eqs_index, idx_vel[j], 0.0);
					SetNodeValue(bc_eqs_index, idx_acc[j], 0.0);
				}
				// Vel
				else if (j >= problem_dimension_dm && j < nv - 1)
				{
					v = GetNodeValue(bc_eqs_index, idx_disp[j]);
					v += bc_value * dt +
					     0.5 * dt * dt *
					         (ARRAY[bc_eqs_index + Shift[j]] +
					          m_num->GetDynamicDamping_beta2() *
					              GetNodeValue(bc_eqs_index, idx_acc0[j]));
					SetNodeValue(bc_eqs_index, idx_disp[j], v);
					SetNodeValue(bc_eqs_index, idx_vel[j], bc_value);
				}
				else if (j == nv - 1)  // Vel
				{                      // p
					SetNodeValue(bc_eqs_index, idx_pre, bc_value);
					SetNodeValue(bc_eqs_index, idx_dpre, 0.0);
				}
			}
		}
	}
	if (BC) return BC;

	// BC
	for (i = 0; i < m_msh->GetNodesNumber(true); i++)
		for (k = 0; k < problem_dimension_dm; k++)
		{
			// If boundary
			if (bc_type[i] & (int)MathLib::fastpow(2, k)) continue;  // u
			//
			v = GetNodeValue(i, idx_disp[k]);
			v += GetNodeValue(i, idx_vel[k]) * dt +
			     0.5 * dt * dt * (ARRAY[i + Shift[k]] +
			                      m_num->GetDynamicDamping_beta2() *
			                          GetNodeValue(i, idx_acc0[k]));
			SetNodeValue(i, idx_disp[k], v);
			if (bc_type[i] & (int)MathLib::fastpow(2, k + problem_dimension_dm))
				continue;
			// v
			v = GetNodeValue(i, idx_vel[k]);
			v += dt * ARRAY[i + Shift[k]] +
			     m_num->GetDynamicDamping_beta1() * dt *
			         GetNodeValue(i, idx_acc0[k]);
			SetNodeValue(i, idx_vel[k], v);
		}

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		if (bc_type[i] & (int)MathLib::fastpow(2, (nv - 1))) continue;
		v = GetNodeValue(i, idx_pre);
		v += ARRAY[i + Shift[problem_dimension_dm]] * dt +
		     m_num->GetDynamicDamping_beta1() * dt * GetNodeValue(i, idx_dpre0);
		SetNodeValue(i, idx_pre, v);
	}

	return BC;
}
