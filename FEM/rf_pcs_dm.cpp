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
#include "MemWatch.h"
#include "display.h"
#include "StringTools.h"

#include "Curve.h"
#include "mathlib.h"
#if defined(NEW_EQS)
#include "equation_class.h"
#endif
#ifdef USE_PETSC
#include "PETSC/PETScLinearSolver.h"
#endif

#include "geo_sfc.h"

#include "msh_elem.h"

#include "ElementValueDM.h"
#include "FEMEnums.h"
#include "fem_ele_std.h"
#include "fem_ele_vec.h"
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_node.h"
#include "rf_msp_new.h"
#include "rf_out_new.h"
#include "rf_pcs.h"
#include "rf_st_new.h"
#include "rf_tim_new.h"
#include "tools.h"


int problem_dimension_dm = 0;

using namespace std;
using namespace FiniteElement;
using namespace SolidProp;
using namespace Math_Group;

CRFProcessDeformation::CRFProcessDeformation()
{
}

CRFProcessDeformation::~CRFProcessDeformation()
{
	// Release memory for element variables
	// alle stationaeren Matrizen etc. berechnen
	// Write Gauss stress // TEST for excavation analysis
	//   if(reload==1)
	// if(reload == 1 || reload == 3)
	if (idata_type == write_all_binary || idata_type == read_write)
	{
		WriteGaussPointStress();
		if (type == 41)  // mono-deformation-liquid
			WriteSolution();
	}

	delete fem_dm;

	for (auto p : ele_value_dm)
		delete p;
	ele_value_dm.clear();
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
	if (msp_vector.empty())
	{
		std::cout << "***ERROR: MSP data not found!\n";
		return;
	}

	//-- NW 25.10.2011
	// this section has to be executed at latest before calling InitGauss()
	// Control for reading and writing solution
	if (reload == 1) idata_type = write_all_binary;
	if (reload == 2) idata_type = read_all_binary;
	if (reload == 3) idata_type = read_write;
	//--
	if ((reload == 2 || reload == 3) && calcDiffFromStress0)
		_isInitialStressNonZero = true;

	//
	const int Axisymm = (m_msh->isAxisymmetry() ? -1 : 1);
	fem_dm = new CFiniteElementVec(this, Axisymm * m_msh->GetCoordinateFlag());
	fem_dm->SetGaussPointNumber(m_num->ele_gauss_points);
	if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
		fem = new CFiniteElementStd(this, Axisymm * m_msh->GetCoordinateFlag());

	pcs_number_deformation = pcs_number;

	// Initialize material transform tensor for tansverse isotropic elasticity
	for (auto msp : msp_vector)
	{
		msp->CalculateTransformMatrixFromNormalVector(problem_dimension_dm);
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
		std::cout << "->Initial stress is given.\n";
		CRFProcess* h_pcs = PCSGet("LIQUID_FLOW");
		if (h_pcs)
		{
			ScreenMessage("->Found LIQUID_FLOW. Store initial liquid pressure.\n");
			int n_nodes = m_msh->GetNodesNumber(false);
			p0.resize(n_nodes);
			const int id_p0 = h_pcs->GetNodeValueIndex("PRESSURE1");
			for (int i = 0; i < n_nodes; i++)
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
	size_t bufferSize(0);
	bool HM_Stagered = false;
	if (getProcessType() == FiniteElement::DEFORMATION)
	{
		bufferSize = GetPrimaryVNumber() * m_msh->GetNodesNumber(true);
		if (H_Process)
			HM_Stagered = true;
	}
	else if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
	{
		bufferSize = (GetPrimaryVNumber() - 1) * m_msh->GetNodesNumber(true) + m_msh->GetNodesNumber(false);
	}

	// Allocate memory for  temporal array
	tempArray.resize(bufferSize);

	// Allocate memory for element variables
	ele_value_dm.reserve(m_msh->ele_vector.size());
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		MeshLib::CElem* elem = m_msh->ele_vector[i];
		ElementValue_DM* ele_val = new ElementValue_DM(elem, m_num->ele_gauss_points, HM_Stagered);
		ele_value_dm.push_back(ele_val);
	}
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
	int Idx_Strain[9] = {-1};

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

	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
		for (int j = 0; j < NS + 1; j++)
			SetNodeValue(i, Idx_Strain[j], 0.0);

	std::vector<CInitialCondition*> stress_ic(6, nullptr);
	for (auto m_ic : ic_vector)
	{
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

	int n_given_ic = 0;
	for (int j = 0; j < NS; j++)
		if (stress_ic[j])
			n_given_ic++;
	if (n_given_ic > 0)
		reload = -1000;

	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		MeshLib::CElem* elem = m_msh->ele_vector[i];
		if (!elem->GetMark())
			continue;

		int MatGroup = elem->GetPatchIndex();
		elem->SetOrder(true);
		fem_dm->ConfigElement(elem);
		ElementValue_DM* eleV_DM = ele_value_dm[i];
		*(eleV_DM->Stress0) = 0.0;
		*(eleV_DM->Stress) = 0.0;

		if (n_given_ic == 0)
			continue;

		int NGS = fem_dm->GetNumGaussPoints();
		for (int gp = 0; gp < NGS; gp++)
		{
			int gp_r = 0, gp_s = 0, gp_t = 0;
			fem_dm->GetGaussData(gp, gp_r, gp_s, gp_t);
			fem_dm->ComputeShapefct(2);
			double xyz[3];
			fem_dm->RealCoordinates(xyz);
			for (int j = 0; j < NS; j++)
			{
				CInitialCondition* m_ic = stress_ic[j];
				if (!m_ic) continue;
				int n_dom = m_ic->GetNumDom();
				for (int k = 0; k < n_dom; k++)
				{
					if (MatGroup != m_ic->GetDomain(k)) continue;
					(*eleV_DM->Stress)(j, gp) =
						m_ic->getLinearFunction()->getValue(m_ic->GetDomain(k), xyz[0], xyz[1], xyz[2]);
					(*eleV_DM->Stress0)(j, gp) = (*eleV_DM->Stress)(j, gp);
				}
			}
		}
		if (eleV_DM->Stress_j)
			(*eleV_DM->Stress_j) = (*eleV_DM->Stress);
		elem->SetOrder(false);
	}

	// Reload the stress results of the previous simulation
	// if(reload >= 2)
	if (idata_type == read_all_binary || idata_type == read_write)
	{
		ReadGaussPointStress();
		if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
			ReadSolution();
	}

	Extropolation_GaussValue();
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
		int local_node_counter = 0;
		for (size_t j = 0; j < m_msh->GetNodesNumber(true); j++)
		{
			if (!m_msh->isNodeLocal(j))
				continue;
			int ish = pcs_number_of_primary_nvals * local_node_counter + i;
			ix[ish] = pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex_Q[j] + i;
			val[ish] = GetNodeValue(j, nidx0 + 1);
			local_node_counter++;
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
	const bool isLinearProblem = FiniteElement::isNewtonKind(m_num->nls_method);

	// setup mesh
	m_msh->SwitchOnQuadraticNodes(true);
	if (hasAnyProcessDeactivatedSubdomains || Deactivated_SubDomain.size() > 0)
		CheckMarkedElement();

#ifdef NEW_EQS
	eqs_new->ConfigNumerics(m_num->ls_precond, m_num->ls_method, m_num->ls_max_iterations, m_num->ls_error_tolerance, m_num->ls_storage_method, m_num->ls_extra_arg);
#endif

	if (this->first_coupling_iteration)
		StoreLastSolution();  // u_n-->temp

	//  Reset stress for each coupling step when partitioned scheme is applied
	//  to HM
	if (H_Process && (type / 10 != 4))
		ResetCouplingStep();

	//
	// Compute the maximum ratio of load increment and
	//   predict the number of load steps
	// ---------------------------------------------------------------
	// Compute the ratio of the current load to initial yield load
	// ---------------------------------------------------------------
	int number_of_load_steps = 1;
	if (type / 10 == 4)  // For monolithic scheme
		number_of_load_steps = 1;
	double damping = 1.0;
	const double Tolerance_global_Newton = m_num->nls_error_tolerance[0];

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
			ScreenMessage("Starting loading step %d/%d.\n", l,
			              number_of_load_steps);
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

#ifndef WIN32
			if (ite_steps == 1)
			{
				BaseLib::MemWatch mem_watch;
				ScreenMessage("\tcurrent mem: %d MB\n", mem_watch.getVirtMemUsage() / (1024 * 1024));
			}
#endif
			// Assemble and solve system equation
			ScreenMessage("Assembling equation system...\n");
			GlobalAssembly();

#ifndef WIN32
			if (ite_steps == 1)
			{
				BaseLib::MemWatch mem_watch;
				ScreenMessage("\tcurrent mem: %d MB\n", mem_watch.getVirtMemUsage() / (1024 * 1024));
			}
#endif

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
			eqs_new->Solver(compress_eqs);
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

	// Recovery the old solution.  Temp --> u_n	for flow proccess
	RecoverSolution();

	// get coupling error
	const double norm_u_k1 = getNormOfDisplacements();  // InitialNormDU_coupling
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
			SetNodeValue(j, p_var_index[i], tempArray[shift + j]);
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
#if defined(NEW_EQS)
	eqs_x = eqs_new->getX();
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
#if defined(USE_PETSC)
	double* eqs_x = eqs_new->GetGlobalSolution();
#elif defined(NEW_EQS)
	double* eqs_x = eqs_new->getX();
#endif

	//
	long shift = 0;
	for (int i = 0; i < problem_dimension_dm; i++)
	{
		long const number_of_nodes = num_nodes_p_var[i];
		//
		int const ColIndex = p_var_index[i] - 1;
		//  Update Newton step: w = w+dw
		if (u_type == 0)
		{
			for (long j = 0; j < number_of_nodes; j++)
			{
#ifdef USE_PETSC
				long k = m_msh->Eqs2Global_NodeIndex_Q[j] * pcs_number_of_primary_nvals + i;
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + eqs_x[k] * damp);
#else
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + eqs_x[j + shift] * damp);
#endif
			}
			shift += number_of_nodes;
		}
		else
		{
			for (long j = 0; j < number_of_nodes; j++)
			{
				SetNodeValue(j, ColIndex + 1, GetNodeValue(j, ColIndex + 1) + GetNodeValue(j, ColIndex));
			}
		}
	}

	// if(type == 42&&m_num->nls_method>0)         //H2M, Newton-Raphson.
	if (type / 10 == 4)  // H2M, HM
	{
		// $p_{n+1}=p_{n+1}+\Delta p$ is already performed when type = 0
		if (u_type == 1) return;

		for (int i = problem_dimension_dm; i < pcs_number_of_primary_nvals; i++)
		{
			long const number_of_nodes = num_nodes_p_var[i];
			//
			int const ColIndex = p_var_index[i];

			for (long j = 0; j < number_of_nodes; j++)
				SetNodeValue(j, ColIndex, GetNodeValue(j, ColIndex) + eqs_x[j + shift] * damp);

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
			tempArray[shift + j] = GetNodeValue(j, p_var_index[i] - ty);
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
	if (getProcessType() == FiniteElement::DEFORMATION_FLOW && ty > 0)
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
				SetNodeValue(j, idx, tempArray[shift + j]);
				if (ty == 1) tempArray[shift + j] = tem;
			}
			else if (ty == 2)
			{
				tem = tempArray[shift + j];
				tempArray[shift + j] = GetNodeValue(j, idx);
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
			fem_dm->ExtropolateGaussStress();
			// TEST        if(!update)
			//           (*eval_DM->Stress) -= (*eval_DM->Stress0);
		}
	}
}


/**************************************************************************
   FEMLib-Method:
   Task: Assemble local matrices and RHS for each element
   Programing:
   02/2005 WW
**************************************************************************/
void CRFProcessDeformation::GlobalAssembly()
{
	GlobalAssembly_DM();

	if (type / 10 == 4)
	{  // p-u monolithic scheme

		// if(!fem_dm->dynamic)   ///
		//  RecoverSolution(1);  // p_i-->p_0
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
#if defined(USE_PETSC)
		eqs_new->EQSV_Viewer("eqs" + number2str(aktueller_zeitschritt) + "a");
#endif
#endif
	}

	if (write_leqs)
	{
		std::string fname = FileName + "_" +
							convertProcessTypeToString(this->getProcessType()) +
							number2str(aktueller_zeitschritt) + "_" +
							number2str(ite_steps) + "_leqs_assembly.txt";
#if defined(NEW_EQS)
		std::ofstream Dum(fname.c_str(), ios::out);
		eqs_new->Write(Dum);
		Dum.close();
#elif defined(USE_PETSC)
		eqs_new->EQSV_Viewer(fname);
#endif
	}


	//   RecoverSolution(2);  // p_i-->p_0


#if 0
	{
		ofstream Dum(std::string("eqs_after_assembly.txt").c_str(), ios::out); // WW
		this->eqs_new->Write(Dum);
		Dum.close();
	}
#endif
	// Apply Neumann BC
	ScreenMessage("-> impose Neumann BC and source/sink terms\n");
	IncorporateSourceTerms();

#if defined(USE_PETSC)
	ScreenMessage2d("assemble PETSc matrix and vectors...\n");
	eqs_new->AssembleUnkowns_PETSc();
	eqs_new->AssembleRHS_PETSc();
	eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);
//		eqs_new->EQSV_Viewer("eqs_after_assembl");
#endif

	// Apply Dirchlete bounday condition
	ScreenMessage("-> impose Dirichlet BC\n");
	IncorporateBoundaryConditions();

	if (write_leqs)
	{
		std::string fname = FileName + "_" +
							convertProcessTypeToString(this->getProcessType()) +
							number2str(aktueller_zeitschritt) + "_" +
							number2str(ite_steps) + "_leqs_assembly_BCST.txt";
#if defined(NEW_EQS)
		std::ofstream Dum(fname.c_str(), ios::out);
		eqs_new->Write(Dum);
		Dum.close();
#elif defined(USE_PETSC)
			eqs_new->EQSV_Viewer(fname);
#endif
	}

#define atest_bin_dump
#ifdef test_bin_dump  // WW
	string fname = FileName + ".eqiation_binary.bin";

	ofstream Dum1(fname.data(), ios::out | ios::binary | ios::trunc);
	if (Dum1.good()) eqs_new->Write_BIN(Dum1);
	Dum1.close();
#endif
}

/*!  \brief Assembe matrix and vectors
      for deformation process

      24.11.2010. WW
 */
void CRFProcessDeformation::GlobalAssembly_DM()
{
	const size_t dn = m_msh->ele_vector.size() / 10;
	const bool print_progress = (dn >= 100);
	if (print_progress)
		ScreenMessage("start local assembly for %d elements...\n",
		              m_msh->ele_vector.size());

	long i;
	MeshLib::CElem* elem = NULL;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		if (print_progress && (i + 1) % dn == 0) ScreenMessage("* ");
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark())  // Marked for use
			continue;

		elem->SetOrder(true);
		fem_dm->ConfigElement(elem);
		fem_dm->LocalAssembly(0);
	}
	if (print_progress)
		ScreenMessage("done\n");
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
