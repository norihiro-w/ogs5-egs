/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_pcs.h"

#ifndef __APPLE__
#include <malloc.h>
#endif

#include <cfloat>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

#if defined(LIS)
#include <lis.h>
#endif

#include "display.h"
#include "FileToolsRF.h"
#include "makros.h"
#include "memory.h"
#include "MemWatch.h"
#include "StringTools.h"


#include "Curve.h"
#if defined(NEW_EQS)
#include "equation_class.h"
#endif
#include "mathlib.h"
#include "InterpolationAlgorithms/InverseDistanceInterpolation.h"
#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#include "PointWithID.h"

#include "msh_faces.h"

#include "DistributionTools.h"
#include "ElementMatrix.h"
#include "eos.h"
#include "fct_mpi.h"
#include "FEMEnums.h"
#include "fem_ele_std.h"
#include "ElementValue.h"
#include "files0.h"
#include "msh_tools.h"
#include "Output.h"
#include "problem.h"
#include "rfmat_cp.h"
#include "rf_ic_new.h"
#include "rf_fct.h"
#include "rf_msp_new.h"
#include "rf_node.h"
#include "rf_pcs_dm.h"
#include "rf_pcs_TH.h"
#ifdef GEM_REACT
// GEMS chemical solver
#include "rf_REACT_GEM.h"
REACT_GEM* m_vec_GEM;
#endif
#ifdef BRNS
// BRNS chemical solver
#include "rf_REACT_BRNS.h"
REACT_BRNS* m_vec_BRNS;
#endif
#include "rf_st_new.h"
#include "tools.h"


using namespace std;
using namespace MeshLib;
using namespace Math_Group;

extern VoidFuncVoid LOPCalcSecondaryVariables_USER;

// Globals, to be checked
VoidXFuncVoidX PCSDestroyELEMatrices[PCS_NUMBER_MAX];
int pcs_no_components = 0;
int pcs_deformation = 0;
int dm_number_of_primary_nvals = 2;
int problem_2d_plane_dm;
int size_eval = 0;
bool hasAnyProcessDeactivatedSubdomains = false;

//--------------------------------------------------------
// Coupling Flag. WW
bool T_Process = false;               // Heat
bool H_Process = false;               // Fluid
bool H2_Process = false;              // Multi-phase
bool H3_Process = false;              // 3-phase
bool M_Process = false;               // Mechanical
bool RD_Process = false;              // Richards
bool MH_Process = false;              // MH monolithic scheme
bool MASS_TRANSPORT_Process = false;  // Mass transport
bool FLUID_MOMENTUM_Process = false;  // Momentum
bool RANDOM_WALK_Process = false;     // RWPT
bool pcs_created = false;
//
int pcs_number_deformation = -1;  // JT2012
int pcs_number_flow = -1;         // JT2012
int pcs_number_heat = -1;         // JT2012
vector<int> pcs_number_mass;      // JT2012

namespace process
{
class CRFProcessDeformation;
}
using process::CRFProcessDeformation;
using MeshLib::CNode;
using MeshLib::CElem;
using FiniteElement::ElementValue;
using Math_Group::vec;

#define noCHECK_EQS
#define noCHECK_ST_GROUP
#define noCHECK_BC_GROUP

extern size_t max_dim;  // OK411 todo

//////////////////////////////////////////////////////////////////////////
// PCS vector
//////////////////////////////////////////////////////////////////////////
// It is better to have space between data type and data name. WW
vector<CRFProcess*> pcs_vector;
vector<double*> ele_val_vector;  // PCH
// vector<string> ele_val_name_vector; // PCH
template <class T>
T* resize(T* array, size_t old_size, size_t new_size);
//////////////////////////////////////////////////////////////////////////
// Construction / destruction
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
   02/2005 WW Local elment assembly (all protected members)
   last modified:
**************************************************************************/
CRFProcess::CRFProcess(void)
    : _problem(NULL),
      p_var_index(NULL),
      num_nodes_p_var(NULL),
      fem(NULL),
      Memory_Type(0),
      Write_Matrix(false),
      matrix_file(NULL),
      WriteSourceNBC_RHS(0),
      ele_val_name_vector(std::vector<std::string>())
{
	iter_lin = 0;
	iter_lin_max = 0;
	iter_nlin = 0;
	iter_nlin_max = 0;
	iter_inner_cpl = 0;
	iter_outer_cpl = 0;
	first_coupling_iteration = false;
	orig_size = 0;
	ite_steps = 0;
	cpl_num_dof_errors = 0;
	continuum_ic = true;
	ExcavDirection = 0;
	ExcavCurve = 0;
	ExcavBeginCoordinate = .0;
	number_of_nvals = 0;
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;
	pcs_number_of_history_values = 0;
	pcs_type_number = 0;
	pcs_number = 0;
	type = -1;
	srand_seed = 0;
	pcs_unknowns_norm = .0;
	size_unknowns = 0;
	temporary_num_dof_errors = 0;
	Tim = NULL;
	saturation_switch = false;

	TempArry = NULL;
	// SB:GS4  pcs_component_number=0; //SB: counter for transport components
	pcs_component_number = pcs_no_components - 1;
	//----------------------------------------------------------------------
	// NUM
	pcs_num_name[0] = NULL;
	pcs_num_name[1] = NULL;
	pcs_sol_name = NULL;
	m_num = NULL;
	cpl_type_name = "PARTITIONED";  // OK
	num_type_name = "FEM";          // OK
	rwpt_app = 0;  // PCH Application types for RWPT such as Cell Dynamics,
	               // Crypto, etc.
	//
	for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
		pcs_number_mass.push_back(
		    -1);        // JT2012 (allow DOF_NUMBER_MAX potential components)
                        //
#if defined(USE_PETSC)  // || defined(using other parallel scheme)//03.3012. WW
	eqs_new = NULL;
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
	MPI_Comm_size(PETSC_COMM_WORLD, &mysize);

#elif NEW_EQS  // WW 07.11.2008
	eqs_new = NULL;
	configured_in_nonlinearloop = false;
#else
	eqs = NULL;                       // WW
#endif
	dof = 1;  // WW
	//
	// ITERATIONS AND COUPLING
	num_notsatisfied = 0;
	num_diverged = 0;
	cpl_max_relative_error = 0.0;
	//
	// PCS designations
	isPCSDeformation = false;
	isPCSFlow = false;
	isPCSMultiFlow = false;
	isPCSHeat = false;
	isPCSMass = false;
	//
	//----------------------------------------------------------------------
	// ELE
	pcs_number_of_evals = 0;
	//----------------------------------------------------------------------
	//
	mobile_nodes_flag = -1;
	//----------------------------------------------------------------------
	// USER
	PCSSetIC_USER = NULL;
	//----------------------------------------------------------------------
	// TIM
	tim_type = FiniteElement::TIM_TRANSIENT;
	time_unit_factor = 1.0;
	timebuffer = 1.0e-5;  // WW
	//_pcs_type_name.empty();
	adaption = false;          // HS 03.2008
	cpl_overlord = NULL;
	cpl_underling = NULL;
	pcs_is_cpl_overlord = false;
	pcs_is_cpl_underling = false;
	//----------------------------------------------------------------------
	// CPL
	for (int i = 0; i < 10; i++)
		Shift[i] = 0;
	selected = true;  // OK
	// MSH OK
	m_msh = NULL;
	// Reload solutions
	reload = -1;
	nwrite_restart = 1;  // kg44 write every timestep is default
	non_linear = false;                   // OK/CMCD
	cal_integration_point_value = false;  // WW
	continuum = 0;
	// adaption = false; JOD removed
	compute_domain_face_normal = false;  // WW
	use_velocities_for_transport = false;
	//
	additioanl2ndvar_print = -1;  // WW
	flow_pcs_type = 0;            // CB default: liquid flow, Sat = 1
	this->simulator = "GEOSYS";   // BG, 09/2009
	this->simulator_model_path =
	    "";  // BG, 09/2009, folder with the Eclipse or DuMux files
	this->simulator_path = "";         // BG, 09/2009, Eclipse or Dumux
	this->simulator_well_path = "";    // KB 02/2011, Eclipse
	this->PrecalculatedFiles = false;  // BG, 01/2011, flag for using already
	                                   // calculated files from Eclipse or DuMux
	this->Phase_Transition_Model = 0;  // BG, 11/2010, flag for using CO2 Phase
	                                   // transition (0...not used, 1...used)
	//----------------------------------------------------------------------
	write_boundary_condition = false;  // 15.01.2008. WW
	OutputMassOfGasInModel = false;    // 05/2012     BG
	WriteProcessed_BC = -1;            // 26.08.2011. WW
	accepted = true;                   // 25.08.2008. WW
	accept_steps = 0;                  // 27.08.1008. WW
	reject_steps = 0;                  // 27.08.1008. WW
	ML_Cap = 0;                        // 23.01.2009 PCH
	PartialPS = 0;                     // 16.02 2009 PCH

#if defined(USE_MPI) || defined(USE_PETSC)  // WW
	cpu_time_assembly = 0;
#endif
// New equation and solver WW
#ifdef NEW_EQS
	eqs_new = NULL;
	configured_in_nonlinearloop = false;
#endif
	flag_couple_GEMS = 0;    // 11.2009 HS
	femFCTmode = false;      // NW
	this->Gl_ML = NULL;      // NW
	this->Gl_Vec = NULL;     // NW
	this->Gl_Vec1 = NULL;    // NW
	this->FCT_AFlux = NULL;  // NW
#ifdef USE_PETSC
	this->FCT_K = NULL;
	this->FCT_d = NULL;
#endif
	ExcavMaterialGroup = -1;  // 01.2010 WX
	PCS_ExcavState = -1;      // WX

	isRSM = false;  // WW
	eqs_x = NULL;
	write_leqs = false;  // NW

	pcs_num_dof_errors = 1;
	for (int i = 0; i < DOF_NUMBER_MAX; i++)
		pcs_absolute_error[i] = .0;

	calcDiffFromStress0 = true;
	resetStrain = true;
	scaleUnknowns = false;
	scaleEQS = false;

	e_n = 1.0;
	e_pre = 1.0;
	e_pre2 = 1.0;
}

void CRFProcess::setProblemObjectPointer(Problem* problem)
{
	_problem = problem;
}

Problem* CRFProcess::getProblemObjectPointer() const
{
	return _problem;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
   01/2005 WW Local assemblier as a member
   last modified:
**************************************************************************/
CRFProcess::~CRFProcess(void)
{
	long i;
	//----------------------------------------------------------------------
	// Finite element
	if (fem) delete fem;  // WW
	fem = NULL;
	//----------------------------------------------------------------------
	// ELE: Element matrices
	ElementMatrix* eleMatrix = NULL;
	ElementValue* gp_ele = NULL;
	if (Ele_Matrices.size() > 0)
	{
		for (i = 0; i < (long)Ele_Matrices.size(); i++)
		{
			eleMatrix = Ele_Matrices[i];
			delete eleMatrix;
			eleMatrix = NULL;
		}
		Ele_Matrices.clear();
	}
	//----------------------------------------------------------------------
	// ELE: Element Gauss point values
	if (ele_gp_value.size() > 0)
	{
		for (i = 0; i < (long)ele_gp_value.size(); i++)
		{
			gp_ele = ele_gp_value[i];
			delete gp_ele;
			gp_ele = NULL;
		}
		ele_gp_value.clear();
	}
	//----------------------------------------------------------------------
	// OUT: Matrix output
	if (matrix_file)
	{
		matrix_file->close();
		delete matrix_file;
	}
	//----------------------------------------------------------------------
	// NOD: Release memory of node values
	for (i = 0; i < (int)nod_val_vector.size(); i++)
	{
		delete[] nod_val_vector[i];  // Add []. WW
		nod_val_vector[i] = NULL;
	}
	nod_val_vector.clear();
	//----------------------------------------------------------------------
	// ST:
	CNodeValue* m_nod_val = NULL;

	// Added &&m_nod_val for RSM model. 15.08.2011. WW
	if (!isRSM)
	{
		for (i = 0; i < (int)st_node_value.size(); i++)
		{
			for (int j = 0; j < (int)st_node_value[i].size(); j++)
			{
				m_nod_val = st_node_value[i][j];
				// OK delete st_node_value[i];
				// OK st_node_value[i] = NULL;
				if (m_nod_val->check_me)  // OK
				{
					m_nod_val->check_me = false;
					delete m_nod_val;
					m_nod_val = NULL;
				}
			}
		}
		st_node_value.clear();
	}
	//----------------------------------------------------------------------
	for (i = 0; i < (int)bc_node_value.size(); i++)
	{
		delete bc_node_value[i];
		bc_node_value[i] = NULL;
	}
	bc_node_value.clear();
	//----------------------------------------------------------------------
	//_pcs_type_name.clear();
	//----------------------------------------------------------------------
	// CON
	continuum_vector.clear();
	// Haibing 13112006------------------------------------------------------
	for (i = 0; i < (long)ele_val_vector.size(); i++)
		delete[] ele_val_vector[i];
	ele_val_vector.clear();
	//----------------------------------------------------------------------
	// 11.08.2010. WW
	DeleteArray(num_nodes_p_var);
	// 20.08.2010. WW
	DeleteArray(p_var_index);
	//----------------------------------------------------------------------
	if (this->m_num && this->m_num->fct_method > 0)  // NW
	{
		delete this->Gl_ML;
		delete this->Gl_Vec;
		delete this->Gl_Vec1;
		delete this->FCT_AFlux;
		this->Gl_ML = NULL;
		this->Gl_Vec = NULL;
		this->Gl_Vec1 = NULL;
		this->FCT_AFlux = NULL;
#ifdef USE_PETSC
		delete this->FCT_K;
		delete this->FCT_d;
		this->FCT_K = NULL;
		this->FCT_d = NULL;
#endif
	}
#if defined(USE_PETSC)  // || defined(other parallel libs)//10.3012. WW
	delete eqs_new;
	eqs_new = NULL;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task:    Gauss point values for CFEMSH
   Programing:
   08/2005 WW Implementation
**************************************************************************/
void CRFProcess::AllocateMemGPoint()
{
	const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
		ele_gp_value.push_back(new ElementValue(this, m_msh->ele_vector[i]));
}

/**************************************************************************
   FEMLib-Method:
   Task:    This function is a part of the monolithic scheme
         and it is used to assign pcs name to IC, ST, BC, TIM and OUT. object
   Programing:
   07/2005 WW Implementation
   10/2010 TF cody style improvements
**************************************************************************/
void CRFProcess::SetOBJNames()
{
	// IC
	const size_t ic_vector_size(ic_vector.size());
	for (size_t i = 0; i < ic_vector_size; i++)
		ic_vector[i]->setProcessType(this->getProcessType());

	// BC
	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
	while (p_bc != bc_list.end())
	{
		(*p_bc)->setProcessType(this->getProcessType());
		++p_bc;
	}

	// ST
	const size_t st_vector_size(st_vector.size());
	for (size_t i = 0; i < st_vector_size; i++)
		st_vector[i]->setProcessType(this->getProcessType());

	// TIM
	const size_t time_vector_size(time_vector.size());
	for (size_t i = 0; i < time_vector_size; i++)
		//		Tim = time_vector[i];
		//		Tim->pcs_type_name = _pcs_type_name;
		time_vector[i]->pcs_type_name =
		    convertProcessTypeToString(this->getProcessType());

	// OUT
	// OK4216
	const size_t out_vector_size(out_vector.size());
	for (size_t i = 0; i < out_vector_size; i++)
		//		m_out = out_vector[i];
		//		m_out->_pcs_type_name = _pcs_type_name;
		out_vector[i]->setProcessType(this->getProcessType());
	// m_out->pcs_pv_name = pcs_primary_function_name[0];//CMCD
	// string temp = pcs_primary_function_name[0];
}

/**************************************************************************
   PCSLib-Method:
   10/2002 OK Implementation
   04/2004 WW Modification for 3D problems
   02/2005 WW New fem calculator
   04/2005 OK MSHCreateNOD2ELERelations
   07/2005 WW Geometry element objects
   02/2006 WW Removed memory leaking
   01/2006 YD MMP for each PCS
   04/2006 WW Unique linear solver for all processes if they share the same mesh
   06/2006 WW Rearrange incorporation of BC and ST. Set BC and ST for domain
decomposition
**************************************************************************/
void CRFProcess::Create()
{
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
#endif
	// we need the string representation of process type at some points
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));

	// Element matrix output. WW
	if (Write_Matrix)
	{
		std::cout << "->Write Matrix" << '\n';
#if defined(USE_MPI) || defined(USE_PETSC)
		char stro[32];
		sprintf(stro, "%d", myrank);
		string m_file_name = FileName + "_" + pcs_type_name + (string)stro +
		                     "_element_matrix.txt";
#else
		std::string m_file_name =
		    FileName + "_" + pcs_type_name + "_element_matrix.txt";
#endif
		matrix_file =
		    new std::fstream(m_file_name.c_str(), ios::trunc | ios::out);
		if (!matrix_file->good())
			std::cout << "Warning in GlobalAssembly: Matrix files are not found"
			          << "\n";
	}
	//----------------------------------------------------------------------------
	if (m_msh)  // OK->MB please shift to Config()

		//		if (_pcs_type_name.compare("GROUNDWATER_FLOW") == 0)
		if (this->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
			MeshLib::MSHDefineMobile(this);
	//----------------------------------------------------------------------------
	int DOF = GetPrimaryVNumber();  // OK should be PCS member variable
	//----------------------------------------------------------------------------
	// MMP - create mmp groups for each process //YD
	ScreenMessage("-> Create MMP\n");
	CMediumPropertiesGroup* m_mmp_group = NULL;
	int continua = 1;  // WW
	if (RD_Process) continua = 2;
	m_mmp_group = MMPGetGroup(pcs_type_name);
	for (int i = 0; i < continua; i++)
		if (!m_mmp_group)
		{
			m_mmp_group = new CMediumPropertiesGroup();
			m_mmp_group->pcs_type_name = pcs_type_name;
			m_mmp_group->Set(this);
			mmp_group_list.push_back(m_mmp_group);
		}
	m_mmp_group = NULL;
	//----------------------------------------------------------------------------
	// NUM_NEW
	ScreenMessage("-> Create NUM\n");
	//	if (pcs_type_name.compare("RANDOM_WALK")) { // PCH RWPT does not need
	// this.
	if (this->getProcessType() !=
	    FiniteElement::RANDOM_WALK)  // PCH RWPT does not need this.
	{
		CNumerics* m_num_tmp = NULL;
		size_t no_numerics(num_vector.size());
		for (size_t i = 0; i < no_numerics; i++)
		{
			m_num_tmp = num_vector[i];

			if ((pcs_type_name.compare(m_num_tmp->pcs_type_name) == 0) ||
			    (m_num_tmp->pcs_type_name.compare(
			         pcs_primary_function_name[0]) == 0))
			{
				m_num = m_num_tmp;
				break;
			}
		}
	}
	if (!m_num)
	{
		ScreenMessage("Warning in CRFProcess::Create() - no NUM data\n");
		m_num = new CNumerics(pcs_type_name);  // OK
		                                       //		m_num = m_num_tmp;
	}
	else
	{
		if (m_num->nls_max_iterations > 1)  // WW
			non_linear = true;
	}
	if (m_num->fct_method > 0)  // NW
	{
// Memory_Type = 1;
#ifdef USE_PETSC
		long gl_size = m_msh->getNumNodesGlobal();
		this->FCT_K = new SparseMatrixDOK(gl_size, gl_size);
		this->FCT_d = new SparseMatrixDOK(gl_size, gl_size);
#else
		long gl_size = m_msh->GetNodesNumber(false);
#endif
		this->FCT_AFlux = new SparseMatrixDOK(gl_size, gl_size);
		this->Gl_ML = new Math_Group::Vector(gl_size);
		this->Gl_Vec = new Math_Group::Vector(gl_size);
		this->Gl_Vec1 = new Math_Group::Vector(gl_size);
	}
	//----------------------------------------------------------------------------
	// EQS - create equation system
	// WW CreateEQS();
	ScreenMessage("-> Create EQS\n");
#if defined(NEW_EQS)
	size_t k;
	for (k = 0; k < fem_msh_vector.size(); k++)
		if (m_msh == fem_msh_vector[k]) break;
	// WW 02.2013. Pardiso
	int eqs_num = 3;

	// if(type==4||type==41)
	//   eqs_new = EQS_Vector[2*k+1];
	if (type == 4 || (type / 10 == 4))  // 03.08.2010. WW
		eqs_new = EQS_Vector[eqs_num * k + 1];
	else
	{
// eqs_new = EQS_Vector[2*k];
		if (getProcessType() == FiniteElement::MULTI_PHASE_FLOW ||
		    getProcessType() == FiniteElement::PS_GLOBAL ||
		    getProcessType() == FiniteElement::TH_MONOLITHIC)
		{
			eqs_new = EQS_Vector[eqs_num * k + 2];
		}
		else
		{
			eqs_new = EQS_Vector[eqs_num * k];
		}
	}  // WW 02.2013. Pardiso
#endif  // If NEW_EQS
	// Set solver properties: EQS<->SOL
	// Internen Speicher allokieren
	// Speicher initialisieren

	//----------------------------------------------------------------------------
	// Time unit factor //WW
	ScreenMessage("-> Create TIM\n");
	// CTimeDiscretization* Tim = TIMGet(_pcs_type_name);
	Tim = TIMGet(pcs_type_name);
	if (!Tim)
	{
		// 21.08.2008. WW
		/* JT->WW: It doesn't seem like a good idea to give a non-existent Tim
		the properties of some specified [0] vector.
		           Why not set default values, and then let other "Tim" control
		the stepping?
		           In other words. If HEAT_TRANSPORT doesn't have time control,
		           we cannot assign a time control type for a FLOW process to a
		HEAT process, this could give incorrect results.
		           THE DEFAULTS ARE NOW SET UP SUCH THAT... if "Tim" doesn't
		exist, this process has no influence on the time step.
		Tim = new CTimeDiscretization(*time_vector[0], pcs_type_name);
		*/
		Tim = new CTimeDiscretization();
		Tim->pcs_type_name = pcs_type_name;
		time_vector.push_back(Tim);  // 21.08.2008. WW
	}
	//	if(Tim->time_control_type == TimeControlType::INVALID &&
	// Tim->time_step_vector.size() > 0)
	//		Tim->time_control_name = "STEPS";
	//
	if (Tim->time_unit.find("MINUTE") != std::string::npos)
		time_unit_factor = 60.0;
	else if (Tim->time_unit.find("HOUR") != std::string::npos)
		time_unit_factor = 3600.0;
	else if (Tim->time_unit.find("DAY") != std::string::npos)
		time_unit_factor = 86400.0;
	else if (Tim->time_unit.find("MONTH") != std::string::npos)
		time_unit_factor = 2592000.0;
	else if (Tim->time_unit.find("YEAR") != std::string::npos)
		time_unit_factor = 31536000;

	//
	if (type == 4 || type / 10 == 4)
		m_msh->SwitchOnQuadraticNodes(true);
	else
		m_msh->SwitchOnQuadraticNodes(false);

	// ELE - config and create element values
	ScreenMessage("-> Config ELE values\n");
	AllocateMemGPoint();
#ifndef WIN32
	ScreenMessaged("\tcurrent mem: %d MB\n",
	               mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif

	// ELE - config element matrices
	// NOD - config and create node values
	ScreenMessage("-> Config NOD values\n");
	double* nod_values = NULL;
	double* ele_values = NULL;  // PCH

	number_of_nvals = 2 * DOF + pcs_number_of_secondary_nvals;
	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		// new time
		nod_val_name_vector.push_back(pcs_primary_function_name[i]);
		// old time //need this MB!
		nod_val_name_vector.push_back(pcs_primary_function_name[i]);
	}
	for (int i = 0; i < pcs_number_of_secondary_nvals; i++)
		// new time
		nod_val_name_vector.push_back(pcs_secondary_function_name[i]);
	//
	long m_msh_nod_vector_size = m_msh->GetNodesNumber(true);
	for (long j = 0; j < number_of_nvals;
	     j++)  // Swap number_of_nvals and mesh size. WW 19.12.2012
	{
		nod_values = new double[m_msh_nod_vector_size];
		for (int i = 0; i < m_msh_nod_vector_size; i++)
			nod_values[i] = 0.0;
		nod_val_vector.push_back(nod_values);
	}
	// Create element values - PCH
	int number_of_evals = 2 * pcs_number_of_evals;  // PCH, increase memory
	if (number_of_evals > 0)  // WW added this "if" condition
	{
		for (int i = 0; i < pcs_number_of_evals; i++)
		{
			// new time
			ele_val_name_vector.push_back(pcs_eval_name[i]);
			// old time
			ele_val_name_vector.push_back(pcs_eval_name[i]);
		}
		size_t m_msh_ele_vector_size(m_msh->ele_vector.size());
		if (ele_val_vector.size() == 0)
			for (size_t j = 0; j < m_msh_ele_vector_size; j++)
			{
				ele_values = new double[number_of_evals];
				size_eval += number_of_evals;  // WW
				for (int i = 0; i < number_of_evals; i++)
					ele_values[i] = 0.0;
				ele_val_vector.push_back(ele_values);
			}
		else
			for (size_t j = 0; j < m_msh_ele_vector_size; j++)
			{
				ele_values = ele_val_vector[j];
				ele_values =
				    resize(ele_values, size_eval, size_eval + number_of_evals);
				size_eval += number_of_evals;
				ele_val_vector[j] = ele_values;
			}
	}
	//
	//--- construct IC
	//-----------------------------------------------------------------
	if (reload >= 2 && ((type != 4 && type / 10 != 4) ||
	                    !resetStrain))  // Modified at 03.08.2010. WW
	{
		// PCH
		ScreenMessage("-> Reloading the primary variables... \n");
		ReadSolution();  // WW
	}

	if (reload < 2)  // PCH: If reload is set, no need to have ICs
	{
		// IC
		ScreenMessage("-> Assign IC\n");
		SetIC();
	}
	else
		// Bypassing IC
		ScreenMessage("-> RELOAD is set to be %d. So bypassing IC's\n", reload);

	if (pcs_type_name_vector.size() &&
	    pcs_type_name_vector[0].find("DYNAMIC") != string::npos)
		setIC_danymic_problems();
	//
	if (pcs_type_name_vector.size() &&
	    pcs_type_name_vector[0].find("DYNAMIC") != string::npos)  // WW
	{
		setBC_danymic_problems();
		setST_danymic_problems();
	}
	else
	{
		// BC - create BC groups for each process
		ScreenMessage("-> Create BC\n");
		CBoundaryConditionsGroup* m_bc_group = NULL;

		// 25.08.2011. WW
		if (WriteProcessed_BC == 2)
			Read_Processed_BC();
		else
		{
			for (int i = 0; i < DOF; i++)
			{
				// OKm_bc_group =
				// BCGetGroup(_pcs_type_name,pcs_primary_function_name[i]);
				// OKif(!m_bc_group){
				BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
				m_bc_group = new CBoundaryConditionsGroup();
				// OK
				m_bc_group->setProcessTypeName(pcs_type_name);
				m_bc_group->setProcessPrimaryVariableName(
				    pcs_primary_function_name[i]);  // OK
				m_bc_group->Set(this, Shift[i]);

				bc_group_list.push_back(
				    m_bc_group);  // Useless, to be removed. WW
				m_bc_group = NULL;
				// OK}
			}
#ifndef USE_PETSC
			if (bc_node_value.size() < 1)  // WW
				cout << "Warning: no boundary conditions specified for "
				     << pcs_type_name << endl;
#endif
			if (WriteProcessed_BC == 1) Write_Processed_BC();
		}
#ifndef WIN32
		ScreenMessaged("\tcurrent mem: %d MB\n",
		               mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
		// ST - create ST groups for each process
		ScreenMessage("-> Create ST\n");
		CSourceTermGroup* m_st_group = NULL;

		if (WriteSourceNBC_RHS == 2)  // Read from file
			ReadRHS_of_ST_NeumannBC();
		else  // WW // Calculate directly
		{
			for (int i = 0; i < DOF; i++)
			{
				// OK m_st_group =
				// m_st_group->Get(pcs_primary_function_name[i]);
				m_st_group =
				    STGetGroup(pcs_type_name, pcs_primary_function_name[i]);
				if (!m_st_group)
				{
					m_st_group = new CSourceTermGroup();
					// OK
					m_st_group->pcs_type_name = pcs_type_name;
					// OK
					m_st_group->pcs_pv_name = pcs_primary_function_name[i];
					m_st_group->Set(this, Shift[i]);
					// Useless, to be removed. WW
					st_group_list.push_back(m_st_group);
				}
			}
			if (WriteSourceNBC_RHS == 1)  // WW
				WriteRHS_of_ST_NeumannBC();
		}
		m_st_group = NULL;
	}
	// Write BC/ST nodes for vsualization.WW
	if (write_boundary_condition && WriteSourceNBC_RHS != 2) WriteBC();

	// Keep all local matrices in the memory
	if (type != 55)  // Not for fluid momentum. WW
	{
		if (Memory_Type != 0) AllocateLocalMatrixMemory();
		if (type == 4 || type / 10 == 4)
		{
			// Set initialization function
			//      CRFProcessDeformation *dm_pcs = (CRFProcessDeformation *)
			//      (this);
			CRFProcessDeformation* dm_pcs =
			    static_cast<CRFProcessDeformation*>(this);
			dm_pcs->Initialization();
		}
		else if (this->getProcessType() == FiniteElement::TH_MONOLITHIC)
		{
			static_cast<CRFProcessTH*>(this)->Initialization();
		}
		else  // Initialize FEM calculator
		{
			int Axisymm = 1;                           // ani-axisymmetry
			if (m_msh->isAxisymmetry()) Axisymm = -1;  // Axisymmetry is true
			fem = new CFiniteElementStd(this,
			                            Axisymm * m_msh->GetCoordinateFlag());
			fem->SetGaussPointNumber(m_num->ele_gauss_points);
		}
	}

	// Initialize the system equations
	if (PCSSetIC_USER) PCSSetIC_USER(pcs_type_number);

	if (compute_domain_face_normal)  // WW
		m_msh->FaceNormal();
	/// Variable index for equation. 20.08.2010. WW
	if (p_var_index)
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
			p_var_index[i] =
			    GetNodeValueIndex(pcs_primary_function_name[i]) + 1;

#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
	size_unknowns = m_msh->GetNodesNumber(true) * pcs_number_of_primary_nvals;
#elif defined(NEW_EQS)
	{
		size_unknowns = eqs_new->getA()->Dim();
	}
#endif

#ifndef WIN32
	ScreenMessage("\tcurrent mem: %d MB\n",
	              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   12/2005 WW
   03/2006 WW Write as acsi
   04/2006 WW
   last modified:
**************************************************************************/
void CRFProcess::WriteRHS_of_ST_NeumannBC()
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_ST_RHS.asc";
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << endl;
		abort();
	}

	os << "$PCS_TYPE  " << endl;

	os << pcs_type_name << endl;
	os << "geo_node_number  ";
	os << "msh_node_number  ";
	os << "CurveIndex ";
	os << "node_value ";
	os << endl;
	os.setf(std::ios::scientific, std::ios::floatfield);
	os.precision(14);
	for (size_t i = 0; i < st_node_value.size(); i++)
	{
		os << st_node_value[i].size() << endl;
		for (size_t j = 0; j < st_node_value[i].size(); j++)
			st_node_value[i][j]->Write(os);
	}
	os.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   03/2006 WW
   last modified: 04/2006
**************************************************************************/
void CRFProcess::ReadRHS_of_ST_NeumannBC()
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_ST_RHS.asc";
	std::ifstream is(m_file_name.c_str(), std::ios::in);
	if (!is.good())
	{
		cout << "File " << m_file_name << " is not found" << endl;
		abort();
	}

	std::string s_buffer;
	getline(is, s_buffer);
	getline(is, s_buffer);
	getline(is, s_buffer);
	st_node_value.resize(st_vector.size());
	for (size_t i = 0; i < st_node_value.size(); i++)
	{
		size_t size;
		is >> size >> ws;
		for (size_t j = 0; j < size; j++)
		{
			CNodeValue* cnodev = new CNodeValue();
			cnodev->Read(is);
			st_node_value[i].push_back(cnodev);
		}
	}
	is.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   08/2011 WW
**************************************************************************/
void CRFProcess::Read_Processed_BC()
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_eqs_BC.asc";
	std::ifstream is(m_file_name.c_str(), std::ios::in);
	if (!is.good())
	{
		cout << "File " << m_file_name << " is not found" << endl;
		abort();
	}

	std::string s_buffer;
	getline(is, s_buffer);
	getline(is, s_buffer);
	getline(is, s_buffer);
	size_t size;
	is >> size >> ws;
	bc_node_value.clear();
	for (size_t i = 0; i < size; i++)
	{
		CBoundaryConditionNode* cnodev = new CBoundaryConditionNode();
		cnodev->Read(is);
		bc_node_value.push_back(cnodev);
	}
	is.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the contribution of ST or Neumann BC to RHS to a file after
      integration
   Programing:
   08/2011 WW
**************************************************************************/
void CRFProcess::Write_Processed_BC()
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_eqs_BC.asc";
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << endl;
		abort();
	}

	os << "$PCS_TYPE  " << endl;

	os << pcs_type_name << endl;
	os << "geo_node_number  ";
	os << "msh_node_number  ";
	os << "CurveIndex ";
	os << "node_value ";
	os << endl;
	os.setf(std::ios::scientific, std::ios::floatfield);
	os.precision(14);
	const size_t bc_node_value_size(bc_node_value.size());
	os << bc_node_value_size << endl;
	for (size_t i = 0; i < bc_node_value_size; i++)
		bc_node_value[i]->Write(os);
	os.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: Get a name for the solution file
   Programing:
   last modified:
**************************************************************************/
std::string CRFProcess::GetSolutionFileName(bool write)
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	std::string m_file_name = FileName + "_" + pcs_type_name + "_" +
	                          pcs_primary_function_name[0] + "_primary_value";
	if (write) m_file_name += "_" + number2str(aktueller_zeitschritt);
#if defined(USE_PETSC)  //|| defined(other parallel libs)//03.3012. WW
	m_file_name += "_rank" + number2str(myrank);
#endif
	m_file_name += ".asc";
	return m_file_name;
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the solution
   Programing:
   04/2006 WW
   last modified:
**************************************************************************/
void CRFProcess::WriteSolution()
{
	if (reload == 2 || reload <= 0) return;
	// kg44 write out only between nwrite_restart timesteps
	if ((aktueller_zeitschritt % nwrite_restart) > 0) return;

	std::string m_file_name = GetSolutionFileName(true);
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		ScreenMessage2("Failure to open file: %s\n", m_file_name.c_str());
		abort();
	}

	os.precision(15);  // 15 digits accuracy seems enough? more fields are
	                   // filled up with random numbers!
	os.setf(std::ios_base::scientific, std::ios_base::floatfield);

	int j;
	std::vector<int> idx(2 * pcs_number_of_primary_nvals);
	for (j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
		idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
	}
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		for (j = 0; j < 2 * pcs_number_of_primary_nvals; j++)
			os << GetNodeValue(i, idx[j]) << "  ";
		os << endl;
	}
	os.close();
	ScreenMessage("Write solutions for timestep %d into file %s\n",
	              aktueller_zeitschritt, m_file_name.c_str());
}

/**************************************************************************
   FEMLib-Method:
   Task: Write the solution
   Programing:
   04/2006 WW
   last modified:
**************************************************************************/
void CRFProcess::ReadSolution()
{
	const std::string m_file_name = GetSolutionFileName(false);
	std::ifstream is(m_file_name.c_str(), ios::in);
	if (!is.good())
	{
		ScreenMessage2("Failure to open file: %s\n", m_file_name.c_str());
		abort();
	}

	std::vector<int> idx(2 * pcs_number_of_primary_nvals);
	std::vector<double> val(2 * pcs_number_of_primary_nvals);

	for (int j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
		idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
	}
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		for (int j = 0; j < 2 * pcs_number_of_primary_nvals; j++)
			is >> val[j];
		is >> ws;
		//		for (j = 0; j < 2 * pcs_number_of_primary_nvals; j++ )
		//			SetNodeValue ( i,idx[j], val[j] );
		// previous and current value should be initially same
		for (int j = 0; j < pcs_number_of_primary_nvals; j++)
		{
			SetNodeValue(i, idx[j], val[j + pcs_number_of_primary_nvals]);
			SetNodeValue(i, idx[j + pcs_number_of_primary_nvals],
			             val[j + pcs_number_of_primary_nvals]);
			//			SetNodeValue (i, idx[j*2], val[j*2+1] );
			//			SetNodeValue (i, idx[j*2+1], val[j*2+1] );
		}
	}
	is.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::setIC_danymic_problems()
{
	const char* function_name[7];
	int i, j, nv;
	nv = 0;
	if (max_dim == 1)  // 2D
	{
		nv = 5;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
	}
	else  // 3D
	{
		nv = 7;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
	}

	CInitialCondition* m_ic = NULL;
	long no_ics = (long)ic_vector.size();
	int nidx;
	for (i = 0; i < nv; i++)
	{
		nidx = GetNodeValueIndex(function_name[i]);
		for (j = 0; j < no_ics; j++)
		{
			m_ic = ic_vector[j];
			if (m_ic->getProcessPrimaryVariable() ==
			    FiniteElement::convertPrimaryVariable(function_name[i]))
				m_ic->Set(nidx);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::setST_danymic_problems()
{
	const char* function_name[7];
	size_t nv = 0;
	if (max_dim == 1)  // 2D
	{
		nv = 5;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
	}  // 3D
	else
	{
		nv = 7;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
	}

	// ST - create ST groups for each process
	CSourceTermGroup* m_st_group = NULL;
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	for (size_t i = 0; i < nv; i++)
	{
		m_st_group = STGetGroup(pcs_type_name, function_name[i]);
		if (!m_st_group)
		{
			m_st_group = new CSourceTermGroup();
			m_st_group->pcs_type_name = pcs_type_name;
			m_st_group->pcs_pv_name = function_name[i];
			m_st_group->Set(this, Shift[i], function_name[i]);
			st_group_list.push_back(m_st_group);  // Useless, to be removed. WW
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::setBC_danymic_problems()
{
	const char* function_name[7];
	size_t nv = 0;
	if (max_dim == 1)  // 2D
	{
		nv = 5;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "VELOCITY_DM_X";
		function_name[3] = "VELOCITY_DM_Y";
		function_name[4] = "PRESSURE1";
	}  // 3D
	else
	{
		nv = 7;
		function_name[0] = "DISPLACEMENT_X1";
		function_name[1] = "DISPLACEMENT_Y1";
		function_name[2] = "DISPLACEMENT_Z1";
		function_name[3] = "VELOCITY_DM_X";
		function_name[4] = "VELOCITY_DM_Y";
		function_name[5] = "VELOCITY_DM_Z";
		function_name[6] = "PRESSURE1";
	}

	cout << "->Create BC" << '\n';
	CBoundaryConditionsGroup* m_bc_group = NULL;
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	for (size_t i = 0; i < nv; i++)
	{
		BCGroupDelete(pcs_type_name, function_name[i]);
		m_bc_group = new CBoundaryConditionsGroup();
		// OK
		m_bc_group->setProcessTypeName(pcs_type_name);
		// OK
		m_bc_group->setProcessPrimaryVariableName(function_name[i]);
		m_bc_group->Set(this, Shift[i], function_name[i]);
		bc_group_list.push_back(m_bc_group);  // Useless, to be removed. WW
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2005 WW Set coupling data
   last modified:
**************************************************************************/
void CRFProcess::ConfigureCouplingForLocalAssemblier()
{
	bool Dyn = false;
	if (pcs_type_name_vector.size() &&
	    pcs_type_name_vector[0].find("DYNAMIC") != string::npos)
		Dyn = true;
	if (fem) fem->ConfigureCoupling(this, Shift, Dyn);
}

/**************************************************************************
   FEMLib-Method:
   06/2003 OK Implementation
        WW 2nd version, PCS_Solver
**************************************************************************/
void PCSDestroyAllProcesses(void)
{

#ifdef NEW_EQS
	for (auto p : EQS_Vector)
		delete p;
	EQS_Vector.clear();
	for (auto p : SparseTable_Vector)
		delete p;
	SparseTable_Vector.clear();
#endif

	//----------------------------------------------------------------------
	// PCS
	CRFProcess* m_process = NULL;
	long i;
	int j;
	for (j = 0; j < (int)pcs_vector.size(); j++)
	{
		m_process = pcs_vector[j];
#ifdef PCS_NOD
		for (i = 0; i < NodeListSize(); i++)
		{
			k = GetNode(i);
			k->values[m_process->pcs_number] =
			    (double*)Free(k->values[m_process->pcs_number]);
		}
#endif
		if (m_process->TempArry)
			m_process->TempArry = (double*)Free(m_process->TempArry);
		delete (m_process);
	}
	pcs_vector.clear();
	//----------------------------------------------------------------------
	// MSH
	for (i = 0; i < (long)fem_msh_vector.size(); i++)
	{
		if (fem_msh_vector[i]) delete fem_msh_vector[i];
		fem_msh_vector[i] = NULL;
	}
	fem_msh_vector.clear();

	//----------------------------------------------------------------------
	// ELE
	for (i = 0; i < (long)ele_val_vector.size(); i++)
		delete ele_val_vector[i];
	ele_val_vector.clear();
	//----------------------------------------------------------------------
	// IC ICDelete()
	for (i = 0; i < (long)ic_vector.size(); i++)
		delete ic_vector[i];
	ic_vector.clear();
	//----------------------------------------------------------------------
	MSPDelete();
	BCDelete();
	ICDelete();
	BCGroupDelete();
	STDelete();
	STGroupsDelete();
	GEOLIB_Clear_GeoLib_Data();
	//......................................................................
	TIMDelete();
	OUTDelete();
	NUMDelete();
	MFPDelete();
	MSPDelete();
	MMPDelete();
	MMPGroupDelete();
	MCPDelete();
	//----------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation
   08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
   last modified:
   10/2010 TF changed process type handling from string to enum
**************************************************************************/
bool PCSRead(std::string file_base_name)
{
	//----------------------------------------------------------------------
	char line[MAX_ZEILE];
	int indexCh1a, indexCh2a;
	std::string CommentK("//");
	std::string line_string;
	std::string pcs_file_name;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
	std::ifstream pcs_file(pcs_file_name.data(), ios::in);
	if (!pcs_file.good())
	{
		cout << "Warning: no PCS data *.pcs file is missing" << endl;
		return false;
	}

	// rewind the file
	pcs_file.clear();
	pcs_file.seekg(0, std::ios::beg);
	//========================================================================
	// Keyword loop
	ScreenMessage("PCSRead ... \n");
	while (!pcs_file.eof())
	{
		pcs_file.getline(line, MAX_ZEILE);
		line_string = line;
		line_string = GetUncommentedLine(line_string);
		if (line_string.find("#STOP") != string::npos) break;
		indexCh1a = (int)line_string.find_first_of(CommentK.c_str());
		indexCh2a = (int)line_string.find("#PROCESS");
		//----------------------------------------------------------------------
		// keyword found
		if (indexCh2a > indexCh1a && (indexCh1a == -1))
		{
			CRFProcess* m_pcs = new CRFProcess();
			m_pcs->file_name_base = file_base_name;  // OK
			position = m_pcs->Read(&pcs_file);
			m_pcs->PCSReadConfigurations();  // JT
			m_pcs->pcs_number = pcs_vector.size();

			// RelocateDeformationProcess(m_pcs);
			//			if (m_pcs->_pcs_type_name.find("DEFORMATION") !=
			// string::npos) { // TF
			if (isDeformationProcess(m_pcs->getProcessType()))
			{
				pcs_vector.push_back(m_pcs->CopyPCStoDM_PCS());
				pcs_vector[pcs_vector.size() - 1]->pcs_number =
				    pcs_vector.size();
				delete m_pcs;
			}
			else if (m_pcs->getProcessType() == FiniteElement::TH_MONOLITHIC)
			{
				pcs_vector.push_back(m_pcs->CopyPCStoTH_PCS());
				pcs_vector[pcs_vector.size() - 1]->pcs_number =
				    pcs_vector.size();
				delete m_pcs;
			}
			else
			{
				pcs_vector.push_back(m_pcs);
			}

			pcs_file.seekg(position, std::ios::beg);
		}  // keyword found
	}      // eof

	ScreenMessage("-> done, read %d processes\n", pcs_vector.size());

	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: Copy data to dm_pcs from PCS read function
   Programing:
   06/2007 OK/WW Implementation
   10/2010 TF many improvements
**************************************************************************/
CRFProcess* CRFProcess::CopyPCStoDM_PCS()
{
	// Numerics
	if (num_type_name.compare("STRONG_DISCONTINUITY") == 0)
		enhanced_strain_dm = 1;

	CRFProcessDeformation* dm_pcs(new CRFProcessDeformation());
	dm_pcs->setProcessType(this->getProcessType());
	dm_pcs->pcs_type_name_vector.push_back(pcs_type_name_vector[0].data());
	dm_pcs->Write_Matrix = Write_Matrix;
	dm_pcs->WriteSourceNBC_RHS = WriteSourceNBC_RHS;
	dm_pcs->num_type_name = num_type_name;
	dm_pcs->Memory_Type = Memory_Type;
	dm_pcs->reload = reload;
	dm_pcs->nwrite_restart = nwrite_restart;
	dm_pcs->isPCSDeformation = true;
	dm_pcs->isPCSFlow = this->isPCSFlow;            // JT
	dm_pcs->isPCSMultiFlow = this->isPCSMultiFlow;  // JT
	// WW
	dm_pcs->write_boundary_condition = write_boundary_condition;
	dm_pcs->Deactivated_SubDomain = Deactivated_SubDomain;
	pcs_deformation = 1;
	// WX:01.2011 for coupled excavation
	if (ExcavMaterialGroup >= 0)
	{
		dm_pcs->ExcavMaterialGroup = ExcavMaterialGroup;
		dm_pcs->ExcavDirection = ExcavDirection;
		dm_pcs->ExcavBeginCoordinate = ExcavBeginCoordinate;
		dm_pcs->ExcavCurve = ExcavCurve;
	}
	dm_pcs->write_leqs = write_leqs;
	dm_pcs->calcDiffFromStress0 = calcDiffFromStress0;
	dm_pcs->resetStrain = resetStrain;
	dm_pcs->scaleUnknowns = scaleUnknowns;
	dm_pcs->tim_type = tim_type;
	//
	return dynamic_cast<CRFProcess*>(dm_pcs);
}

CRFProcess* CRFProcess::CopyPCStoTH_PCS()
{
	CRFProcessTH* dm_pcs(new CRFProcessTH());
	dm_pcs->setProcessType(this->getProcessType());
	dm_pcs->pcs_type_name_vector.push_back(pcs_type_name_vector[0].data());
	dm_pcs->Write_Matrix = Write_Matrix;
	dm_pcs->WriteSourceNBC_RHS = WriteSourceNBC_RHS;
	dm_pcs->num_type_name = num_type_name;
	dm_pcs->Memory_Type = Memory_Type;
	dm_pcs->reload = reload;
	dm_pcs->nwrite_restart = nwrite_restart;
	dm_pcs->isPCSDeformation = false;
	dm_pcs->isPCSFlow = this->isPCSFlow;            // JT
	dm_pcs->isPCSMultiFlow = this->isPCSMultiFlow;  // JT
	// WW
	dm_pcs->write_boundary_condition = write_boundary_condition;
	dm_pcs->Deactivated_SubDomain = Deactivated_SubDomain;
	pcs_deformation = 1;
	// WX:01.2011 for coupled excavation
	if (ExcavMaterialGroup >= 0)
	{
		dm_pcs->ExcavMaterialGroup = ExcavMaterialGroup;
		dm_pcs->ExcavDirection = ExcavDirection;
		dm_pcs->ExcavBeginCoordinate = ExcavBeginCoordinate;
		dm_pcs->ExcavCurve = ExcavCurve;
	}
	dm_pcs->write_leqs = write_leqs;
	dm_pcs->scaleUnknowns = scaleUnknowns;
	dm_pcs->vec_scale_dofs = vec_scale_dofs;
	dm_pcs->scaleEQS = scaleEQS;
	dm_pcs->vec_scale_eqs = vec_scale_eqs;
	dm_pcs->tim_type = tim_type;
	//
	return dynamic_cast<CRFProcess*>(dm_pcs);
}

/**************************************************************************
   FEMLib-Method:
   Task: Initial needed configurations following a read of PCS
   Programing:
   03/2012 JT
**************************************************************************/
void CRFProcess::PCSReadConfigurations()
{
	if (pcs_type_name_vector.size() > 1)
	{
		string pname = pcs_type_name_vector[0] + pcs_type_name_vector[1];
		pcs_type_name_vector.pop_back();
		if (pname.find("FLOW") != string::npos &&
		    pname.find("DEFORMATION") != string::npos)
		{
			setProcessType(FiniteElement::DEFORMATION_FLOW);
			MH_Process = true;  // MH monolithic scheme
			if (pname.find("DYNAMIC") != string::npos)
				pcs_type_name_vector[0] = "DYNAMIC";
		}
	}
	else if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
	{
		// NW
		std::cout << "***Error: DEFORMATION_FLOW is vague definition."
		          << "\n";
		exit(0);
	}

	if (isFlowProcess(getProcessType()))
	{
		this->isPCSFlow = true;
		H_Process = true;
	}
	if (isMultiFlowProcess(getProcessType()))
	{
		this->isPCSMultiFlow = true;
		H2_Process = true;
	}
	if (isDeformationProcess(getProcessType()))
	{
		this->isPCSDeformation = true;
		M_Process = true;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: PCS read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Read deformation process
   11/2004 OK file streaming
   12/2005 OK MSH_TYPE
   01/2006 OK GEO_TYPE
**************************************************************************/
std::ios::pos_type CRFProcess::Read(std::ifstream* pcs_file)
{
	char line[MAX_ZEILE];
	string line_string;
	string CommentK("//");
	string hash("#");
	bool new_keyword = false;
	bool new_subkeyword = false;
	ios::pos_type position;
	ios::pos_type position_subkeyword;
	std::stringstream line_stream;
	saturation_switch = false;  // JOD for Richards
	//----------------------------------------------------------------------
	while (!new_keyword)
	{
		position = pcs_file->tellg();
		pcs_file->getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			break;
		}
		//....................................................................
		// WW Comment line
		if (line_string.find_first_of(CommentK.c_str()) != string::npos)
			return position;
		// SB check for comment sign ;
		line_string = GetUncommentedLine(line_string);
		//....................................................................
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
		{
			while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
			{
				position = pcs_file->tellg();
				line_string = GetLineFromFile1(pcs_file);
				if (line_string.find("#") != string::npos) return position;
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				line_stream.str(line_string);
				std::string pcs_type_name;
				line_stream >> pcs_type_name;
				ScreenMessage("-> found %s\n", pcs_type_name.c_str());
				pcs_type_name_vector.push_back(pcs_type_name);
				this->setProcessType(
				    FiniteElement::convertProcessType(pcs_type_name));
				line_stream.clear();

				if (isFlowProcess(this->getProcessType()))
				{
					H_Process = true;
					this->isPCSFlow = true;               // JT2012
					pcs_number_flow = pcs_vector.size();  // JT2012
					if (this->getProcessType() == FiniteElement::PS_GLOBAL ||
					    this->getProcessType() ==
					        FiniteElement::MULTI_PHASE_FLOW ||
					    pcs_type_name.find("H2") != string::npos)
					{
						this->isPCSMultiFlow = true;
					}
				}
				if (isDeformationProcess(this->getProcessType()))
				{
					M_Process = true;
					this->isPCSDeformation = true;  // JT2012
					// JT: "pcs_number_deformation" is set in
					// CRFProcessDeformation::Initialization()
				}
				if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
				{
					H_Process = true;
					MASS_TRANSPORT_Process = true;
					this->isPCSMass = true;  // JT2012
					pcs_number_mass[pcs_no_components] =
					    pcs_vector.size();  // JT2012
					pcs_no_components++;
					this->setProcessPrimaryVariable(
					    FiniteElement::CONCENTRATION);
				}
				if (this->getProcessType() == FiniteElement::HEAT_TRANSPORT)
				{
					T_Process = true;
					this->isPCSHeat = true;               // JT2012
					pcs_number_heat = pcs_vector.size();  // JT2012
				}
				if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM)
				{
					FLUID_MOMENTUM_Process = true;
				}
				if (this->getProcessType() == FiniteElement::RANDOM_WALK)
				{
					RANDOM_WALK_Process = true;
				}
			}
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$NUM_TYPE") != string::npos)
		{
			*pcs_file >> num_type_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$CPL_TYPE") != string::npos)
		{
			*pcs_file >> cpl_type_name;
			if (cpl_type_name.compare("MONOLITHIC") == 0)
			{
				pcs_deformation = 11;
			}
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIM_TYPE") != string::npos)
		{
			std::string tim_type_name;
			*pcs_file >> tim_type_name;
			this->tim_type = FiniteElement::convertTimType(tim_type_name);
			pcs_file->ignore(MAX_ZEILE, '\n');
			ScreenMessage("-> $TIM_TYPE = %s\n", tim_type_name.c_str());
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$APP_TYPE") != string::npos)
		{
			*pcs_file >> rwpt_app;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PRIMARY_VARIABLE") != string::npos)
		{
			*pcs_file >> primary_variable_name;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ELEMENT_MATRIX_OUTPUT") != string::npos)
		{
			*pcs_file >> Write_Matrix;  // WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// WW
		if (line_string.find("$BOUNDARY_CONDITION_OUTPUT") != string::npos)
		{
			write_boundary_condition = true;
			continue;
		}
		//....................................................................
		// BG 05/2012
		if (line_string.find("$OutputMassOfGasInModel") != string::npos)
		{
			OutputMassOfGasInModel = true;
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ST_RHS") != string::npos)
		{
			*pcs_file >> WriteSourceNBC_RHS;  // WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		if (line_string.find("$PROCESSED_BC") != string::npos)  // 25.08.2011.
		                                                        // WW
		{
			*pcs_file >> WriteProcessed_BC;
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}

		//....................................................................
		// subkeyword found
		if (line_string.find("$MEMORY_TYPE") != string::npos)
		{
			*pcs_file >> Memory_Type;  // WW
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$RELOAD") != string::npos)
		{
			*pcs_file >> reload;  // WW
			if (reload == 1 || reload == 3)
				*pcs_file >> nwrite_restart;  // kg44 read number of timesteps
			                                  // between writing restart files
			pcs_file->ignore(MAX_ZEILE, '\n');
			continue;
		}
		// subkeyword found
		if (line_string.find("$DEACTIVATED_SUBDOMAIN") != string::npos)
		{
			int NumDeactivated_SubDomains = 0;
			*pcs_file >> NumDeactivated_SubDomains;
			pcs_file->ignore(MAX_ZEILE, '\n');

			if (NumDeactivated_SubDomains > 0) {
				Deactivated_SubDomain.resize(NumDeactivated_SubDomains);
				std::stringstream ss;
				for (int i = 0; i < NumDeactivated_SubDomains; i++) {
					*pcs_file >> Deactivated_SubDomain[i] >> ws;
					ss << Deactivated_SubDomain[i] << " ";
				}
				ScreenMessage("-> Deactivate subdomain(s) : %s\n", ss.str().c_str());
			}

			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$MSH_TYPE") != string::npos)
		{
			*pcs_file >> msh_type_name >> ws;
			continue;
		}
		//....................................................................
		//		if (line_string.find("$GEO_TYPE") != string::npos) { //OK
		//			*pcs_file >> geo_type >> geo_type_name >> ws;
		//			continue;
		//		}
		//
		//....................................................................
		// subkeyword found
		if (line_string.find("$MEDIUM_TYPE") != string::npos)
		{
			while ((!new_keyword) || (!new_subkeyword) || (!pcs_file->eof()))
			{
				position_subkeyword = pcs_file->tellg();
				*pcs_file >> line_string;
				if (line_string.size() == 0) break;
				if (line_string.find("#") != string::npos)
				{
					new_keyword = true;
					break;
				}
				if (line_string.find("$") != string::npos)
				{
					new_subkeyword = true;
					break;
				}
				if (line_string.find("CONTINUUM") != string::npos)
				{
					*pcs_file >> line_string;
					// WW
					double w_m = strtod(line_string.data(), NULL);
					continuum_vector.push_back(w_m);
					// WW
					continuum_vector.push_back(1.0 - w_m);
					break;  // WW
				}
				pcs_file->ignore(MAX_ZEILE, '\n');
			}
			continue;
		}
		// OK
		if (line_string.find("$SATURATION_SWITCH") != string::npos)
		{
			saturation_switch = true;
			continue;
		}
		// SB4900
		if (line_string.find("$USE_VELOCITIES_FOR_TRANSPORT") != string::npos)
		{
			// Only for fluid momentum process
			if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM)
				use_velocities_for_transport = true;
			continue;
		}
		// Interface to Eclipse and Dumux, BG, 09/2010
		//	if(line_string.find("$SIMULATOR")!=string::npos) { //OK
		if (line_string.compare("$SIMULATOR") ==
		    0)  // BG, 09/2010, coupling to Eclipse and DuMux
		{
			*pcs_file >> this->simulator;
			continue;
		}
		if (line_string.find("$SIMULATOR_PATH") ==
		    0)  // BG, 09/2010, coupling to Eclipse and DuMux
		{
			*pcs_file >> this->simulator_path;
			continue;
		}
		// BG, 09/2010, coupling to Eclipse and DuMux
		if (line_string.find("$SIMULATOR_MODEL_PATH") == 0)
		{
			*pcs_file >> this->simulator_model_path;
			continue;
		}
		// BG, 09/2010, coupling to Eclipse and DuMux
		if (line_string.find("$USE_PRECALCULATED_FILES") == 0)
		{
			this->PrecalculatedFiles = true;
			continue;
		}
		// KB, 02/2011, coupling to Eclipse and DuMux
		if (line_string.find("$SIMULATOR_WELL_PATH") == 0)
		{
			*pcs_file >> this->simulator_well_path;
			continue;
		}
		// BG, NB 11/2010, calculating phase transition for CO2
		if (line_string.find("$PHASE_TRANSITION") == 0)
		{
			string tempstring;
			*pcs_file >> tempstring;
			if (tempstring == "CO2_H2O_NaCl") this->Phase_Transition_Model = 1;
			continue;
		}
		// WX:07.2011
		if (line_string.find("$TIME_CONTROLLED_EXCAVATION") == 0)
		{
			*pcs_file >> ExcavMaterialGroup >> ExcavDirection >>
			    ExcavBeginCoordinate >> ExcavCurve;
			continue;
		}
		if (line_string.find("$LEQS_OUTPUT") == 0)
		{
			write_leqs = true;
			continue;
		}
		if (line_string.find("$CALC_DIFF_FROM_STRESS0") == 0)
		{
			// a flag to calculate deformation in total stress or differential
			// stress from reference state
			// default is true
			int dummy = 0;
			*pcs_file >> dummy;
			calcDiffFromStress0 = (dummy != 0);
			ScreenMessage("-> $CALC_DIFF_FROM_STRESS0 = %d\n", dummy);
			continue;
		}
		if (line_string.find("$RESET_STRAIN") == 0)
		{
			int dummy = 0;
			*pcs_file >> dummy;
			resetStrain = (dummy != 0);
			ScreenMessage("-> $RESET_STRAIN = %d\n", dummy);
			if (!resetStrain)
				ScreenMessage("-> $RESET_STRAIN=0 is not supported yet\n",
				              dummy);
			continue;
		}
		if (line_string.find("$TOTAL_STRESS_COUPLING") == 0)
		{
			int dummy = 0;
			*pcs_file >> dummy;
			use_total_stress_coupling = (dummy != 0);
			if (use_total_stress_coupling)
				ScreenMessage("-> coupling with deformation will be done via total stress instead of strain\n", dummy);
			continue;
		}
		if (line_string.find("$SCALE_DOF") == 0)
		{
			int n = 0;
			*pcs_file >> n;
			vec_scale_dofs.resize(n);
			for (int i = 0; i < n; i++)
				*pcs_file >> vec_scale_dofs[i];
			scaleUnknowns = true;
			std::stringstream ss;
			ss << "-> scale DOFs:";
			for (int i = 0; i < n; i++)
				ss << "[" << i << "] " << vec_scale_dofs[i] << " ";
			ScreenMessage("%s\n", ss.str().c_str());
			continue;
		}
		if (line_string.find("$SCALE_EQS") == 0)
		{
			int n = 0;
			*pcs_file >> n;
			vec_scale_eqs.resize(n);
			for (int i = 0; i < n; i++)
				*pcs_file >> vec_scale_eqs[i];
			scaleEQS = true;
			std::stringstream ss;
			ss << "-> scale EQSs:";
			for (int i = 0; i < n; i++)
				ss << "[" << i << "] " << vec_scale_eqs[i] << " ";
			ScreenMessage("%s\n", ss.str().c_str());
			continue;
		}
		//....................................................................
	}
	//----------------------------------------------------------------------
	return position;
}

/**************************************************************************
   FEMLib-Method:
   01/2004 OK Implementation
   08/2004 WW Read the deformation process
           Check the comment key '//' in .pcs
   06/2009 OK Write only if existing
**************************************************************************/
void PCSWrite(string file_base_name)
{
	if ((int)pcs_vector.size() < 1) return;
	//----------------------------------------------------------------------
	// File handling
	string pcs_file_name = file_base_name + PCS_FILE_EXTENSION;
	fstream pcs_file(pcs_file_name.data(), ios::trunc | ios::out);
	pcs_file.clear();
	//----------------------------------------------------------------------
	// PCS loop
	cout << "PCSWrite" << endl;
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		pcs_file << "#PROCESS" << endl;
		m_pcs->Write(&pcs_file);
	}
	//----------------------------------------------------------------------
	pcs_file << "#STOP" << endl;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2004 OK Implementation
   12/2005 OK MSH_TYPE
   last modified:
**************************************************************************/
void CRFProcess::Write(std::fstream* pcs_file)
{
	*pcs_file << " $PCS_TYPE" << endl;
	*pcs_file << "  " << convertProcessTypeToString(this->getProcessType())
	          << endl;

	*pcs_file << " $NUM_TYPE" << endl;
	*pcs_file << "  " << num_type_name << endl;

	*pcs_file << " $CPL_TYPE" << endl;
	*pcs_file << "  " << cpl_type_name << endl;

	*pcs_file << " $TIM_TYPE" << endl;
	*pcs_file << "  " << FiniteElement::convertTimTypeToString(tim_type)
	          << endl;

	if (msh_type_name.size() > 0)
	{
		*pcs_file << " $MSH_TYPE" << endl;
		*pcs_file << "  " << msh_type_name << endl;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2004 OK Implementation
   last modified:
   10/2010 changed access to process type
**************************************************************************/
CRFProcess* PCSGet(const std::string& pcs_type_name)
{
	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(pcs_type_name));
	for (size_t i = 0; i < pcs_vector.size(); i++)
		//		m_pcs = pcs_vector[i];
		//		if(m_pcs->pcs_type_name.compare(pcs_type_name)==0) { TF
		if (pcs_vector[i]->getProcessType() == pcs_type) return pcs_vector[i];

	return NULL;
}

CRFProcess* PCSGet(FiniteElement::ProcessType pcs_type)
{
	for (size_t i = 0; i < pcs_vector.size(); i++)
		if (pcs_vector[i]->getProcessType() == pcs_type) return pcs_vector[i];

	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2008 TK New Version with Primary Variable Comparision
   last modified:
   10/2010 TF changed access to process type
**************************************************************************/
CRFProcess* PCSGetNew(const string& pcs_type_name,
                      const string& primary_variable_name)
{
	CRFProcess* m_pcs_return = NULL;

	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(pcs_type_name));

	int matches = 0;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		CRFProcess* pcs = pcs_vector[i];
		//		if (pcs->pcs_type_name.compare(pcs_type_name) == 0) { TF
		if (pcs->getProcessType() == pcs_type)
		{
			for (size_t j = 0; j < pcs->GetPrimaryVNumber(); j++)
				if (primary_variable_name.compare(pcs->GetPrimaryVName(j)) == 0)
				{
					m_pcs_return = pcs;
					matches++;
					if (matches > 1) return NULL;
				}
		}
	}
	if (matches == 0)
		return NULL;
	else
		return m_pcs_return;
}

//////////////////////////////////////////////////////////////////////////
// Access
//////////////////////////////////////////////////////////////////////////

// OK->SB please try Get function
CRFProcess* CRFProcess::GetProcessByFunctionName(char* name)
{
	CRFProcess* m_process = NULL;
	/* Tests */
	if (!name) return m_process;
	int i;
	int no_processes = (int)pcs_vector.size();
	for (i = 0; i < no_processes; i++)
	{
		m_process = pcs_vector[i];
		if (strcmp(StrUp(m_process->pcs_primary_function_name[0]),
		           StrUp(name)) == 0)
			break;
	}
	return m_process;
}

// SB: new 3912
CRFProcess* CRFProcess::GetProcessByNumber(int number)
{
	CRFProcess* m_process = NULL;
	/* Tests */
	if (number < 1) return m_process;
	int i;
	int no_processes = (int)pcs_vector.size();
	for (i = 0; i < no_processes; i++)
	{
		m_process = pcs_vector[i];
		if (m_process->pcs_number == number) break;
	}
	return m_process;
}

//////////////////////////////////////////////////////////////////////////
// Configuration
//////////////////////////////////////////////////////////////////////////

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2003 OK Implementation
   08/2003 WW Modified to fit monolithic scheme
   02/2005 OK Unsaturated flow (Richards model)
   02/2005 MB string
   05/2005 WW/DL Dymanic problem
   01/2006 YD Dual Richards
   OKToDo switch to char
   03/2009 PCH PS_GLOBAL
   get rid of type
**************************************************************************/
void CRFProcess::Config(void)
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	// Set mesh pointer to corresponding mesh
	m_msh = MeshLib::FEMGet(pcs_type_name);
	if (!m_msh)
	{
		cout << "Error in CRFProcess::Config - no MSH data" << endl;
		return;
	}
	CheckMarkedElement();

	if (continuum_vector.empty())
		continuum_vector.push_back(1.0);

	if (this->getProcessType() == FiniteElement::LIQUID_FLOW)
	{
		type = 1;
		ConfigLiquidFlow();
	}
	if (this->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
	{
		type = 1;
		ConfigGroundwaterFlow();
	}
	if (this->getProcessType() == FiniteElement::RICHARDS_FLOW)
	{
		if (continuum_vector.size() > 1)
		{
			RD_Process = true;
			type = 22;
		}
		else
			type = 14;
		ConfigUnsaturatedFlow();
	}
	if (this->getProcessType() == FiniteElement::AIR_FLOW)
	{
		type = 5;
		ConfigGasFlow();
	}
	if (this->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
	{
		type = 12;
		ConfigMultiphaseFlow();
	}
	if (this->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	{
		type = 3;
		ConfigHeatTransport();
	}
	if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		type = 2;
		ConfigMassTransport();
	}
	if (isDeformationProcess(getProcessType())) ConfigDeformation();
	if (this->getProcessType() == FiniteElement::FLUID_MOMENTUM)
	{
		type = 55;
		ConfigFluidMomentum();
	}
	if (this->getProcessType() == FiniteElement::RANDOM_WALK)
	{
		type = 55;
		ConfigRandomWalk();
	}
	if (this->getProcessType() ==
		FiniteElement::MULTI_PHASE_FLOW)
	{
		type = 1212;
		ConfigMultiPhaseFlow();
	}
	if (this->getProcessType() == FiniteElement::PS_GLOBAL)
	{
		type = 1313;
		ConfigPS_Global();
	}
	if (this->getProcessType() == FiniteElement::TH_MONOLITHIC)
	{
		type = 1414;
		ConfigTH();
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
   02/2005 MB Pressure version for LIQUID Flow
**************************************************************************/
void CRFProcess::ConfigLiquidFlow()
{
	// pcs_num_name[0] = "PRESSURE0";
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;
	pcs_number_of_evals = 0;
	Def_Variable_LiquidFlow();  // NW

	// Output material parameters
	configMaterialParameters();
}

/**************************************************************************
   FEMLib-Method:
   03/2003 OK Implementation
        WW Splitted for processes
   02/2005 MB head version for GroundwaterFlow
   08/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigGroundwaterFlow()
{
	pcs_num_name[0] = "HEAD";
	// NOD values
	pcs_number_of_primary_nvals = 1;
	pcs_primary_function_name[0] = "HEAD";
	pcs_primary_function_unit[0] = "m";
	// ELE values
	pcs_number_of_evals = 6;
	pcs_eval_name[0] = "VOLUME";
	pcs_eval_unit[0] = "m3";
	pcs_eval_name[1] = "VELOCITY1_X";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Y";
	pcs_eval_unit[2] = "m/s";
	pcs_eval_name[3] = "VELOCITY1_Z";
	pcs_eval_unit[3] = "m/s";
	pcs_eval_name[4] = "PERMEABILITY";
	pcs_eval_unit[4] = "m^2";
	pcs_eval_name[5] = "POROSITY";
	pcs_eval_unit[5] = "-";
	//----------------------------------------------------------------------
	// Secondary variables
	pcs_number_of_secondary_nvals = 5;
	pcs_secondary_function_name[0] = "FLUX";
	pcs_secondary_function_unit[0] = "m3/s";
	pcs_secondary_function_timelevel[0] = 1;
	pcs_secondary_function_name[1] = "WDEPTH";
	pcs_secondary_function_unit[1] = "m";
	pcs_secondary_function_timelevel[1] = 1;
	pcs_secondary_function_name[2] = "COUPLING";  // JOD
	pcs_secondary_function_unit[2] = "m/s";
	pcs_secondary_function_timelevel[2] = 0;
	pcs_secondary_function_name[3] = "COUPLING";  // JOD
	pcs_secondary_function_unit[3] = "m/s";
	pcs_secondary_function_timelevel[3] = 1;
	pcs_secondary_function_name[4] = "STORE";  // JOD  for subtiming 4.7.10
	pcs_secondary_function_unit[4] = "m";
	pcs_secondary_function_timelevel[4] = 1;

	pcs_number_of_secondary_nvals = 5;  // WW
	// WW
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;  // WW
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;  // WW
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;  // WW
	//----------------------------------------------------------------------
	// WW / TF
	// Output material parameters
	configMaterialParameters();

	if (m_msh) DefineMobileNodes(m_msh, this);
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   10/2004 OK Implementation
   last modified:
**************************************************************************/
void CRFProcess::ConfigGasFlow()
{
	//----------------------------------------------------------------------
	// Primary variables - NOD values
	pcs_number_of_primary_nvals = 1;
	pcs_number_of_secondary_nvals = 0;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	//----------------------------------------------------------------------
	// Secondary variables - NOD values
	pcs_number_of_secondary_nvals = 1;
	pcs_secondary_function_name[0] = "NOD_MASS_FLUX";
	pcs_secondary_function_unit[0] = "kg/s";
	//----------------------------------------------------------------------
	// ELE values
	pcs_number_of_evals = 3;
	pcs_eval_name[0] = "VELOCITY1_X";
	pcs_eval_unit[0] = "m/s";
	pcs_eval_name[1] = "VELOCITY1_Y";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Z";
	pcs_eval_unit[2] = "m/s";
	//----------------------------------------------------------------------
	// NUM
	pcs_num_name[0] = "PRESSURE0";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
void CRFProcess::ConfigMultiphaseFlow()
{
	switch (pcs_type_number)
	{
		case 0:
			pcs_num_name[0] = "PRESSURE0";
			break;
		case 1:
			pcs_num_name[0] = "SATURATION0";
			break;
	}
	//----------------------------------------------------------------------
	// NOD Primary variables
	pcs_number_of_primary_nvals = 1;
	switch (pcs_type_number)
	{
		case 0:
			pcs_primary_function_name[0] = "PRESSURE1";
			pcs_primary_function_unit[0] = "Pa";
			break;
		case 1:
			pcs_primary_function_name[0] = "SATURATION2";
			pcs_primary_function_unit[0] = "m3/m3";
			break;
	}
	//----------------------------------------------------------------------
	// NOD Secondary variables
	pcs_number_of_secondary_nvals = 12;  // BG
	switch (pcs_type_number)
	{
		case 0:
			pcs_secondary_function_name[0] = "PRESSURE2";
			pcs_secondary_function_unit[0] = "Pa";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_name[1] = "PRESSURE2";
			pcs_secondary_function_unit[1] = "Pa";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_name[2] = "PRESSURE_CAP";
			pcs_secondary_function_unit[2] = "Pa";
			pcs_secondary_function_timelevel[2] = 0;
			pcs_secondary_function_name[3] = "FLUX";
			pcs_secondary_function_unit[3] = "m3/s";
			pcs_secondary_function_timelevel[3] = 0;
			pcs_secondary_function_name[4] = "DENSITY1";
			pcs_secondary_function_unit[4] = "kg/m3";
			pcs_secondary_function_timelevel[4] = 1;
			pcs_secondary_function_name[5] = "VISCOSITY1";
			pcs_secondary_function_unit[5] = "Pa s";
			pcs_secondary_function_timelevel[5] = 1;
			// BG
			pcs_secondary_function_name[6] = "VELOCITY_X1";
			pcs_secondary_function_unit[6] = "m/s";
			pcs_secondary_function_timelevel[6] = 1;
			pcs_secondary_function_name[7] = "VELOCITY_Y1";
			pcs_secondary_function_unit[7] = "m/s";
			pcs_secondary_function_timelevel[7] = 1;
			pcs_secondary_function_name[8] = "VELOCITY_Z1";
			pcs_secondary_function_unit[8] = "m/s";
			pcs_secondary_function_timelevel[8] = 1;
			pcs_secondary_function_name[9] = "VELOCITY_X2";
			pcs_secondary_function_unit[9] = "m/s";
			pcs_secondary_function_timelevel[9] = 1;
			pcs_secondary_function_name[10] = "VELOCITY_Y2";
			pcs_secondary_function_unit[10] = "m/s";
			pcs_secondary_function_timelevel[10] = 1;
			pcs_secondary_function_name[11] = "VELOCITY_Z2";
			pcs_secondary_function_unit[11] = "m/s";
			pcs_secondary_function_timelevel[11] = 1;
			break;
		case 1:
			pcs_secondary_function_name[0] = "SATURATION1";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_unit[0] = "m3/m3";
			pcs_secondary_function_name[1] = "SATURATION1";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_unit[1] = "m3/m3";
			pcs_secondary_function_name[2] = "PRESSURE_CAP";
			pcs_secondary_function_unit[2] = "Pa";
			pcs_secondary_function_timelevel[2] = 1;
			pcs_secondary_function_name[3] = "FLUX";
			pcs_secondary_function_unit[3] = "m3/s";
			pcs_secondary_function_timelevel[3] = 1;
			pcs_secondary_function_name[4] = "DENSITY2";
			pcs_secondary_function_unit[4] = "kg/m3";
			pcs_secondary_function_timelevel[4] = 1;
			pcs_secondary_function_name[5] = "VISCOSITY2";
			pcs_secondary_function_unit[5] = "Pa s";
			pcs_secondary_function_timelevel[5] = 1;
			// BG
			pcs_secondary_function_name[6] = "VELOCITY_X1";
			pcs_secondary_function_unit[6] = "m/s";
			pcs_secondary_function_timelevel[6] = 1;
			pcs_secondary_function_name[7] = "VELOCITY_Y1";
			pcs_secondary_function_unit[7] = "m/s";
			pcs_secondary_function_timelevel[7] = 1;
			pcs_secondary_function_name[8] = "VELOCITY_Z1";
			pcs_secondary_function_unit[8] = "m/s";
			pcs_secondary_function_timelevel[8] = 1;
			pcs_secondary_function_name[9] = "VELOCITY_X2";
			pcs_secondary_function_unit[9] = "m/s";
			pcs_secondary_function_timelevel[9] = 1;
			pcs_secondary_function_name[10] = "VELOCITY_Y2";
			pcs_secondary_function_unit[10] = "m/s";
			pcs_secondary_function_timelevel[10] = 1;
			pcs_secondary_function_name[11] = "VELOCITY_Z2";
			pcs_secondary_function_unit[11] = "m/s";
			pcs_secondary_function_timelevel[11] = 1;
			break;
	}
	//----------------------------------------------------------------------
	// ELE values
	pcs_number_of_evals = 7;
	switch (pcs_type_number)
	{
		case 0:
			pcs_eval_name[0] = "VELOCITY1_X";
			pcs_eval_unit[0] = "m/s";
			pcs_eval_name[1] = "VELOCITY1_Y";
			pcs_eval_unit[1] = "m/s";
			pcs_eval_name[2] = "VELOCITY1_Z";
			pcs_eval_unit[2] = "m/s";
			pcs_eval_name[3] = "POROSITY1";  // MX 03.2005
			pcs_eval_unit[3] = "-";
			pcs_eval_name[4] = "POROSITY1_IL";  // MX 03.2005
			pcs_eval_unit[4] = "-";
			pcs_eval_name[5] = "PERMEABILITY1";  // MX 03.2005
			pcs_eval_unit[5] = "-";
			pcs_eval_name[6] = "POROSITY1_SW";  // MX 03.2005
			pcs_eval_unit[6] = "-";
			break;
		case 1:
			pcs_eval_name[0] = "VELOCITY2_X";
			pcs_eval_unit[0] = "m/s";
			pcs_eval_name[1] = "VELOCITY2_Y";
			pcs_eval_unit[1] = "m/s";
			pcs_eval_name[2] = "VELOCITY2_Z";
			pcs_eval_unit[2] = "m/s";
			pcs_eval_name[3] = "POROSITY";  // MX 03.2005
			pcs_eval_unit[3] = "-";
			pcs_eval_name[4] = "POROSITY_IL";  // MX 03.2005
			pcs_eval_unit[4] = "-";
			pcs_eval_name[5] = "PERMEABILITY";  // MX 03.2005
			pcs_eval_unit[5] = "-";
			pcs_eval_name[6] = "POROSITY_SW";  // MX 03.2005
			pcs_eval_unit[6] = "-";
			break;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
void CRFProcess::ConfigNonIsothermalFlow()
{
	//----------------------------------------------------------------------
	// Primary variables
	pcs_number_of_primary_nvals = 1;
	switch (pcs_type_number)
	{
		case 0:
			pcs_num_name[0] = "PRESSURE0";
			pcs_primary_function_name[0] = "PRESSURE1";
			pcs_primary_function_unit[0] = "Pa";
			break;
		case 1:
			pcs_num_name[0] = "SATURATION0";
			pcs_primary_function_name[0] = "SATURATION2";
			pcs_primary_function_unit[0] = "m3/m3";
			break;
	}
	//----------------------------------------------------------------------
	// Secondary variables
	pcs_number_of_secondary_nvals = 6;
	switch (pcs_type_number)
	{
		case 0:
			pcs_secondary_function_name[0] = "PRESSURE2";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_unit[0] = "Pa";
			pcs_secondary_function_name[1] = "PRESSURE2";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_unit[1] = "Pa";
			pcs_secondary_function_name[2] = "MASS_FRACTION1";
			pcs_secondary_function_timelevel[2] = 0;
			pcs_secondary_function_unit[2] = "kg/kg";
			pcs_secondary_function_name[3] = "MASS_FRACTION1";
			pcs_secondary_function_timelevel[3] = 1;
			pcs_secondary_function_unit[3] = "kg/kg";
			pcs_secondary_function_name[4] = "PRESSURE_CAP";
			pcs_secondary_function_timelevel[4] = 0;
			pcs_secondary_function_unit[4] = "Pa";
			pcs_secondary_function_name[5] = "DENSITY1";
			pcs_secondary_function_timelevel[5] = 1;
			pcs_secondary_function_unit[5] = "kg/m3";
			break;
		case 1:
			pcs_secondary_function_name[0] = "SATURATION1";
			pcs_secondary_function_timelevel[0] = 0;
			pcs_secondary_function_unit[0] = "m3/m3";
			pcs_secondary_function_name[1] = "SATURATION1";
			pcs_secondary_function_timelevel[1] = 1;
			pcs_secondary_function_unit[1] = "m3/m3";
			pcs_secondary_function_name[2] = "MASS_FRACTION2";
			pcs_secondary_function_timelevel[2] = 0;
			pcs_secondary_function_unit[2] = "kg/kg";
			pcs_secondary_function_name[3] = "MASS_FRACTION2";
			pcs_secondary_function_timelevel[3] = 1;
			pcs_secondary_function_unit[3] = "kg/kg";
			pcs_secondary_function_name[4] = "PRESSURE_CAP";
			pcs_secondary_function_timelevel[4] = 1;
			pcs_secondary_function_unit[4] = "Pa";
			pcs_secondary_function_name[5] = "DENSITY2";
			pcs_secondary_function_timelevel[5] = 1;
			pcs_secondary_function_unit[5] = "kg/m3";
			break;
	}
	// Node
	pcs_number_of_primary_nvals = 1;
	// ELE values
	pcs_number_of_evals = 14;
	pcs_eval_name[0] = "COMP_FLUX";
	pcs_eval_name[1] = "POROSITY";
	pcs_eval_name[2] = "PERMEABILITY";
	pcs_eval_name[3] = "VELOCITY1_X";
	pcs_eval_name[4] = "VELOCITY1_Y";
	pcs_eval_name[5] = "VELOCITY1_Z";
	pcs_eval_name[6] = "VELOCITY2_X";
	pcs_eval_name[7] = "VELOCITY2_Y";
	pcs_eval_name[8] = "VELOCITY2_Z";
	pcs_eval_name[9] = "POROSITY_IL";
	pcs_eval_name[10] = "VoidRatio";
	pcs_eval_name[11] = "PorosityChange";
	pcs_eval_name[12] = "n_sw_Rate";
	pcs_eval_name[13] = "POROSITY_SW";
	pcs_eval_unit[0] = "kg/s";
	pcs_eval_unit[1] = "m3/m3";
	pcs_eval_unit[2] = "m2";
	pcs_eval_unit[3] = "m/s";
	pcs_eval_unit[4] = "m/s";
	pcs_eval_unit[5] = "m/s";
	pcs_eval_unit[6] = "m/s";
	pcs_eval_unit[7] = "m/s";
	pcs_eval_unit[8] = "m/s";
	pcs_eval_unit[9] = "-";
	pcs_eval_unit[10] = "-";
	pcs_eval_unit[11] = "-";
	pcs_eval_unit[12] = "-";
	pcs_eval_unit[13] = "-";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
// void CRFProcess::ConfigNonIsothermalFlowRichards()

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2004 SB Implementation
        WW Splitted for processes
   01/2006 OK Tests
   08/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigMassTransport()
{
	long comp = 1;
	/* count transport processes */
	pcs_component_number++;
	comp = pcs_component_number;
	// 1 NOD values
	// 1.1 primary variables
	pcs_number_of_primary_nvals = 1;
	pcs_primary_function_name[0] = new char[80];
	//  sprintf(pcs_primary_function_name[0], "%s%li","CONCENTRATION",comp);
	//----------------------------------------------------------------------
	// Tests
	// WW int size;
	// WW  size = (int)cp_vec.size();
	// int comb;                                      //OK411
	// comb = pcs_component_number;

	if ((int)cp_vec.size() < pcs_component_number + 1)
	{
		cout << "Error in CRFProcess::ConfigMassTransport - not enough MCP data"
		     << endl;
		return;
	}
	//----------------------------------------------------------------------
	pcs_primary_function_name[0] =
	    cp_vec[pcs_component_number]->compname.c_str();
	// sprintf(pcs_primary_function_name[0], "%s",
	// cp_vec[pcs_component_number]->compname.c_str());
	pcs_primary_function_unit[0] = "kg/m3";              // SB
	/* SB: Eintrag component name in Ausgabestruktur */  // SB:todo : just one
	                                                     // phase todo:name
	                                                     /*
	                                                           pcs_primary_function_name[0] =
	                                                        GetTracerCompName(0,this->pcs_component_number-1);
	                                                           name_initial_condition_tracer_component =
	                                                        pcs_primary_function_name[0];
	                                                           pcs_ic_name_mass = pcs_primary_function_name[0];
	                                                      */
	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 2;  // SB3909
	pcs_secondary_function_name[0] = new char[80];
	char pcs_secondary_function_name_tmp[80];
	sprintf(pcs_secondary_function_name_tmp, "%s%li", "MASS_FLUX_", comp);
	pcs_secondary_function_name[0] = pcs_secondary_function_name_tmp;
	//      pcs_secondary_function_name[0] = "MASS_FLUX1";
	pcs_secondary_function_unit[0] = "kg/m3/s";
	pcs_secondary_function_timelevel[0] = 0;
	pcs_secondary_function_name[1] = new char[80];
	sprintf(pcs_secondary_function_name_tmp, "%s%li", "MASS_FLUX_", comp);
	pcs_secondary_function_name[1] = pcs_secondary_function_name_tmp;
	pcs_secondary_function_unit[1] = "kg/m3/s";
	pcs_secondary_function_timelevel[1] = 1;
	// KG44 added secondary function for adaptive time stepping
	if (adaption)
	{
		pcs_number_of_secondary_nvals = 3;
		pcs_secondary_function_name[2] = new char[80];
		sprintf(pcs_secondary_function_name_tmp, "%s%li", "CONC_BACK_", comp);
		pcs_secondary_function_name[2] = pcs_secondary_function_name_tmp;
		pcs_secondary_function_unit[2] = "kg/m3";
		pcs_secondary_function_timelevel[2] = 0;
	}
	// OK  LOPCalcSecondaryVariables_USER = MTM2CalcSecondaryVariables;
	// //SB:todo
	// 2 ELE values
	pcs_number_of_evals = 0;
//	  pcs_eval_name[0] = "Darcy velocity";
#ifdef REACTION_ELEMENT
	pcs_number_of_evals = 1;
	pcs_eval_name[0] = pcs_primary_function_name[0];
	pcs_eval_unit[0] = "mol/kgH2O";
#endif
	// 3 ELE matrices
	// NUM
	pcs_num_name[0] = "CONCENTRATION0";
	/* SB: immer solver properties der ersten Komponente nehmen */
	// SB ??
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 OK Implementation
        WW Splitted for processes
   last modified:
**************************************************************************/
void CRFProcess::ConfigHeatTransport()
{
	pcs_num_name[0] = "TEMPERATURE0";
	// NOD
	if ((int)continuum_vector.size() == 1)
	{
		pcs_number_of_primary_nvals = 1;
		pcs_primary_function_name[0] = "TEMPERATURE1";
		pcs_primary_function_unit[0] = "K";
		pcs_number_of_secondary_nvals = 0;
#ifdef REACTION_ELEMENT
		pcs_number_of_evals = 1;  // MX
		pcs_eval_name[0] = "TEMPERATURE1";
		pcs_eval_unit[0] = "K";
#endif
	}
	if ((int)continuum_vector.size() == 2)
	{
		pcs_number_of_primary_nvals = 2;
		pcs_primary_function_name[0] = "TEMPERATURE1";
		pcs_primary_function_unit[0] = "K";
		pcs_number_of_primary_nvals = 2;
		pcs_primary_function_name[1] = "TEMPERATURE2";
		pcs_primary_function_unit[1] = "K";
		pcs_number_of_secondary_nvals = 0;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2003 WW Implementation
   last modified:
**************************************************************************/
void CRFProcess::ConfigDeformation()
{
#ifndef USE_PETSC
	// Generate high order nodes for all elements.
	ScreenMessage("-> generate higher order nodes\n");
	m_msh->GenerateHighOrderNodes();
	ScreenMessage("-> %d quadratic nodes are created. the number of total nodes becomes %d\n",
				  m_msh->GetNodesNumber(true) - m_msh->GetNodesNumber(false), m_msh->GetNodesNumber(true));
#endif
	type = 4;
	if (getProcessType() == FiniteElement::DEFORMATION_FLOW ||
	    getProcessType() == FiniteElement::DEFORMATION_H2)
	{
		type = 41;
		if (getProcessType() == FiniteElement::DEFORMATION_H2) type = 42;
		cpl_type_name = "MONOLITHIC";
		pcs_deformation = 11;
	}

	CNumerics* num = NULL;

	for (size_t ii = 0; ii < num_vector.size(); ii++)
	{
		num = num_vector[ii];
		if (num->pcs_type_name.find("DEFORMATION") != string::npos)
		{
			num->pcs_type_name = FiniteElement::convertProcessTypeToString(
			    this->getProcessType());
			if (FiniteElement::isNewtonKind(num->nls_method))  // Newton-Raphson
			{
				pcs_deformation = 101;
				if (type / 10 == 4) pcs_deformation = 110;
			}
			break;
		}
	}

	// Prepare for restart
	// RFConfigRenumber();

	// NUM
	pcs_num_name[0] = "DISPLACEMENT0";
	pcs_num_name[1] = "PRESSURE0";
	if (pcs_type_name_vector[0].find("DYNAMIC") != string::npos)
		VariableDynamics();
	else
		VariableStaticProblem();
	// OBJ names are set to PCS name
	// Geometry dimension
	problem_dimension_dm =
	    m_msh->GetMaxElementDim();  // m_msh->GetCoordinateFlag() / 10;
	problem_2d_plane_dm = 1;

	// Coupling
	int i;
	for (i = 0; i < problem_dimension_dm; i++)
		Shift[i] = i * m_msh->GetNodesNumber(true);

	/// 11-20.08.2010 WW
	long nn_H = (long)m_msh->GetNodesNumber(true);
	if (getProcessType() == FiniteElement::DEFORMATION)
	{
		num_nodes_p_var = new long[problem_dimension_dm];
		p_var_index = new int[problem_dimension_dm];
		for (int i = 0; i < problem_dimension_dm; i++)
			num_nodes_p_var[i] = nn_H;
	}
	else if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
	{
		num_nodes_p_var = new long[problem_dimension_dm + 1];
		p_var_index = new int[problem_dimension_dm + 1];
		for (i = 0; i < problem_dimension_dm; i++)
			num_nodes_p_var[i] = nn_H;
		num_nodes_p_var[problem_dimension_dm] =
		    (long)m_msh->GetNodesNumber(false);
		Shift[problem_dimension_dm] =
		    problem_dimension_dm * m_msh->GetNodesNumber(true);
	}
	else if (type == 42)
	{
		num_nodes_p_var = new long[problem_dimension_dm + 2];
		p_var_index = new int[problem_dimension_dm + 2];

		for (i = 0; i < problem_dimension_dm; i++)
			num_nodes_p_var[i] = nn_H;
		for (i = problem_dimension_dm; i < problem_dimension_dm + 2; i++)
			num_nodes_p_var[i] = (long)m_msh->GetNodesNumber(false);
		Shift[problem_dimension_dm] =
		    problem_dimension_dm * m_msh->GetNodesNumber(true);
		Shift[problem_dimension_dm + 1] =
		    Shift[problem_dimension_dm] + m_msh->GetNodesNumber(false);
	}

	if (type / 10 == 4)
		SetOBJNames();  // if(type==41) SetOBJNames(); //OK->WW please put to
	                    // Config()
}

/**************************************************************************
   FEMLib-Method: Static problems
   Task:
   Programing:
   05/2005 WW Implementation
   last modified:
**************************************************************************/
void CRFProcess::VariableStaticProblem()
{
	//----------------------------------------------------------------------
	// NOD Primary functions
	pcs_number_of_primary_nvals =
	    2;  // OK distinguish 2/3D problems, problem_dimension_dm;
	dm_number_of_primary_nvals = 2;
	pcs_number_of_evals = 0;
	pcs_primary_function_name[0] = "DISPLACEMENT_X1";
	pcs_primary_function_name[1] = "DISPLACEMENT_Y1";
	pcs_primary_function_unit[0] = "m";
	pcs_primary_function_unit[1] = "m";
	if (m_msh->GetMaxElementDim() == 3)
	{
		pcs_number_of_primary_nvals = 3;
		dm_number_of_primary_nvals = 3;
		pcs_primary_function_name[2] = "DISPLACEMENT_Z1";
		pcs_primary_function_unit[2] = "m";
	}
	//----------------------------------------------------------------------
	// NOD Secondary functions
	pcs_number_of_secondary_nvals = 0;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_PLS";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	//  pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
	//  "POROPRESSURE0";
	//  pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	//  pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	//  pcs_number_of_secondary_nvals++;

	if (m_msh->GetMaxElementDim() == 3)  // 3D
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRESS_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRESS_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRAIN_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "--";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRAIN_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "--";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}

	if (getProcessType() == FiniteElement::DEFORMATION_FLOW)
	{
		Def_Variable_LiquidFlow();
		configMaterialParameters();
	}
	else if (type == 42)  // Monolithic scheme H2M. 03.08.2010. WW
	{
		Def_Variable_MultiPhaseFlow();

		// Output material parameters
		configMaterialParameters();
	}
}

/**************************************************************************
   FEMLib-Method: Dynamic problems
   Task:
   Programing:
   05/2005 WW/LD Implementation
   last modified:
**************************************************************************/
void CRFProcess::VariableDynamics()
{
	//----------------------------------------------------------------------
	// NOD Primary functions
	pcs_number_of_primary_nvals = 2;
	dm_number_of_primary_nvals = 2;
	pcs_primary_function_name[0] = "ACCELERATION_X1";
	pcs_primary_function_name[1] = "ACCELERATION_Y1";
	pcs_primary_function_unit[0] = "m/s^2";
	pcs_primary_function_unit[1] = "m/s^2";
	if (max_dim == 2)
	{
		pcs_number_of_primary_nvals = 3;
		dm_number_of_primary_nvals = 3;
		pcs_primary_function_name[2] = "ACCELERATION_Z1";
		pcs_primary_function_unit[2] = "m/s^2";
	}
	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE_RATE1";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa/s";
	pcs_number_of_primary_nvals++;

	//----------------------------------------------------------------------
	// NOD Secondary functions
	pcs_number_of_secondary_nvals = 0;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRESS_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	//  pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
	//  "POROPRESSURE0";
	//  pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	//  pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	//  pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XX";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_XY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_YY";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_ZZ";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "STRAIN_PLS";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
	    "DISPLACEMENT_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
	    "DISPLACEMENT_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
	    "VELOCITY_DM_X";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
	    "VELOCITY_DM_Y";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// 3D
	if (max_dim == 2)  // 3D
	{
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRESS_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRESS_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRAIN_XZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "STRAIN_YZ";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "DISPLACEMENT_Z1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_DM_Z";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2005 OK Implementation
   02/2006 OK FLUX
**************************************************************************/
void CRFProcess::ConfigUnsaturatedFlow()
{
	if ((int)continuum_vector.size() == 1)
	{
		// 1.1 primary variables
		pcs_number_of_primary_nvals = 1;
		pcs_primary_function_name[0] = "PRESSURE1";
		pcs_primary_function_unit[0] = "Pa";
		// 1.2 secondary variables
		// OK LOPCalcSecondaryVariables_USER =
		// MMPCalcSecondaryVariablesRichards; // p_c and S^l
		pcs_number_of_secondary_nvals = 0;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "PRESSURE_CAP1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		// MB
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "FLUX";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		// MB
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "FLUX";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		// Nodal velocity. WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_X1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_Y1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_Z1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
		// JOD
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "COUPLING";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;
		// JOD
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "COUPLING";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;
	}
	else if ((int)continuum_vector.size() == 2)
	{
		dof = 2;  // WW
		// 1.1 primary variables
		pcs_number_of_primary_nvals = 2;  // YD
		pcs_primary_function_name[0] = "PRESSURE1";
		pcs_primary_function_unit[0] = "Pa";
		pcs_primary_function_name[1] = "PRESSURE2";
		pcs_primary_function_unit[1] = "Pa";
		// 1.2 secondary variables
		// OK LOPCalcSecondaryVariables_USER =
		// MMPCalcSecondaryVariablesRichards; // p_c and S^l
		pcs_number_of_secondary_nvals = 0;
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "SATURATION1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "SATURATION2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "SATURATION2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "PRESSURE_CAP1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "PRESSURE_CAP2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
		pcs_number_of_secondary_nvals++;  // WW
		// Nodal velocity. WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_X1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_Y1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_Z1";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_X2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_Y2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
		    "VELOCITY_Z2";
		pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
		pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
		pcs_number_of_secondary_nvals++;  // WW
		// 03.03.2008. WW
		for (size_t i = 0; i < GetPrimaryVNumber(); i++)
			Shift[i] = i * m_msh->GetNodesNumber(true);
	}

	// Output material parameters
	// WW // TF
	configMaterialParameters();

	// 2 ELE values
	pcs_number_of_evals = 8;
	pcs_eval_name[0] = "VELOCITY1_X";
	pcs_eval_unit[0] = "m/s";
	pcs_eval_name[1] = "VELOCITY1_Y";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Z";
	pcs_eval_unit[2] = "m/s";
	pcs_eval_name[3] = "POROSITY";  // MX 11.2005
	pcs_eval_unit[3] = "-";
	pcs_eval_name[4] = "POROSITY_IL";  // MX 11.2005
	pcs_eval_unit[4] = "-";
	pcs_eval_name[5] = "PERMEABILITY";  // MX 11.2005
	pcs_eval_unit[5] = "-";
	pcs_eval_name[6] = "n_sw";  // MX 11.2005
	pcs_eval_unit[6] = "-";
	pcs_eval_name[7] = "n_sw_rate";  // MX 11.2005
	pcs_eval_unit[7] = "-";
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 PCH Implementation
   last modified:
**************************************************************************/
void CRFProcess::ConfigFluidMomentum()
{
	// pcs_num_name[0] = "VELOCITY1_X";
	// Nothing added in terms of matrix solver.
	// Just linear solver is good enough.
	// NOD values
	pcs_number_of_primary_nvals = 3;
	pcs_primary_function_name[0] = "VELOCITY1_X";
	pcs_primary_function_unit[0] = "m/s";
	pcs_primary_function_name[1] = "VELOCITY1_Y";
	pcs_primary_function_unit[1] = "m/s";
	pcs_primary_function_name[2] = "VELOCITY1_Z";
	pcs_primary_function_unit[2] = "m/s";

	// I'm adding this to initialize for Fluid Momentum process
	pcs_number_of_secondary_nvals = 0;
	pcs_number_of_evals = 3;

	pcs_eval_name[0] = "VELOCITY1_X";
	pcs_eval_unit[0] = "m/s";
	pcs_eval_name[1] = "VELOCITY1_Y";
	pcs_eval_unit[1] = "m/s";
	pcs_eval_name[2] = "VELOCITY1_Z";
	pcs_eval_unit[2] = "m/s";
}

/**************************************************************************/
void CRFProcess::ConfigRandomWalk()
{
	// Nothing added in terms of matrix solver.
	// Just linear solver is good enough.

	// NOD values
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;

	// 2 ELE values
	pcs_number_of_evals = 1;
	pcs_eval_name[0] = "CONCENTRATION0";
	pcs_eval_unit[0] = "kg/m3";
}

////////////////////////////////////////////////////////////////////////////
//
///  Define variables of multi-phase flow model  (WW 08.2010)
//
////////////////////////////////////////////////////////////////////////////
void CRFProcess::Def_Variable_MultiPhaseFlow()
{
	// 1.1 primary variables
	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE1";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
	pcs_number_of_primary_nvals++;

	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE2";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
	pcs_number_of_primary_nvals++;

	// 1.2 secondary variables
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_W";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	// 1.3 elemental variables								// BG, 04/2012
	// pcs_number_of_evals = 0;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY2_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY2_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY2_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
}

////////////////////////////////////////////////////////////////////////////
//
///  Define variables of liquid flow model  (NW 09.2011)
//
////////////////////////////////////////////////////////////////////////////
void CRFProcess::Def_Variable_LiquidFlow()
{
	// 1.1 primary variables
	pcs_primary_function_name[pcs_number_of_primary_nvals] = "PRESSURE1";
	pcs_primary_function_unit[pcs_number_of_primary_nvals] = "Pa";
	pcs_number_of_primary_nvals++;

	// 1.2 secondary variables
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "HEAD";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	// 1.3 elemental variables
	// pcs_number_of_evals = 0;
	pcs_eval_name[pcs_number_of_evals] = "VOLUME";
	pcs_eval_unit[pcs_number_of_evals] = "m3";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_X";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Y";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] = "VELOCITY1_Z";
	pcs_eval_unit[pcs_number_of_evals] = "m/s";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] =
	    "POROSITY";  // MX, test for n=n(c), 04.2005
	pcs_eval_unit[pcs_number_of_evals] = "-";
	pcs_number_of_evals++;
	pcs_eval_name[pcs_number_of_evals] =
	    "PERMEABILITY";  // JT 2010 -- need this for index call of heterogeneous
	                     // permeability
	pcs_eval_unit[pcs_number_of_evals] = "m2";
	pcs_number_of_evals++;
}

/**************************************************************************
   FEMLib-Method: For non-isothermal multi-phase flow
   Task:
   Programing:
   02/2007 WW Implementation
   04/2011 WW Apdate for H2M
**************************************************************************/
void CRFProcess::ConfigMultiPhaseFlow()
{
	dof = 2;
	pcs_number_of_primary_nvals = 0;
	pcs_number_of_secondary_nvals = 0;

	Def_Variable_MultiPhaseFlow();

	// Output material parameters
	configMaterialParameters();

	// 11.08.2010. WW
	long nn = m_msh->GetNodesNumber(false);
	//
	for (size_t i = 0; i < GetPrimaryVNumber(); i++)  // 03.03.2008. WW
		Shift[i] = i * nn;

	num_nodes_p_var = new long[2];
	num_nodes_p_var[0] = num_nodes_p_var[1] = nn;
}

/**************************************************************************
   FEMLib-Method: For PS model for multiphase flow
   Task:
   Programing:
   03/2009 PCH Implementation
**************************************************************************/
void CRFProcess::ConfigPS_Global()
{
	dof = 2;
	// 1.1 primary variables
	pcs_number_of_primary_nvals = 2;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	pcs_primary_function_name[1] = "SATURATION2";
	pcs_primary_function_unit[1] = "m3/m3";
	// 1.2 secondary variables
	pcs_number_of_secondary_nvals = 0;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "PRESSURE_CAP";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "Pa";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 0;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "SATURATION1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m3/m3";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	// Nodal velocity.
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z2";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	//
	for (size_t i = 0; i < GetPrimaryVNumber(); i++)  // 03.03.2008. WW
		Shift[i] = i * m_msh->GetNodesNumber(true);
}

void CRFProcess::ConfigTH()
{
	dof = 2;
	pcs_number_of_primary_nvals = 2;
	pcs_primary_function_name[0] = "PRESSURE1";
	pcs_primary_function_unit[0] = "Pa";
	pcs_primary_function_name[1] = "TEMPERATURE1";
	pcs_primary_function_unit[1] = "K";
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_X1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Y1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;
	pcs_secondary_function_name[pcs_number_of_secondary_nvals] = "VELOCITY_Z1";
	pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "m/s";
	pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
	pcs_number_of_secondary_nvals++;

	//
	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
		Shift[i] = i * m_msh->GetNodesNumber(false);

	num_nodes_p_var = new long[2];
	num_nodes_p_var[0] = num_nodes_p_var[1] = m_msh->GetNodesNumber(false);
	p_var_index = new int[2];
}


/**************************************************************************
   FEMLib-Method:
   Task:  Activate or deactivate elements specified in .pcs file
   Programing:
   05/2005 WW Implementation
**************************************************************************/
void CRFProcess::CheckMarkedElement()
{
	for (CNode* node : m_msh->nod_vector)
		node->SetMark(false);

	for (CElem* elem : m_msh->ele_vector)
	{
		bool done = false;
		for (auto mat_id : Deactivated_SubDomain)
		{
			if (elem->GetPatchIndex() == (unsigned)mat_id)
			{
				elem->SetMark(false);
				//elem->MarkingAll(false);
				done = true;
				break;
			}
		}
		if (!done)
			elem->MarkingAll(true);
	}

	for (auto node : m_msh->nod_vector)
		node->getConnectedElementIDs().clear();

	for (size_t l = 0; l < m_msh->ele_vector.size(); l++)
	{
		CElem* elem = m_msh->ele_vector[l];
		if (!elem->GetMark()) continue;
		for (size_t i = 0; i < elem->GetNodesNumber(m_msh->getOrder()); i++)
		{
			bool done = false;
			auto node = elem->GetNode(i);
			for (size_t j = 0; j < node->getConnectedElementIDs().size(); j++)
			{
				if (l == node->getConnectedElementIDs()[j])
				{
					done = true;
					break;
				}
			}
			if (!done)
				node->getConnectedElementIDs().push_back(l);
		}
	}
}


/**************************************************************************
   FEMLib-Method:
   Task:  check the excavation state of each aktive element
   Programing:
   01/2011 WX Implementation
**************************************************************************/
void CRFProcess::CheckExcavedElement()
{
#ifndef OGS_ONLY_TH
	int valid;
	long l;
	// bool done;
	CElem* elem = NULL;
	// CNode *node = NULL;
	for (l = 0; l < (long)m_msh->ele_vector.size(); l++)
	{
		elem = m_msh->ele_vector[l];
		if (elem->GetPatchIndex() == static_cast<size_t>(ExcavMaterialGroup) &&
		    elem->GetMark())
		{
			double const* ele_center(elem->GetGravityCenter());
			if ((GetCurveValue(ExcavCurve, 0, aktuelle_zeit, &valid) +
			     ExcavBeginCoordinate) > (ele_center[ExcavDirection]) &&
			    (ele_center[ExcavDirection] - ExcavBeginCoordinate) > -0.001)
				elem->SetExcavState(1);
		}
	}
#endif
}

//////////////////////////////////////////////////////////////////////////
// PCS Execution
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   02/2003 OK Implementation
   04/2003 OK Storing time step results
   08/2004 OK PCS2
   01/2005 WW new ELE concept
   03/2005 OK MultiMSH concept
   06/2005 MB NEWTON error calculation
   05/2007 WW DOF>1, unerrelaxation for Picard, and removing the old or spurious
 stuff
   12/2007 WW Classes of sparse matrix (jagged diagonal storage) and linear
 solver
   and parellelisation of them
   07/2011 WW Add Newton-Rahson and reduce #ifdef
   3/2012  JT Clean, and correct the error obtainment for NLS and CPL
   last modified:
 **************************************************************************/
double CRFProcess::Execute()
{
#ifndef USE_PETSC
	int ii;
	double val_n;
	long nshift;
#endif
	int nidx1;
	double pcs_error, nl_theta;
	long j, k, g_nnodes;  // 07.01.07 WW
	double* eqs_x = NULL;

	pcs_error = DBL_MAX;
	g_nnodes = m_msh->GetNodesNumber(false);

#ifdef NEW_EQS
	double implicit_lim = 1.0 - DBL_EPSILON;
	eqs_x = eqs_new->getX();
#endif


#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
	eqs_new->Initialize();
#elif defined(NEW_EQS)
	if (!configured_in_nonlinearloop)
	// Also allocate temporary memory for linear solver. WW
	{
		//_new 02/2010. WW
		eqs_new->SetDOF(pcs_number_of_primary_nvals);
		eqs_new->ConfigNumerics(m_num->ls_precond, m_num->ls_method, m_num->ls_max_iterations, m_num->ls_error_tolerance, m_num->ls_storage_method, m_num->ls_extra_arg);
	}
	eqs_new->Initialize();
#endif
	/*
	   //TEST_MPI
	   string test = "rank";
	   char stro[1028];
	   sprintf(stro, "%d",myrank);
	   string test1 = test+(string)stro+"Assemble.txt";
	   ofstream Dum(test1.c_str(), ios::out); // WW
	   MPI_Finalize();
	   exit(1);
	 */
	//----------------------------------------------------------------------
	nl_theta = 1.0 - m_num->nls_relaxation;
	if (nl_theta < DBL_EPSILON) nl_theta = 1.0;

	// NW. should mark active elements if any process uses deactivation
	// if(NumDeactivated_SubDomains>0)
	// TODO if it's nonlinear, CheckMarkedElement() has been already called
	if (hasAnyProcessDeactivatedSubdomains)
#ifdef NEW_EQS  // WW
		if (!configured_in_nonlinearloop)
#endif
			CheckMarkedElement();
	m_msh->SwitchOnQuadraticNodes(false);

	// If not Newton-Raphson method. 20.07.2011. WW
	if (!FiniteElement::isNewtonKind(m_num->nls_method))
	{
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
		InitializeRHS_with_u0();
#else
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			long const ish = ii * g_nnodes;
			for (j = 0; j < g_nnodes; j++)  // WW
				eqs_x[j + ish] =
				    GetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1);
		}
#endif
	}

#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	ScreenMessage("\tcurrent mem: %d MB\n", mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
//---------------------------------------------------------------------
// Assembly
#if defined(USE_MPI) || defined(USE_PETSC)  // WW
	clock_t cpu_time = 0;                   // WW
	cpu_time = -clock();
	if (myrank == 0)
#endif
		cout << "Assembling equation system..." << endl;
	GlobalAssembly();
#ifndef WIN32
	ScreenMessage("\tcurrent mem: %d MB\n", mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
#if defined(USE_MPI) || defined(USE_PETSC)  // WW
	cpu_time += clock();
	cpu_time_assembly += cpu_time;
	if (myrank == 0)
#endif
		cout << "Calling linear solver..." << endl;
//----------------------------------------------------------------------
// Execute linear solver
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
	if (write_leqs)
	{
		std::string eqs_name_base =
		    convertProcessTypeToString(this->getProcessType());
		eqs_name_base += "_" + number2str<size_t>(aktueller_zeitschritt);
		eqs_new->EQSV_Viewer(eqs_name_base, true);
	}
#if 0
	if (aktueller_zeitschritt==19 && this->getProcessType()==FiniteElement::HEAT_TRANSPORT) {
		eqs_new->CheckIfMatrixIsSame("/home/localadmin/tasks/20131217_LiWah/petsc/2units2faults/HEAT_TRANSPORT_19_eqs_A.dat");
	}
#endif
	iter_lin = eqs_new->Solver();
	// TEST 	double x_norm = eqs_new->GetVecNormX();
	eqs_new->MappingSolution();
#elif defined(NEW_EQS)  // WW
	bool compress_eqs = (type / 10 == 4 || this->Deactivated_SubDomain.size() > 0);
	iter_lin = eqs_new->Solver(compress_eqs);
#endif
#ifndef WIN32
	ScreenMessage("\tcurrent mem: %d MB\n", mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif

	if (iter_lin == -1)
	{
		ScreenMessage("*** Linear solve failed\n");
// abort
#ifdef NEW_EQS  // WW
		if (!configured_in_nonlinearloop)
			// Also allocate temporary memory for linear solver. WW
			eqs_new->Clean();
#endif
		return -1;
	}
	iter_lin_max = std::max(iter_lin_max, iter_lin);

	//----------------------------------------------------------------------
	// Linearized Flux corrected transport (FCT) by Kuzmin 2009
	//----------------------------------------------------------------------
	if (m_num->fct_method > 0)  // NW
	{
#if defined(USE_PETSC)
		eqs_x = eqs_new->GetGlobalSolution();
//		pcs_error = CalcIterationNODError(1);
#endif
		//#else
		pcs_error = CalcIterationNODError(m_num->getNonLinearErrorMethod(),
		                                  true, false);  // JT
//#endif

#if defined(USE_MPI) || defined(USE_PETSC)
		if (myrank == 0)
		{
#endif
			cout << "    Relative PCS error: " << pcs_error << "\n";
			cout << "    Start FCT calculation"
			     << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
		}
#endif
		// Set u^H: use the solution as the higher-order solution
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC)
				k = m_msh->Eqs2Global_NodeIndex[j] *
				        pcs_number_of_primary_nvals +
				    ii;
				SetNodeValue(j, nidx1, eqs_x[k]);
#else
				k = m_msh->Eqs2Global_NodeIndex[j];
				SetNodeValue(k, nidx1, eqs_x[j + ii * g_nnodes]);
#endif
			}
		}

// Initialize the algebra system
#if defined(USE_PETSC)
		eqs_new->Initialize();
#else
		if (!configured_in_nonlinearloop)
			eqs_new->Initialize();
#endif

// Set initial guess
#if !defined(USE_PETSC)
		for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
				k = m_msh->Eqs2Global_NodeIndex[j];
				eqs_x[j + ii * g_nnodes] = GetNodeValue(k, nidx1);
			}
		}
#endif

// Assembly
#if defined(USE_MPI) || defined(USE_PETSC)
		clock_t cpu_time = 0;  // WW
		cpu_time = -clock();
#endif
		femFCTmode = true;
		GlobalAssembly();
		femFCTmode = false;
#if defined(USE_MPI) || defined(USE_PETSC)
		cpu_time += clock();
		cpu_time_assembly += cpu_time;
#endif

// Solve the algebra
#if defined(USE_PETSC)
		//		std::string eqs_output_file = FileName +
		// number2str(aktueller_zeitschritt);
		//		eqs_new->EQSV_Viewer(eqs_output_file);
		eqs_new->Solver();
		eqs_new->MappingSolution();
#else
#ifdef NEW_EQS  // WW
		eqs_new->Solver(this->m_num);  // NW
#endif
#endif  // USE_PETSC
	}
	//----------------------------------------------------------------------
	// END OF FLUX CORRECTED TRANSPORT
	//----------------------------------------------------------------------

	// PCSDumpModelNodeValues();
	//----------------------------------------------------------------------
	// ERROR CALCULATION
	//----------------------------------------------------------------------

	//
	// Save the solution of the prevoius iteration of the nonlinear step for
	// the automatic time stepping. Modified for PETsc solver. 03.07.2012. WW
	if (Tim->GetPITimeStepCrtlType() > 0)
	{
		double* x_k = NULL;
		bool get_buffer_u_k = true;
		x_k = _problem->GetBufferArray(get_buffer_u_k);
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)
			const long ish = i * g_nnodes;
#endif
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC)  // || defined(other parallel libs)
				x_k[j * pcs_number_of_primary_nvals + i] =
				    GetNodeValue(j, nidx1);
#else
				x_k[j + ish] =
				    GetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1);
#endif
			}
		}
	}
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
	// PICARD
	//----------------------------------------------------------------------
	// Error calculation
	//----------------------------------------------------------------------
	if (m_num->nls_method_name.find("PICARD") != string::npos)
	{
		eqs_x = eqs_new->GetGlobalSolution();
		//......................................................................
		//		pcs_error = CalcIterationNODError(1); //OK4105//WW4117
		if (iter_nlin > 0)
		{  // Just getting NL error
			pcs_error =
			    CalcIterationNODError(m_num->getNonLinearErrorMethod(), true,
			                          false);  // OK4105//WW4117//JT
		}
		else
		{  // Getting NL and CPL error
			pcs_error = CalcIterationNODError(m_num->getCouplingErrorMethod(),
			                                  true, true);  // JT2012
			if (m_num->getNonLinearErrorMethod() !=
			    m_num->getCouplingErrorMethod())  // JT: If CPL error method is
				                                  // different, must call
				                                  // separately
				pcs_error = CalcIterationNODError(
				    m_num->getNonLinearErrorMethod(), true,
				    false);  // JT2012 // get the NLS error. CPL was obtained
			                 // before.
		}
		//		ScreenMessage("PCS error: %g\n", pcs_error);

		//--------------------------------------------------------------------
		// 7 Store solution vector in model node values table
		//....................................................................
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
			for (j = 0; j < g_nnodes; j++)
			{
				k = m_msh->Eqs2Global_NodeIndex[j] *
				        pcs_number_of_primary_nvals +
				    i;
				SetNodeValue(j, nidx1,
				             (1. - nl_theta) * GetNodeValue(j, nidx1) +
				                 nl_theta * eqs_x[k]);
			}
		}

	}  // END PICARD
#else
	// JT: Coupling error was wrong. Now ok.
	if (iter_nlin > 0)
	{  // Just getting NL error
		pcs_error = CalcIterationNODError(m_num->getNonLinearErrorMethod(),
		                                  true, false);  // OK4105//WW4117//JT
	}
	else
	{  // Getting NL and CPL error
		pcs_error = CalcIterationNODError(m_num->getCouplingErrorMethod(), true,
		                                  true);  // JT2012
		if (m_num->getNonLinearErrorMethod() !=
		    m_num->getCouplingErrorMethod())  // JT: If CPL error method is
			                                  // different, must call separately
			pcs_error =
			    CalcIterationNODError(m_num->getNonLinearErrorMethod(), true,
			                          false);  // JT2012 // get the NLS error.
		                                       // CPL was obtained before.
	}

	//----------------------------------------------------------------------
	// PICARD
	//----------------------------------------------------------------------
	if (m_num->nls_method_name.find("PICARD") != string::npos)
	{
		if (pcs_error <
		    1.0)  // JT: Then the solution has converged, take the final value
			nl_theta = 1.0;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
		eqs_x = eqs_new->GetGlobalSolution();
#endif
		//
		if (nl_theta > implicit_lim)  // This is most common. So go for the
		                              // lesser calculations.
		{
			for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
				const long nshift = ii * g_nnodes;
				for (j = 0; j < g_nnodes; j++)
				{
					k = m_msh->Eqs2Global_NodeIndex[j];
					const double val_n =
					    GetNodeValue(k, nidx1);  // 03.04.2009. WW
					SetNodeValue(k, nidx1, eqs_x[j + nshift]);
					eqs_x[j + nshift] =
					    val_n;  // Used for time stepping. 03.04.2009. WW
				}
			}
		}
		else
		{
			for (int ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
				const long nshift = ii * g_nnodes;
				for (j = 0; j < g_nnodes; j++)
				{
					k = m_msh->Eqs2Global_NodeIndex[j];
					const double val_n =
					    GetNodeValue(k, nidx1);  // 03.04.2009. WW
					SetNodeValue(k, nidx1, (1.0 - nl_theta) * val_n +
					                           nl_theta * eqs_x[j + nshift]);
					eqs_x[j + nshift] =
					    val_n;  // Used for time stepping. 03.04.2009. WW
				}
			}
		}
	}
//----------------------------------------------------------------------
// END OF PICARD
//----------------------------------------------------------------------
#endif

#ifdef NEW_EQS  // WW
	if (!configured_in_nonlinearloop)
		// Also allocate temporary memory for linear solver. WW
		eqs_new->Clean();
#endif

	if (this->getProcessType() == FiniteElement::LIQUID_FLOW &&
	    Tim->usePIDControlInSelfAdaptive)
		this->CalIntegrationPointValue();

	return pcs_error;
}

/*************************************************************************
   GEOSYS - Function:
   Task:
   Programming:
   08/2008 WW Implementation
   11/2008 WW Update
   last modified:
 **************************************************************************/
void CRFProcess::CopyU_n()
{
	int i, nidx1;
	long g_nnodes, j;

	double* temp_v = _problem->GetBufferArray();

	for (i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		/// H2M with monilithic scheme. 02.2011. WW
		if (type == 42)
		{
			nidx1 = p_var_index[i] + 1;
			g_nnodes = num_nodes_p_var[i];
		}
		else
		{
			nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
			g_nnodes = m_msh->GetNodesNumber(false);  // DOF>1, WW
		}
#if !defined(USE_PETSC)  // && !defined(other parallel libs)
		const long ish = i * g_nnodes;
#endif
		for (j = 0; j < g_nnodes; j++)
		{
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
			temp_v[j * pcs_number_of_primary_nvals + i] =
			    GetNodeValue(j, nidx1);
#else
			temp_v[j + ish] =
			    GetNodeValue(m_msh->Eqs2Global_NodeIndex[j], nidx1);
#endif
		}
	}
}


/*************************************************************************
   GeoSys-Function:
   Task: Algebraic operation for the flux corrected transport (FCT)
   Programming:
   04/2010 NW Implementation
   last modified:
   05/2013 NW Support PETSc parallelization
 **************************************************************************/
void CRFProcess::AddFCT_CorrectionVector()
{
	int idx0 = 0;
	int idx1 = idx0 + 1;
	const double theta = this->m_num->ls_theta;
	const size_t node_size = m_msh->GetNodesNumber(false);
	SparseMatrixDOK::mat_t& fct_f = this->FCT_AFlux->GetRawData();
	SparseMatrixDOK::col_t* col;
	SparseMatrixDOK::mat_t::const_iterator ii;
	SparseMatrixDOK::col_t::const_iterator jj;
	Math_Group::Vector* ML = this->Gl_ML;
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;  // WW
	A = this->eqs_new->getA();
#endif

#ifdef USE_PETSC
	// gather K
	FCT_MPI::gatherK(FCT_MPI::ct, *FCT_K);
	// compute D
	FCT_MPI::computeD(m_msh, *FCT_K, *FCT_d);
#endif

	// List of Dirichlet nodes
	std::set<long> list_bc_nodes;
	for (size_t i = 0; i < bc_node_value.size(); i++)
	{
		CBoundaryConditionNode* bc_node = bc_node_value[i];
		long nod_id = bc_node->geo_node_number;
		list_bc_nodes.insert(nod_id);
	}

	//----------------------------------------------------------------------
	// Construct global matrices: antidiffusive flux(f_ij), positivity matrix(L)
	// - f_ij =
	// 1/dt*m_ij*(DeltaU_ij^H-DeltaU_ij^n)-theta*d_ij^H*DeltaU_ij^H-(1-theta)*d_ij^n*DeltaU_ij^n
	// - L = K + D
	// - D_ij = min(0, -K_ij, -K_ji)
	// * f_ji = -f_ij
	// * K:original coefficient matrix, D:artificial diffusion operator
	// Implementation memo:
	// - K is stored in A matrix in the element assembly.
	// - the first part of the antidiffusive flux is done in the element
	// assembly.
	//   -> f_ij = m_ij
	//----------------------------------------------------------------------
	// f_ij*=1/dt*(DeltaU_ij^H-DeltaU_ij^n)  for i!=j
	for (size_t i = 0; i < node_size; i++)
	{
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			if (i > j) continue;  // symmetric part, off-diagonal
			double diff_uH =
			    this->GetNodeValue(i, idx1) - this->GetNodeValue(j, idx1);
			double diff_u0 =
			    this->GetNodeValue(i, idx0) - this->GetNodeValue(j, idx0);
			double v = 1.0 / dt * (diff_uH - diff_u0);
			(*FCT_AFlux)(i, j) *= v;  // MC is already done in local ele
			                          // assembly
			(*FCT_AFlux)(j, i) *=
			    -v;  // MC is already done in local ele assembly
		}
	}

	// Complete f, L
	// Remark: Using iteration is only possible after the sparse table has been
	// constructed.
	for (size_t i = 0; i < node_size; i++)
	{
#ifdef USE_PETSC
		const size_t i_global = FCT_GLOB_ADDRESS(i);
#endif
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
#ifdef USE_PETSC
			const size_t j_global = FCT_GLOB_ADDRESS(j);
#endif
			if (i > j || i == j)
				continue;  // do below only for upper triangle due to symmetric

// Get artificial diffusion operator D
#ifdef USE_PETSC
			double d1 = (*FCT_d)(i_global, j_global);
#else
#if defined(NEW_EQS)
			double K_ij = (*A)(i, j);
			double K_ji = (*A)(j, i);
#else
			double K_ij = MXGet(i, j);
			double K_ji = MXGet(j, i);
#endif
			if (K_ij == 0.0 && K_ji == 0.0) continue;
			double d1 = GetFCTADiff(K_ij, K_ji);
#endif
			if (d1 == 0.0) continue;
			double d0 =
			    d1;  // TODO should use AuxMatrix at the previous time step
			// if (list_bc_nodes.find(i)!=list_bc_nodes.end() ||
			// list_bc_nodes.find(j)!=list_bc_nodes.end()) {
			//  d1 = d0 = 0.0;
			//}

			// Complete antidiffusive flux: f_ij += -theta*d_ij^H*DeltaU_ij^H -
			// (1-theta)*d_ij^n*DeltaU_ij^n
			double diff_uH =
			    this->GetNodeValue(i, idx1) - this->GetNodeValue(j, idx1);
			double diff_u0 =
			    this->GetNodeValue(i, idx0) - this->GetNodeValue(j, idx0);
			double v = -(theta * d1 * diff_uH + (1.0 - theta) * d0 * diff_u0);
			(*FCT_AFlux)(i, j) += v;

			// prelimiting f
			v = (*FCT_AFlux)(i, j);
			if (this->m_num->fct_prelimiter_type == 0)
			{
				if (v * (-diff_uH) > 0.0) v = 0.0;
			}
			else if (this->m_num->fct_prelimiter_type == 1)
				v = MinMod(v, -d1 * diff_uH);
			else if (this->m_num->fct_prelimiter_type == 2)
				v = SuperBee(v, -d1 * diff_uH);
			(*FCT_AFlux)(i, j) = v;
#ifdef USE_PETSC
			(*FCT_AFlux)(j, i) = -v;
#else
			(*FCT_AFlux)(j, i) = v;
#endif

#ifdef USE_PETSC
			// A += theta * D
			if (i < (size_t)m_msh->getNumNodesLocal())
			{
				eqs_new->addMatrixEntry(i_global, i_global, -d1 * theta);
				eqs_new->addMatrixEntry(i_global, j_global, d1 * theta);
			}
			if (j < (size_t)m_msh->getNumNodesLocal())
			{
				eqs_new->addMatrixEntry(j_global, i_global, d1 * theta);
				eqs_new->addMatrixEntry(j_global, j_global, -d1 * theta);
			}
#else
// L = K + D
#if defined(NEW_EQS)
			(*A)(i, i) += -d1;
			(*A)(i, j) += d1;
			(*A)(j, i) += d1;
			(*A)(j, j) += -d1;
#endif
#endif
		}
	}

	//----------------------------------------------------------------------
	// Assemble RHS: b_i += [- (1-theta) * L_ij] u_j^n
	//----------------------------------------------------------------------
	Math_Group::Vector* V1 = this->Gl_Vec1;
	Math_Group::Vector* V = this->Gl_Vec;
	(*V1) = 0.0;
	(*V) = 0.0;
#ifdef NEW_EQS
	double* eqs_rhs;
	eqs_rhs = eqs_new->getRHS();
#endif
	// b = [-(1-theta) * L] u^n
	if (1.0 - theta > .0)
	{
		// u^n
		for (size_t i = 0; i < node_size; i++)
			(*V1)(i) = this->GetNodeValue(i, idx0);
		// L*u^n
		for (size_t i = 0; i < node_size; i++)
		{
#ifdef USE_PETSC
			const size_t i_global = FCT_GLOB_ADDRESS(i);
#endif
			for (size_t j = 0; j < node_size; j++)
			{
#ifdef USE_PETSC
				const size_t j_global = FCT_GLOB_ADDRESS(j);
				// b+=-(1-theta)*D*u^n
				(*V)(i) += (*FCT_d)(i_global, j_global) * (*V1)(j);
#else
				(*V)(i) += (*A)(i, j) * (*V1)(j);
#endif
			}
		}
		for (size_t i = 0; i < node_size; i++)
		{
#if defined(USE_PETSC)
			if (i < (size_t)m_msh->getNumNodesLocal())
			{
				const size_t i_global = FCT_GLOB_ADDRESS(i);
				eqs_new->add_bVectorEntry(i_global, -(1.0 - theta) * (*V)(i),
				                          ADD_VALUES);
			}
#else
			eqs_rhs[i] -= (1.0 - theta) * (*V)(i);
//(*RHS)(i+LocalShift) +=  NodalVal[i];
#endif
		}
	}

#ifndef USE_PETSC
	//----------------------------------------------------------------------
	// Assemble A matrix: 1/dt*ML + theta * L
	//----------------------------------------------------------------------
	// A matrix: theta * L
	if (theta == 0.0)
	{
#ifdef NEW_EQS
		(*A) = 0.0;
#endif
	}
	else if (theta == 1.0)
	{
		// keep
	}
	else
	{
#ifdef NEW_EQS
		(*A) *= theta;
#endif
	}
	// A matrix: += 1/dt * ML
	for (size_t i = 0; i < node_size; i++)
	{
		double v = 1.0 / dt * (*ML)(i);
#ifdef NEW_EQS
		(*A)(i, i) += v;
#endif
	}
#endif

	//----------------------------------------------------------------------
	// Assemble RHS: b += alpha * f
	//----------------------------------------------------------------------
	// Calculate R+, R-
	Math_Group::Vector* R_plus = this->Gl_Vec1;
	Math_Group::Vector* R_min = this->Gl_Vec;
	(*R_plus) = 0.0;
	(*R_min) = 0.0;
	for (size_t i = 0; i < node_size; i++)
	{
		const size_t i_global = FCT_GLOB_ADDRESS(i);
		double P_plus, P_min;
		double Q_plus, Q_min;
		P_plus = P_min = 0.0;
		Q_plus = Q_min = 0.0;
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			if (i == j) continue;
			double f = (*jj).second;  // double f = (*FCT_AFlux)(i,j);
#ifndef USE_PETSC
			if (i > j) f *= -1.0;
#endif
			double diff_uH =
			    this->GetNodeValue(j, idx1) - this->GetNodeValue(i, idx1);

			P_plus += max(0.0, f);
			P_min += min(0.0, f);
			Q_plus = max(Q_plus, diff_uH);
			Q_min = min(Q_min, diff_uH);
		}
		double ml = (*ML)(i_global);

		if (P_plus == 0.0)
			(*R_plus)(i_global) = 0.0;
		else
			(*R_plus)(i_global) = min(1.0, ml * Q_plus / (dt * P_plus));
		if (P_min == 0.0)
			(*R_min)(i_global) = 0.0;
		else
			(*R_min)(i_global) = min(1.0, ml * Q_min / (dt * P_min));
	}

#ifdef USE_PETSC
	FCT_MPI::gatherR(FCT_MPI::ct, *R_plus, *R_min);
#endif

	// for Dirichlet nodes
	for (size_t i = 0; i < bc_node_value.size(); i++)
	{
		CBoundaryConditionNode* bc_node = bc_node_value[i];
		long nod_id = bc_node->geo_node_number;
		const size_t i_global = FCT_GLOB_ADDRESS(nod_id);
		(*R_plus)(i_global) = 1.0;
		(*R_min)(i_global) = 1.0;
	}

	// b_i += alpha_i * f_ij
	for (size_t i = 0; i < node_size; i++)
	{
		const size_t i_global = FCT_GLOB_ADDRESS(i);
		col = &fct_f[i];
		for (jj = col->begin(); jj != col->end(); jj++)
		{
			const size_t j = (*jj).first;
			const size_t j_global = FCT_GLOB_ADDRESS(j);
			if (i == j) continue;

			double f = (*jj).second;  // double f = (*FCT_AFlux)(i,j);
#ifndef USE_PETSC
			if (i > j) f *= -1;  // symmetric
#endif
			double alpha = 1.0;
			if (f > 0)
				alpha = min((*R_plus)(i_global), (*R_min)(j_global));
			else
				alpha = min((*R_plus)(j_global), (*R_min)(i_global));

			double val = .0;
			if (this->m_num->fct_const_alpha < 0.0)
				val = alpha * f;
			else
				val = this->m_num->fct_const_alpha * f;

#ifdef USE_PETSC
			if (i < (size_t)m_msh->getNumNodesLocal())
				eqs_new->add_bVectorEntry(i_global, val, ADD_VALUES);
#else
			eqs_rhs[i] += val;
#endif

			// Note: Galerkin FEM is recovered if alpha = 1 as below,
			// eqs_rhs[i] += 1.0*f;
		}
	}
}

/*************************************************************************
   GeoSys-Function:
   Task: Assemble the global system equation
   Programming:
   01/2005 WW/OK Implementation
   04/2005 OK MSH
   07/2005 WW Change due to the geometry element objects applied
   10/2005 OK DDC
   11/2005 YD time step control
   01/2006 OK/TK Tests
   12/2007 WW Spase matrix class and condensation sequential and parallel loop
   10/2010 TF changed access to process type
   06/2012 WW Node based decompostion
 **************************************************************************/
void CRFProcess::GlobalAssembly()
{
	// Tests
	if (!Tim) Tim = TIMGet(convertProcessTypeToString(this->getProcessType()));
	if (!Tim)
	{
		cout << "Error in CRFProcess::GlobalAssembly() - no TIM data" << endl;
		return;
	}
	if (Write_Matrix)
	{
		std::string pcs_type_name(
		    convertProcessTypeToString(this->getProcessType()));
#if defined(USE_MPI)
		char stro[32];
		sprintf(stro, "%d", myrank);
		string m_file_name = FileName + "_" + pcs_type_name + (string)stro +
		                     "_element_matrix.txt";
#else
		std::string m_file_name =
		    FileName + "_" + pcs_type_name + "_element_matrix.txt";
#endif
#if 0
		if (matrix_file) matrix_file->close();
		delete matrix_file;
		matrix_file = new std::fstream(m_file_name.c_str(), ios::trunc | ios::out);
		if (!matrix_file->good())
			std::cout << "Warning in GlobalAssembly: Matrix files are not found"
			          << "\n";
#endif
	}

	assert(fem);
#if 0
	if (!fem)
		// Which process needs this?
		// Only one instance of CFiniteElementStd is required for each process
		// Use "new" in such way will cause memory problem.
		// Please move this declaration to pcs configuration.     WW
		if (m_msh)
			fem = new CFiniteElementStd(this, m_msh->GetCoordinateFlag());
#endif

	CElem* elem = NULL;

	bool Check2D3D;
	Check2D3D = false;
	if (this->femFCTmode)  // NW
	{
		(*this->FCT_AFlux) = 0.0;
		(*this->Gl_ML) = 0.0;
		(*this->Gl_Vec) = 0.0;
		(*this->Gl_Vec1) = 0.0;
#ifdef USE_PETSC
		(*this->FCT_K) = 0.0;
		(*this->FCT_d) = .0;
#endif
	}

	{  // STD
	const size_t dn = m_msh->ele_vector.size() / 10;
	const bool print_progress = (dn >= 100);
	if (print_progress)
		ScreenMessage("start local assembly for %d elements...\n",
		              m_msh->ele_vector.size());
	// const long n_eles = (long)m_msh->ele_vector.size();

#ifdef USE_PETSC
	for (size_t i = 0; i < eqs_new->vec_subRHS.size(); i++)
	{
		VecGetSubVector(eqs_new->b, eqs_new->vec_isg[i],
		                &eqs_new->vec_subRHS[i]);
	}
#endif

	// YDTEST. Changed to DOF 15.02.2007 WW
	for (size_t ii = 0; ii < continuum_vector.size(); ii++)
	{
		continuum = ii;
		for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
		{
			if (print_progress && (i + 1) % dn == 0) ScreenMessage("* ");
			// ScreenMessage("%d \%\n", ((i+1)*100/n_eles));
			elem = m_msh->ele_vector[i];
// Marked for use //WX: modified for coupled excavation
#ifndef OGS_ONLY_TH
			if (elem->GetMark() && elem->GetExcavState() == -1)
#else
				if (elem->GetMark())
#endif
			{
				elem->SetOrder(false);
				fem->ConfigElement(elem, Check2D3D);
				fem->Assembly();
			}
		}
	}
	ScreenMessage("done\n");

#ifdef USE_PETSC
	for (size_t i = 0; i < eqs_new->vec_subRHS.size(); i++)
	{
		VecAssemblyBegin(eqs_new->vec_subRHS[i]);
		VecAssemblyEnd(eqs_new->vec_subRHS[i]);
		VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[i],
		                    &eqs_new->vec_subRHS[i]);
	}
	eqs_new->AssembleRHS_PETSc();
#endif

	if (femFCTmode)
		AddFCT_CorrectionVector();

	if (write_leqs)
	{
		std::string fname = FileName + "_" +
		                    convertProcessTypeToString(this->getProcessType()) +
		                    number2str(aktueller_zeitschritt) + "_" +
		                    number2str(iter_nlin) + "_leqs_assembly.txt";
#if defined(NEW_EQS)
		std::ofstream Dum(fname.c_str(), ios::out);
		eqs_new->Write(Dum);
		Dum.close();
#elif defined(USE_PETSC)
		eqs_new->EQSV_Viewer(fname);
#endif
	}

	// eqs_new->Write();
	ScreenMessage("-> impose Neumann BC and source/sink terms\n");
	IncorporateSourceTerms();
	if (false && write_leqs)
	{
		std::string fname = FileName + "_" +
		                    convertProcessTypeToString(this->getProcessType()) +
		                    number2str(aktueller_zeitschritt) + "_" +
		                    number2str(iter_nlin) + "_leqs_st.txt";
#if defined(NEW_EQS)
		std::ofstream Dum(fname.c_str(), ios::out);
		eqs_new->Write(Dum);
		Dum.close();
#elif defined(USE_PETSC)
			eqs_new->EQSV_Viewer(fname);
#endif
	}

#ifdef GEM_REACT
	if (getProcessType() == FiniteElement::MASS_TRANSPORT &&
	    aktueller_zeitschritt > 1)
	{  // JT->KG. New coupling system.
		if (_problem->GetCPLMaxIterations() >
		        1 ||  // Then there is a coupling on all processes
		    (this->pcs_is_cpl_overlord &&
		     this->m_num->cpl_max_iterations >
		         1) ||  // This process (the overlord) controls coupling for
		                // another (the underling) process.
		    (this->pcs_is_cpl_underling &&
		     this->cpl_overlord->m_num->cpl_max_iterations >
		         1))  // This process (the underling) has a coupling with a
		              // controlling (the overlord) process
		{
			// IncorporateSourceTerms_GEMS();
		}
	}
#endif

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012.
	eqs_new->AssembleRHS_PETSc(false);
	//		ScreenMessage("-> assemble a global matrix\n");
	eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY);
#endif
	ScreenMessage("-> impose Dirichlet BC\n");
	IncorporateBoundaryConditions();

	// ofstream Dum("rf_pcs.txt", ios::out); // WW
	// eqs_new->Write(Dum);   Dum.close();
	if (write_leqs)
	{
		std::string fname = FileName + "_" +
		                    convertProcessTypeToString(this->getProcessType()) +
		                    number2str(aktueller_zeitschritt) + "_" +
		                    number2str(iter_nlin) + "_leqs_bc.txt";
#if defined(NEW_EQS)
		std::ofstream Dum(fname.c_str(), ios::out);
		eqs_new->Write(Dum);
		Dum.close();
#elif defined(USE_PETSC)
		eqs_new->EQSV_Viewer(fname);
#endif
	}

#define nOUTPUT_EQS_BIN
#ifdef OUTPUT_EQS_BIN
	string fname = FileName + "_equation.bin";
	ofstream dum_bin(fname.c_str(), ios::out | ios::binary | ios::trunc);
	if (dum_bin.good())
	{
		eqs_new->Write_BIN(dum_bin);
		dum_bin.close();
		exit(1);
	}
#endif
//
//

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012.
	MPI_Barrier(MPI_COMM_WORLD);
//	eqs_new->AssembleRHS_PETSc();
// eqs_new->AssembleMatrixPETSc(MAT_FINAL_ASSEMBLY );
#endif
}
}

//--------------------------------------------------------------------
/*! \brief Assmble eqiations
     for all PDEs excluding deformation;

     24.11.2010. WW
 */
void CRFProcess::GlobalAssembly_std(bool is_quad, bool Check2D3D)
{
	long i;
	CElem* elem = NULL;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark())  // Marked for use
			continue;          // For OpenMP. WW

		elem->SetOrder(is_quad);
		fem->ConfigElement(elem, Check2D3D);
		fem->Assembly();
	}
}

/*************************************************************************
   GeoSys-Function:
   Task: Integration
   Programming:
   05/2009 WW Implementation
 **************************************************************************/
void CRFProcess::Integration(vector<double>& node_velue)
{
	//----------------------------------------------------------------------
	size_t k;
	long i;
	CElem* elem = NULL;
	bool Check2D3D;
	Check2D3D = false;
	double n_val[8];

	vector<double> buffer((long)node_velue.size());
	for (i = 0; i < (long)buffer.size(); i++)
		buffer[i] = 0.;

	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (!elem->GetMark()) continue;

		for (k = 0; k < elem->GetNodesNumber(false); k++)
			n_val[k] = node_velue[elem->GetNodeIndex(k)];

		elem->SetOrder(false);
		fem->ConfigElement(elem, Check2D3D);
		fem->FaceIntegration(n_val);

		for (k = 0; k < elem->GetNodesNumber(false); k++)
			buffer[elem->GetNodeIndex(k)] += n_val[k];
	}
	//----------------------------------------------------------------------
}

/*************************************************************************
   GeoSys-Function:
   Task: Calculate integration point velocity
   Programming:
   08/2005 WW Implementation
   last modified:
 **************************************************************************/
void CRFProcess::CalIntegrationPointValue()
{
	CElem* elem = NULL;
	cal_integration_point_value = false;
	continuum = 0;  // 15.02.2007/
	if (getProcessType() == FiniteElement::LIQUID_FLOW ||
	    getProcessType() == FiniteElement::RICHARDS_FLOW ||
	    getProcessType() == FiniteElement::MULTI_PHASE_FLOW ||
	    getProcessType() == FiniteElement::GROUNDWATER_FLOW ||
	    getProcessType() == FiniteElement::TWO_PHASE_FLOW ||
	    getProcessType() == FiniteElement::DEFORMATION_H2 ||
	    getProcessType() == FiniteElement::AIR_FLOW ||
	    getProcessType() == FiniteElement::PS_GLOBAL ||
	    getProcessType() == FiniteElement::DEFORMATION_FLOW  ||
	    getProcessType() == FiniteElement::TH_MONOLITHIC)
		cal_integration_point_value = true;
	if (!cal_integration_point_value) return;
	const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
	for (size_t i = 0; i < mesh_ele_vector_size; i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())
		{
			fem->ConfigElement(elem);
			fem->Config();
			fem->Cal_Velocity();
		}
	}
	if (getProcessType() == FiniteElement::TWO_PHASE_FLOW)
		cal_integration_point_value = false;
}

/*************************************************************************
   GeoSys-Function: CalGPVelocitiesfromFluidMomentum
   Task: Calculate gauss point velocities from fluid momentum solution
      extrapolate velocitiues from nodes to gauss points
   Programming:
   09/2009 SB Implementation

 **************************************************************************/
void CRFProcess::CalGPVelocitiesfromFluidMomentum()
{
	long i;
	MeshLib::CElem* elem = NULL;
	int i_ind[3];

	cout << "      CalGPVelocitiesfromFluidMomentum()" << endl;

	// Get fluid_momentum process
	CRFProcess* m_pcs_fm = PCSGet(FiniteElement::FLUID_MOMENTUM);

	//  check all possibilities for grid orientation (ccord_flag)
	int coordinateflag = this->m_msh->GetCoordinateFlag();
	// get index of velocity
	i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_X") + 1;
	if (coordinateflag == 11)
		// get index of velocity
		i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Y") + 1;
	if (coordinateflag == 12)
		// get index of velocity
		i_ind[0] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z") + 1;
	// get index of velocity
	i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Y") + 1;
	if (coordinateflag == 22)
		// get index of velocity
		i_ind[1] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z") + 1;
	// get index of velocity
	i_ind[2] = m_pcs_fm->GetNodeValueIndex("VELOCITY1_Z") + 1;

	if ((i_ind[0] < 0) || (i_ind[1] < 0) || (i_ind[2] < 0))
		cout << " Error - wrong index in Cal_GP_Velocity_FM " << endl;

	// Loop over all elements
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];  // get element
		if (elem->GetMark())          // Marked for use
		{
			fem->ConfigElement(elem);
			fem->Cal_GP_Velocity_FM(i_ind);
		}
	}  // end element loop
}

/*************************************************************************
   ROCKFLOW - Function: AllocateLocalMatrixMemory
   Task: As the function name
   Programming:
   01/2005 WW/OK Implementation
   06/2005 OK MSH
   last modified:
 **************************************************************************/
void CRFProcess::AllocateLocalMatrixMemory()
{
	long i;
	//----------------------------------------------------------------------
	int up_type = 0;
	if (!M_Process) up_type = 0;
	if (H_Process && M_Process)
	{
		if (type != 4 && type != 41)
			up_type = 1;
		else
		{
			if (type == 4) up_type = 3;
			if (type == 41) up_type = 4;
		}
	}
	if (!H_Process) up_type = 2;
	// SB for steady state element matrices in transport
	if (MASS_TRANSPORT_Process || T_Process) up_type = 5;
	//----------------------------------------------------------------------
	ElementMatrix* eleMatrix = NULL;
	CElem* elem = NULL;
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			eleMatrix = new ElementMatrix();
			eleMatrix->AllocateMemory(elem, up_type);
			Ele_Matrices.push_back(eleMatrix);
		}
	}
}

/**************************************************************************
   FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
   Task: set PCS boundary conditions
   Programing:
   02/2004 OK Implementation
   ????    WW  and time curve
   04/2005 OK MSH
   05/2005 OK conditional BCs
   04/2006 WW Changes due to the geometry object applied
   04/2006 OK Conditions by PCS coupling OK
   05/2006 WW Re-implement
   05/2006 WW DDC
   10/2007 WW Changes for the new classes of sparse matrix and linear solver
   last modification:
**************************************************************************/
void CRFProcess::IncorporateBoundaryConditions(const int rank, bool updateA,
                                               bool updateRHS, bool isResidual)
{
	static long i;
	static double bc_value, fac = 1.0, time_fac = 1.0;
	long bc_msh_node = -1;
#ifndef USE_PETSC
	long bc_eqs_index;
#endif
	long shift;
	int interp_method = 0;
	int curve, valid = 0;
	int ii, idx0 = -1;
	CBoundaryConditionNode* m_bc_node;  // WW
	CBoundaryCondition* m_bc;           // WW
	CFunction* m_fct = NULL;            // OK
	bool is_valid = false;              // OK
	// WW bool onExBoundary = false;                     //WX
	bool excavated = false;  // WX
#if defined(USE_PETSC)       // || defined(other parallel libs)//03~04.3012. WW
	vector<int> bc_eqs_id;
	vector<double> bc_eqs_value;
	vector<vector<int> > dof_node_id(this->GetPrimaryVNumber());
	vector<vector<double> > dof_node_value(this->GetPrimaryVNumber());
#else
	double* eqs_rhs = NULL;
#endif
//
#ifdef NEW_EQS
	Linear_EQS* eqs_p = NULL;
#endif
//------------------------------------------------------------WW

	// WW
	double Scaling = 1.0;
	bool quadr = false;
	if (type == 4 || type / 10 == 4)
	{
		fac = Scaling;
		quadr = true;
	}
	long begin = 0;
	long end = 0;
	long gindex = 0;
	if (rank == -1)
	{
		begin = 0;
		end = (long)bc_node_value.size();
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
// TODO
#ifdef NEW_EQS  // WW
		eqs_p = eqs_new;
		eqs_rhs = eqs_new->getRHS();  // 27.11.2007 WW
#else
			eqs_rhs = eqs->b;
#endif
#endif
	}
	size_t count_constrained_excluded = 0;

	for (i = begin; i < end; i++)
	{
		gindex = i;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03.3012. WW
		if (rank > -1) gindex = bc_node_value_in_dom[i];
#endif
		m_bc_node = bc_node_value[gindex];
		m_bc = bc_node[gindex];
		//
		// WX: check if bc is aktive, when Time_Controlled_Aktive for this bc is
		// defined
		if (m_bc->getTimeContrCurve() > 0)
			if (GetCurveValue(m_bc->getTimeContrCurve(), 0, aktuelle_zeit,
			                  &valid) < MKleinsteZahl)
				continue;

		// WX: 01.2011. for excavation bc, check if excavated and if on boundary
		if (m_bc->getExcav() > 0)
		{
			CNode* node;
			CElem* elem;
			// unsigned int counter;	//void warning
			// WW onExBoundary = true;                     //WX:01.2011
			excavated = false;
			node = m_msh->nod_vector[m_bc_node->geo_node_number];
			double const* node_coordinate(
			    node->getData());  // Coordinates(node_coordinate);
			                       // tmp_counter3++;
			                       // counter = 0;
			                       /*if((node_coordinate[ExcavDirection]-(GetCurveValue(ExcavCurve,0,aktuelle_zeit,&valid)-ExcavBeginCoordinate)<0.001
			                             &&(node_coordinate[ExcavDirection]-ExcavBeginCoordinate)>-0.001))*/
			// used with deactive subdomain
			if ((node_coordinate[ExcavDirection] -
			     (GetCurveValue(ExcavCurve, 0, aktuelle_zeit, &valid) -
			      ExcavBeginCoordinate)) < 0.001)
			{
				excavated = true;
				for (unsigned int j = 0;
				     j < node->getConnectedElementIDs().size();
				     j++)
				{
					elem = m_msh->ele_vector[node->getConnectedElementIDs()[j]];
					double const* tmp_ele_coor(elem->GetGravityCenter());
					// if(elem->GetPatchIndex()!=ExcavMaterialGroup){
					// if(elem->GetExcavState()==-1)
					if (elem->GetPatchIndex() !=
					    static_cast<size_t>(ExcavMaterialGroup))
						continue;
					else if (tmp_ele_coor[ExcavDirection] -
					             (GetCurveValue(ExcavCurve, 0, aktuelle_zeit,
					                            &valid) -
					              ExcavBeginCoordinate) <
					         0.001)
						// WW onExBoundary = false;
						// tmp_counter1++;
						break;
				}
			}
			// if(fabs(node_coordinate[ExcavDirection]-GetCurveValue(ExcavCurve,0,aktuelle_zeit,&valid)-ExcavBeginCoordinate)<0.001)
			//{
			// excavated = true;
			// onExBoundary = true;
			//}
			/*for(unsigned int j=0; j<node->connected_elements.size(); j++)
			   {
			   elem = m_msh->ele_vector[node->connected_elements[j]];
			   if(elem->GetExcavState()>0)
			      counter ++;
			   }
			   if ( counter!=0 && counter<node->connected_elements.size())
			   {
			   onExBoundary = true;
			   excavated = true;
			   }
			   else
			   if(counter==0)
			   excavated = false;
			   else
			   if(counter==node->connected_elements.size())
			   {
			   excavated = true;
			   onExBoundary = false;
			   }
			 */
		}

		if ((m_bc->getExcav() > 0) &&
		    !excavated)  // WX:01.2011. excav bc but is not excavated jet
			continue;
//
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
		bc_msh_node = m_bc_node->geo_node_number;
		// Check whether the node is in this subdomain
		if (quadr)
		{
			if (!m_msh->isNodeLocal(bc_msh_node)) continue;
		}
		else
		{
			if (!m_msh->isNodeLocal(bc_msh_node)) continue;
		}

		int dof_per_node = 0;
		if (m_msh->GetNodesNumber(false) == m_msh->GetNodesNumber(true))
		{
			dof_per_node = pcs_number_of_primary_nvals;
			shift = m_bc_node->msh_node_number / m_msh->GetNodesNumber(false);
		}
		else
		{
			if (bc_msh_node < static_cast<long>(m_msh->GetNodesNumber(false)))
				dof_per_node = pcs_number_of_primary_nvals;
			else
				dof_per_node = m_msh->GetMaxElementDim();
			shift = m_bc_node->msh_node_number / m_msh->GetNodesNumber(true);
		}

#else
		shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
		bc_msh_node = m_bc_node->geo_node_number;
#endif  // END: if defined(USE_PETSC) // || defined(other parallel libs
		//------------------------------------------------------------WW
		if (m_msh)  // OK
			//	    if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
			//          continue;
			time_fac = 1.0;
		if (bc_msh_node >= 0)
		{
			//................................................................
			// Constrain condition
			if (m_bc->has_constrain)
			{
				CBoundaryCondition* bc = m_bc_node->_bc;
				if (bc->constrain_var_id < 0)
				{
					bc->constrain_var_id =
					    GetNodeValueIndex(bc->constrain_var_name) + 1;
				}
				double val = GetNodeValue(bc_msh_node, bc->constrain_var_id);
				if (!FiniteElement::compare(val, bc->constrain_value,
				                            bc->constrain_operator))
				{
					count_constrained_excluded++;
					continue;  // skip this bc node
				}
			}

			//................................................................
			// Time dependencies - CURVE
			curve = m_bc_node->CurveIndex;
			if (curve > 0)
			{
				interp_method = m_bc->TimeInterpolation;

				if (curve > 10000000)  /// 16.08.2010. WW
					time_fac = GetCurveValue(curve - 10000000, interp_method,
					                         aktuelle_zeit, &valid);
				else
					time_fac = GetCurveValue(curve, interp_method,
					                         aktuelle_zeit, &valid);
				if (!valid) continue;
			}
			else
				time_fac = 1.0;
			//................................................................
			// Time dependencies - FCT
			if (!m_bc_node->fct_name.empty())
			{
				m_fct = FCTGet(m_bc_node->fct_name);
				if (m_fct) time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
				// if(!valid) continue;
				else
					cout << "Warning in "
					        "CRFProcess::IncorporateBoundaryConditions - no "
					        "FCT data" << endl;
			}
			//................................................................
			// Conditions
			if (m_bc_node->conditional)
			{
				int idx_1 = -1;               // 28.2.2007 WW
				for (ii = 0; ii < dof; ii++)  // 28.2.2007 WW

					if (convertPrimaryVariableToString(
					        m_bc->getProcessPrimaryVariable())
					        .find(pcs_primary_function_name[ii]) !=
					    string::npos)
					{
						idx_1 =
						    GetNodeValueIndex(pcs_primary_function_name[ii]) +
						    1;
						break;
					}
				bc_value =
				    time_fac * fac *
				    GetNodeValue(m_bc_node->msh_node_number_subst, idx_1);
			}
			else
				// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
				bc_value = time_fac * fac * m_bc_node->node_value;
			//----------------------------------------------------------------
			// MSH
			/// 16.08.2010. WW
			if (curve > 10000000 && fabs(time_fac) > DBL_EPSILON)
				bc_value =
				    bc_value / time_fac + time_fac;  /// bc_value +time_fac;

			if (m_bc->isPeriodic())  // JOD
				bc_value *= sin(2 * 3.14159 * aktuelle_zeit /
				                    m_bc->getPeriodeTimeLength() +
				                m_bc->getPeriodePhaseShift());

//----------------------------------------------------------------
#ifndef USE_PETSC
			if (rank > -1)
				bc_eqs_index = bc_msh_node;
			else
				// WW#
				bc_eqs_index =
					m_msh->nod_vector[bc_msh_node]->GetEquationIndex(quadr);
#endif
			//..............................................................
			if (FiniteElement::isNewtonKind(m_num->nls_method)
				|| type == 4 || type / 10 == 4)
			{  // Solution is in the manner of increment !
				idx0 = GetNodeValueIndex(convertPrimaryVariableToString(
				    m_bc->getProcessPrimaryVariable()));
				if ((type == 4 || type / 10 == 4) &&
				    m_bc_node->pcs_pv_name.find("DISPLACEMENT") != string::npos)
				{
					bc_value -=
					    GetNodeValue(m_bc_node->geo_node_number, idx0) +
					    GetNodeValue(m_bc_node->geo_node_number, idx0 + 1);
				}
				else
				{
					/// p_{n+1} = p_b,
					//  SetNodeValue(m_bc_node->geo_node_number, idx0++,
					//  bc_value);
					/// dp = u_b-u_n
					if (isResidual)
					{
						// SetNodeValue(m_bc_node->geo_node_number, idx0+1,
						// bc_value);
						bc_value -=
						    GetNodeValue(m_bc_node->geo_node_number, ++idx0);
						// bc_value *= -1;
					}
					else
						bc_value -=
						    GetNodeValue(m_bc_node->geo_node_number, ++idx0);
					// SetNodeValue(m_bc_node->geo_node_number, idx0, bc_value);
					// SetNodeValue(m_bc_node->geo_node_number, idx0+1,
					// bc_value);
					// bc_value = 0.;
				}
			}
#if !defined(USE_PETSC)  // && !defined(other parallel solver). //WW 04.2012. WW
			bc_eqs_index += shift;
#endif

				if (this->scaleUnknowns)
				{
					if (m_bc->getProcessPrimaryVariable() ==
					    FiniteElement::PRESSURE)
						bc_value *= vec_scale_dofs[0];
					else if (m_bc->getProcessPrimaryVariable() ==
					         FiniteElement::TEMPERATURE)
						bc_value *= vec_scale_dofs[1];
				}
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
				bc_eqs_id.push_back(static_cast<int>(
				    m_msh->nod_vector[bc_msh_node]->GetEquationIndex() *
				        dof_per_node +
				    shift));
				bc_eqs_value.push_back(bc_value);
				if (m_num->petsc_split_fields)
				{
					int dof_id = m_bc->getProcessPrimaryVariable() ==
					                     FiniteElement::PRESSURE
					                 ? 0
					                 : 1;  // TODO
					dof_node_id[dof_id].push_back(static_cast<int>(
					    m_msh->nod_vector[bc_msh_node]->GetEquationIndex()));
					dof_node_value[dof_id].push_back(bc_value);
				}

#elif defined(NEW_EQS)  // WW
			eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#else
			MXRandbed(bc_eqs_index, bc_value, eqs_rhs);
#endif
		}
	}
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
	if (m_num->petsc_split_fields)
	{
		if (updateRHS)
		{
			for (unsigned i = 0; i < dof_node_id.size(); i++)
			{
				VecGetSubVector(eqs_new->b, eqs_new->vec_isg[i],
				                &eqs_new->vec_subRHS[i]);
				//					for (unsigned j=0; j<dof_node_id[i].size();
				// j++)
				//						ScreenMessage2("-> bc: dof=%d, node
				//id=%d,
				// val=%g\n", i, dof_node_id[i][j], dof_node_value[i][j]);
				//					VecView(eqs_new->vec_subRHS[i],
				// PETSC_VIEWER_STDOUT_WORLD);
				int nrow = dof_node_id[i].size();
				if (nrow > 0)
					VecSetValues(eqs_new->vec_subRHS[i], nrow,
					             &dof_node_id[i][0], &dof_node_value[i][0],
					             INSERT_VALUES);
				VecAssemblyBegin(eqs_new->vec_subRHS[i]);
				VecAssemblyEnd(eqs_new->vec_subRHS[i]);
				//					ScreenMessage2("-> after bc\n");
				//					VecView(eqs_new->vec_subRHS[i],
				// PETSC_VIEWER_STDOUT_WORLD);
				VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[i],
				                    &eqs_new->vec_subRHS[i]);
			}
			eqs_new->AssembleRHS_PETSc();
		}

		if (updateA)
		{
			for (unsigned i = 0; i < dof_node_id.size(); i++)
			{
				int nrow = dof_node_id[i].size();
				for (unsigned j = 0; j < dof_node_id.size(); j++)
				{
					PetscScalar diag_val = 1.0;
					if (i != j) diag_val = .0;
					if (nrow > 0)
						MatZeroRows(
						    eqs_new->vec_subA[i * dof_node_id.size() + j], nrow,
						    &dof_node_id[i][0], diag_val, PETSC_NULL,
						    PETSC_NULL);
					else
						MatZeroRows(
						    eqs_new->vec_subA[i * dof_node_id.size() + j], 0,
						    PETSC_NULL, diag_val, PETSC_NULL, PETSC_NULL);
				}
			}
			eqs_new->AssembleMatrixPETSc();
		}
	}
	else
	{
		int nbc = static_cast<int>(bc_eqs_id.size());
		if (updateRHS)
		{
			if (nbc > 0)
			{
				eqs_new->setArrayValues(0, nbc, &bc_eqs_id[0], &bc_eqs_value[0],
				                        INSERT_VALUES);
				eqs_new->setArrayValues(1, nbc, &bc_eqs_id[0], &bc_eqs_value[0],
				                        INSERT_VALUES);
			}

#ifdef petsc_zero_row_test
			// We have do the following collection because MatZeroR must be
			// called by all processes
			const int mpi_size = eqs_new->getMPI_Size();
			vector<int> r_cnt(mpi_size);
			vector<int> r_disp(mpi_size);
			vector<int> r_vec(mpi_size);
			int k;
			for (k = 0; k < mpi_size; k++)
			{
				r_cnt[k] = 1;
				r_disp[k] = k;
			}
			// Get nbc
			MPI_Allgatherv(&nbc, 1, MPI_INT, &r_vec[0], &r_cnt[0], &r_disp[0],
			               MPI_INT, PETSC_COMM_WORLD);
			int v_disp = 0;
			for (k = 0; k < mpi_size; k++)
			{
				r_disp[k] = v_disp;
				v_disp += r_vec[k];
			}
			r_cnt.resize(v_disp);
			r_vec.resize(v_disp);
			for (k = 0; k < v_disp; k++)
			{
				r_cnt[k] = 0;
			}
			const int v_shift = r_disp[eqs_new->getMPI_Rank()];
			for (k = 0; k < nbc; k++)
			{
				r_cnt[v_shift + k] = bc_eqs_id[k];
			}
			MPI_Allreduce(&r_cnt[0], &r_vec[0], v_disp, MPI_INT, MPI_SUM,
			              PETSC_COMM_WORLD);

			MPI_Barrier(MPI_COMM_WORLD);
			eqs_new->zeroRows_in_Matrix(v_disp, &r_vec[0]);
#endif  //  petsc_zero_row_test
			eqs_new->AssembleUnkowns_PETSc();
			eqs_new->AssembleRHS_PETSc();
		}

		// TEST
		// PetscViewer viewer;
		// eqs_new->EQSV_Viewer(FileName, viewer);

		if (updateA)
		{
			eqs_new->zeroRows_in_Matrix(nbc, &bc_eqs_id[0]);
			eqs_new->AssembleMatrixPETSc();
		}
	}

#endif

	if (count_constrained_excluded > 0)
		std::cout
		    << "-> " << count_constrained_excluded
		    << " nodes are excluded from BC because of constrained conditions"
		    << "\n";
}

/**************************************************************************
   FEMLib-Method: CRFProcess::IncorporateBoundaryConditions
   Task: set PCS boundary conditions for FLUID_MOMENTUM depending on axis
   Programing:
   01/2007 PCH Implementation
   last modification:
**************************************************************************/
void CRFProcess::IncorporateBoundaryConditions(const int rank, const int axis)
{
	static long i;
	static double bc_value, fac = 1.0, time_fac = 1.0;
	long bc_msh_node;
	long bc_eqs_index, shift;
	int interp_method = 0;
	int curve, valid = 0;
	int idx0, idx1;
	CBoundaryConditionNode* m_bc_node;  // WW
	CBoundaryCondition* m_bc;           // WW

	CFunction* m_fct = NULL;  // OK
	bool is_valid = false;    // OK
#if defined(USE_PETSC)        // || defined(other parallel libs)//03~04.3012. WW
	vector<int> bc_eqs_id;
	vector<double> bc_eqs_value;
#endif
#ifdef NEW_EQS
	Linear_EQS* eqs_p = NULL;
#endif
	//------------------------------------------------------------WW
	// WW
	double Scaling = 1.0;
	if (type == 4 || type == 41) fac = Scaling;

	long begin = 0;
	long end = 0;
	long gindex = 0;

	// WW CBoundaryConditionsGroup *m_bc_group = NULL;
	//  m_bc_group =
	//  BCGetGroup(this->_pcs_type_name,this->pcs_primary_function_name[axis]);
	// TF
	// WW m_bc_group =
	// BCGetGroup(convertProcessTypeToString(this->getProcessType()),this->pcs_primary_function_name[axis]);

	if (rank == -1)
	{
		begin = 0;
		end = (long)bc_node_value.size();
#ifdef NEW_EQS           // WW
		eqs_p = eqs_new;
#endif
	}
	for (i = begin; i < end; i++)
	{
		gindex = i;
		m_bc_node = bc_node_value[gindex];

		// PCH
		if (axis == 0 &&
		    m_bc_node->pcs_pv_name.find("VELOCITY1_X") != string::npos)
		{
			m_bc = bc_node[gindex];

#if defined(USE_PETSC)  // ||defined(other parallel libs)//03~04.3012. WW
			bc_msh_node = m_bc_node->geo_node_number;

#else
			shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
			bc_msh_node = m_bc_node->geo_node_number;
#endif
			//------------------------------------------------------------WW
			if (m_msh)  // OK
				//			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
				//				continue;
				time_fac = 1.0;
			if (bc_msh_node >= 0)
			{
				//................................................................
				// Time dependencies - CURVE
				curve = m_bc_node->CurveIndex;
				if (curve > 0)
				{
					time_fac = GetCurveValue(curve, interp_method,
					                         aktuelle_zeit, &valid);
					if (!valid) continue;
				}
				else
					time_fac = 1.0;
				//................................................................
				// Time dependencies - FCT
				if (m_bc_node->fct_name.length() > 0)
				{
					m_fct = FCTGet(m_bc_node->fct_name);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					// if(!valid) continue;
					else
						cout << "Warning in "
						        "CRFProcess::IncorporateBoundaryConditions - "
						        "no FCT data" << endl;
				}
				//................................................................
				// Conditions
				if (m_bc_node->conditional)
					bc_value =
					    time_fac * fac *
					    GetNodeValue(
					        m_bc_node->msh_node_number_subst,
					        // WW  bc_value = time_fac*fac*
					        // GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1);
					        // // YD-----TEST---
					        GetNodeValueIndex(pcs_primary_function_name[0]) +
					            1);
				else
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				//----------------------------------------------------------------
				// MSH
				if (rank > -1)
					bc_eqs_index = bc_msh_node;
				else
					// WW#
					bc_eqs_index =
					    m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
				//..............................................................
				// NEWTON WW
				if (FiniteElement::isNewtonKind(
				        m_num->nls_method)  //_name.find("NEWTON")!=string::npos
				    ||
				    type == 4 ||
				    type == 41)  // Solution is in the manner of increment !
				{
					idx0 = GetNodeValueIndex(convertPrimaryVariableToString(
					    m_bc->getProcessPrimaryVariable()));
					if (type == 4 || type == 41)
					{
						// 30.08.2010. WW
						if (m_bc_node->pcs_pv_name.find("PRESSURE") ==
						    string::npos)
						{
							idx1 = idx0 + 1;
							bc_value -=
							    GetNodeValue(m_bc_node->geo_node_number, idx0) +
							    GetNodeValue(m_bc_node->geo_node_number, idx1);
						}
					}
					else
						bc_value =
						    bc_value -
						    GetNodeValue(m_bc_node->geo_node_number, idx0);
				}
				//----------------------------------------------------------------
				bc_eqs_index += shift;
				if ((int)continuum_vector.size() > 1)
					// YD/WW
					if (m_bc_node->pcs_pv_name.find(
					        pcs_primary_function_name[continuum]) ==
					    string::npos)
						continue;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif NEW_EQS  // WW
				eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#endif
			}
		}
		// PCH
		else if (axis == 1 &&
		         m_bc_node->pcs_pv_name.find("VELOCITY1_Y") != string::npos)
		{
			m_bc = bc_node[gindex];
			shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
			bc_msh_node = m_bc_node->geo_node_number;
			//------------------------------------------------------------WW
			if (m_msh)  // OK
				//			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
				//				continue;
				time_fac = 1.0;
			if (bc_msh_node >= 0)
			{
				//................................................................
				// Time dependencies - CURVE
				curve = m_bc_node->CurveIndex;
				if (curve > 0)
				{
					time_fac = GetCurveValue(curve, interp_method,
					                         aktuelle_zeit, &valid);
					if (!valid) continue;
				}
				else
					time_fac = 1.0;
				//................................................................
				// Time dependencies - FCT
				if (m_bc_node->fct_name.length() > 0)
				{
					m_fct = FCTGet(m_bc_node->fct_name);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					// if(!valid) continue;
					else
						cout << "Warning in "
						        "CRFProcess::IncorporateBoundaryConditions - "
						        "no FCT data" << endl;
				}
				//................................................................
				// Conditions
				if (m_bc_node->conditional)
					bc_value =
					    time_fac * fac *
					    GetNodeValue(
					        m_bc_node->msh_node_number_subst,
					        // WW  bc_value = time_fac*fac*
					        // GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1);
					        // // YD-----TEST---
					        GetNodeValueIndex(pcs_primary_function_name[0]) +
					            1);
				else
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				//----------------------------------------------------------------
				// MSH
				if (rank > -1)
					bc_eqs_index = bc_msh_node;
				else
					// WW#
					bc_eqs_index =
					    m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
				//..............................................................
				// NEWTON WW
				if (FiniteElement::isNewtonKind(
				        m_num->nls_method)  //_name.find("NEWTON")!=string::npos
				    ||
				    type == 4 ||
				    type == 41)  // Solution is in the manner of increment !
				{
					idx0 = GetNodeValueIndex(
					    convertPrimaryVariableToString(
					        m_bc->getProcessPrimaryVariable()).c_str());
					if (type == 4 || type == 41)
					{
						idx1 = idx0 + 1;
						bc_value -=
						    GetNodeValue(m_bc_node->geo_node_number, idx0) +
						    GetNodeValue(m_bc_node->geo_node_number, idx1);
					}
					else
						bc_value =
						    bc_value -
						    GetNodeValue(m_bc_node->geo_node_number, idx0);
				}
				//----------------------------------------------------------------
				bc_eqs_index += shift;
				if ((int)continuum_vector.size() > 1)
					// YD/WW
					if (m_bc_node->pcs_pv_name.find(
					        pcs_primary_function_name[continuum]) ==
					    string::npos)
						continue;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif NEW_EQS  // WW
				eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#endif
			}
		}
		// PCH
		else if (axis == 2 &&
		         m_bc_node->pcs_pv_name.find("VELOCITY1_Z") != string::npos)
		{
			m_bc = bc_node[gindex];
			shift = m_bc_node->msh_node_number - m_bc_node->geo_node_number;
			bc_msh_node = m_bc_node->geo_node_number;
			//------------------------------------------------------------WW
			if (m_msh)  // OK
				//			if(!m_msh->nod_vector[bc_msh_node]->GetMark()) //WW
				//				continue;
				time_fac = 1.0;
			if (bc_msh_node >= 0)
			{
				//................................................................
				// Time dependencies - CURVE
				curve = m_bc_node->CurveIndex;
				if (curve > 0)
				{
					time_fac = GetCurveValue(curve, interp_method,
					                         aktuelle_zeit, &valid);
					if (!valid) continue;
				}
				else
					time_fac = 1.0;
				//................................................................
				// Time dependencies - FCT
				if (m_bc_node->fct_name.length() > 0)
				{
					m_fct = FCTGet(m_bc_node->fct_name);
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					// if(!valid) continue;
					else
						cout << "Warning in "
						        "CRFProcess::IncorporateBoundaryConditions - "
						        "no FCT data" << endl;
				}
				//................................................................
				// Conditions
				if (m_bc_node->conditional)
					bc_value =
					    time_fac * fac *
					    GetNodeValue(
					        m_bc_node->msh_node_number_subst,
					        // WW  bc_value = time_fac*fac*
					        // GetNodeVal(bc_msh_node+1,GetNODValueIndex(pcs_primary_function_name[0])+1);
					        // // YD-----TEST---
					        GetNodeValueIndex(pcs_primary_function_name[0]) +
					            1);
				else
					// time_fac*fac*PCSGetNODValue(bc_msh_node,"PRESSURE1",0);
					bc_value = time_fac * fac * m_bc_node->node_value;
				//----------------------------------------------------------------
				// MSH
				if (rank > -1)
					bc_eqs_index = bc_msh_node;
				else
					// WW#
					bc_eqs_index =
					    m_msh->nod_vector[bc_msh_node]->GetEquationIndex();
				//..............................................................
				// NEWTON WW
				if (FiniteElement::isNewtonKind(
				        m_num
				            ->nls_method)  //_name.find("NEWTON")!=std::string::npos
				    ||
				    type == 4 ||
				    type == 41)  // Solution is in the manner of increment !
				{
					idx0 = GetNodeValueIndex(
					    convertPrimaryVariableToString(
					        m_bc->getProcessPrimaryVariable()).c_str());
					if (type == 4 || type == 41)
					{
						idx1 = idx0 + 1;
						bc_value -=
						    GetNodeValue(m_bc_node->geo_node_number, idx0) +
						    GetNodeValue(m_bc_node->geo_node_number, idx1);
					}
					else
						bc_value =
						    bc_value -
						    GetNodeValue(m_bc_node->geo_node_number, idx0);
				}
				//----------------------------------------------------------------
				bc_eqs_index += shift;
				if ((int)continuum_vector.size() > 1)
					// YD/WW
					if (m_bc_node->pcs_pv_name.find(
					        pcs_primary_function_name[continuum]) ==
					    std::string::npos)
						continue;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
// TODO
#elif NEW_EQS  // WW
				eqs_p->SetKnownX_i(bc_eqs_index, bc_value);
#endif
			}
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: PCS source terms into EQS
   Programing:
   04/2004 OK Implementation
   08/2004 WW Extension for monolithic PCS and time curve
   last modification:
   02/2005 MB River Condition and CriticalDepth
   05/2005 WW Dynamic problems
   07/2005 WW Changes due to the geometry object applied
   03/2006 WW Re-arrange
   04/2006 OK CPL
   05/2006 WW DDC
   08/2006 YD FCT use
**************************************************************************/
void CRFProcess::IncorporateSourceTerms(const int rank)
{
	double value = 0, fac = 1.0, time_fac;
	int interp_method = 0;
	int curve, valid = 0;
	long msh_node, shift;
	MshElemType::type EleType;  // ii
	double q_face = 0.0;
	CElem* elem = NULL;
	CElem* face = NULL;
	ElementValue* gp_ele = NULL;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
	vector<int> st_eqs_id;
	vector<double> st_eqs_value;
	vector<vector<int> > dof_node_id(this->GetPrimaryVNumber());
	vector<vector<double> > dof_node_value(this->GetPrimaryVNumber());
#else
	double* eqs_rhs = NULL;
	long bc_eqs_index = -1;
	int dim_space = 0;  // kg44 better define here and not in a loop!
#endif
	double vel[3];
	bool is_valid;            // YD
	CFunction* m_fct = NULL;  // YD
	long i;                   //, group_vector_length;

	double Scaling = 1.0;
	bool quadr = false;
	if (type == 4 || type / 10 == 4)
	{
		fac = Scaling;
		quadr = true;
	}

	CNodeValue* cnodev = NULL;
	CSourceTerm* m_st = NULL;
	//
	long begin = 0;
	long end = 0;
	long gindex = 0;

	//====================================================================
	// Look for active boundary elements if constrain is given
	for (unsigned i = 0; i < st_vector.size(); i++)
	{
		CSourceTerm* st = st_vector[i];
		if (!st->has_constrain) continue;

		CSourceTermGroup m_st_group;
		m_st_group.pcs_type_name =
		    FiniteElement::convertProcessTypeToString(st->getProcessType());
		m_st_group.pcs_pv_name = FiniteElement::convertPrimaryVariableToString(
		    st->getProcessPrimaryVariable());
		int idx = GetNodeValueIndex(m_st_group.pcs_pv_name) / 2;
		m_st_group.Set(this, Shift[idx]);
	}

	for (size_t is = 0; is < st_node_value.size(); is++)
	{
		m_st = NULL;
		if (is < st_vector.size()) m_st = st_vector[is];
		if (rank == -1)
		{
			begin = 0;
			end = (long)st_node_value[is].size();
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
#ifdef NEW_EQS           // WW
			eqs_rhs = eqs_new->getRHS();  // 27.11.2007 WW
#else
				eqs_rhs = eqs->b;
#endif
#endif
		}
		std::vector<bool> active_elements;

		// constrain
		if (m_st->has_constrain && m_st->getSTType() == FiniteElement::NEUMANN)
		{
			std::cout << "-> update constrained ST " << is << "\n";
			// get distributed values
			DistributionData distData;
			setDistributionData(m_st, distData);
			std::vector<long> nodes_vector(st_node_value[is].size());
			for (size_t ii = 0; ii < st_node_value[is].size(); ii++)
			{
				nodes_vector[ii] = st_node_value[is][ii]->msh_node_number;
			}
			std::vector<double> node_value(nodes_vector.size());
			setDistribution(distData, *m_msh, nodes_vector, node_value);

			// boundary integration
			if (m_st->constrain_var_id < 0)
				m_st->constrain_var_id =
				    GetNodeValueIndex(m_st->constrain_var_name) + 1;
			active_elements.resize(m_st->st_boundary_elements.size(), true);
			size_t count_active = 0;
			if (m_st->getGeoType() == GEOLIB::POINT)
			{
				double val =
				    GetNodeValue(m_st->geo_node_number, m_st->constrain_var_id);
				if (!FiniteElement::compare(val, m_st->constrain_value,
				                            m_st->constrain_operator))
				{
					node_value[0] = .0;
					active_elements[0] = false;
				}
				else
				{
					count_active++;
				}
				if (m_st->is_transfer_bc)
				{
					node_value[0] *= m_st->transfer_h_values[0];
				}
			}
			else
			{
				for (size_t in = 0; in < m_st->st_boundary_elements.size();
				     in++)
				{
					MeshLib::CElem* face = m_st->st_boundary_elements[in];
					const unsigned nen =
					    face->GetNodesNumber(m_msh->getOrder());
					double avg_val = .0;
					for (unsigned k = 0; k < nen; k++)
						avg_val += GetNodeValue(face->GetNode(k)->GetIndex(),
						                        m_st->constrain_var_id);
					avg_val /= (double)nen;
					active_elements[in] =
					    FiniteElement::compare(avg_val, m_st->constrain_value,
					                           m_st->constrain_operator);
					if (active_elements[in]) count_active++;
				}
				if (m_msh->GetMaxElementDim() == 2 &&
				    m_st->getGeoType() == GEOLIB::POLYLINE)
					m_st->EdgeIntegration(m_msh, nodes_vector, node_value,
					                      &active_elements);
				else if (m_msh->GetMaxElementDim() == 3 &&
				         m_st->getGeoType() == GEOLIB::SURFACE)
					m_st->FaceIntegration(m_msh, nodes_vector, node_value,
					                      &active_elements);
			}

			// update ST values
			for (size_t ii = 0; ii < st_node_value[is].size(); ii++)
			{
				st_node_value[is][ii]->node_value = node_value[ii];
			}
			std::cout << "-> " << count_active
			          << " nodes/elements are active in total "
			          << active_elements.size() << " nodes/elements"
			          << "\n";
		}

		// exchange condition needs to update a coefficient matrix
		if (m_st->is_transfer_bc)
		{
			// only Neumann BC
			if (m_st->getSTType() != FiniteElement::NEUMANN) continue;

			double mass[100] = {};
			if (m_st->getGeoType() == GEOLIB::POINT)
			{
				if (m_st->has_constrain && !active_elements[0]) continue;
				cnodev = st_node_value[is][0];
				const int k_eqs_id = m_msh->nod_vector[cnodev->geo_node_number]
				                         ->GetEquationIndex();
#if defined(USE_PETSC)
				eqs_new->addMatrixEntry(k_eqs_id, k_eqs_id,
				                        m_st->transfer_h_values[0]);
#elif defined(NEW_EQS)
				(*eqs_new->getA())(k_eqs_id, k_eqs_id) += m_st->transfer_h_values[0];
#endif
			}
			else if (m_st->getGeoType() == GEOLIB::SURFACE ||
			         m_st->getGeoType() == GEOLIB::POLYLINE)
			{
				for (size_t in = 0; in < m_st->st_boundary_elements.size();
				     in++)
				{
					if (m_st->has_constrain && !active_elements[in]) continue;
					MeshLib::CElem* face = m_st->st_boundary_elements[in];
					const unsigned nen = face->GetNodesNumber(false);
					fem->setOrder(m_msh->getOrder() + 1);
					fem->ConfigElement(face, true);
					for (unsigned k = 0; k < nen; k++)
						for (unsigned l = 0; l < nen; l++)
							mass[k * nen + l] = .0;
					fem->CalcFaceMass(mass);
					const double h =
					    m_st->transfer_h_values[face->GetPatchIndex()];
					for (unsigned k = 0; k < nen; k++)
					{
						const int k_eqs_id =
						    face->GetNode(k)->GetEquationIndex();
						for (unsigned l = 0; l < nen; l++)
						{
							const int l_eqs_id =
							    face->GetNode(l)->GetEquationIndex();
#if defined(USE_PETSC)
							eqs_new->addMatrixEntry(k_eqs_id, l_eqs_id,
							                        mass[k * nen + l] * h);
#elif defined(NEW_EQS)
							(*eqs_new->getA())(k_eqs_id, l_eqs_id) +=
							    mass[k * nen + l] * h;
#endif
						}
					}
				}
			}
		}

		//====================================================================
		// Add ST to RHS
		for (long i = begin; i < end; i++)
		{
			gindex = i;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03.3012. WW
			if (rank > -1) gindex = st_node_value_in_dom[i];
#endif

			cnodev = st_node_value[is][gindex];

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
			msh_node = cnodev->geo_node_number;
			// Check whether the node is in this subdomain
			if (quadr)
			{
				if (!m_msh->isNodeLocal(msh_node)) continue;
			}
			else
			{
				if (!m_msh->isNodeLocal(msh_node)) continue;
			}

			int dof_per_node = 0;
			if (m_msh->GetNodesNumber(false) == m_msh->GetNodesNumber(true))
			{
				dof_per_node = pcs_number_of_primary_nvals;
				shift = cnodev->msh_node_number / m_msh->GetNodesNumber(false);
			}
			else
			{
				if (msh_node < static_cast<long>(m_msh->GetNodesNumber(false)))
					dof_per_node = pcs_number_of_primary_nvals;
				else
					dof_per_node = m_msh->GetMaxElementDim();
				shift = cnodev->msh_node_number / m_msh->GetNodesNumber(true);
			}

#else
			shift = cnodev->msh_node_number - cnodev->geo_node_number;
			msh_node = cnodev->msh_node_number;
			msh_node -= shift;
#endif
			value = cnodev->node_value;
			//--------------------------------------------------------------------
			// Tests
			if (msh_node < 0) continue;
			//--------------------------------------------------------------------
			// CPL
			// if(m_st->_pcs_type_name_cond.size()>0) continue; // this is a CPL
			// source term, JOD removed
			//--------------------------------------------------------------------
			// system dependent YD
			if (cnodev->getProcessDistributionType() ==
			    FiniteElement::SYSTEM_DEPENDENT)
			{
				long no_st_ele = (long)m_st->element_st_vector.size();
				for (long i_st = 0; i_st < no_st_ele; i_st++)
				{
					long ele_index = m_st->element_st_vector[i_st];
					elem = m_msh->ele_vector[ele_index];
					if (elem->GetMark())
					{
						fem->ConfigElement(elem);
						fem->Cal_Velocity();
					}
					gp_ele = ele_gp_value[ele_index];
					gp_ele->GetEleVelocity(vel);
					EleType = elem->GetElementType();
					if (EleType == MshElemType::LINE)  // Line
						cnodev->node_value += vel[0];
					// Traingle & Qua
					if (EleType == MshElemType::TRIANGLE ||
					    EleType == MshElemType::QUAD)
					{
						for (size_t i_face = 0;
						     i_face < m_msh->face_vector.size();
						     i_face++)
						{
							face = m_msh->face_vector[i_face];
							if ((size_t)m_st->element_st_vector[i_st] ==
							    face->GetOwner()->GetIndex())
								//
								q_face = PointProduction(
								             vel, m_msh->face_normal[i_face]) *
								         face->GetVolume();
							// for(i_node)
						}
						cnodev->node_value = +q_face / 2;
					}
				}
				//--------------------------------------------------------------------
				// MB
				// if(m_st->conditional && !m_st->river)
				//{

				GetNODValue(value, cnodev, m_st);
			}  // st_node.size()>0&&(long)st_node.size()>i
			//----------------------------------------------------------------------------------------
			//--------------------------------------------------------------------
			// Please do not move the this section
			curve = cnodev->CurveIndex;
			if (curve > 0)
			{
				// Reading Time interpolation method; BG
				if (m_st != NULL)  // in some cases the m_st is not defined ->
					               // interp_method is not changed for this
					               // cases
					if (interp_method != m_st->TimeInterpolation)
						interp_method = m_st->TimeInterpolation;

				time_fac =
				    GetCurveValue(curve, interp_method, aktuelle_zeit, &valid);
				// cout << "step: " << this->Tim->step_current << " Time: " <<
				// aktuelle_zeit << " Laenge: " << this->Tim->this_stepsize << "
				// Beginn: " << this->Tim->time_start << " Ende " <<
				// this->Tim->time_end << " Faktor: " << time_fac << endl;
				if (!valid)
				{
					cout << "\n!!! Time dependent curve is not found. Results "
					        "are not guaranteed " << endl;
					cout << " in void CRFProcess::IncorporateSourceTerms(const "
					        "double Scaling)" << endl;
					time_fac = 1.0;
				}
			}
			else
				time_fac = 1.0;

			// Time dependencies
			if (m_st)
			{
				if (m_msh && !m_msh->geo_name.empty() &&
				    m_msh->geo_name.find("LOCAL") != string::npos)
				{
					if (m_st->getFunctionName().length() > 0)
					{
						m_fct = FCTGet(pcs_number);
						if (m_fct)
							time_fac = m_fct->GetValue(
							    aktuelle_zeit,
							    &is_valid,
							    m_st->getFunctionMethod());
						else
							cout << "Warning in "
							        "CRFProcess::IncorporateSourceTerms - no "
							        "FCT data" << endl;
					}
				}
				else if (!m_st->getFunctionName().empty())
				{
					m_fct = FCTGet(m_st->getFunctionName());
					if (m_fct)
						time_fac = m_fct->GetValue(aktuelle_zeit, &is_valid);
					else
						cout << "Warning in CRFProcess::IncorporateSourceTerms "
						        "- no FCT data" << endl;
				}
			}
			//----------------------------------------------------------------------------------------
			value *= time_fac * fac;
			if (scaleEQS)
			{
				if (m_st->getProcessPrimaryVariable() == FiniteElement::PRESSURE)
				{
					value *= vec_scale_eqs[0];
				}
				else if (m_st->getProcessPrimaryVariable() == FiniteElement::TEMPERATURE)
				{
					value *= vec_scale_eqs[1];
				}
			}
//------------------------------------------------------------------
// EQS->RHS
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW

			st_eqs_id.push_back(static_cast<int>(
			    m_msh->nod_vector[msh_node]->GetEquationIndex() * dof_per_node +
			    shift));
			st_eqs_value.push_back(value);
			if (m_num->petsc_split_fields)
			{
				int dof_id =
				    m_st->getProcessPrimaryVariable() == FiniteElement::PRESSURE
				        ? 0
				        : 1;  // TODO
				dof_node_id[dof_id].push_back(static_cast<int>(
				    m_msh->nod_vector[msh_node]->GetEquationIndex()));
				dof_node_value[dof_id].push_back(value);
			}

#else
			if (rank > -1)
				bc_eqs_index = msh_node + shift;
			else
				bc_eqs_index =
					m_msh->nod_vector[msh_node]->GetEquationIndex(quadr) + shift;
			eqs_rhs[bc_eqs_index] += value;
#endif
		}
	}

	//====================================================================

	// if coupling to GEMS
	// exist----------------------------------------------------
	// HS, added 11.2008
	// KG44 03/03/2010 modified to hopefully soon work with parallel solvers
	long gem_node_index = -1, glocalindex = -1;
	if (flag_couple_GEMS == 1 && aktueller_zeitschritt > 1)
	{
		begin = 0;
		if (rank == -1)  // serial version and also Version for PETSC!!

			end = (long)Water_ST_vec.size();
		else  // parallel version
		{
			end = 0;
			if (rank_stgem_node_value_in_dom.size() > 0)
				end = rank_stgem_node_value_in_dom[0];
		}
		// only when switch is on and not in the first time step
		// loop over the Water_ST_vec vector,
		// add the excess water to the right-hand-side of the equation
		for (long i = begin; i < end; i++)
		{
			if (rank > -1)  // parallel version: stgem_node_value_in_dom and
			                // stgem_local_index_in_dom contain only values for
			                // the corresponding domain  == rank
			{
				//				cout << "rank " << rank ;
				gindex = stgem_node_value_in_dom[i];  // contains indexes to
				                                      // water-st_vec
				//				cout << " gindex " << gindex << " i " << i <<
				//endl
				//;
				// contains index to node
				glocalindex = stgem_local_index_in_dom[i];
//				cout << " gem_node_index " << gem_node_index << "\n";
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW

				st_eqs_id.push_back(static_cast<int>(
				    m_msh->nod_vector[glocalindex]->GetEquationIndex()));
				st_eqs_value.push_back(Water_ST_vec[gindex].water_st_value);
#else
				eqs_rhs[glocalindex] += Water_ST_vec[gindex].water_st_value;
#endif
			}
			else  // serial version
			{
				gem_node_index = Water_ST_vec[i].index_node;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
				st_eqs_id.push_back(static_cast<int>(
				    m_msh->nod_vector[gem_node_index]->GetEquationIndex()));
				st_eqs_value.push_back(Water_ST_vec[i].water_st_value);

#else
				eqs_rhs[gem_node_index] += Water_ST_vec[i].water_st_value;
#endif
			}
		}
		// after finished adding to RHS, clear the vector
		Water_ST_vec.clear();
		if (rank > -1)
		{
			stgem_node_value_in_dom.clear();
			stgem_local_index_in_dom.clear();
			rank_stgem_node_value_in_dom.clear();
		}
	}
#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
	if (m_num->petsc_split_fields)
	{
		//				for (unsigned i=0; i<eqs_new->vec_subRHS.size(); i++)
		//					VecGetSubVector(eqs_new->b, eqs_new->vec_isg[i],
		//&eqs_new->vec_subRHS[i]);
		for (unsigned i = 0; i < dof_node_id.size(); i++)
		{
			VecGetSubVector(eqs_new->b, eqs_new->vec_isg[i],
			                &eqs_new->vec_subRHS[i]);
			if (!st_eqs_id.empty())
			{
				int nrow = dof_node_id[i].size();
				if (nrow > 0)
					VecSetValues(eqs_new->vec_subRHS[i], nrow,
					             &dof_node_id[i][0], &dof_node_value[i][0],
					             ADD_VALUES);
			}
			VecAssemblyBegin(eqs_new->vec_subRHS[i]);
			VecAssemblyEnd(eqs_new->vec_subRHS[i]);
			VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[i],
			                    &eqs_new->vec_subRHS[i]);
		}
		//				for (size_t i=0; i<eqs_new->vec_subRHS.size(); i++)
		//					VecRestoreSubVector(eqs_new->b, eqs_new->vec_isg[i],
		//&eqs_new->vec_subRHS[i]);
		// eqs_new->AssembleRHS_PETSc();
	}
	else
	{
		if (st_eqs_id.size() > 0)
			eqs_new->setArrayValues(1, static_cast<int>(st_eqs_id.size()),
			                        &st_eqs_id[0], &st_eqs_value[0]);
		// eqs_new->AssembleRHS_PETSc();
	}
#endif
}


/*-------------------------------------------------------------------------
   ROCKFLOW - Function: PCSRestart
   Task: Insert process to list
   Programming:
   06/2003 OK Implementation
   11/2004 OK file_name_base
   last modified:
   -------------------------------------------------------------------------*/
void PCSRestart()
{
	/*OK411
	   int j;
	   CRFProcess *m_pcs = NULL;
	   int nidx0,nidx1;
	   int i;
	   int no_processes =(int)pcs_vector.size();
	   if(no_processes==0)
	    return; //OK41
	   int ok = 0;
	   //----------------------------------------------------------------------
	   string file_name_base = pcs_vector[0]->file_name_base;
	   //OK  ok = ReadRFRRestartData(file_name_base);
	   if(ok==0){
	   cout << "RFR: no restart data" << endl;
	   return;
	   }
	   //----------------------------------------------------------------------
	   for(i=0;i<no_processes;i++){
	   m_pcs = pcs_vector[i];
	   for(j=0;j<m_pcs->GetPrimaryVNumber();j++) {
	   // timelevel=0;
	   nidx0 = m_pcs->GetNodeValueIndex(m_pcs->GetPrimaryVName(j));
	   // timelevel= 1;
	   nidx1 = nidx0+1;
	   CopyNodeVals(nidx1,nidx0);
	   }
	   }
	 */
}

/**************************************************************************
   FEMLib-Method:
   Task: Relocate Deformation process
   Programing:
   09/2004 WW Implementation
   10/2010 TF changes due to conversion from std::string to enum for process
type
**************************************************************************/
void RelocateDeformationProcess(CRFProcess* m_pcs)
{
	//   string pcs_name_dm = m_pcs->_pcs_type_name;
	FiniteElement::ProcessType pcs_name_dm(m_pcs->getProcessType());

	string num_type_name_dm;
	// Numerics
	if (m_pcs->num_type_name.compare("STRONG_DISCONTINUITY") == 0)
	{
		num_type_name_dm = m_pcs->num_type_name;
		enhanced_strain_dm = 1;
	}

	delete m_pcs;
	m_pcs = dynamic_cast<CRFProcess*>(new CRFProcessDeformation());
	m_pcs->setProcessType(pcs_name_dm);

	if (enhanced_strain_dm == 1) m_pcs->num_type_name = num_type_name_dm;
	pcs_deformation = 1;
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::PCSMoveNOD
   Task:
   Programming:
   08/2004 MB/OK Implementation
   last modified:
 **************************************************************************/
void CRFProcess::PCSMoveNOD(void)
{
	switch (this->type)
	{
		case 1:
			MeshLib::MSHMoveNODUcFlow(this);
			break;
		default:
			DisplayMsgLn("PCSMoveNOD: no valid process");
			abort();
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   10/2004 OK 2nd version
**************************************************************************/
std::string PCSProblemType()
{
	std::string pcs_problem_type;
	size_t no_processes(pcs_vector.size());

	for (size_t i = 0; i < no_processes; i++)
	{
		switch (pcs_vector[i]->getProcessType())
		{
			case FiniteElement::LIQUID_FLOW:
				pcs_problem_type = "LIQUID_FLOW";
				break;
			case FiniteElement::GROUNDWATER_FLOW:
				pcs_problem_type = "GROUNDWATER_FLOW";
				break;
			case FiniteElement::TWO_PHASE_FLOW:
				pcs_problem_type = "TWO_PHASE_FLOW";
				break;
			case FiniteElement::RICHARDS_FLOW:  // MX test 04.2005
				pcs_problem_type = "RICHARDS_FLOW";
				break;
			case FiniteElement::DEFORMATION:
				if (pcs_problem_type.empty())
					pcs_problem_type = "DEFORMATION";
				else
					pcs_problem_type += "+DEFORMATION";
				break;
			case FiniteElement::DEFORMATION_FLOW:
				if (pcs_problem_type.empty())
					pcs_problem_type = "DEFORMATION";
				else
					pcs_problem_type += "+DEFORMATION";
				break;
			case FiniteElement::HEAT_TRANSPORT:
				if (pcs_problem_type.empty())
					pcs_problem_type = "HEAT_TRANSPORT";
				else
					pcs_problem_type += "+HEAT_TRANSPORT";
				break;
			case FiniteElement::MASS_TRANSPORT:
				if (pcs_problem_type.empty())
					pcs_problem_type = "MASS_TRANSPORT";
				else
					pcs_problem_type += "+MASS_TRANSPORT";
				break;
			case FiniteElement::FLUID_MOMENTUM:
				if (pcs_problem_type.empty())
					pcs_problem_type = "FLUID_MOMENTUM";
				else
					pcs_problem_type += "+FLUID_MOMENTUM";
				break;
			case FiniteElement::RANDOM_WALK:
				if (pcs_problem_type.empty())
					pcs_problem_type = "RANDOM_WALK";
				else
					pcs_problem_type += "+RANDOM_WALK";
				break;
			default:
				pcs_problem_type = "";
		}
	}

	return pcs_problem_type;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2004 OK Implementation
**************************************************************************/
// void CRFProcess::CalcELEMassFluxes(void)
//{
/*OK411
   int i;
   double e_value = -1.0;
   string e_value_name;
   double geo_factor, density;
   double velocity = 0.0;
   int e_idx;
   int phase = 0;
   long e;
   CMediumProperties* m_mmp = NULL;
   CFluidProperties* m_mfp = NULL;
   m_mfp = mfp_vector[phase]; //OK ToDo
   //======================================================================
   for(e=0;e<ElementListLength;e++){
   m_mmp = mmp_vector[ElGetElementGroupNumber(e)];
   geo_factor = m_mmp->geo_area;
   density = m_mfp->Density();
   for(i=0;i<pcs_number_of_evals;i++){
   e_value_name = pcs_eval_data[i].name;
   e_idx = PCSGetELEValueIndex(pcs_eval_data[i].name);
   if(e_value_name.find("MASS_FLUX1_X")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_X"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   if(e_value_name.find("MASS_FLUX1_Y")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_Y"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   if(e_value_name.find("MASS_FLUX1_Z")!=string::npos){
   velocity = ElGetElementVal(e,PCSGetELEValueIndex("VELOCITY1_Z"));
   e_value = geo_factor * density * velocity;
   ElSetElementVal(e,e_idx,e_value);
   }
   }
   }
   //======================================================================
 */
//}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2005 OK Implementation
   02/2006 YD Dual Richards
   02/2007 WW General function for all unsaturated flow
**************************************************************************/
void CRFProcess::CalcSecondaryVariables(bool initial)
{
	//  char pcsT;
	//  pcsT = _pcs_type_name[0];
	//  if(type==1212) pcsT = 'V'; //WW
	//  switch(pcsT){
	//    case 'L':
	//      break;
	//    case 'U':
	//      break;
	//    case 'G':
	//      break;
	//    case 'T':
	//      break;
	//    case 'C':
	//      break;
	//    case 'R': // Richards flow
	//	  if(_pcs_type_name[1] == 'I')	// PCH To make a distinction with RANDOM
	// WALK.
	//        CalcSecondaryVariablesUnsaturatedFlow(initial); // WW
	//      break;
	//    case 'D':
	//      break;
	//    case 'V':
	//      CalcSecondaryVariablesUnsaturatedFlow(initial); //WW
	//      break;
	//	case 'P':
	//      CalcSecondaryVariablesPSGLOBAL(); //WW
	//      break;
	//  }

	switch (getProcessType())
	{
		case FiniteElement::LIQUID_FLOW:
			break;
		case FiniteElement::GROUNDWATER_FLOW:
			break;
		case FiniteElement::TWO_PHASE_FLOW:
			break;
		case FiniteElement::RICHARDS_FLOW:  // Richards flow
			// WW
			CalcSecondaryVariablesUnsaturatedFlow(initial);
			break;
		case FiniteElement::DEFORMATION || FiniteElement::DEFORMATION_FLOW ||
		    FiniteElement::DEFORMATION_DYNAMIC:
			if (type == 42)  // H2M //WW
				CalcSecondaryVariablesUnsaturatedFlow(initial);

			break;
		case FiniteElement::PS_GLOBAL:
			CalcSecondaryVariablesPSGLOBAL();  // WW
			break;
		default:
			if (type == 1212)
				// WW
				CalcSecondaryVariablesUnsaturatedFlow(initial);
			break;
	}
}

//////////////////////////////////////////////////////////////////////////
// ReMove site
//////////////////////////////////////////////////////////////////////////

/*************************************************************************
   ROCKFLOW - Function: GetCompNamehelp
   Task: Namepatch, until primary function names are finally sorted out
 //SB:todo
   Programming:	08/2003 SB Implementation
   last modified:
   superseded by GetPFNamebyCPName() but left here, as not all files are already
 in the new concept
 **************************************************************************/
/* SB: namepatch
   Repariert kurzfristig die Ausgabe
   input: datafield_n[j].name
   wenn dar Name "CONCENTRATIONx" ist, wird er durch den enstprechenden
   Komponentennamen ersetzt, sonst bleibts */

char* GetCompNamehelp(char* inname)
{
	int comp;  // WW, phase;
	char* outname, help[MAX_ZEILE];
	CRFProcess* m_pcs = NULL;
	outname = inname;
	// WW phase = 0;
	for (comp = 0; comp < GetRFProcessNumComponents(); comp++)
	{
		sprintf(help, "%s%d", "CONCENTRATION", comp + 1);
		/*  help has to be a part of inname (strstr) and also have the same
		 * length (strcmp) */
		if (strstr(inname, help) && (strcmp(inname, help) == 0))
		{
			m_pcs = m_pcs->GetProcessByFunctionName(help);
			if (m_pcs == NULL) break;
			//		outname =
			// GetTracerCompName(phase,m_pcs->GetProcessComponentNumber()-1);
			//		outname =
			// cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname;
			outname = (char*)cp_vec[m_pcs->GetProcessComponentNumber() - 1]
			              ->compname.data();
			return outname;
		}
	}
	return outname;
}  // SB:namepatch

/*************************************************************************
   ROCKFLOW - Function: GetCPNamebyPFName(string )
   Task: Replaces CP Name by Primary function name for output input
   Input:	component property name
   Output: primary function name
   Programming:	10/2004 SB Implementation
   10/2010 TF changed access to process type
 **************************************************************************/
string GetPFNamebyCPName(string inname)
{
	int i, j;  // WW, k;
	int pcs_vector_size = (int)pcs_vector.size();
	string outname;
	char help[MAX_ZEILE];
	CRFProcess* m_pcs = NULL;
	outname = "dummy";
	if (pcs_vector_size > 0)
		for (i = 0; i < pcs_vector_size; i++)
		{
			m_pcs = pcs_vector[i];
			//	if(m_pcs->_pcs_type_name.compare("MASS_TRANSPORT") == 0){ // if
			// this is mass transport // TF
			if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				j = m_pcs->GetProcessComponentNumber();
				// WW k = cp_vec[j]->transport_phase;
				outname = cp_vec[m_pcs->GetProcessComponentNumber()]->compname;
				if (outname == inname)  // right process found
				{
					sprintf(help, "%s%d", "CONCENTRATION", j);
					outname = help;
					return outname;
				}
			}
		}
	// for(i=0;comp<GetRFProcessNumComponents();i++) {
	//	sprintf(help,"%s%d","CONCENTRATION",i);
	/*  help has to be a part of inname (strstr) and also have the same length
	 * (strcmp) */
	//	if(strstr(inname, help) && (strcmp(inname,help) == 0)){
	//		m_pcs = m_pcs->GetProcessByFunctionName(help);
	//		if(m_pcs == NULL) break;
	//		outname = cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname;
	//		outname = (char *)
	// cp_vec[m_pcs->GetProcessComponentNumber()-1]->compname.data();
	//		if(outname.compare(inname) == 0)
	//			return outname;
	//	};
	// };
	// Inname is not from a mass transport process, therefore return inname
	return inname;
}  // SB:namepatch

//========================================================================

int GetRFProcessChemicalModel(void)
{
	cout << "GetRFProcessChemicalModel - to be removed" << endl;
	return 0;
}


int GetRFProcessHeatReactModel(void)
{
	cout << "GetRFProcessHeatReactModel - to be removed" << endl;
	return 0;
}

int GetRFProcessNumPhases(void)
{
	// DisplayMsgLn("GetRFProcessNumPhases - to be removed");
	int no_phases = (int)mfp_vector.size();
	return no_phases;
}

int GetRFProcessProcessing(char* rfpp_type)
{
	bool pcs_flow = false;
	bool pcs_deform = false;
	CRFProcess* m_pcs = NULL;
	size_t no_processes = pcs_vector.size();
	for (size_t i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->_pcs_type_name.find("DEFORMATION") != string::npos)
		if (isDeformationProcess(m_pcs->getProcessType())) pcs_deform = true;
		//		if (m_pcs->_pcs_type_name.find("FLOW") != string::npos)
		if (isFlowProcess(m_pcs->getProcessType())) pcs_flow = true;
	}

	if (strcmp(rfpp_type, "SD") == 0)
	{
		if (pcs_flow && pcs_deform) return 1;
	}
	else
		cout << "GetRFProcessProcessing - to be removed" << endl;
	return 0;
}

int GetRFProcessProcessingAndActivation(const char*)
{
	ScreenMessage("GetRFProcessProcessingAndActivation - to be removed\n");
	return 0;
}

long GetRFProcessNumComponents(void)
{
	// DisplayMsgLn("GetRFProcessNumComponents - to be removed");
	int no_components = (int)cp_vec.size();
	return no_components;
}

int GetRFControlModex(void)
{
	cout << "GetRFControlModex - to be removed" << endl;
	return 0;
}

int GetRFProcessNumContinua(void)
{
	cout << "GetRFProcessNumContinua - to be removed" << endl;
	return 0;
}

int GetRFProcessNumElectricFields(void)
{
	cout << "GetRFProcessNumElectricFields - to be removed" << endl;
	return 0;
}

int GetRFProcessNumTemperatures(void)
{
	cout << "GetRFProcessNumTemperatures - to be removed" << endl;
	return -1;
}

int GetRFProcessSimulation(void)
{
	cout << "GetRFProcessSimulation - to be removed" << endl;
	return -1;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void PCSDelete()
{
	for (int i = 0; i < (int)pcs_vector.size(); i++)
		delete pcs_vector[i];
	pcs_vector.clear();
	pcs_no_components = 0;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modified:
**************************************************************************/
void CRFProcess::SetNodeValue(long n, int nidx, double value)
{
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in  CRFProcess::SetNodeValue() " << endl;
		abort();
	}
#endif
	// WW 11.12.2012 	nod_val_vector[n][nidx] = value;
	nod_val_vector[nidx][n] = value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 PCH Implementation
   last modified:
**************************************************************************/
void CRFProcess::SetElementValue(long n, int nidx, double value)
{
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in  CRFProcess::SetElementValue() " << endl;
		abort();
	}
#endif
	ele_val_vector[n][nidx] = value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modified:
**************************************************************************/
double CRFProcess::GetNodeValue(size_t n, int nidx)
{
	double value;
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in  CRFProcess::GetNodeValue() " << endl;
		abort();
	}
#endif
	// WW 11.12.2012		value = nod_val_vector[n][nidx];
	value = nod_val_vector[nidx][n];
	return value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 PCH Implementation
   last modified:
**************************************************************************/
double CRFProcess::GetElementValue(size_t n, int nidx)
{
	double value;
#ifdef gDEBUG
	if (nidx < 0)
	{
		cout << " Fatal error in CRFProcess::GetElementValue() " << endl;
		abort();
	}
#endif
	value = ele_val_vector[n][nidx];
	return value;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   01/2006 OK Test
   03/2012 JT: Checks are unnecessary. Allow reverse order.
**************************************************************************/
int CRFProcess::GetNodeValueIndex(
    const std::string& var_name, bool reverse_order)  // JT: Allow reverse order
{
	if (!reverse_order)
	{
		for (size_t i = 0; i < nod_val_name_vector.size(); i++)
		{
			if (nod_val_name_vector[i].compare(var_name) == 0) return i;
		}
	}
	else
	{
		int nvals = ((int)nod_val_name_vector.size()) - 1;
		for (int j = nvals; j > -1; j--)
		{
			if (nod_val_name_vector[j].compare(var_name) == 0) return j;
		}
	}
	//
	return -2;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   09/2005 PCH Implementation
   last modified:
**************************************************************************/
int CRFProcess::GetElementValueIndex(const string& var_name, bool reverse_order)
{
	if (!reverse_order)
	{
		for (size_t i = 0; i < ele_val_name_vector.size(); i++)
		{
			if (ele_val_name_vector[i].compare(var_name) == 0) return i;
		}
	}
	else
	{
		int nvals = ((int)ele_val_name_vector.size()) - 1;
		for (int j = nvals; j > -1; j--)
		{
			if (ele_val_name_vector[j].compare(var_name) == 0) return j;
		}
	}
	//
	return -2;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   05/2005 OK pcs_pv_name,
   12/2005 OK RESTART
   07/2006 OK/MX MSH
**************************************************************************/
void CRFProcess::SetIC()
{
	CInitialCondition* m_ic = NULL;
	// HS, for MASS_TRANSPORT PCS,
	// it is not necessary to use PrimaryVarible as second check.
	// nidx will give the proper IC pointer.
	if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			int nidx = GetNodeValueIndex(pcs_primary_function_name[i]);

			for (size_t j = 0; j < ic_vector.size(); j++)
			{
				m_ic = ic_vector[j];
				m_ic->m_msh = m_msh;

				if (m_ic->getProcess() == this)
				{
					ScreenMessage("* %s on %s\n", pcs_primary_function_name[i], m_ic->getGeoName().data());
					m_ic->Set(nidx);
					m_ic->Set(nidx + 1);
				}
			}
		}
	}
	else  // otherwise PrimaryVariable check is still performed.
	{
		for (int i = 0; i < pcs_number_of_primary_nvals; i++)
		{
			int nidx = GetNodeValueIndex(pcs_primary_function_name[i]);
			FiniteElement::PrimaryVariable pv_i(
			    FiniteElement::convertPrimaryVariable(
			        pcs_primary_function_name[i]));
			for (size_t j = 0; j < ic_vector.size(); j++)
			{
				m_ic = ic_vector[j];
				m_ic->m_msh = m_msh;

				if (m_ic->getProcessType() != this->getProcessType()) continue;

				m_ic->setProcess(this);
				if (m_ic->getProcessPrimaryVariable() == pv_i)
				{
					ScreenMessage("* %s on %s\n", pcs_primary_function_name[i], m_ic->getGeoName().data());
					m_ic->Set(nidx);
					m_ic->Set(nidx + 1);
				}  // end of if
			}      // end of for j
		}          // end of for i
	}

	// end of if-else
}


/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   last modified:
**************************************************************************/
double CRFProcess::CalcNodeValueChanges(int ii)
{
	const int nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]);

#ifdef USE_PETSC
	long g_nnodes = m_msh->getNumNodesLocal();
#else
	long g_nnodes = m_msh->GetNodesNumber(false);
#endif

	double diff_norm = 0.0, current_norm = 0.0;
	//
	for (long i = 0; i < g_nnodes; i++)
	{
		double u0 = GetNodeValue(i, nidx1);
		double u1 = GetNodeValue(i, nidx1 + 1);
		diff_norm += (u1 - u0) * (u1 - u0);
		current_norm += u1 * u1;
	}

#ifdef USE_PETSC
	double diff_norm_global = 0;
	MPI_Allreduce(&diff_norm, &diff_norm_global, 1, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	diff_norm = std::sqrt(diff_norm_global);
	double current_norm_global = 0;
	MPI_Allreduce(&current_norm, &current_norm_global, 1, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	current_norm = std::sqrt(current_norm_global);
#else
	diff_norm = std::sqrt(diff_norm);
	current_norm = std::sqrt(current_norm);
#endif

	return diff_norm / (current_norm + DBL_EPSILON);
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   last modified:
**************************************************************************/
double CRFProcess::CalcVelocityChanges()
{
	assert(getProcessType() == FiniteElement::LIQUID_FLOW);

	double diff_norm = 0.0, current_norm = 0.0;
	//
	for (size_t i = 0; i < ele_gp_value.size(); i++)
	{
		const Matrix& v0 = ele_gp_value[i]->Velocity0;
		const Matrix& v1 = ele_gp_value[i]->Velocity;
		for (unsigned gp = 0; gp < v0.Cols(); gp++)
		{
			for (unsigned k = 0; k < v0.Rows(); k++)
			{
				diff_norm += std::pow(v1(k, gp) - v0(k, gp), 2);
				current_norm += std::pow(v1(k, gp), 2);
			}
		}
	}
#ifdef USE_PETSC
	double diff_norm_global = 0;
	MPI_Allreduce(&diff_norm, &diff_norm_global, 1, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	diff_norm = std::sqrt(diff_norm_global);
	double current_norm_global = 0;
	MPI_Allreduce(&current_norm, &current_norm_global, 1, MPI_DOUBLE, MPI_SUM,
	              MPI_COMM_WORLD);
	current_norm = std::sqrt(current_norm_global);
#else
	diff_norm = std::sqrt(diff_norm);
	current_norm = std::sqrt(current_norm);
#endif

	return diff_norm / (current_norm + DBL_EPSILON);
}

/**************************************************************************
   FEMLib-Method:
   Task: Ermittelt den Fehler bei Iterationen
   new_iteration     : Vektor des neuen Iterationsschritts
   old_iteration_ndx : Knotenindex fuer Werte des alten Iterationsschritts
   reference_ndx     : Knotenindex fuer Werte des alten Zeitschritts (als
Referenz)
   method            : Methode der Fehlerermittlung
   Programing:
   01/2005 OK NUM implementation
   05/2005 OK MSH
   08/2005 WW Re-implememtation based on NUMCalcIterationError
   01/2007 WW For DOF>1
   11/2007 WW Changes for the new classes of sparse and linear solver
   3/2012  JT Clean, add newton, add CPL vs. NLS, go to enum system
   last modification:
**************************************************************************/
//#if !defined(USE_PETSC) // && !defined(other parallel libs)//02.3013. WW
double CRFProcess::CalcIterationNODError(FiniteElement::ErrorMethod method,
                                         bool nls_error, bool cpl_error)
{
	long i, k, g_nnodes;
	double error, error_g, val1, val2, value;
	int nidx1, ii;
#ifndef USE_PETSC
	double* eqs_x = NULL;  // 11.2007. WW
#endif
	int num_dof_errors = pcs_number_of_primary_nvals;
	double unknowns_norm = 0.0;
	double absolute_error[DOF_NUMBER_MAX];
#ifdef USE_PETSC
	g_nnodes = m_msh->getNumNodesLocal();
#else
	g_nnodes = m_msh->GetNodesNumber(false);
#endif

#if defined(NEW_EQS)
	eqs_x = eqs_new->getX();
#elif !defined(USE_PETSC)
	eqs_x = eqs->x;
#endif

	switch (method)
	{
		//
		// --> ENORM:	|x1-x0|
		//     Norm of the solution vector delta (absolute error).
		//     Norm taken over entire solution vector (all primary variables)
		//     and checked against a single tolerance.
		//
		case FiniteElement::ENORM:
		{
			if (FiniteElement::isNewtonKind(m_num->nls_method))
			{  // NEWTON-RAPHSON
#ifndef USE_PETSC
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					for (i = 0; i < g_nnodes; i++)
					{
						val1 = eqs_x[i + ii * g_nnodes];
						unknowns_norm += val1 * val1;
					}
				}
#endif
			}
			else
			{  // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 =
					    GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
						k = m_msh->Eqs2Global_NodeIndex[i];
						val1 =
						    GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
						unknowns_norm += val1 * val1;
					}
				}
			}
			num_dof_errors = 1;
			unknowns_norm = sqrt(unknowns_norm);
			absolute_error[0] = unknowns_norm;
			break;
		}
		//
		// --> ERNORM:	|(x1-x0)/x0)|
		//     Norm of the solution vector delta divided by the solution vector
		//     (relative error).
		//     A single tolerance applied to all primary variables.
		//
		case FiniteElement::ERNORM:
		{
			value = 0.0;
			if (FiniteElement::isNewtonKind(m_num->nls_method))
			{  // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 =
					    GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
						k = m_msh->Eqs2Global_NodeIndex[i];
						val1 = eqs_x[i + ii * g_nnodes];
						val2 = GetNodeValue(k, nidx1);
						//
						unknowns_norm += val1 * val1;
						value += val2 * val2;
					}
				}
			}
			else
			{  // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 =
					    GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
						k = m_msh->Eqs2Global_NodeIndex[i];
						val1 = GetNodeValue(k, nidx1);
						val2 = val1 - eqs_x[i + ii * g_nnodes];
						//
						unknowns_norm += val2 * val2;
						value += val1 * val1;
					}
				}
			}
			num_dof_errors = 1;
			unknowns_norm = sqrt(unknowns_norm);
			absolute_error[0] = unknowns_norm / (sqrt(value) + DBL_EPSILON);
			break;
		}
		//
		// --> EVNORM:	|x1-x0|
		//     Norm of the solution vector delta (absolute error).
		//     Norm taken over solution vector of each primary variable, checked
		//     againes a tolerence specific to each variable.
		//
		case FiniteElement::EVNORM:
		{
			if (FiniteElement::isNewtonKind(m_num->nls_method))
			{  // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
						val1 = eqs_x[i + ii * g_nnodes];
						error += val1 * val1;
					}
					unknowns_norm += error;
					absolute_error[ii] = sqrt(error);
				}
			}
			else
			{  // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 =
					    GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
						k = m_msh->Eqs2Global_NodeIndex[i];
						val1 =
						    GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
						error += val1 * val1;
					}
					unknowns_norm += error;
					absolute_error[ii] = sqrt(error);
				}
			}
			unknowns_norm = sqrt(unknowns_norm);
			break;
		}
		//
		// --> BNORM: Get norm of solution vector, same as ENORM. RHS norm will
		// be calculated later.
		//     Norm of the solution vector delta (absolute error).
		//     Norm taken over entire solution vector (all primary variables)
		//     and checked against a single tolerance.
		//
		case FiniteElement::BNORM:
		{
			if (FiniteElement::isNewtonKind(m_num->nls_method))
			{  // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					for (i = 0; i < g_nnodes; i++)
					{
						val1 = eqs_x[i + ii * g_nnodes];
						unknowns_norm += val1 * val1;
					}
				}
			}
			else
			{  // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					nidx1 =
					    GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
						k = m_msh->Eqs2Global_NodeIndex[i];
						val1 =
						    GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
						unknowns_norm += val1 * val1;
					}
				}
			}
			num_dof_errors = 1;
			unknowns_norm = sqrt(unknowns_norm);
			absolute_error[0] = unknowns_norm;
			break;
		}
		//
		// --> LMAX:	max(x1-x0)
		//     Local max error (across all elements) of solution vector delta
		//     (absolute error).
		//     Tolerance required for each primary variable.
		//
		case FiniteElement::LMAX:
		{
			if (FiniteElement::isNewtonKind(m_num->nls_method))
			{  // NEWTON-RAPHSON
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					error = 0.0;
					for (i = 0; i < g_nnodes; i++)
					{
						val1 = eqs_x[i + ii * g_nnodes];
						unknowns_norm += val1 * val1;
						val1 = fabs(val1);
						if (val1 > error) error = val1;
					}
					absolute_error[ii] = error;
				}
			}
			else
			{  // PICARD
				for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
				{
					error = 0.0;
					nidx1 =
					    GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
					//
					for (i = 0; i < g_nnodes; i++)
					{
#ifdef USE_PETSC
						val1 = GetNodeValue(i, nidx1);
						val1 -= eqs_x[pcs_number_of_primary_nvals *
						                  m_msh->Eqs2Global_NodeIndex[i] +
						              ii];
#else
						k = m_msh->Eqs2Global_NodeIndex[i];
						val1 =
						    GetNodeValue(k, nidx1) - eqs_x[i + ii * g_nnodes];
#endif
						unknowns_norm += val1 * val1;
						error = std::max(error, fabs(val1));
					}
#ifdef USE_PETSC
					double error_l = error;
					MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX,
					              MPI_COMM_WORLD);
#endif
					absolute_error[ii] = error;
				}
			}
#ifdef USE_PETSC
			double unknowns_norm_l = unknowns_norm;
			MPI_Allreduce(&unknowns_norm_l, &unknowns_norm, 1, MPI_DOUBLE,
			              MPI_SUM, MPI_COMM_WORLD);
#endif
			unknowns_norm = sqrt(unknowns_norm);
			break;
		}
		//
		default:
			ScreenMessage(
			    "ERROR: Invalid error method for Iteration or Coupling Node "
			    "error.\n");
			return 0.0;
			//
			/*
			-----------------------------------------------------------------------------------------------
			ALTERNATIVE METHODS NOT YET IMPLEMENTED. MODIFY THEM AND ADD THEIR
			ENUM VALUES IF YOU WANT THEM.
			-----------------------------------------------------------------------------------------------
			// METHOD 4
			case 4:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = max_c = 0.0;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) +
			1;
			        //
			        for (i = 0l; i < g_nnodes; i++){
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] -
			GetNodeValue(k, nidx1)));
			           max_c = MMax(MMax(max_c,
			fabs(fabs(eqs_x[i+ii*g_nnodes]))),fabs(GetNodeValue(k, nidx1)));
			        }
			        pcs_absolute_error[ii] = error / (max_c + MKleinsteZahl);
			    }
			    break;
			//
			// METHOD 5
			case 5:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = max_c = 0.0;
			        min_c = 1.e99;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) +
			1;
			        //
			        for (i = 0l; i < g_nnodes; i++){
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] -
			GetNodeValue(k, nidx1)));
			           min_c = MMin(min_c, fabs(eqs_x[i+ii*g_nnodes]));
			           max_c = MMax(max_c, fabs(eqs_x[i+ii*g_nnodes]));
			        }
			        pcs_absolute_error[ii] = error / (max_c - min_c +
			MKleinsteZahl) ;
			    }
			    break;
			//
			// METHOD 6
			case 6:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = 0.0;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) +
			1;
			        //
			        for (i = 0l; i < g_nnodes; i++) {
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] -
			GetNodeValue(k, nidx1))
			              / (fabs(eqs_x[i+ii*g_nnodes] - GetNodeValue(k,
			nidx1-1)) + MKleinsteZahl));
			        }
			        pcs_absolute_error[ii] = error;
			    }
			    break;
			//
			// METHOD 7
			case 7:
			    for(ii=0;ii<pcs_number_of_primary_nvals;ii++)
			    {
			        error = change = max_c = 0.0;
			        min_c = 1.e99;
			        nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) +
			1;
			        //
			        for (i = 0l; i < g_nnodes; i++){
			           k = m_msh->Eqs2Global_NodeIndex[i];
			           error = MMax(error, fabs(eqs_x[i+ii*g_nnodes] -
			GetNodeValue(k, nidx1)));
			           change = MMax(change, fabs(eqs_x[i+ii*g_nnodes] -
			GetNodeValue(k, nidx1-1)));
			        }
			        pcs_absolute_error[ii] = error / (change + MKleinsteZahl);
			    }
			    break;
			//
			*/
	}
	//
	// Store the error (JT)
	// JT: now returning RELATIVE error. NECESSARY BECAUSE DOF MAY BE > 1 AND
	// EACH DOF MAY HAVE DIFFERENT CHARACTER
	if (cpl_error)
	{  // Return coupling error
		error_g = 0.0;
		for (ii = 0; ii < num_dof_errors; ii++)
		{
			cpl_absolute_error[ii] = absolute_error[ii];
			error = absolute_error[ii] / m_num->cpl_error_tolerance[ii];
			error_g = std::max(
			    error_g, error);  // Coupling error just stores the maximum
		}
		cpl_max_relative_error = error_g;
		cpl_num_dof_errors = num_dof_errors;
	}
	if (nls_error)
	{  // Return Non-Linear iteration error
		error_g = 0.0;
		for (ii = 0; ii < num_dof_errors; ii++)
		{
			pcs_absolute_error[ii] = absolute_error[ii];
			pcs_relative_error[ii] =
			    absolute_error[ii] / m_num->nls_error_tolerance[ii];
			error_g = std::max(error_g, pcs_relative_error[ii]);
		}
		pcs_num_dof_errors = num_dof_errors;
		pcs_unknowns_norm = unknowns_norm;
	}
	if (!nls_error && !cpl_error)
	{   // Then this routine called from somewhere else to get the error (i.e.
		// time control). Store it in a temporary vector for access.
		for (ii = 0; ii < num_dof_errors; ii++)
		{
			temporary_absolute_error[ii] = absolute_error[ii];
		}
		temporary_num_dof_errors = num_dof_errors;
	}
	//
	return error_g;  // Always returns the maximum relative error
}
//#endif // #if !defined(USE_PETSC)  WW

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation
   04/2006 YD Add contiuum coupling OK???Why ere
   04/2007 WW Remove the spurious stuff
   08/2008 WW Time step size control (First)
   12/2008 WW Time step size control (Update)
   07/2011 WW Newton-Raphson method
   3/2012  JT Clean, correct error obtainment, modify Newton convergence
criteria
**************************************************************************/
double CRFProcess::ExecuteNonLinear(int loop_process_number, bool print_pcs)
{
	double nl_itr_err = 1.0;
	double nl_theta, damping, norm_x0, norm_b0, norm_x, norm_b, val;
	double error_x1, error_x2, error_b1, error_b2 = .0, error, nl_itr_err_pre,
	                                     percent_difference;
	// double* eqs_x = NULL;     //
	double* eqs_b = NULL;
	bool converged;
	int ii, nidx1, num_fail = 0;
	size_t j, g_nnodes;

	string delim = " ";
	damping = 1.0;
	norm_x0 = norm_b0 = norm_x = norm_b = 0.;
	error = 1.;
	error_x2 = DBL_MAX;
	nl_theta = 1.0 - m_num->nls_relaxation;  // JT
	if (nl_theta < DBL_EPSILON) nl_theta = 1.0;
	g_nnodes = m_msh->GetNodesNumber(false);

#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
	eqs_x = eqs_new->GetGlobalSolution();
#endif
#ifdef NEW_EQS
	int k;
	eqs_x = eqs_new->getX();
	eqs_b = eqs_new->getRHS();
	configured_in_nonlinearloop = true;
// Also allocate temporary memory for linear solver. WW
//
	eqs_new->SetDOF(pcs_number_of_primary_nvals);  //_new 02/2010. WW
	eqs_new->ConfigNumerics(m_num->ls_precond, m_num->ls_method, m_num->ls_max_iterations, m_num->ls_error_tolerance, m_num->ls_storage_method, m_num->ls_extra_arg);
#endif
	//..................................................................
	// PI time step size control. 29.08.2008. WW
	if (Tim->GetPITimeStepCrtlType() > 0) CopyU_n();
	if (hasAnyProcessDeactivatedSubdomains) this->CheckMarkedElement();  // NW
	Tim->last_dt_accepted = true;  // JT2012

#if defined(USE_PETSC) || \
    defined(USE_MPI)  // || defined(other parallel libs)//01.3013. WW
	if (myrank == 0)
	{
#endif
		if (print_pcs)
		{  // JT: need check because of Regional Richards
			std::cout << "\n================================================"
			          << "\n";
			if (getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				std::cout << "->Process   " << loop_process_number << ": "
				          << convertProcessTypeToString(getProcessType())
				          << "\n";
				std::cout << "->Component " << pcs_component_number << ": "
				          << pcs_primary_function_name[0] << "\n";
			}
			else
			{
				std::cout << "->Process " << loop_process_number << ": "
				          << convertProcessTypeToString(getProcessType())
				          << "\n";
			}
			std::cout << "================================================"
			          << "\n";
		}
#if defined(USE_PETSC) || \
    defined(USE_MPI)  // || defined(other parallel libs)//01.3013. WW#ifdef
		              // USE_MPI
	}
#endif

	// ------------------------------------------------------------
	// NON-LINEAR ITERATIONS (OR SINGLE ITERATION IF LINEAR)
	// ------------------------------------------------------------
	diverged = false;
	converged = false;
	accepted = true;
	nl_itr_err_pre = 1.0;
	num_fail = 0;
	for (iter_nlin = 0; iter_nlin < m_num->nls_max_iterations; iter_nlin++)
	{
		nl_itr_err = Execute();
		if (nl_itr_err == -1.)  // linear solve broken
		{
			ScreenMessage("*** Nonlinear solve failed\n");
			accepted = false;
			Tim->last_dt_accepted = false;
			break;
		}
		//
		// ---------------------------------------------------
		// LINEAR SOLUTION
		// ---------------------------------------------------
        if (m_num->nls_method == FiniteElement::NL_LINEAR)
		{
			PrintStandardIterationInformation(true, nl_itr_err);
			converged = true;
		}
		else
		{
			// ---------------------------------------------------
			// NON-LINEAR SOLUTION
			// ---------------------------------------------------
			//
			damping = nl_theta;
			switch (m_num->getNonLinearErrorMethod())
			{
				// For most error methods (also works for Newton)
				default:
					PrintStandardIterationInformation(true, nl_itr_err);
					if (iter_nlin == 0)
						percent_difference = .0;
					else
						percent_difference =
						    100 *
						    ((nl_itr_err_pre - nl_itr_err) / nl_itr_err_pre);
					//						conv_rate = std::max(conv_rate, nl_itr_err
					///
					// nl_itr_err_pre);
					//						if (conv_rate)
					ScreenMessage(
					    "\tNonlinear error: %e, Conv. rate: %2.1f\%\n",
					    nl_itr_err, percent_difference);
					ScreenMessage(
					    "------------------------------------------------\n");
					//
					if (nl_itr_err <= 1.0)
					{
						converged = true;
					}
					else
					{  // Check for stagnation
						if (iter_nlin > 0)
						{
							// std::cout << "-> error_k=" << last_error << ",
							// error_k1=" << nonlinear_iteration_error << ",
							// convergence rate=" << percent_difference << "%"<<
							// "\n";
							if (percent_difference <
							    1.0)  // less than 1% difference (or an error
								      // increase) from previous error
								num_fail++;
							else
								num_fail = 0;
							//
							if (num_fail > 1)
								diverged =
								    true;  // require 2 consecutive failures
						}
						else
						{
							num_fail = 0;
						}
						nl_itr_err_pre = nl_itr_err;
					}
					break;

				// For (OGS) classic Newton error control
				case FiniteElement::BNORM:
					PrintStandardIterationInformation(false, nl_itr_err);
//
#if defined(USE_PETSC)  // || defined(other parallel libs)//06.3012. WW
					norm_x = eqs_new->GetVecNormX();
					norm_b = eqs_new->GetVecNormRHS();
#else
					norm_x = pcs_unknowns_norm;  // JT: this is already obtained
					                             // in CalcIterationNodeError.
					norm_b = 0.0;                // must calculate this
					for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
					{
						for (j = 0; j < g_nnodes; j++)
						{
							val = eqs_b[j + ii * g_nnodes];
							norm_b += val * val;
						}
					}
					norm_b = sqrt(norm_b);
#endif
					//
					if (iter_nlin == 0)
					{
						norm_x0 = norm_x;
						norm_b0 = norm_b;
						error_x2 = error_b2 = DBL_MAX;
					}
					else
					{
						error_x1 = norm_x / norm_x0;
						error_b1 = norm_b / norm_b0;

						if (norm_x < m_num->nls_error_tolerance[0] &&
						    error_x1 > norm_x)
							error_x1 = norm_x;
						if (norm_b < m_num->nls_error_tolerance[0] &&
						    error_b1 > norm_b)
							error_b1 = norm_b;
						if (error_x1 / error_x2 > 0.1 ||
						    error_b1 / error_b2 > 0.1)
							damping *= 0.5;  // take 1/2 of original theta
						//
						error = max(error_x1, error_b1);
						error_x2 = error_x1;
						error_b2 = error_b1;
						//
						// Check for divergence
						if (error > 10.0 && iter_nlin > 1)
						{
							diverged = true;
							if (Tim->GetPITimeStepCrtlType() > 0)
							{  // if PI automatic time control
								accepted = false;
								PI_TimeStepSize();
							}
						}
						//
						// Check convergence
						if (norm_x0 < m_num->nls_error_tolerance[0])
						{
							error = norm_x0;
							converged = true;
						}
						if (norm_b0 < 10 * m_num->nls_error_tolerance[0])
						{
							error = norm_b0;
							converged = true;
						}
						if (norm_b < 0.001 * norm_b0)
						{
							error = norm_b;
							converged = true;
						}
						if (error <= m_num->nls_error_tolerance[0])
						{
							converged = true;
						}
					}
// Newton information printout.
#if defined(USE_MPI) || defined(USE_PETSC)
					if (myrank == 0)
					{
#endif
						cout.width(10);
						cout.precision(3);
						cout.setf(ios::scientific);
						cout << "         NR-Error  |"
						     << "    RHS Norm|"
						     << "  Unknowns Norm|"
						     << " Damping\n";
						cout << "         " << setw(10) << error << "|  "
						     << setw(9) << norm_b << "| ";
						cout << setw(14) << norm_x << "| " << setw(9) << damping
						     << "\n";
#if defined(USE_MPI) || defined(USE_PETSC)
					}
#endif
					break;
			}
		}

		// FOR NEWTON: COPY DAMPED CHANGES TO NEW TIME
		// ---------------------------------------------------
		if (FiniteElement::isNewtonKind(m_num->nls_method))
		{
			if (converged)
				damping = 1.0;  // Solution has converged. Take newest values.
			//
			for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
				for (j = 0; j < g_nnodes; j++)
				{
					SetNodeValue(
					    j, nidx1,
					    GetNodeValue(j, nidx1) +
					        damping * eqs_x[m_msh->Eqs2Global_NodeIndex[j] *
					                            pcs_number_of_primary_nvals +
					                        ii]);
				}
#else
				const long ish = ii * g_nnodes;
				for (j = 0; j < g_nnodes; j++)
				{
					k = m_msh->Eqs2Global_NodeIndex[j];
					val = GetNodeValue(k, nidx1) + damping * eqs_x[j + ish];
					SetNodeValue(k, nidx1, val);
				}
#endif
			}
		}

		// OTHER CONFIGURATIONS AT THE END OF THIS NON-LINEAR ITERATION
		// ---------------------------------------------------
		if (mobile_nodes_flag == 1)
		{
			PCSMoveNOD();
		}

		// CHECK FOR TIME STEP FAILURE
		// ---------------------------------------------------
		if (diverged || !accepted || Tim->isDynamicTimeFailureSuggested(this))
		{
			accepted = false;
			Tim->last_dt_accepted = false;
			break;
		}

		/* JT: I don't know if this time control method is used anymore. But it
		   relies on a single error
		       produced from CalcIterationNodeError(), but this now depends on
		   the type of error to use.
		       Therefore, I simply provide the error of the first dof, and not
		   depending on the error type. If
		       this time step is still used, someone will need to find another
		   way to calculate the error it uses.
		*/
		Tim->repeat = true;
		if (converged)
		{
			Tim->repeat = false;
			Tim->nonlinear_iteration_error = pcs_absolute_error[0];
		}

		// BREAK CRITERIA
		if (converged || diverged)
		{
			break;
		}
	}
	if ((!converged || diverged) && m_num->nls_max_iterations > 1)
		accepted = false;
	iter_nlin_max = std::max(iter_nlin_max, iter_nlin);
	// ------------------------------------------------------------
	// NON-LINEAR ITERATIONS COMPLETE
	// ------------------------------------------------------------
	// PI time step size control. 27.08.2008. WW
	if (accepted && Tim->GetPITimeStepCrtlType() > 0)
	{
		PI_TimeStepSize();  // might also set accepted to false here.
	}
	//
	if (m_num->nls_max_iterations > 1)  // only for non-linear iterations
	{
		if (diverged)
		{
			if (accepted)
			{  // only increment if not fixed by a failed time step.
				num_diverged++;
				ScreenMessage("\nNon-linear iteration stabilized.\n");
			}
			else
			{
				ScreenMessage("\nNon-linear iteration diverged.\n");
			}
		}
		else if (!converged)
		{
			if (accepted)  // only increment if not fixed by a failed time step.
				num_notsatisfied++;
			if (Tim->GetPITimeStepCrtlType() < 1)  // PI has the intrinsic
				                                   // property of doing this. So
				                                   // don't print it.
				ScreenMessage(
				    "\nMax number of non-linear iterations reached.\n");
		}
	}
	//
	// Calculate secondary variables
	if (accepted)
	{
		CalcSecondaryVariables();
	}
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03.3012. WW
                         // Release temporary memory of linear solver. WW
#ifdef NEW_EQS           // WW
	eqs_new->Clean();  // Release buffer momery WW
	configured_in_nonlinearloop = false;
#endif
#endif
	return nl_itr_err;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   3/2012  JT
**************************************************************************/
void CRFProcess::PrintStandardIterationInformation(bool write_std_errors,
                                                   double nl_error)
{
	int ii;
	//
	// LINEAR SOLUTION
    if (m_num->nls_method == FiniteElement::NL_LINEAR)
	{
		ScreenMessage("      -->LINEAR solution complete. \n");
		if (write_std_errors)
		{
			for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				ScreenMessage("         PCS error DOF[%d]: %g\n", ii,
				              pcs_absolute_error[ii]);
			}
		}
		return;
	}
	//
	// NON-LINEAR METHODS
	if (m_num->nls_method == FiniteElement::NL_PICARD)
		ScreenMessage("-->End of PICARD iteration: %d/%d, error=%g \n",
		              iter_nlin, m_num->nls_max_iterations, nl_error);
	else
		ScreenMessage("-->End of NEWTON-RAPHSON iteration: %d/%d \n", iter_nlin,
		              m_num->nls_max_iterations);
	//
	// Errors
	// --------------------------------------------------
	if (write_std_errors)
	{
		if (pcs_num_dof_errors == 1)
		{
			ScreenMessage("\tPCS error: %g\n", pcs_absolute_error[0]);
		}
		else
		{
			for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
			{
				ScreenMessage("\tPCS error DOF[%d]: %e\n", ii,
				              pcs_absolute_error[ii]);
			}
		}
		ScreenMessage("\tEuclidian norm of unknowns: %e\n", pcs_unknowns_norm);
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   02/2007 WW Implementation
**************************************************************************/
void CRFProcess::Extropolation_GaussValue()
{
	int k, NS;
	long i = 0;
	int idx[3];
	// const long LowOrderNodes= m_msh->GetNodesNumber(false);
	MeshLib::CElem* elem = NULL;

	//
	NS = m_msh->GetCoordinateFlag() / 10;
	idx[0] = GetNodeValueIndex("VELOCITY_X1");
	idx[1] = GetNodeValueIndex("VELOCITY_Y1");
	idx[2] = GetNodeValueIndex("VELOCITY_Z1");
	for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		for (k = 0; k < NS; k++)
			SetNodeValue(i, idx[k], 0.0);
	if (type == 1212 || type == 1313)  // Multi-phase flow
	{
		idx[0] = GetNodeValueIndex("VELOCITY_X2");
		idx[1] = GetNodeValueIndex("VELOCITY_Y2");
		idx[2] = GetNodeValueIndex("VELOCITY_Z2");
		for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
			for (k = 0; k < NS; k++)
				SetNodeValue(i, idx[k], 0.0);
	}
	//
	for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			fem->ConfigElement(elem);
			for (k = 0; k < NS; k++)
				fem->ExtropolateGauss(this, k);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:    Calculate the material values at Gauss points and extropolate them
         to element node
   Programing:
   04/2007 WW Implementation
**************************************************************************/
void CRFProcess::Extropolation_MatValue()
{
	//	if (_pcs_type_name.find("FLOW") == string::npos)
	if (!isFlowProcess(this->getProcessType())) return;
	if (additioanl2ndvar_print < 0) return;

	//
	int NS = m_msh->GetMaxElementDim();
	//
	if ((additioanl2ndvar_print > 0) && (additioanl2ndvar_print < 3))
	{
		int idx[3];
		idx[0] = GetNodeValueIndex("PERMEABILITY_X1");
		idx[1] = GetNodeValueIndex("PERMEABILITY_Y1");
		if (NS > 2) idx[2] = GetNodeValueIndex("PERMEABILITY_Z1");
		for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			for (int k = 0; k < NS; k++)
				SetNodeValue(i, idx[k], 0.0);
	}
	if (additioanl2ndvar_print > 1)
	{
		int idxp = GetNodeValueIndex("POROSITY");
		for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
			SetNodeValue(i, idxp, 0.0);
	}
	//
	continuum = 0;
	if (continuum_vector.size() == 2) continuum = 1;

	MeshLib::CElem* elem = NULL;
	for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
	{
		elem = m_msh->ele_vector[i];
		if (elem->GetMark())  // Marked for use
		{
			fem->ConfigElement(elem);
			fem->CalcNodeMatParatemer();
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
**************************************************************************/
void PCSDelete(const std::string& m_pcs_type_name)
{
	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(m_pcs_type_name));
	CRFProcess* m_pcs = NULL;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->_pcs_type_name.compare(m_pcs_type_name) == 0) { TF
		if (m_pcs->getProcessType() == pcs_type)
		{
			delete m_pcs;
			pcs_vector.erase(pcs_vector.begin() + i);
		}
	}
}

/**************************************************************************
   GeoSys - Function: Reallocation

   Aufgabe:
        Reallocte memory by new operator
   09/2005   WW    Erste Version

**************************************************************************/
template <class T>
T* resize(T* array, size_t old_size, size_t new_size)
{
	T* temp = new T[new_size];
	for (size_t i = 0; i < old_size; i++)
		temp[i] = array[i];
	for (size_t i = old_size; i < new_size; i++)
		temp[i] = 0;
	delete[] array;
	array = temp;
	return temp;
}


/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 MB implementation
**************************************************************************/
void CRFProcess::CopyCouplingNODValues()
{
	// Carefull if cpl_variable = primary variable -> need extra column in
	// NodeValueTable !
	int nidx0 = GetNodeValueIndex(m_num->cpl_variable_JOD);
	for (size_t l = 0; l < m_msh->GetNodesNumber(false); l++)
		SetNodeValue(l, nidx0, GetNodeValue(l, nidx0 + 1));
	//	if (_pcs_type_name.find("RICHARDS") != string::npos) { //WW
	if (this->getProcessType() == FiniteElement::RICHARDS_FLOW)  // WW
	{
		nidx0 = GetNodeValueIndex("SATURATION1");
		for (size_t l = 0; l < m_msh->GetNodesNumber(false); l++)
			SetNodeValue(l, nidx0, GetNodeValue(l, nidx0 + 1));
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2005 MB implementation
   02/2006 WW Modified for the cases of high order element and saturation
   08/2008 WW Make it twofold copy: forward and backward
**************************************************************************/
void CRFProcess::CopyTimestepNODValues(bool forward)
{
	bool Quadr = false;  // WW
	if (type == 4 || type == 41) Quadr = true;

	for (int j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		int nidx0 = GetNodeValueIndex(pcs_primary_function_name[j]);
		int nidx1 = nidx0 + 1;
		if (!forward)  // 08.2008. WW
		{
			nidx0++;
			nidx1--;
		}
		for (size_t l = 0; l < m_msh->GetNodesNumber(Quadr); l++)
			SetNodeValue(l, nidx0, GetNodeValue(l, nidx1));
		// WW
		//		if (_pcs_type_name.find("RICHARDS") != string::npos || type ==
		// 1212) { //Multiphase. WW
		// Multiphase. WW
		if (this->getProcessType() == FiniteElement::RICHARDS_FLOW ||
		    getProcessType() == FiniteElement::MULTI_PHASE_FLOW || type == 42)
		{
			if (j == 1 && (getProcessType() == FiniteElement::MULTI_PHASE_FLOW || type == 42))  // Multiphase. WW
				continue;
			if (j == 0)
				nidx0 = GetNodeValueIndex("SATURATION1");
			else if (j == 1)
				nidx0 = GetNodeValueIndex("SATURATION2");
			nidx1 = nidx0 + 1;
			if (!forward)  // 27.08.2008. WW
			{
				nidx0++;
				nidx1--;
			}
			//
			for (size_t l = 0; l < m_msh->GetNodesNumber(false); l++)
				SetNodeValue(l, nidx0, GetNodeValue(l, nidx1));
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   3/2012 JT. Based on the nodal version
**************************************************************************/
void CRFProcess::CopyTimestepELEValues(bool forward)
{
	size_t j, nvals, nidx0, nidx1;
	long iel, num_ele;
	bool copy_porosity = true;
	nvals = ele_val_name_vector.size();
	if (nvals < 2) return;  // then we don't have 2 time levels of anything
	num_ele = (long)m_msh->ele_vector.size();
//
#ifdef GEM_REACT
	// do nothing as porosity update is handled by REACT_GEMS after!! flow and
	// transport solution
	copy_porosity = false;
#endif
	//
	for (j = 0; j < nvals - 1; j++)
	{
		if (ele_val_name_vector[j].compare(ele_val_name_vector[j + 1]) !=
		    0)  // If not the same, then we only have a single time slot
			continue;
		if (ele_val_name_vector[j].find("POROSITY") != string::npos &&
		    !copy_porosity)  // Porosity
			continue;
		//
		nidx0 = j;
		nidx1 = j + 1;
		if (!forward)
		{
			nidx0++;
			nidx1--;
		}
		for (iel = 0; iel < num_ele; iel++)
		{
			SetElementValue(iel, nidx0, GetElementValue(iel, nidx1));
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2006 SB Implementation
   last modified:
**************************************************************************/
int PCSGetPCSIndex(const string& pcs_type_name, const string& comp_name)
{
	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(pcs_type_name));

	CRFProcess* m_pcs = NULL;
	int i, pcs_no;
	int no_processes = (int)pcs_vector.size();
	string testname;
	pcs_no = -1;
	for (i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		//		if (m_pcs->pcs_type_name.compare(pcs_type_name) == 0) {
		if (m_pcs->getProcessType() == pcs_type)
		{
			testname = m_pcs->pcs_primary_function_name[0];
			if (testname.compare(comp_name) == 0)
			{
				//        cout << " Found in PCSGetbyTypeandCompName for
				//        PCSType/Compname " << pcs_type_name << ", " <<
				//        comp_name;
				//        cout << " Process number " << m_pcs->pcs_number << ",
				//        compnumber " << m_pcs->pcs_component_number << endl;
				pcs_no = i;
				return pcs_no;
			}
		}
	}
	return pcs_no;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2005 SB Implementation
   last modified:
   10/2010 TF restructured function
**************************************************************************/
CRFProcess* PCSGet(const std::string& pcs_type_name,
                   const std::string& comp_name)
{
	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(pcs_type_name));
	size_t no_processes(pcs_vector.size());
	for (size_t i = 0; i < no_processes; i++)
		//		if (pcs_vector[i]->pcs_type_name.compare(pcs_type_name) == 0) {
		////
		// TF
		if (pcs_vector[i]->getProcessType() == pcs_type)
			if (comp_name.compare(
			        pcs_vector[i]->pcs_primary_function_name[0]) == 0)
				return pcs_vector[i];

	return NULL;
}

CRFProcess* PCSGet(FiniteElement::ProcessType pcs_type,
                   const std::string& comp_name)
{
	size_t no_processes(pcs_vector.size());
	for (size_t i = 0; i < no_processes; i++)
		if (pcs_vector[i]->getProcessType() == pcs_type)
			if (comp_name.compare(
			        pcs_vector[i]->pcs_primary_function_name[0]) == 0)
				return pcs_vector[i];

	return NULL;
}

/**************************************************************************
   PCSLib-Method:
   12/2005 OK Implementation
**************************************************************************/
CRFProcess* PCSGet(const string& var_name, bool bdummy)
{
	string pcs_var_name;
	CRFProcess* m_pcs = NULL;
	bdummy = bdummy;  // WW
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		for (size_t j = 0; j < m_pcs->GetPrimaryVNumber(); j++)
		{
			pcs_var_name = m_pcs->pcs_primary_function_name[j];
			if (pcs_var_name.compare(var_name) == 0) return m_pcs;
		}
		for (size_t j = 0; j < m_pcs->GetSecondaryVNumber(); j++)
		{
			pcs_var_name = m_pcs->pcs_secondary_function_name[j];
			if (pcs_var_name.compare(var_name) == 0) return m_pcs;
		}
	}
	return NULL;
}

/**************************************************************************
   PCSLib-Method:
   05/2006 CMCD Implementation
**************************************************************************/
CRFProcess* PCSGetFluxProcess()
{
	CRFProcess* m_pcs = NULL;
	bool found = false;
	const size_t no_processes(pcs_vector.size());

	for (size_t i = 0; i < no_processes; i++)
	{
		//		if (pcs_vector[i]->_pcs_type_name == "LIQUID_FLOW") { // TF
		if (pcs_vector[i]->getProcessType() == FiniteElement::LIQUID_FLOW)
		{
			m_pcs = pcs_vector[i];
			found = true;
		}
		//		if (pcs_vector[i]->_pcs_type_name == "GROUNDWATER_FLOW") {
		if (pcs_vector[i]->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
		{
			m_pcs = pcs_vector[i];
			found = true;
		}
		if (found) return m_pcs;
	}
	return NULL;
}

/*************************************************************************
   GeoSys-FEM Function:
   Task:
   For Richards, dual-Richards and multi-phase flow.
   Programming:
   02/2007 WW Implementation
 **************************************************************************/
void CRFProcess::CalcSecondaryVariablesUnsaturatedFlow(bool initial)
{
	int jj;
	long i;
	double p_cap;
	int idxp, idxcp, idxS;
	//----------------------------------------------
	vector<string> secondSNames;
	vector<string> secondCPNames;
	secondSNames.push_back("SATURATION1");
	if (getProcessType() == FiniteElement::MULTI_PHASE_FLOW || type == 42)  // Multiphase flow
		secondCPNames.push_back("PRESSURE1");
	else
		secondCPNames.push_back("PRESSURE_CAP1");
	if (continuum_vector.size() == 2)
	{
		secondSNames.push_back("SATURATION2");
		secondCPNames.push_back("PRESSURE_CAP2");
	}

	//
	CElem* elem = NULL;
	CFiniteElementStd* fem = GetAssembler();
	//----------------------------------------------------------------------
	for (int ii = 0; ii < (int)secondSNames.size(); ii++)
	{
		idxS = GetNodeValueIndex(secondSNames[ii].c_str()) + 1;
		if (getProcessType() == FiniteElement::MULTI_PHASE_FLOW || type == 42)  // Multiphase flow
		{
			jj = ii;
			if (type == 42) jj += problem_dimension_dm;
			idxcp = GetNodeValueIndex(pcs_primary_function_name[jj]) + 1;
			idxp = GetNodeValueIndex("PRESSURE_W");
		}
		else
		{
			idxp = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
			idxcp = GetNodeValueIndex(secondCPNames[ii].c_str());
		}
		//----------------------------------------------------------------------
		// Capillary pressure
		for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		{
			// Copy the new saturation to the old level
			//
			if (getProcessType() == FiniteElement::MULTI_PHASE_FLOW || type == 42)  // Multiphase flow
			{
				// p_w = p_g-p_c
				p_cap = GetNodeValue(i, idxcp + 2) - GetNodeValue(i, idxcp);
				SetNodeValue(i, idxp, p_cap);
			}
			else
			{
				p_cap = -GetNodeValue(i, idxp);
				SetNodeValue(i, idxcp, p_cap);
			}
			SetNodeValue(i, idxS, 0.0);
		}
	}
	// Cal S
	for (int ii = 0; ii < (int)secondSNames.size(); ii++)
	{
		continuum = ii;
		for (i = 0; i < (long)m_msh->ele_vector.size(); i++)
		{
			elem = m_msh->ele_vector[i];
			if (elem->GetMark())  // Marked for use
			{
				elem->SetOrder(false);
				fem->ConfigElement(elem, false);
				fem->CalcSatution();
			}
		}
	}
	//
	if (!initial) return;
	//----------
	for (int ii = 0; ii < (int)secondSNames.size(); ii++)
	{
		idxS = GetNodeValueIndex(secondSNames[ii].c_str()) + 1;
		for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
			SetNodeValue(i, idxS - 1, GetNodeValue(i, idxS));
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   Task: Updating secondary variables for Multiphase flow in PS_GLOBAL

   Programming:
   03/2009 PCH Implementation
 **************************************************************************/
void CRFProcess::CalcSecondaryVariablesPSGLOBAL()
{
	long i;
	int ndx_pressure1, ndx_p_cap, ndx_pressure2,
	    ndx_s_wetting;  // WW, ndx_s_nonwetting;

	// The primary variables
	ndx_pressure1 = GetNodeValueIndex("PRESSURE1");
	// WW ndx_s_nonwetting = GetNodeValueIndex("SATURATION2");

	// The secondary variables
	ndx_pressure2 = GetNodeValueIndex("PRESSURE2");
	ndx_p_cap = GetNodeValueIndex("PRESSURE_CAP");
	ndx_s_wetting = GetNodeValueIndex("SATURATION1");

	double pressure1, pressure2, p_cap, s_wetting;
	for (i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
	{
		pressure1 = GetNodeValue(i, ndx_pressure1 + 1);  // New
		pressure2 = GetNodeValue(i, ndx_pressure1);      // Old

		// Let's get capillary pressure before updating pressure2
		// by accessing the primary variable of the saturation equation
		// not the secondary variable of it.
		int ndx_sat2 = GetNodeValueIndex("SATURATION2");
		double sat2 = GetNodeValue(i, ndx_sat2 + 1);
		// Due to the iterative solution scheme in solving Snw with no
		// explicit boundary condition for non-zero flux condition,
		// Snw may become negative particularly the density difference
		// between two fluids is big. To prevent negative Snw, the
		// saturation restriction added.
		CMediumProperties* mmp = NULL;
		if (mmp_vector.size() > 1)
		{
			double sum = 0.0;
			CNode* thisNode = m_msh->nod_vector[i];
			int NumOfNeighborElements =
			    (int)thisNode->getConnectedElementIDs().size();
			// Harmonic mean
			for (int i = 0; i < NumOfNeighborElements; ++i)
			{
				// Mount neighboring elemenets and get the corresponding
				// material group one by one.
				size_t eleIdx = thisNode->getConnectedElementIDs()[i];
				CElem* thisEle = m_msh->ele_vector[eleIdx];
				size_t matgrp = thisEle->GetPatchIndex();
				mmp = mmp_vector[matgrp];
				mmp->mode = 2;
				sum += 1.0 / sat2;
			}
			sat2 = (double)NumOfNeighborElements / sum;
		}
		else
			mmp = mmp_vector[0];
		s_wetting = 1.0 - sat2;
		// Assigning the secondary variable, Sw
		SetNodeValue(i, ndx_s_wetting, s_wetting);
		// Assigning the primary variable Snw here one more time
		// to completely bound the range of saturation
		//	SetNodeValue(i,ndx_s_nonwetting,sat2);
		//	SetNodeValue(i,ndx_s_nonwetting+1,sat2);

		// Assigning the secondary variable, Pc
		if (mmp_vector.size() > 1)
			p_cap = GetCapillaryPressureOnNodeByNeighobringElementPatches(
			    i, 2, 1.0 - sat2);
		else
			p_cap = mmp->CapillaryPressureFunction(1.0 - sat2);

		SetNodeValue(i, ndx_p_cap, p_cap);

		pressure2 = pressure1 + p_cap;
		// Assigning the secondary variables, Pnw
		SetNodeValue(i, ndx_pressure2, pressure2);
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate saturation on node by averaging the patches of the
         neighboring elements in three means

   0: Arithmatic mean
   1: Geomtric mean
   2: Harmonic mean

   Programing:
   03/2009 PCH Implementation
   last modification:
 *************************************************************************/
double CRFProcess::GetCapillaryPressureOnNodeByNeighobringElementPatches(
    int nodeIdx, int meanOption, double Sw)
{
	double p_cap = 0.0, sum = 0.0;

	CNode* thisNode = m_msh->nod_vector[nodeIdx];
	int NumOfNeighborElements = (int)thisNode->getConnectedElementIDs().size();

	switch (meanOption)
	{
		case 0:
			break;
		case 1:
			break;
		case 2:  // Harmonic mean
			for (int i = 0; i < NumOfNeighborElements; ++i)
			{
				// Mount neighboring elemenets and get the corresponding
				// material group one by one.
				int eleIdx = thisNode->getConnectedElementIDs()[i];
				CElem* thisEle = m_msh->ele_vector[eleIdx];
				int matgrp = thisEle->GetPatchIndex();
				CMediumProperties* mmp = mmp_vector[matgrp];
				sum += 1.0 / mmp->CapillaryPressureFunction(Sw);
			}
			p_cap = (double)NumOfNeighborElements / sum;
			break;

		default:
			cout << "Please define the option for various means!" << endl;
			cout << "The code stops at "
			        "GetCapillaryPressureOnNodeByNeighobringElementPatches "
			        "function!" << endl;
			abort();
			break;
	}

	return p_cap;
}

/*************************************************************************
   GeoSys-FEM Function:
   Task: Calculates saturation for richards flow,
      alternative to CRFProcess::CalcSecondaryVariablesRichards
      uses pressure of only one node
     evoked by switch case in pcs-file: SATURATION_SWITCH = true
   Programming:
   06/2007 JOD Implementation
 **************************************************************************/
void CRFProcess::CalcSaturationRichards(int timelevel, bool update)
{
	double p_cap, saturation, volume_sum;
	int idxp, idxcp, idxS, idx_tS = -1;
	const size_t number_continuum(continuum_vector.size());
	size_t i_s, i_e;

	CMediumProperties* m_mmp = NULL;
	CElem* elem = NULL;
	// WW  CFiniteElementStd* fem = GetAssembler();

	if (continuum_ic)  // Create IC: for both continua
	{
		i_s = 0;
		i_e = number_continuum;
	}
	else
	{
		i_s = continuum;
		i_e = continuum + 1;
	}

	for (size_t i_pv = i_s; i_pv < i_e; i_pv++)
	{
		idxp = GetNodeValueIndex(pcs_primary_function_name[i_pv]) + timelevel;
		idxS = GetNodeValueIndex(
		           pcs_secondary_function_name[i_pv * number_continuum]) +
		       timelevel;
		idxcp = GetNodeValueIndex(
		            pcs_secondary_function_name[i_pv * number_continuum +
		                                        number_continuum * 2]) +
		        timelevel;
		if (continuum_vector.size() > 1)
			idx_tS = GetNodeValueIndex("TOTAL_SATURATION") + timelevel;

		for (long i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		{
			// Capillary pressure
			p_cap = -GetNodeValue(i, idxp);
			if (timelevel == 1 && update)
				SetNodeValue(i, idxcp - 1, GetNodeValue(i, idxcp));
			SetNodeValue(i, idxcp, p_cap);
			if (timelevel == 1 && update)
				SetNodeValue(i, idxS - 1, GetNodeValue(i, idxS));

			// Liquid saturation
			if (continuum_vector.size() > 1) SetNodeValue(i, idx_tS, 0.0);
			//

			saturation = 0., volume_sum = 0.;
			size_t elemsCnode =
			    m_msh->nod_vector[i]->getConnectedElementIDs().size();

			for (size_t j = 0; j < elemsCnode; j++)
			{
				elem = m_msh->ele_vector[m_msh->nod_vector[i]
				                             ->getConnectedElementIDs()[j]];
				m_mmp = mmp_vector[elem->GetPatchIndex()];
				volume_sum += elem->GetVolume();
				saturation +=
				    m_mmp->SaturationCapillaryPressureFunction(p_cap) *
					elem->GetVolume();
			}
			saturation /= volume_sum;
			SetNodeValue(i, idxS, saturation);
		}
	}

	if (continuum > 0)
		for (long i = 0; i < (long)m_msh->GetNodesNumber(false); i++)
		{
			double total_S = 0;
			for (size_t j = 0; j < continuum_vector.size(); j++)
			{
				idxS = GetNodeValueIndex(
				           pcs_secondary_function_name[j * number_continuum]) +
				       timelevel;
				total_S += GetNodeValue(i, idxS) * continuum_vector[j];
			}
			SetNodeValue(i, idx_tS, total_S);
		}
}

/**************************************************************************
   GeoSys - Function: Get mean element value for element index from secondary
node values
                  of process pcs_name and for variable var_name; old and new
timelevel
    01/2006   SB    Implementation
    02/2008   CB    generalization
**************************************************************************/
double PCSGetEleMeanNodeSecondary_2(long index,
                                    int pcsT,
                                    const string& var_name,
                                    int timelevel)
{
	double val =
	    1.0;  // As this returns saturation, default is fully saturated = 1.0;
	int idx = 0, j;  // OK411
	long enode;
	CRFProcess* m_pcs = NULL;
	CRFProcess* cplpcs = NULL;
	CElem* elem = NULL;

	// Get index of secondary node value
	switch (pcsT)
	{
		case 0:  // Liquid_Flow
			break;
		case 1:  // Groundwater Flow
			break;
		case 5:  // Air Flow
			break;
		case 11:  // Componental Flow
			break;
		case 1212:                               // Multiphase Flow
			m_pcs = PCSGet("MULTI_PHASE_FLOW");  // SB, BG
			if (m_pcs)
			{
				idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
				cplpcs = m_pcs;
			}
			break;
		case 12:  // Two_phase_Flow
			m_pcs = PCSGet(FiniteElement::TWO_PHASE_FLOW);
			if (m_pcs)
			{
				if (m_pcs->pcs_type_number == 0)
					cplpcs = pcs_vector[m_pcs->pcs_number + 1];
				else if (m_pcs->pcs_type_number == 1)
					cplpcs = pcs_vector[m_pcs->pcs_number - 1];
				idx = cplpcs->GetNodeValueIndex(var_name) + timelevel;
			}
			break;
		case 22:  // Richards flow
			m_pcs = PCSGet(FiniteElement::RICHARDS_FLOW);
			if (m_pcs)
			{
				idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
				cplpcs = m_pcs;
			}
			break;
		default:
			break;
	}

	if (m_pcs)
	{
		// Get element with index index
		elem = m_pcs->m_msh->ele_vector[index];
		val = 0.0;
		for (j = 0; j < elem->GetVertexNumber();
		     j++)  // average all adjoining nodes
		{
			enode = elem->GetNodeIndex(j);
			val += cplpcs->GetNodeValue(enode, idx);
		}
		val = val / ((double)elem->GetVertexNumber());
	}
	return val;
}

/**************************************************************************
   GeoSys - Function: Get mean element value for element index from secondary
node values
                  of process pcs_name and for variable var_name; old and new
timelevel
    01/2006   SB    Implementation
**************************************************************************/
double PCSGetEleMeanNodeSecondary(long index, const string& pcs_name,
                                  const string& var_name, int timelevel)
{
	double val =
	    1.0;  // As this returns saturation, default is fully saturated = 1.0;
	int idx, j;
	long enode;
	CRFProcess* m_pcs = NULL;

	// Get process by process name
	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(pcs_name));
	m_pcs = PCSGet(pcs_type);
	if (m_pcs)
	{
		// Get index of secondary node value
		idx = m_pcs->GetNodeValueIndex(var_name) + timelevel;
		// Get element with index index
		CElem* elem = NULL;
		elem = m_pcs->m_msh->ele_vector[index];
		val = 0.0;
		// average all adjoining nodes
		for (j = 0; j < elem->GetVertexNumber(); j++)
		{
			enode = elem->GetNodeIndex(j);
			val += m_pcs->GetNodeValue(enode, idx);
		}
		val = val / ((double)elem->GetVertexNumber());
	}
	return val;
}


/*************************************************************************
   GeoSys-FEM Function:
   06/2006 YD Implementation
   02/2008 JOD removed
   03/2008 HS/KG activated for adaptive time step
   Reload primary variable
 **************************************************************************/
void CRFProcess::PrimaryVariableReload()
{
	//  char pcsT;
	//  pcsT = _pcs_type_name[0];
	//  switch(pcsT){
	//    case 'L':
	//      break;
	//    case 'U':
	//      break;
	//    case 'G':
	//      break;
	//    case 'T':
	//      break;
	//    case 'C':
	//      break;
	//    case 'M':
	//      PrimaryVariableReloadTransport();
	//      break;
	//    case 'R': // Richards flow
	//      PrimaryVariableReloadRichards();
	//      break;
	//  }

	switch (this->getProcessType())
	{
		case FiniteElement::MASS_TRANSPORT:
			PrimaryVariableReloadTransport();
			break;
		case FiniteElement::RICHARDS_FLOW:
			PrimaryVariableReloadRichards();
			break;
		default:
			break;
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   06/2006 YD Implementation
   02/2008 JOD removed
   Reload primary variable of Richards Flow
 **************************************************************************/
void CRFProcess::PrimaryVariableReloadRichards()
{
	size_t i;
	int idxp, idx_storage;
	double storage_p;

	idxp = GetNodeValueIndex(pcs_primary_function_name[0]);
	idx_storage = GetNodeValueIndex("STORAGE_P");
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		storage_p = GetNodeValue(i, idx_storage);
		SetNodeValue(i, idxp, storage_p);
		SetNodeValue(i, idxp + 1, storage_p);
	}
	CalcSecondaryVariables(0);
	CalcSecondaryVariables(1);
}

/*************************************************************************
   GeoSys-FEM Function:
   12/2007 kg44 Implementation
   Reload primary variable for Transport
 **************************************************************************/
void CRFProcess::PrimaryVariableReloadTransport()
{
	size_t i;
	int idxp, idx_storage;
	double conc_back;
	char* mcomp_name;

	// kg44 test
	idxp = GetNodeValueIndex(pcs_primary_function_name[0]);

	mcomp_name = new char[80];
	sprintf(mcomp_name, "%s%i", "CONC_BACK_", pcs_component_number);

	idx_storage = GetNodeValueIndex(mcomp_name);

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		conc_back = GetNodeValue(i, idx_storage);
		SetNodeValue(i, idxp, conc_back);
		SetNodeValue(i, idxp + 1, conc_back);
	}

	CalcSecondaryVariables(0);
	CalcSecondaryVariables(1);
}

/*************************************************************************
   GeoSys-FEM Function:
   12/2007 kg44 Implementation
   Reload primary variable for Transport
 **************************************************************************/
void CRFProcess::PrimaryVariableStorageTransport()
{
	size_t i;
	int idxp, idx_storage;
	double concentration;
	char* mcomp_name;

	idxp = GetNodeValueIndex(pcs_primary_function_name[0]);
	mcomp_name = new char[80];
	sprintf(mcomp_name, "%s%i", "CONC_BACK_", pcs_component_number);
	//   cout << "mcomp_name"<< mcomp_name<<endl;
	idx_storage = GetNodeValueIndex(mcomp_name);

	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		concentration = GetNodeValue(i, idxp);
		SetNodeValue(i, idx_storage, concentration);
		//    SetNodeValue(i,idx_storage+1,concentration);
	}
}

/*************************************************************************
   GeoSys-FEM Function:
   06/2006 YD Implementation
   Reload primary variable of Richards Flow
 **************************************************************************/
void CRFProcess::PrimaryVariableStorageRichards()
{
	size_t i;
	int idxp, idx_storage;
	double pressure;

	idxp = GetNodeValueIndex(pcs_primary_function_name[0]) + 1;
	idx_storage = GetNodeValueIndex("STORAGE_P");
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		pressure = GetNodeValue(i, idxp);
		SetNodeValue(i, idx_storage, pressure);
		SetNodeValue(i, idx_storage + 1, pressure);
	}
}

//*************************************************************************
// GeoSys-FEM Function:
// 12/2007 kg44 Implementation
// check change of concentration and set new time step factor
//**************************************************************************/
#ifdef kg44  // WW
double CRFProcess::GetNewTimeStepSizeTransport(double mchange)
{
	size_t i;
	long mnode = -1;
	int comp;
	int idxn, idxo;
	double conc_new, conc_old, /*time_coeff,*/ max_change = 1.0e-10,
	                                           tchange = 1.0;
	char* mcomp_name;

	idxo = GetNodeValueIndex(pcs_primary_function_name[0]);
	comp = pcs_component_number;  // get component number
	//   cout << "comp number "<<comp<<endl;
	mcomp_name = new char[80];
	sprintf(mcomp_name, "%s%li", "CONC_BACK_", comp);
	//   cout << "mcomp_name"<< mcomp_name<<endl;
	idxn = GetNodeValueIndex(mcomp_name);
	for (i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		conc_old = abs(GetNodeValue(i, idxo));
		conc_new = abs(GetNodeValue(i, idxn));
		if (((conc_old) > MKleinsteZahl) && ((conc_new) > MKleinsteZahl))
		{
			max_change = MMax(max_change, conc_new / conc_old);
			mnode = i;
		}
	}
	tchange = mchange / max_change;
	if (tchange > 2.0) tchange = 2.0;
	cout << "Transport: max change of " << max_change << " at node " << mnode
	     << " factor " << tchange << endl;
	return tchange;
}
#endif


/**************************************************************************
   PCSLib-Method:
   04/2006 OK Implementation
**************************************************************************/
void CRFProcess::CreateBCGroup()
{
	cout << "->Create BC" << '\n';
	CBoundaryConditionsGroup* m_bc_group = NULL;

	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
	{
		BCGroupDelete(pcs_type_name, pcs_primary_function_name[i]);
		m_bc_group = new CBoundaryConditionsGroup();
		// OK
		m_bc_group->setProcessTypeName(pcs_type_name);
		// OK
		m_bc_group->setProcessPrimaryVariableName(pcs_primary_function_name[i]);
		m_bc_group->Set(this, Shift[i]);
		bc_group_list.push_back(m_bc_group);
	}
}

/**************************************************************************
   PCSLib-Method:
   04/2006 OK Implementation
**************************************************************************/
void CRFProcess::CreateSTGroup()
{
	cout << "->Create ST" << '\n';
	CSourceTermGroup* m_st_group = NULL;
	// WW
	std::ifstream* iSourceNBC_RHS_file = NULL;
	std::ofstream* oSourceNBC_RHS_file = NULL;

	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));

	if (WriteSourceNBC_RHS == 1)
	{
		string m_file_name =
		    FileName + "_" + pcs_type_name + "_source_Neumann_RHS.bin";
		oSourceNBC_RHS_file = new ofstream(m_file_name.c_str(), ios::binary);
	}
	else if (WriteSourceNBC_RHS == 2)
	{
		string m_file_name =
		    FileName + "_" + pcs_type_name + "_source_Neumann_RHS.bin";
		iSourceNBC_RHS_file = new ifstream(m_file_name.c_str(), ios::binary);
		if (!iSourceNBC_RHS_file->good())
			cout << "_source_Neumann_RHS file is not found" << endl;
	}

	for (size_t i = 0; i < GetPrimaryVNumber(); i++)
	{
		// OK m_st_group = m_st_group->Get(pcs_primary_function_name[i]);
		m_st_group = STGetGroup(pcs_type_name, pcs_primary_function_name[i]);
		if (!m_st_group)
		{
			m_st_group = new CSourceTermGroup();
			// OK
			m_st_group->pcs_type_name = pcs_type_name;
			// OK
			m_st_group->pcs_pv_name = pcs_primary_function_name[i];
			//      if(iSourceNBC_RHS_file)  // Read from data. WW
			//        m_st_group->Read(*iSourceNBC_RHS_file);
			//	  else
			m_st_group->Set(this, Shift[i]);
			st_group_list.push_back(m_st_group);
		}
	}
	if (oSourceNBC_RHS_file)  // WW
		//    WriteRHS_of_ST_NeumannBC(*oSourceNBC_RHS_file);

		if (iSourceNBC_RHS_file)  // WW
		{
			iSourceNBC_RHS_file->close();
			delete iSourceNBC_RHS_file;
			iSourceNBC_RHS_file = NULL;
		}
	if (oSourceNBC_RHS_file)  // WW
	{
		oSourceNBC_RHS_file->close();
		delete oSourceNBC_RHS_file;
		oSourceNBC_RHS_file = NULL;
	}
}

/**************************************************************************
    PCSLib-Method:
    08/2006 OK Implementation
    04/2012 BG Extension to 2 Phases
**************************************************************************/
void CRFProcess::CalcELEFluxes(const GEOLIB::Polyline* const ply,
                               double* result)
{
	int coordinateflag, dimension = 0, axis = 0;
	bool Edge_already_used;
	// bool Node_already_used;
	bool Use_Element;
	vector<size_t> vecConsideredEdges;
	vector<CNode*> vec_nodes_edge;
	bool Edge_on_Geo = false, Point_on_Geo = false;
	std::vector<size_t> nod_vector_at_geo;
	CElem* m_ele = NULL;
	CEdge* m_edg = NULL;
	vec<CEdge*> ele_edges_vector(15);
	double edg_normal_vector[3];
	double vn;
	double edg_length = 0.0;
	double vn_vec[3];
	double edge_vector[3];
	double f_n_sum = 0.0;
	vec<long> element_nodes(20);
	double f[3], v[3];
	double v_2[3], f_2[3];
	double f_n_sum_2 = 0.0;

	result[0] = result[1] = 0;
	FiniteElement::ProcessType pcs_type(getProcessType());

	CRFProcess* m_pcs_flow = NULL;
	//	if (_pcs_type_name.find("FLOW") != string::npos) {
	if (isFlowProcess(this->getProcessType()))
		m_pcs_flow = this;
	else
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);

	// calculates element velocity based on 1 GP
	// CalcELEVelocities();

	int v_eidx[3];
	int v_eidx_2[3];
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");

	if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
	{
		v_eidx_2[0] = m_pcs_flow->GetElementValueIndex("VELOCITY2_X");
		v_eidx_2[1] = m_pcs_flow->GetElementValueIndex("VELOCITY2_Y");
		v_eidx_2[2] = m_pcs_flow->GetElementValueIndex("VELOCITY2_Z");
	}

	for (size_t i = 0; i < 3; i++)
		if (v_eidx[i] < 0)
		{
			std::cout << i << " " << v_eidx[i]
			          << "Velocity output is not specified"
			          << "\n";
			exit(0);  // return 0.0;
		}
	if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
		for (size_t i = 0; i < 3; i++)
			if (v_eidx_2[i] < 0)
			{
				std::cout << i << " " << v_eidx[i]
				          << "Velocity output is not specified"
				          << "\n";
				exit(0);  // return 0.0;
			}

	// determine the dimension and the orientation of the mesh
	coordinateflag = m_msh->GetCoordinateFlag();
	if (coordinateflag == 10)
	{
		dimension = 1;
		axis = 0;
	}  // x only
	else if (coordinateflag == 11)
	{
		dimension = 1;
		axis = 1;
	}  // y only
	else if (coordinateflag == 12)
	{
		dimension = 1;
		axis = 2;
	}  // z only
	else if (coordinateflag == 21)
	{
		dimension = 2;
		axis = 1;
	}  // x, y only
	else if (coordinateflag == 22)
	{
		dimension = 2;
		axis = 2;
	}  // x, z only
	else if (coordinateflag == 32)
	{
		dimension = 3;
		axis = 3;
	}  // x, y, z only

	// Get elements at GEO
	std::vector<size_t> ele_vector_at_geo;
	m_msh->GetELEOnPLY(ply, ele_vector_at_geo, true);
	// BG: 04/2011 nodes are needed to provide the correct edge
	m_msh->GetNODOnPLY(ply, nod_vector_at_geo);

	for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
	{
		m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
		m_ele->SetNormalVector();
		m_ele->GetNodeIndeces(element_nodes);
		Use_Element = true;

		// Configure Element for interpolation of node velocities to GP
		// velocities
		fem->ConfigElement(m_ele);
		// velocity vector
		for (size_t j = 0; j < 3; j++)
		{
			// v[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[j]);
			// Calculate Element velocity
			v[j] =
			    fem->Get_Element_Velocity(m_ele->GetIndex(), m_pcs_flow, 0, j);
		}
		// Test mit Knotengeschwindigkeiten
		// double temp_v[3];
		// temp_v[0] = temp_v[1] = temp_v[2] = 0.0;
		// int variable_index[3];
		// variable_index[0] = m_pcs_flow->GetNodeValueIndex("VELOCITY_X1");
		// variable_index[1] = m_pcs_flow->GetNodeValueIndex("VELOCITY_Y1");
		// variable_index[2] = m_pcs_flow->GetNodeValueIndex("VELOCITY_Z1");
		//
		// for (size_t j = 0; j < 3; j++)
		//{
		//	for (size_t k = 0; k < m_ele->GetNodesNumber(false); k++)
		//	{
		//		temp_v[j] += m_pcs_flow->GetNodeValue(element_nodes[k],
		// variable_index[j]);
		//	}
		//	temp_v[j] /=  m_ele->GetNodesNumber(false);
		//	v[j] = temp_v[j];
		//}

		// BG 04/2011: MassFlux Calculation is not working correctly if it is a
		// 2D mesh in x-z-direction
		// z-velocity is stored a y-position -> this is corrected here
		if ((dimension == 2) && (axis == 2))
		{
			v[2] = v[1];
			v[1] = 0;
		}
		// edge projection // edge marked
		m_ele->GetEdges(ele_edges_vector);
		// cout << "Element: " << endl;
		edg_length = 0;
		// loop over the edges of the element to find the edge at the polyline
		for (size_t j = 0; j < static_cast<size_t>(m_ele->GetEdgesNumber());
		     j++)
		{
			m_edg = ele_edges_vector[j];
			// check if edge was already used
			Edge_already_used = false;
			for (size_t k = 0; k < vecConsideredEdges.size(); k++)
				if (m_edg->GetIndex() == vecConsideredEdges[k])
					Edge_already_used = true;
			if (Edge_already_used == true)
			{
				Edge_on_Geo = false;
				Use_Element = false;
			}
			else
			{
				vec_nodes_edge.clear();
				// BG 04/2011: check if edge is completely on the polyline
				Edge_on_Geo = true;
				vec_nodes_edge.push_back(m_edg->GetNode(0));
				vec_nodes_edge.push_back(m_edg->GetNode(1));
				// loop over the nodes of the edge to check if all of them are
				// on the polyline
				for (int k = 0; k < int(vec_nodes_edge.size()); k++)
				{
					Point_on_Geo = false;
					for (int l = 0; l < int(nod_vector_at_geo.size()); l++)
					{
						if (vec_nodes_edge[k]->GetIndex() ==
						    nod_vector_at_geo[l])
						{
							Point_on_Geo = true;
							l = nod_vector_at_geo.size();
						}
					}
					// cout << "     Node: " << vec_nodes_edge[k]->GetIndex() <<
					// " " << Point_on_Geo << endl;
					if (Point_on_Geo == false) Edge_on_Geo = false;
				}
				if (Edge_on_Geo == true) j = m_ele->GetEdgesNumber();
			}
		}
		if ((m_edg->GetMark()) && (Use_Element == true) &&
		    (Edge_on_Geo == true))
		{
			vecConsideredEdges.push_back(m_edg->GetIndex());  // all edges that
			                                                  // were already
			                                                  // used are stored
			m_edg->SetNormalVector(m_ele->normal_vector, edg_normal_vector);
			edg_length = m_edg->getLength();
			// cout << "Element: " << m_ele->GetIndex() << " LÃ¤nge: " <<
			// edg_length << " Normalvektor: x=" << edg_normal_vector[0] << "
			// y=" << edg_normal_vector[1] << " z=" << edg_normal_vector[2] <<
			// endl;
			m_edg->GetEdgeVector(edge_vector);

			vn = MSkalarprodukt(v, edg_normal_vector, 3);
			for (size_t j = 0; j < 3; j++)
				vn_vec[j] = vn * edg_normal_vector[j];

			for (size_t j = 0; j < 3; j++)
				f[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea();

			f_n_sum += MBtrgVec(f, 3);

			if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)  // BG, 04/2012
			{
				// velocity vector
				for (size_t j = 0; j < 3; j++)
				{
					// v_2[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(),
					// v_eidx_2[j]);
					// Calculate Element velocity
					v_2[j] = fem->Get_Element_Velocity(m_ele->GetIndex(),
					                                   m_pcs_flow, 1, j);
				}

				// variable_index[0] =
				// m_pcs_flow->GetNodeValueIndex("VELOCITY_X2");
				// variable_index[1] =
				// m_pcs_flow->GetNodeValueIndex("VELOCITY_Y2");
				// variable_index[2] =
				// m_pcs_flow->GetNodeValueIndex("VELOCITY_Z2");
				//
				// for (size_t j = 0; j < 3; j++)
				//{
				//	for (size_t k = 0; k < m_ele->GetNodesNumber(false); k++)
				//	{
				//		temp_v[j] += m_pcs_flow->GetNodeValue(element_nodes[k],
				// variable_index[j]);
				//	}
				//	temp_v[j] /=  m_ele->GetNodesNumber(false);
				//	v_2[j] = temp_v[j];
				//}

				// BG 04/2011: MassFlux Calculation is not working correctly if
				// it is a 2D mesh in x-z-direction
				// z-velocity is stored a y-position -> this is corrected here
				if ((dimension == 2) && (axis == 2))
				{
					v_2[2] = v_2[1];
					v_2[1] = 0;
				}

				vn = MSkalarprodukt(v_2, edg_normal_vector, 3);
				for (size_t j = 0; j < 3; j++)
					vn_vec[j] = vn * edg_normal_vector[j];

				for (size_t j = 0; j < 3; j++)
					f_2[j] = vn_vec[j] * edg_length * m_ele->GetFluxArea();

				f_n_sum_2 += MBtrgVec(f_2, 3);
			}
		}
	}
	result[0] = f_n_sum;
	result[1] = f_n_sum_2;
	ele_vector_at_geo.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   11/2004 OK Implementation
   08/2006 OK new
**************************************************************************/
void CRFProcess::CalcELEVelocities(void)
{
	int eidx[3];

	// If not FLUID_MOMENTUM,
	//	if (_pcs_type_name.compare("RANDOM_WALK") != 0) {
	FiniteElement::ProcessType pcs_type(getProcessType());  // BG

	if (this->getProcessType() != FiniteElement::RANDOM_WALK)
	{
		int eidx[3];
		eidx[0] = GetElementValueIndex("VELOCITY1_X");
		eidx[1] = GetElementValueIndex("VELOCITY1_Y");
		eidx[2] = GetElementValueIndex("VELOCITY1_Z");
		for (size_t i = 0; i < 3; i++)
			if (eidx[i] < 0)
				cout << "Fatal error in CRFProcess::CalcELEVelocities - abort"
				     << endl;
		// abort();	// PCH commented abort() out for FM.

		FiniteElement::ElementValue* gp_ele = NULL;
		double vx, vy, vz;
		const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		for (size_t i = 0; i < mesh_ele_vector_size; i++)
		{
			gp_ele = ele_gp_value[i];
			vx = gp_ele->Velocity(0, 0);
			SetElementValue(i, eidx[0], vx);
			SetElementValue(i, eidx[0] + 1, vx);
			vy = gp_ele->Velocity(1, 0);
			SetElementValue(i, eidx[1], vy);
			SetElementValue(i, eidx[1] + 1, vy);
			vz = gp_ele->Velocity(2, 0);
			SetElementValue(i, eidx[2], vz);
			SetElementValue(i, eidx[2] + 1, vz);
		}
	}
	if (pcs_type == FiniteElement::MULTI_PHASE_FLOW)
	{
		eidx[0] = GetElementValueIndex("VELOCITY2_X");
		eidx[1] = GetElementValueIndex("VELOCITY2_Y");
		eidx[2] = GetElementValueIndex("VELOCITY2_Z");
		for (size_t i = 0; i < 3; i++)
			if (eidx[i] < 0)
				cout << "Fatal error in CRFProcess::CalcELEVelocities - abort"
				     << endl;
		// abort();	// PCH commented abort() out for FM.

		FiniteElement::ElementValue* gp_ele = NULL;
		double vx, vy, vz;
		const size_t mesh_ele_vector_size(m_msh->ele_vector.size());
		for (size_t i = 0; i < mesh_ele_vector_size; i++)
		{
			gp_ele = ele_gp_value[i];
			vx = gp_ele->Velocity_g(0, 0);
			SetElementValue(i, eidx[0], vx);
			SetElementValue(i, eidx[0] + 1, vx);
			vy = gp_ele->Velocity_g(1, 0);
			SetElementValue(i, eidx[1], vy);
			SetElementValue(i, eidx[1] + 1, vy);
			vz = gp_ele->Velocity_g(2, 0);
			SetElementValue(i, eidx[2], vz);
			SetElementValue(i, eidx[2] + 1, vz);
		}
	}
}

/**************************************************************************
GeoSys - Function: CalcELEMassFluxes
Task: Calculate the Mass Flux for Elements at Polylines
Return: MassFlux
Programming: 05/2011 BG
Modification:
 **************************************************************************/
void CRFProcess::CalcELEMassFluxes(const GEOLIB::Polyline* const ply,
                                   std::string const& NameofPolyline,
                                   double* result)
{
	CRFProcess* m_pcs_flow = NULL;
	std::vector<size_t> ele_vector_at_geo;
	std::vector<size_t> nod_vector_at_geo;
	vector<CNode*> vec_nodes_edge;
	vec<CEdge*> ele_edges_vector(15);
	vec<long> element_nodes(20);
	vector<size_t> vecConsideredEdges;
	vector<size_t> vecConsideredNodes;
	CElem* m_ele = NULL;
	CEdge* m_edg = NULL;
	CompProperties* m_cp = NULL;
	CMediumProperties* m_mat_mp = NULL;
	CNode* m_node = NULL;
	bool Edge_on_Geo, Point_on_Geo;
	double edg_normal_vector[3];
	double vn;
	double edg_length = 0.0;
	double vn_vec[3];
	double edge_vector[3];
	//	double f_n_sum = 0.0; // TF 2012-08 // removed warning unused variable
	double J_adv_sum = 0.0, J_disp_sum = 0.0, J_diff_sum = 0.0, J_sum = 0.0;
	double C_ele = 0.0;  // C1 = 0, C2 = 0; // TF 2012-08 // removed warning
	                     // unused variable
	//	int number_nodes_at_edge = 0; // TF 2012-08 // removed warning unused
	// variable
	double J_adv[3], J_disp[3], J_diff[3], v[3], temp_j, j_diff[3], j_disp[3],
	    norm_v;
	double J_adv_temp, J_disp_temp, J_diff_temp, J_temp;
	int coordinateflag, dimension = 0, axis = 0;
	bool Edge_already_used, Node_already_used;
	bool Use_Element;
	double ElementConcentration[4];
	double Dm, Dm_eff, tortuosity, porosity, alpha_l, alpha_t, Disp_xx, Disp_yy,
	    Disp_xy;
	// double Disp_zz, Disp_xz, Disp_yz;
	int group;
	//	double totalmass; // 2012-08 TF not used
	//	double totalmassflux; // 2012-08 TF / unused
	double* ConcentrationGradient(new double[3]);
	int numberPolyline = 0;

	if (this->Tim->step_current == 0)
	{
		this->PolylinesforOutput.push_back(NameofPolyline);
		numberPolyline = this->PolylinesforOutput.size() - 1;
		this->TotalMass[numberPolyline] = 0;
	}
	else
	{
		for (int i = 0; i < int(this->PolylinesforOutput.size()); i++)
			if (NameofPolyline == this->PolylinesforOutput[i])
				numberPolyline = i;
	}

	m_msh->GetELEOnPLY(ply, ele_vector_at_geo, true);
	// BG: 04/2011 nodes are needed to provide the correct edge
	m_msh->GetNODOnPLY(ply, nod_vector_at_geo);

	// get the current flow and transport process
	if (isFlowProcess(this->getProcessType()))
	{
		m_pcs_flow = this;
	}
	else
	{
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);
	}

	// get the indices of the velocity
	int v_eidx[3];
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	// check if velocity exists
	for (size_t i = 0; i < 3; i++)
	{
		if (v_eidx[i] < 0)
		{
			std::cout << i << " " << v_eidx[i]
			          << "Velocity output is not specified"
			          << "\n";
			exit(0);
		}
	}
	// determine the dimension and the orientation of the mesh
	coordinateflag = m_msh->GetCoordinateFlag();
	if (coordinateflag == 10)
	{
		dimension = 1;
		axis = 0;
	}  // x only
	else if (coordinateflag == 11)
	{
		dimension = 1;
		axis = 1;
	}  // y only
	else if (coordinateflag == 12)
	{
		dimension = 1;
		axis = 2;
	}  // z only
	else if (coordinateflag == 21)
	{
		dimension = 2;
		axis = 1;
	}  // x, y only
	else if (coordinateflag == 22)
	{
		dimension = 2;
		axis = 2;
	}  // x, z only
	else if (coordinateflag == 32)
	{
		dimension = 3;
		axis = 3;
	}  // x, y, z only

	// loop over all elements at the polyline
	for (size_t i = 0; i < ele_vector_at_geo.size(); i++)
	{
		m_ele = m_msh->ele_vector[ele_vector_at_geo[i]];
		if (m_ele->GetIndex() == 4421) cout << i << endl;
		m_ele->SetNormalVector();
		m_ele->GetNodeIndeces(element_nodes);
		Use_Element = true;
		// velocity vector
		for (size_t j = 0; j < 3; j++)
			v[j] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[j]);

		// BG 04/2011: MassFlux Calculation is not working correctly if it is a
		// 2D mesh in x-z-direction
		// z-velocity is stored a y-position -> this is corrected here
		if ((dimension == 2) && (axis == 2))
		{
			v[2] = v[1];
			v[1] = 0;
		}
		//--------------------------------------------------------------------
		// get properties
		// get porosity of the element
		group = m_ele->GetPatchIndex();
		m_mat_mp = mmp_vector[group];
		porosity = m_mat_mp->Porosity(
		    m_ele->GetIndex(),
		    1);  // CB Now provides also heterogeneous porosity, model 11
		// get diffusion coefficient
		m_cp = cp_vec[this->pcs_component_number];
		Dm = m_cp->CalcDiffusionCoefficientCP(0, 0, this);
		tortuosity = m_mat_mp->TortuosityFunction(m_ele->GetIndex(), 0, 0);
		Dm_eff = Dm * tortuosity * porosity;
		// get element concentration
		C_ele = 0.0;
		for (int j = 0; j < int(m_ele->GetNodesNumber(false)); j++)
		{
			// element concentration = average over all points
			size_t nidx = GetNodeValueIndex(pcs_primary_function_name[0]) + 1;
			ElementConcentration[j] = GetNodeValue(element_nodes[j], nidx);
			// if (ElementConcentration[j] > 0)
			C_ele += ElementConcentration[j];  // BG: 04/2011 sum was not
			                                   // calculated correctly
		}
		C_ele /= (double)m_ele->GetNodesNumber(false);
		// calculate the norm (Betrag) of the velocity vector
		norm_v = MBtrgVec(v, 3);
		for (size_t j = 0; j < 3; j++)
		{
			J_adv[j] = 0;
			J_disp[j] = 0;
			J_diff[j] = 0;
			j_diff[j] = 0;
			j_disp[j] = 0;
		}

		//--------------------------------------------------------------------
		// Calculate Flux depending on element dimension
		if (m_ele->GetDimension() == 2)
		{
			// edge projection // edge marked
			m_ele->GetEdges(ele_edges_vector);
			// cout << "Element: " << endl;
			edg_length = 0;
			// loop over the edges of the element to find the edge at the
			// polyline
			for (size_t j = 0; j < static_cast<size_t>(m_ele->GetEdgesNumber());
			     j++)
			{
				m_edg = ele_edges_vector[j];
				// check if edge was already used
				Edge_already_used = false;
				for (int k = 0; k < int(vecConsideredEdges.size()); k++)
					if (m_edg->GetIndex() == vecConsideredEdges[k])
						Edge_already_used = true;
				if (Edge_already_used == true)
				{
					Edge_on_Geo = false;
					Use_Element = false;
				}
				else
				{
					vec_nodes_edge.clear();
					// BG 04/2011: check if edge is completely on the polyline
					Edge_on_Geo = true;
					vec_nodes_edge.push_back(m_edg->GetNode(0));
					vec_nodes_edge.push_back(m_edg->GetNode(1));
					// loop over the nodes of the edge to check if all of them
					// are on the polyline
					for (int k = 0; k < int(vec_nodes_edge.size()); k++)
					{
						Point_on_Geo = false;
						for (int l = 0; l < int(nod_vector_at_geo.size()); l++)
						{
							if (vec_nodes_edge[k]->GetIndex() ==
							    nod_vector_at_geo[l])
							{
								Point_on_Geo = true;
								l = nod_vector_at_geo.size();
							}
						}
						// cout << "     Node: " <<
						// vec_nodes_edge[k]->GetIndex() << " " << Point_on_Geo
						// << endl;
						if (Point_on_Geo == false) Edge_on_Geo = false;
					}
					if (Edge_on_Geo == true) j = m_ele->GetEdgesNumber();
				}
				if ((m_edg->GetMark()) && (Use_Element == true) &&
				    (Edge_on_Geo == true))
				{
					vecConsideredEdges.push_back(
					    m_edg->GetIndex());  // all edges that were already used
					                         // are stored
					m_edg->SetNormalVector(m_ele->normal_vector,
					                       edg_normal_vector);
					edg_length = m_edg->getLength();
					// cout << "Element: " << m_ele->GetIndex() << " LÃ¤nge: "
					// << edg_length << " Normalvektor: x=" <<
					// edg_normal_vector[0] << " y=" << edg_normal_vector[1] <<
					// " z=" << edg_normal_vector[2] << endl;
					m_edg->GetEdgeVector(edge_vector);

					// cout << i << " " << m_ele->GetIndex() << endl;
					// calculate the velocity vector perpendicular to the edge
					vn = MSkalarprodukt(
					    v, edg_normal_vector,
					    3);  // vn = MSkalarprodukt(v, edg_normal_vector, 3);
					for (size_t j = 0; j < 3; j++)
					{
						vn_vec[j] = vn * edg_normal_vector[j];
					}

					//------------------------------------------------------------------------------
					// calculate advective mass flux = v_n * l^e * z^e * C^e
					for (size_t j = 0; j < 3;
					     j++)  // for (size_t j = 0; j < 3; j++)
					{
						// Darcy velocity times concentration times
						// cross-section area
						J_adv[j] = vn_vec[j] * edg_length *
						           m_ele->GetFluxArea() *
						           C_ele;  // unit: [mol/s]
					}
					// cout << " element: " << m_ele->GetIndex() << " vx: " <<
					// vn_vec[0] << " vy: " << vn_vec[1] << " vz: " << vn_vec[2]
					// << " c1: " << ElementConcentration[0] << " c2: " <<
					// ElementConcentration[1] << " c3: " <<
					// ElementConcentration[2] << " c: " << C_ele << endl;
					// cout << "      J_adv_x: " << J_adv[0] << " J_adv_y: " <<
					// J_adv[1] << " J_adv_z: " << J_adv[2] << endl;

					//------------------------------------------------------------------------------
					// diffusive mass flux
					// calculate concentration gradient for the element
					Calc2DElementGradient(m_ele, ElementConcentration,
					                      ConcentrationGradient);
					// calculate diffusive mass flux
					for (size_t j = 0; j < 3; j++)
						j_diff[j] =
						    Dm_eff * ConcentrationGradient[j];  // [mol/m2.s]
					// calculated the mass flux perpendicular to the edge
					temp_j = MSkalarprodukt(j_diff, edg_normal_vector, 3);
					for (size_t j = 0; j < 3; j++)
					{
						j_diff[j] = temp_j * edg_normal_vector[j];
						J_diff[j] = j_diff[j] * edg_length *
						            m_ele->GetFluxArea();  // [mol/s]
					}
					// cout << "      Dm: " << Dm_eff << " dc_x: " <<
					// ConcentrationGradient[0] << " dc_y: " <<
					// ConcentrationGradient[1] << " J_diff_x: " << J_diff[0] <<
					// " J_diff_y: " << J_diff[1] << " J_diff_z: " << J_diff[2]
					// << endl;

					//------------------------------------------------------------------------------
					// dispersive mass flux
					// get dispersivities
					alpha_l = m_mat_mp->mass_dispersion_longitudinal;
					alpha_t = m_mat_mp->mass_dispersion_transverse;

					if (norm_v != 0)
					{
						Disp_xx = alpha_l * pow(v[0], 2) / norm_v +
						          alpha_t * pow(v[1], 2) /
						              norm_v;  // ToDo: z-Dimension
						Disp_yy = alpha_l * pow(v[1], 2) / norm_v +
						          alpha_t * pow(v[0], 2) /
						              norm_v;  // ToDo: z-Dimension
						Disp_xy = (alpha_l - alpha_t) * v[0] * v[1] /
						          norm_v;  // ToDo: z-Dimension
						// calculate dispersive mass flux
						j_disp[0] =
						    porosity *
						    (Disp_xx * ConcentrationGradient[0] +
						     Disp_xy * ConcentrationGradient[1]);  // [mol/m2.s]
						j_disp[1] =
						    porosity *
						    (Disp_yy * ConcentrationGradient[1] +
						     Disp_xy * ConcentrationGradient[0]);  // [mol/m2.s]
						j_disp[2] = 0;                             // [mol/m2.s]
					}
					else
						for (size_t j = 0; j < 3; j++)
							j_disp[j] = 0;
					// calculated the mass flux perpendicular to the edge
					temp_j = MSkalarprodukt(j_disp, edg_normal_vector, 3);
					for (size_t j = 0; j < 3; j++)
					{
						j_disp[j] = temp_j * edg_normal_vector[j];
						J_disp[j] = j_disp[j] * edg_length *
						            m_ele->GetFluxArea();  // [mol/s]
					}
				}
			}
		}
		if (m_ele->GetDimension() == 1)
		{
			// loop over the nodes of the element to check if the node at the
			// polyline was already used, BG 11/2011
			Use_Element = true;
			for (size_t j = 0;
			     j < static_cast<size_t>(m_ele->GetNodesNumber(false));
			     j++)
			{
				m_node = m_ele->GetNode(j);
				// check if edge was already used
				Node_already_used = false;
				for (int k = 0; k < int(vecConsideredNodes.size()); k++)
					if (m_node->GetIndex() == vecConsideredNodes[k])
						Node_already_used = true;
				if (Node_already_used == true)
				{
					Use_Element = false;
				}
			}

			if (Use_Element == true)
			{
				for (int l = 0; l < int(nod_vector_at_geo.size()); l++)
				{
					if (static_cast<size_t>(m_ele->GetNodeIndex(0)) ==
					    nod_vector_at_geo[l])
					{
						vecConsideredNodes.push_back(m_ele->GetNodeIndex(0));
						l = int(nod_vector_at_geo.size());
					}
					if (static_cast<size_t>(m_ele->GetNodeIndex(1)) ==
					    nod_vector_at_geo[l])
					{
						vecConsideredNodes.push_back(m_ele->GetNodeIndex(1));
						l = int(nod_vector_at_geo.size());
					}
				}
				// calculate advective mass flux = v_n * l^e * z^e * C^e
				// Darcy velocity times concentration times cross-section area
				J_adv[0] =
				    norm_v * m_ele->GetFluxArea() * C_ele;  // unit: [mol/s]
				// cout << " element: " << m_ele->GetIndex() << " vx: " <<
				// vn_vec[0] << " vy: " << vn_vec[1] << " vz: " << vn_vec[2] <<
				// " c1: " << ElementConcentration[0] << " c2: " <<
				// ElementConcentration[1] << " c3: " << ElementConcentration[2]
				// << " c: " << C_ele << endl;
				// cout << "      J_adv_x: " << J_adv[0] << " J_adv_y: " <<
				// J_adv[1] << " J_adv_z: " << J_adv[2] << endl;

				//------------------------------------------------------------------------------
				// diffusive mass flux
				// calculate concentration gradient for the element
				Calc2DElementGradient(m_ele, ElementConcentration,
				                      ConcentrationGradient);
				// calculate diffusive mass flux
				j_diff[0] = Dm_eff * ConcentrationGradient[0];  // [mol/m2.s]
				J_diff[0] = j_diff[0] * m_ele->GetFluxArea();   // [mol/s]
				// cout << "      Dm: " << Dm_eff << " dc_x: " <<
				// ConcentrationGradient[0] << " dc_y: " <<
				// ConcentrationGradient[1] << " J_diff_x: " << J_diff[0] << "
				// J_diff_y: " << J_diff[1] << " J_diff_z: " << J_diff[2] <<
				// endl;

				//------------------------------------------------------------------------------
				// dispersive mass flux
				// get dispersivity along the element direction
				alpha_l = m_mat_mp->mass_dispersion_longitudinal;

				// calculate dispersion coefficients based on velocity of the
				// element
				if (norm_v != 0)
				{
					Disp_xx = alpha_l * norm_v;
					// calculate dispersive mass flux
					j_disp[0] = porosity * Disp_xx *
					            ConcentrationGradient[0];  // [mol/m2.s]
				}
				else
					j_disp[0] = 0;
				// calculated the mass flux perpendicular to the edge
				J_disp[0] = j_disp[0] * m_ele->GetFluxArea();  // [mol/s]
			}
		}
		//--------------------------------------------------------------------
		// sum the mass flux over all considered elements
		J_adv_temp = MBtrgVec(J_adv, 3);
		J_disp_temp = MBtrgVec(J_disp, 3);
		J_diff_temp = MBtrgVec(J_diff, 3);
		J_temp = J_adv_temp + J_disp_temp + J_diff_temp;
		// mass check
		//		totalmassflux = 0.0; // 2012-08 TF not used
		//		totalmass = C_ele * porosity * m_ele->GetVolume(); // 2012-08 TF
		// not used

		// 2012-08 TF, since totalmassflux is not used we can comment out the
		// following two lines
		//		if (this->Tim->time_step_length > 0)
		//			totalmassflux = totalmass / this->Tim->time_step_length;
		// cout << "Element: " << i << " Masse: " << totalmass << " Abstrom_adv:
		// " << this->Tim->this_stepsize * MBtrgVec(J_adv, 3)  << "
		// Abstrom_diff: " << this->Tim->this_stepsize * MBtrgVec(J_diff, 3) <<
		// " Abstrom_disp: " << this->Tim->this_stepsize * MBtrgVec(J_disp, 3);
		// cout << endl;
		// correct mass flux
		// if (J_temp > totalmassflux)
		//{
		//	J_adv_temp += J_adv_temp / J_temp * (totalmassflux - J_temp);
		//	J_disp_temp += J_disp_temp / J_temp * (totalmassflux - J_temp);
		//	J_diff_temp += J_diff_temp / J_temp * (totalmassflux - J_temp);
		//	J_temp = totalmassflux;
		//}
		J_adv_sum += J_adv_temp;
		J_disp_sum += J_disp_temp;
		J_diff_sum += J_diff_temp;
		J_sum += J_temp;
		if (this->Tim->time_step_length > 0)
			this->TotalMass[numberPolyline] +=
			    J_temp * this->Tim->time_step_length;
	}
	result[0] = J_adv_sum;
	result[1] = J_disp_sum;
	result[2] = J_diff_sum;
	result[3] = J_sum;
	result[4] = this->TotalMass[numberPolyline];
	ele_vector_at_geo.clear();

	delete[] ConcentrationGradient;
}

/**************************************************************************
GeoSys - Function: Calc2DElementGradient
Task: Calculate the Gradient for an 2D Element
Return: nothing
Programming: 05/2011 BG
Modification:
 **************************************************************************/
void CRFProcess::Calc2DElementGradient(MeshLib::CElem* m_ele,
                                       double ElementConcentration[4],
                                       double* grad)
{
	double coord_Point1[3];
	double coord_Point2[3];
	double coord_Point3[3];
	double coord_Point4[3];
	double* normal_vector;

	// get the nodes of the 2D Element
	vector<CNode*> vecElementNodes;
	m_ele->GetNodes(vecElementNodes);
	CPlaneEquation* PlaneEquation;
	PlaneEquation = new CPlaneEquation();

	if (m_ele->GetDimension() == 1)
	{
		if (m_ele->GetElementType() == MshElemType::LINE)
		{
			// calculate gradient based on the concentration difference and the
			// length of the element
			grad[0] = (ElementConcentration[1] - ElementConcentration[0]) /
			          (m_ele->GetVolume() / m_ele->GetFluxArea());
			grad[1] = 0;
			grad[2] = 0;
		}
		else
			cout << "This element option is not yet considered for calculating "
			        "the concentration gradient!" << endl;
	}

	else if (m_ele->GetDimension() == 2)
	{
		if (m_ele->GetElementType() == MshElemType::QUAD)  // BG 04/2011:
		                                                   // calculation of the
		                                                   // normal vector of a
		                                                   // quad element
		{
			// define the points used for the plane equation
			// order of points: 1, 2, 3, 4 against clock direction -> points are
			// given in the order 1, 2, 4 to get a positive normal vector
			double const* const p0(vecElementNodes[0]->getData());
			double const* const p1(vecElementNodes[1]->getData());
			double const* const p2(vecElementNodes[2]->getData());
			double const* const p3(vecElementNodes[3]->getData());

			coord_Point1[0] = p0[0];
			coord_Point1[1] = p0[1];
			coord_Point1[2] = ElementConcentration[0];
			coord_Point2[0] = p1[0];
			coord_Point2[1] = p1[1];
			coord_Point2[2] = ElementConcentration[1];
			coord_Point3[0] = p2[0];
			coord_Point3[1] = p2[1];
			coord_Point3[2] = ElementConcentration[2];
			coord_Point4[0] = p3[0];
			coord_Point4[1] = p3[1];
			coord_Point4[2] = ElementConcentration[3];

			// Calculate the plane equation
			PlaneEquation->CalculatePlaneEquationFrom3Points(
			    coord_Point1, coord_Point2, coord_Point4);
			// check if 4. point lies within the plane
			if (PlaneEquation->CheckIfPointInPlane(coord_Point3) == false)
				return;
		}

		if (m_ele->GetElementType() == MshElemType::TRIANGLE)
		{
			// define the points used for the plane equation
			// order of points: 1, 2, 3 against clock direction -> points are
			// given in the order 2, 3, 1 to get a positive normal vector
			double const* const p0(vecElementNodes[0]->getData());
			double const* const p1(vecElementNodes[1]->getData());
			double const* const p2(vecElementNodes[2]->getData());

			coord_Point1[0] = p0[0];
			coord_Point1[1] = p0[1];
			coord_Point1[2] = ElementConcentration[0];
			coord_Point2[0] = p1[0];
			coord_Point2[1] = p1[1];
			coord_Point2[2] = ElementConcentration[1];
			coord_Point3[0] = p2[0];
			coord_Point3[1] = p2[1];
			coord_Point3[2] = ElementConcentration[2];

			// Calculate the plane equation
			PlaneEquation->CalculatePlaneEquationFrom3Points(
			    coord_Point2, coord_Point3, coord_Point1);
		}
		else
			cout << "This element option is not yet considered for calculating "
			        "the concentration gradient!" << endl;

		// calculate gradient
		normal_vector = PlaneEquation->GetNormalVector();
		grad[0] = normal_vector[0] / normal_vector[2];
		grad[1] = normal_vector[1] / normal_vector[2];
		grad[2] = 0;  // grad[2] = normal_vector[2] / normal_vector[2];
	}

	else
		cout << "This element option is not yet considered for calculating the "
		        "concentration gradient!" << endl;
}

/*************************************************************************
   GeoSys-FEM Function:
   08/2006 OK Implementation
 **************************************************************************/
#if !defined(NEW_EQS) && !defined(USE_PETSC)  // WW. 07.11.2008. 04.2012
//(vector<long>&ele_number_vector)
void CRFProcess::AssembleParabolicEquationRHSVector(CNode* m_nod)
{
	// cout << "CRFProcess::AssembleParabolicEquationRHSVector" << endl;
	// int i;
	// WW long ldummy;
	// WW double ddummy;
	//----------------------------------------------------------------------
	// Init
	for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
		eqs->b[m_nod->getConnectedElementIDs()[i]] = 0.0;
	//----------------------------------------------------------------------
	CElem* m_ele = NULL;
	CEdge* m_edg = NULL;
	double edg_normal_vector[3];
	double edge_mid_point[3];
	vec<CEdge*> ele_edges_vector(15);
	int j;
	double aux_vector[3];
	double check_sign;
	//----------------------------------------------------------------------
	// Element velocity
	int v_eidx[3];
	// kg44  v_eidx[0] = GetElementValueIndex("VELOCITY1_X");
	// kg44  v_eidx[1] = GetElementValueIndex("VELOCITY1_Y");
	// kg44  v_eidx[2] = GetElementValueIndex("VELOCITY1_Z");
	CRFProcess* m_pcs_flow = NULL;
	//  if(_pcs_type_name.find("FLOW")!=string::npos) {
	if (isFlowProcess(this->getProcessType()))
		m_pcs_flow = this;
	else
		m_pcs_flow = PCSGet(FiniteElement::GROUNDWATER_FLOW);
	v_eidx[0] = m_pcs_flow->GetElementValueIndex("VELOCITY1_X");
	v_eidx[1] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Y");
	v_eidx[2] = m_pcs_flow->GetElementValueIndex("VELOCITY1_Z");
	for (size_t i = 0; i < 3; i++)
		if (v_eidx[i] < 0)
		{
			cout << v_eidx[i] << i << " Warning in "
			                          "CRFProcess::"
			                          "AssembleParabolicEquationRHSVector - no "
			                          "PCS-VEL data" << endl;
			return;
		}
	double v[3];
	//======================================================================
	// Topology
	for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
	{
		m_ele = m_msh->ele_vector[m_nod->getConnectedElementIDs()[i]];
		m_ele->SetNormalVector();  // OK_BUGFIX
		v[0] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[0]);
		v[1] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[1]);
		v[2] = m_pcs_flow->GetElementValue(m_ele->GetIndex(), v_eidx[2]);
		m_ele->SetMark(false);
		switch (m_ele->GetElementType())
		{
			//------------------------------------------------------------------
			// line elements
			case MshElemType::LINE:
			{
				v[1] = GetElementValue(m_ele->GetIndex(), v_eidx[0]);
				v[0] = GetElementValue(m_ele->GetIndex(), v_eidx[1]);
				if (m_nod->getConnectedElementIDs().size() == 1)
				{
					m_ele->SetMark(true);
					break;
				}
				double const* gravity_center(m_ele->GetGravityCenter());
				double const* const pnt(m_nod->getData());
				aux_vector[0] = gravity_center[0] - pnt[0];
				aux_vector[1] = gravity_center[1] - pnt[1];
				aux_vector[2] = gravity_center[2] - pnt[2];
				check_sign = MSkalarprodukt(v, aux_vector, 3);
				if (check_sign < 0.0) m_ele->SetMark(true);
				break;
			}
			//------------------------------------------------------------------
			// tri elements
			case MshElemType::TRIANGLE:
				m_ele->GetEdges(ele_edges_vector);
				for (j = 0; j < (int)m_ele->GetEdgesNumber(); j++)
				{
					m_edg = ele_edges_vector[j];
					if (m_edg->GetMark())
					{
						m_edg->SetNormalVector(m_ele->normal_vector,
						                       edg_normal_vector);
						break;
						/*
						            m_edg->GetEdgeMidPoint(edge_mid_point);
						           gravity_center = m_ele->GetGravityCenter();
						           aux_vector[0] = gravity_center[0] -
						   edge_mid_point[0];
						           aux_vector[1] = gravity_center[1] -
						   edge_mid_point[1];
						           aux_vector[2] = gravity_center[2] -
						   edge_mid_point[2];
						           check_sign =
						   MSkalarprodukt(edg_normal_vector,aux_vector,3);
						           if(check_sign<0.0) break;
						 */
					}
				}
				if (m_edg->GetMark()) break;
			//----------------------------------------------------------------
			// ToDo
			default:
				cout << "Warning in "
				        "CRFProcess::AssembleParabolicEquationRHSVector - not "
				        "implemented for this element type" << endl;
				break;
		}  // switch
	}
	//======================================================================
	for (size_t i = 0; i < m_nod->getConnectedElementIDs().size(); i++)
	{
		m_ele = m_msh->ele_vector[m_nod->getConnectedElementIDs()[i]];
		switch (m_ele->GetElementType())
		{
			//------------------------------------------------------------------
			// line elements
			case MshElemType::LINE:
				if (m_ele->GetMark())
				{
					cout << m_ele->GetIndex() << endl;
					// WW ldummy = m_nod->GetIndex();
					// WW ddummy = eqs->b[m_nod->GetIndex()];
					fem->ConfigElement(m_ele, false);
					fem->AssembleParabolicEquationRHSVector();
					// WW ddummy = eqs->b[m_nod->GetIndex()];
				}
				break;
			//------------------------------------------------------------------
			// tri elements
			case MshElemType::TRIANGLE:
			{
				m_edg->GetEdgeMidPoint(edge_mid_point);
				double const* gravity_center(m_ele->GetGravityCenter());
				aux_vector[0] = gravity_center[0] - edge_mid_point[0];
				aux_vector[1] = gravity_center[1] - edge_mid_point[1];
				aux_vector[2] = gravity_center[2] - edge_mid_point[2];
				check_sign = MSkalarprodukt(edg_normal_vector, aux_vector, 3);
				if (check_sign < 0.0) continue;
				{
					// cout << m_ele->GetIndex() << endl;
					fem->ConfigElement(m_ele, false);
					fem->AssembleParabolicEquationRHSVector();
				}
				break;
			}
			//----------------------------------------------------------------
			// ToDo
			default:
				cout << "Warning in "
				        "CRFProcess::AssembleParabolicEquationRHSVector - not "
				        "implemented for this element type" << endl;
				break;
		}  // switch
	}
	//======================================================================
}
#endif

/**************************************************************************
PCSLib-Method:
08/2006 OK Implementation
compare with CMCDs PCSGetFluxProcess
03/2012 JT No more loops.
**************************************************************************/
CRFProcess* PCSGetFlow()
{
	if (pcs_number_flow >= 0) return pcs_vector[pcs_number_flow];
	return NULL;
}

/**************************************************************************
PCSLib-Method:
03/2012 JT Implementation
**************************************************************************/
CRFProcess* PCSGetHeat()
{
	if (pcs_number_heat >= 0) return pcs_vector[pcs_number_heat];
	return NULL;
}

/**************************************************************************
PCSLib-Method:
03/2012 JT Implementation
**************************************************************************/
CRFProcess* PCSGetDeformation()
{
	if (pcs_number_deformation >= 0) return pcs_vector[pcs_number_deformation-1];
	return NULL;
}

/**************************************************************************
PCSLib-Method:
03/2012 JT Implementation
**************************************************************************/
CRFProcess* PCSGetMass(size_t component_number)
{
	if (component_number < DOF_NUMBER_MAX)
	{  // don't exceed array dimensions
		if (pcs_number_mass[component_number] >= 0)
		{
			return pcs_vector[pcs_number_mass[component_number]];
		}
	}
	return NULL;
}

/**************************************************************************
   PCSLib-Method:
   based on MMPCalcSecondaryVariables
   01/2007 OK Implementation
**************************************************************************/
void CRFProcess::SetBC()
{
	// WW CBoundaryCondition *m_bc = NULL;
	CBoundaryConditionNode* m_node = NULL;
	int nidx = GetNodeValueIndex(pcs_primary_function_name[0]);
	for (long i = 0; i < (long)bc_node_value.size(); i++)
	{
		m_node = bc_node_value[i];
		// WW m_bc = bc_node[i];
		// old time
		SetNodeValue(m_node->msh_node_number, nidx, m_node->node_value);
		// new time
		SetNodeValue(m_node->msh_node_number, nidx + 1, m_node->node_value);
	}
}

/**************************************************************************
   Task: Postprocessing function calculates the NAPL saturation after
      Flow, Transport and kkinetic NAPL dissolution

   Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void CalcNewNAPLSat(CRFProcess* m_pcs)
{
	long i, j, k, l;
	int idx0, idx1, idx2, idxC;
	long nnodes, nNAPLcomps;
	double conc, rho_N_new, rho_N_old, rho_N_fluid;
	double mass, volume;
	double satu_N_new, satu_N_old;

	vector<int> pcs_napl_comps_vector;
	vector<double> molar_weights_vector;
	vector<double> molar_densities_vector;

	string var_name;
	int no_processes;

	CRFProcess* n_pcs = NULL;
	CFluidProperties* m_mfp = NULL;

	nnodes = (long)fem_msh_vector[0]->nod_vector.size();

	// make sure m_pcs is NAPL-phase
	if (m_pcs->pcs_type_number == 0)
		m_pcs = pcs_vector[m_pcs->pcs_number + 1];  // this is the NAPL phase

	m_mfp = mfp_vector[m_pcs->pcs_type_number];
	// m_mfp->mode = 1; // CB ??

	// Get indices of node value variables for phase 2
	idx0 = m_pcs->GetNodeValueIndex("SATURATION2");  // old timelevel
	idx1 = m_pcs->GetNodeValueIndex("DENSITY2");
	idx2 = m_pcs->GetNodeValueIndex("SATURATION1");  // old timelevel

	i = j = k = l = 0;
	no_processes = (int)pcs_vector.size();

	for (i = 0; i < no_processes; i++)
	{
		n_pcs = pcs_vector[i];
		//     if(n_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0){ // TF
		if (n_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
		{
			j = n_pcs->GetProcessComponentNumber();
			if (cp_vec[j]->transport_phase == 3)  // is in napl
			{
				pcs_napl_comps_vector.push_back(
				    i);  // store processes that are NAPL dissolution components
				         // get the corresponding molar weight
				molar_weights_vector.push_back(cp_vec[j]->molar_weight);
				// get the corresponding densities
				molar_densities_vector.push_back(cp_vec[j]->molar_density);
				l++;
			}
		}
	}
	nNAPLcomps = l;

	for (i = 0; i < nnodes; i++)
	{
		conc = rho_N_new = rho_N_old = rho_N_fluid = satu_N_new = satu_N_old =
		    mass = volume = 0;

		// deteremine the old NAPL DENSITY after flow / transport step
		for (j = 0; j < nNAPLcomps; j++)
		{
			l = pcs_napl_comps_vector[j];
			idxC = pcs_vector[l]->GetNodeValueIndex(
			    pcs_vector[l]->pcs_primary_function_name[0]);
			// old timelevel
			conc = pcs_vector[l]->GetNodeValue(i, idxC);
			mass += conc * molar_weights_vector[j];
		}
		if (mass > 0)
			rho_N_old = mass;  // [kg/m�³REV]
		else
			rho_N_old = 0;
		// determine the new NAPL density rho_N_neu at current node
		mass = 0;
		for (j = 0; j < nNAPLcomps; j++)
		{
			l = pcs_napl_comps_vector[j];
			idxC = pcs_vector[l]->GetNodeValueIndex(
			    pcs_vector[l]->pcs_primary_function_name[0]);
			// +1 new timelevel
			conc = pcs_vector[l]->GetNodeValue(i, idxC + 1);
			if (fabs(conc) < 1e-19) conc = 0.0;
			mass += conc * molar_weights_vector[j];
			// this is required calculating the napl fluid density
			volume += conc / molar_densities_vector[j];
		}
		if (mass > 0)
			rho_N_new = mass;  // [kg/m�³REV]
		else
			rho_N_new = 0;

		// get the old SATURATION2 of NAPL after flow / transport step
		satu_N_old = m_pcs->GetNodeValue(i, idx0);
		// calculate new NAPL Saturation: dSatu = Satu(t+dt)-Satu(t) =
		// Satu(t)*(1-rho(t)/rho(t+dt))
		if (satu_N_old * rho_N_new * rho_N_old > 0)
			// fast
			satu_N_new = satu_N_old + satu_N_old * (1 - rho_N_old / rho_N_new);
		// satu_N_new = satu_N_old + satu_N_old * (rho_N_new / rho_N_old - 1);
		// //slow
		else
			satu_N_new = satu_N_old;
		if (satu_N_new < 0)
		{
			cout << " Warning in fct CalcNewNAPLSat: NAPL-Sat: " << satu_N_new
			     << endl;
			satu_N_new = MRange(0.0, satu_N_new, 1.0);
		}
		// set new SATURATION2
		// if(satu_N_new > 0.00001)  satu_N_new = 1-0.95; // 0.985 CB 14.01.09
		// Hansen+Kueper BM

		// m_pcs->SetNodeValue(i, idx0, satu_N_new);  // idx0 for timelevel 0 ??
		m_pcs->SetNodeValue(i, idx0 + 1,
		                    satu_N_new);  // idx0+1 for timelevel 1 ??
		                                  // set new SATURATION1
		// m_pcs->SetNodeValue(i, idx2, 1-satu_N_new);     // idx2 for timelevel
		// 0 ??
		// idx2+1 for timelevel 1 ??
		m_pcs->SetNodeValue(i, idx2 + 1, 1 - satu_N_new);

		// finally determine the new napl fluid density
		if (mass * volume > 0)
			rho_N_fluid =
			    mass / volume;  // [kg/m�³N] = [kg/m�³REV] / [m�³N/m�³REV]
		else
			rho_N_fluid = m_mfp->Density();  // use NAPL phase fluid density as
		                                     // defined in .mfp
		// set new DENSITY2
		m_pcs->SetNodeValue(i, idx1, rho_N_fluid);
	}

	pcs_napl_comps_vector.clear();
	molar_weights_vector.clear();
	molar_densities_vector.clear();

	// m_pcs->WriteAllVariables();
}

/**************************************************************************
   Task: Calculates the NAPL and Water phase saturations
      in case of NAPL dissolution
   Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
double CalcNAPLDens(CRFProcess* m_pcs, int node)
{
	long i, j, k, l;
	int idxC;  // WWidx1, idxC;
	int nNAPLcomps;
	double conc, rho_N_new, mass, volume;
	string var_name;

	vector<int> pcs_napl_comps_vector;
	vector<double> molar_weights_vector;
	vector<double> molar_densities_vector;

	int no_processes;
	CRFProcess* n_pcs = NULL;
	CFluidProperties* m_mfp = NULL;

	// Get indices of node value variable DENSITY2 for NAPL phase
	if (m_pcs->pcs_type_number == 0) m_pcs = pcs_vector[m_pcs->pcs_number + 1];
	// WW idx1 = m_pcs->GetNodeValueIndex("DENSITY2");

	m_mfp = mfp_vector[m_pcs->pcs_type_number];
	// m_mfp->mode = 1; // CB ??

	i = j = k = l = 0;
	no_processes = (int)pcs_vector.size();

	// collect parameters and concentrations for NAPL-components
	for (i = 0; i < no_processes; i++)
	{
		n_pcs = pcs_vector[i];
		//	  if(n_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0){
		if (n_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
		{
			j = n_pcs->GetProcessComponentNumber();
			if (cp_vec[j]->transport_phase == 3)  // is in NAPL component
			{
				pcs_napl_comps_vector.push_back(
				    i);  // store processes that are NAPL dissolution components
				         // get the corresponding molar weight
				molar_weights_vector.push_back(cp_vec[j]->molar_weight);
				molar_densities_vector.push_back(
				    cp_vec[j]
				        ->molar_density);  // get the corresponding densities
				l++;
			}
		}
	}
	nNAPLcomps = l;

	// determine the NAPL density rho_N_neu at current node
	conc = rho_N_new = mass = volume = 0;
	for (j = 0; j < nNAPLcomps; j++)
	{
		l = pcs_napl_comps_vector[j];
		idxC = pcs_vector[l]->GetNodeValueIndex(
		    pcs_vector[l]->pcs_primary_function_name[0]);
		conc = pcs_vector[l]->GetNodeValue(node, idxC);
		if (fabs(conc) < 1e-19) conc = 0.0;
		mass +=
		    conc *
		    molar_weights_vector[j];  // [kg/m³REV] = [molN/m³REV] * [kg/molN]
		volume += conc / molar_densities_vector[j];  // [m³N/m³REV] =
		                                             // [molN/m³REV] /
		                                             // [molN/m³REV]
	}
	if (mass * volume > 0)
		rho_N_new = mass / volume;  // [kg/m³N] = [kg/m³REV] / [m³N/m³REV]
	else
		rho_N_new =
		    m_mfp
		        ->Density();  // use NAPL phase fluid density as defined in .mfp

	pcs_napl_comps_vector.clear();
	molar_weights_vector.clear();
	molar_densities_vector.clear();

	return rho_N_new;
}

/**************************************************************************
   Task: Preprocessing function set flag pcs->flow_pcs_type for
      calculation of velocity and Saturation in mass transport element
      matrices

   Programing:
   01/2008   CB   Implementation                                          */
/**************************************************************************/
void SetFlowProcessType()
{
	int i;
	int no_processes, flowtype;
	CRFProcess* m_pcs = NULL;

	m_pcs = PCSGetFlow();
	flowtype = m_pcs->type;
	no_processes = (int)pcs_vector.size();

	for (i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		// if(m_pcs->_pcs_type_name.compare("MASS_TRANSPORT")==0)
		// m_pcs->twophaseflow=true;
		m_pcs->flow_pcs_type = flowtype;
	}
}

/**************************************************************************
   Task: Postprocessing function copies the new time step node values of
      secondary variables PRESSURE1 and SATURATION1 to the old time level

   Programing:
   13/2008   CB   Implementation                                          */
/**************************************************************************/
void CopyTimestepNODValuesSVTPhF()
{
	long i, j;
	int idx0, idx1;
	long nnodes;
	double val = 0;
	CRFProcess* m_pcs = NULL;

	nnodes = (long)fem_msh_vector[0]->nod_vector.size();

	for (j = 0; j < 2; j++)  // pcs 1 and 2
	{
		m_pcs = pcs_vector[j];
		//		if (m_pcs->_pcs_type_name.compare("TWO_PHASE_FLOW") != 0)
		if (m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW) break;
		if (j == 0)
			// old timelevel
			idx0 = m_pcs->GetNodeValueIndex("PRESSURE2");
		else
			// old timelevel
			idx0 = m_pcs->GetNodeValueIndex("SATURATION1");
		idx1 = idx0 + 1;
		for (i = 0; i < nnodes; i++)
		{
			val = m_pcs->GetNodeValue(i, idx1);
			m_pcs->SetNodeValue(i, idx0, val);
		}
		// m_pcs->WriteAllVariables();
	}
}

/**************************************************************************
   PCSLib-Method:
   based on WriteSolution by WW
   01/2007 OK Implementation
**************************************************************************/
void CRFProcess::WriteAllVariables()
{
	std::string pcs_type_name(
	    convertProcessTypeToString(this->getProcessType()));
	string m_file_name = FileName + "_" + pcs_type_name + "_" +
	                     pcs_primary_function_name[0] + ".asc";
	ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << endl;
		abort();
	}
	//
	int j;
	int idx[20];
	for (j = 0; j < pcs_number_of_primary_nvals; j++)
	{
		os << pcs_primary_function_name[j] << " ";
		idx[j] = GetNodeValueIndex(pcs_primary_function_name[j]);
		os << pcs_primary_function_name[j] << " ";
		idx[j + pcs_number_of_primary_nvals] = idx[j] + 1;
	}
	if (this->type == 12)  // output of old & new time step for secondary
	                       // variables PRESSURE2 and SATURATION1
	{
		os << pcs_secondary_function_name[0] << " ";
		idx[2 * pcs_number_of_primary_nvals + 0] =
		    GetNodeValueIndex(pcs_secondary_function_name[0]);
		os << pcs_secondary_function_name[0] << " ";
		idx[2 * pcs_number_of_primary_nvals + 1] =
		    GetNodeValueIndex(pcs_secondary_function_name[0]) + 1;
		// other secondary variables
		for (j = 2; j < pcs_number_of_secondary_nvals; j++)
		{
			os << pcs_secondary_function_name[j] << " ";
			idx[2 * pcs_number_of_primary_nvals + j] =
			    GetNodeValueIndex(pcs_secondary_function_name[j]);
		}
	}
	else
		for (j = 0; j < pcs_number_of_secondary_nvals; j++)
		{
			os << pcs_secondary_function_name[j] << " ";
			idx[2 * pcs_number_of_primary_nvals + j] =
			    GetNodeValueIndex(pcs_secondary_function_name[j]);
		}
	os << endl;
	for (size_t i = 0; i < m_msh->GetNodesNumber(false); i++)
	{
		for (j = 0; j < 2 * pcs_number_of_primary_nvals; j++)
			os << GetNodeValue(i, idx[j]) << "  ";
		for (j = 0; j < pcs_number_of_secondary_nvals; j++)
			os << GetNodeValue(i, idx[2 * pcs_number_of_primary_nvals + j])
			   << "  ";
		os << endl;
	}
	os.close();
}

/**************************************************************************
   PCSLib-Method:
   based on MMPCalcSecondaryVariables
   01/2007 OK Implementation
   08/2008 CB NAPLdissolution
**************************************************************************/
void MMPCalcSecondaryVariablesNew(CRFProcess* m_pcs, bool NAPLdiss)
{
	long i;

	//----------------------------------------------------------------------
	int ndx_density_phase;
	int ndx_viscosity_phase;
	ndx_density_phase = -1;         // WW
	ndx_viscosity_phase = -1;       // WW
	CFEMesh* m_msh = m_pcs->m_msh;  // PCH
	//----------------------------------------------------------------------
	m_pcs->SetBC();

	// For accessing the other process
	CRFProcess* cpl_pcs = NULL;
	if (m_pcs->pcs_type_number == 0)
		cpl_pcs = pcs_vector[m_pcs->pcs_number + 1];
	else if (m_pcs->pcs_type_number == 1)
		cpl_pcs = pcs_vector[m_pcs->pcs_number - 1];

	int ndx_pressure1, ndx_p_cap, ndx_pressure2, ndx_s_wetting,
	    ndx_s_nonwetting;
	//======================================================================
	//----------------------------------------------------------------------
	// Capillary pressure - p_c (S) <- This is always the secondary variable
	// in both phase1 and phase2	// PCH
	CMediumProperties* mmp = NULL;

	//======================================================================
	switch (m_pcs->pcs_type_number)
	{
		case 0:
			//..................................................................
			//..................................................................
			// PCH
			// The primary variable is PRESSURE1
			// From PRESSURE1, we are assigning PRESSURE2 which is
			// the secondary variables of PRESSURE1.
			ndx_pressure1 = m_pcs->GetNodeValueIndex("PRESSURE1");
			ndx_pressure2 = m_pcs->GetNodeValueIndex("PRESSURE2");
			ndx_p_cap = m_pcs->GetNodeValueIndex("PRESSURE_CAP");
			double pressure1, pressure2, p_cap;
			for (i = 0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
			{
				// New
				pressure1 = m_pcs->GetNodeValue(i, ndx_pressure1 + 1);
				// Old
				pressure2 = m_pcs->GetNodeValue(i, ndx_pressure1);

				// Let's get capillary pressure before updating pressure2
				// by accessing the primary variable of the saturation equation
				// not the secondary variable of it.
				int cpl_ndx_sat2 = cpl_pcs->GetNodeValueIndex("SATURATION2");
				double cpl_sat2 = cpl_pcs->GetNodeValue(i, cpl_ndx_sat2 + 1);

				if (mmp_vector.size() > 1)
				{
					double sum = 0.0;
					CNode* thisNode = m_msh->nod_vector[i];
					int NumOfNeighborElements =
					    (int)thisNode->getConnectedElementIDs().size();
					// Harmonic mean
					for (int p = 0; p < NumOfNeighborElements; ++p)
					{
						// Mount neighboring elemenets and get the corresponding
						// material group one by one.
						int eleIdx = thisNode->getConnectedElementIDs()[p];
						CElem* thisEle = m_msh->ele_vector[eleIdx];
						int matgrp = thisEle->GetPatchIndex();
						mmp = mmp_vector[matgrp];
						mmp->mode = 2;
						sum += 1.0 / cpl_sat2;
					}
					cpl_sat2 = (double)NumOfNeighborElements / sum;
				}
				// Assigning the secondary variable, Pc
				if (mmp_vector.size() > 1)
					p_cap =
					    m_pcs
					        ->GetCapillaryPressureOnNodeByNeighobringElementPatches(
					            i, 2, 1.0 - cpl_sat2);
				else
					p_cap = mmp->CapillaryPressureFunction(1.0 - cpl_sat2);

				m_pcs->SetNodeValue(i, ndx_p_cap, p_cap);
				m_pcs->SetNodeValue(i, ndx_p_cap + 1, p_cap);

				pressure2 = pressure1 + p_cap;
				// Assigning the secondary variables
				// Previous
				m_pcs->SetNodeValue(i, ndx_pressure2, pressure2);
				// Now
				m_pcs->SetNodeValue(i, ndx_pressure2 + 1, pressure2);
			}
			//......................................................................
			ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY1");
			ndx_viscosity_phase = m_pcs->GetNodeValueIndex("VISCOSITY1");
			printf(
			    "Pressure2 from the known Pressure1 is updated for Process "
			    "0\n");
			break;

		case 1:
			//..................................................................
			// PCH
			// Calc secondary variable saturation Snonwetting = 1-Swetting
			// Don't forget here the primary variable is SATURATION2
			// From SATURATION2, we are assigning SATURATION1 which is
			// the secondary variables of SATURATION2.
			ndx_s_wetting = m_pcs->GetNodeValueIndex("SATURATION1");
			ndx_s_nonwetting = m_pcs->GetNodeValueIndex("SATURATION2");
			ndx_p_cap = cpl_pcs->GetNodeValueIndex("PRESSURE_CAP");

			double s_wetting, s_nonwetting;
			for (i = 0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
			{
				s_nonwetting = m_pcs->GetNodeValue(i, ndx_s_nonwetting + 1);
				// Due to the iterative solution scheme in solving Snw with no
				// explicit boundary condition for non-zero flux condition,
				// Snw may become negative particularly the density difference
				// between two fluids is big. To prevent negative Snw, the
				// saturation restriction added.
				if (mmp_vector.size() > 1)
				{
					double sum = 0.0;
					CNode* thisNode = m_msh->nod_vector[i];
					int NumOfNeighborElements =
					    (int)thisNode->getConnectedElementIDs().size();
					// Harmonic mean
					for (int p = 0; p < NumOfNeighborElements; ++p)
					{
						// Mount neighboring elemenets and get the corresponding
						// material group one by one.
						int eleIdx = thisNode->getConnectedElementIDs()[p];
						CElem* thisEle = m_msh->ele_vector[eleIdx];
						int matgrp = thisEle->GetPatchIndex();
						mmp = mmp_vector[matgrp];
						mmp->mode = 2;
						sum += 1.0 / s_nonwetting;
					}
					s_nonwetting = (double)NumOfNeighborElements / sum;
				}
				// Assigning the secondary variable, Pc
				if (mmp_vector.size() > 1)
					p_cap =
					    m_pcs
					        ->GetCapillaryPressureOnNodeByNeighobringElementPatches(
					            i, 2, 1.0 - s_nonwetting);
				else
					p_cap = mmp->CapillaryPressureFunction(1.0 - s_nonwetting);

				m_pcs->SetNodeValue(i, ndx_s_nonwetting, s_nonwetting);
				m_pcs->SetNodeValue(i, ndx_s_nonwetting + 1, s_nonwetting);
				s_wetting = 1.0 - s_nonwetting;

				// Assigning the secondary variables
				// Previous
				m_pcs->SetNodeValue(i, ndx_s_wetting, s_wetting);
				// Now
				m_pcs->SetNodeValue(i, ndx_s_wetting + 1, s_wetting);

				cpl_pcs->SetNodeValue(i, ndx_p_cap, p_cap);
				cpl_pcs->SetNodeValue(i, ndx_p_cap + 1, p_cap);
			}
			//......................................................................
			ndx_density_phase = m_pcs->GetNodeValueIndex("DENSITY2");
			ndx_viscosity_phase = m_pcs->GetNodeValueIndex("VISCOSITY2");
			printf(
			    "Saturation1 from the known Saturation2 is updated for Process "
			    "1\n");
			break;
	}

	//----------------------------------------------------------------------
	// Fluid properties
	double density;
	double viscosity;
	CFluidProperties* m_mfp = NULL;
	m_mfp = mfp_vector[m_pcs->pcs_type_number];
	m_mfp->mode = 1;
	for (i = 0; i < (long)m_pcs->m_msh->nod_vector.size(); i++)
	{
		// CB NAPL dissolution reqiuires update of Density based on new
		// composition of NAPL phase
		// CB phase 2
		if (NAPLdiss == true && m_pcs->pcs_type_number == 1)
			density = CalcNAPLDens(m_pcs, i);
		else
			density = m_mfp->Density();
		m_pcs->SetNodeValue(i, ndx_density_phase, density);
		viscosity = m_mfp->Viscosity();
		m_pcs->SetNodeValue(i, ndx_viscosity_phase, viscosity);
	}
	m_mfp->mode = 0;
	//----------------------------------------------------------------------
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool CRFProcess::Check()
{
	// MMP
	MeshLib::MSHTestMATGroups();
	return true;
}

/**************************************************************************
   PCSLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool PCSCheck()
{
	if ((int)pcs_vector.size() == 0) return false;
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		// if(m_pcs->m_bCheck)
		if (!m_pcs->Check()) return false;
	}

	return true;
}



/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   09/2007 WW Implementation
 **************************************************************************/
#if defined(NEW_EQS)  // 1.09.2007 WW

void CreateEQS_LinearSolver()
{
	size_t i;
	// CB_merge_0513
	// int dof_DM = 1;
	int dof_DM = 0;    // WW 02.2023. Pardiso
	int DM_type = -1;  // 03.08.2010. WW
	CRFProcess* m_pcs = NULL;
	CFEMesh* a_msh = NULL;
	SparseTable* sp = NULL;
	SparseTable* spH = NULL;
	Linear_EQS* eqs = NULL;
	Linear_EQS* eqs_dof = NULL;  // WW 02.2023. Pardiso
	Linear_EQS* eqsH = NULL;

	bool need_eqs = false;      // WW 02.2023. Pardiso
	bool need_eqs_dof = false;  // WW 02.2023. Pardiso
	int dof = 1;
	//
	// size_t dof_nonDM (1);     //WW 02.2023. Pardiso
	size_t dof_nonDM(0);

	for (i = 0; i < pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		if (m_pcs->type ==
		    1212)  // Important for parallel computing. 24.1.2011 WW
		{
			dof_nonDM = m_pcs->GetPrimaryVNumber();
			dof = dof_nonDM;
		}
		if (m_pcs->type == 4 || m_pcs->type / 10 == 4)  // Deformation
		{
			dof_DM = m_pcs->GetPrimaryVNumber();
			dof = dof_DM;
			DM_type = m_pcs->type;  // 03.08.2010. WW
			if (m_pcs->type == 42) dof = m_pcs->m_msh->GetMaxElementDim();
		}
		else  // Monolithic scheme for the process with linear elements
		{
			// CB_merge_0513
			// if(dof_nonDM < m_pcs->GetPrimaryVNumber()) //WW 02.2023. Pardiso
			//{
			//   dof_nonDM = m_pcs->GetPrimaryVNumber();
			//   // PCH: DOF Handling for FLUID_MOMENTUM in case that the LIS
			//   and PARDISO solvers
			//   // are chosen.
			//   //
			//   if(m_pcs->_pcs_type_name.compare("FLUID_MOMENTUM")==0)
			//   if(m_pcs->getProcessType() == FLUID_MOMENTUM)
			//      dof_nonDM = 1;
			//} //WW 02.2023. Pardiso

			// 02.2013. WW //WW 02.2023. Pardiso
			// Assume that the system with linear element only have one equation
			// with DOF >1;
			if (m_pcs->GetPrimaryVNumber() > 1)
			{
				dof_nonDM = m_pcs->GetPrimaryVNumber();
				dof = dof_nonDM;
				need_eqs_dof = true;
			}
			else
			{
				dof = 1;
				need_eqs = true;
			}  // WW 02.2023. Pardiso
		}
	}

	//
	for (i = 0; i < fem_msh_vector.size(); i++)
	{
		a_msh = fem_msh_vector[i];
		SparseTable* sp = nullptr, *spH = nullptr;
		CreateSparseTable(a_msh, sp, spH);
		if (sp)
			SparseTable_Vector.push_back(sp);
		if (spH)
			SparseTable_Vector.push_back(spH);
		//
		eqs = NULL;
		eqsH = NULL;
		// CB_merge_0513
		eqs_dof = NULL;  // WW 02.2023. Pardiso
		if (sp)          // WW 02.2023. Pardiso
		{
			if (need_eqs)  // 02.2013. WW
			{
				// eqs = new Linear_EQS(*sp, dof_nonDM);//WW 02.2023. Pardiso
				eqs = new Linear_EQS(*sp, 1);
			}
			if (need_eqs_dof)
			{
				eqs_dof = new Linear_EQS(*sp, dof_nonDM);
			}
		}  // WW 02.2023. Pardiso
		if (spH) eqsH = new Linear_EQS(*spH, dof_DM);
		EQS_Vector.push_back(eqs);
		EQS_Vector.push_back(eqsH);
		EQS_Vector.push_back(eqs_dof);  // WW 02.2023. Pardiso
	}
}

#endif  // ifdef NEW_QES

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   01/2008 WW Implementation
 **************************************************************************/
void CRFProcess::WriteBC()
{
	const size_t size_bc(bc_node_value.size());
	const size_t size_st(st_node_value.size());

	if (size_bc == 0 && size_st == 0) return;

	std::string m_file_name =
	    FileName + "_" + convertProcessTypeToString(this->getProcessType()) + "_BC_ST";
#ifdef USE_PETSC
	m_file_name += "_rank" + std::to_string(myrank) + ".asc";
#else
	m_file_name += ".asc";
#endif
	std::ofstream os(m_file_name.c_str(), ios::trunc | ios::out);
	if (!os.good())
	{
		cout << "Failure to open file: " << m_file_name << endl;
		abort();
	}
	os.setf(ios::scientific, ios::floatfield);
	os.precision(12);
	long nindex = 0;
	if (size_bc > 0)
	{
		os << "#Dirchilet BC  (from " << m_file_name << ".bc file) " << endl;
		os << "#Total BC nodes  " << size_bc << endl;
		os << "#Node index, name, x, y, z,   value: " << endl;
		for (size_t i = 0; i < size_bc; i++)
		{
			nindex = bc_node_value[i]->geo_node_number;
			//         anode = m_msh->nod_vector[nindex];
			//         os << nindex << "  " << bc_node_value[i]->pcs_pv_name <<
			//         " "
			//            << std::setw(14) << anode->X() << " " << std::setw(14)
			//            << anode->Y()
			//            << " " << std::setw(14) << anode->Z() << " " <<
			//            std::setw(14)
			//            << bc_node_value[i]->node_value << endl;
			double const* const pnt(m_msh->nod_vector[nindex]->getData());
			os << m_msh->nod_vector[nindex]->GetGlobalIndex();
			os << "  " << bc_node_value[i]->pcs_pv_name << " "
			   << std::setw(14) << pnt[0] << " " << std::setw(14) << pnt[1]
			   << " " << std::setw(14) << pnt[2] << " " << std::setw(14)
			   << bc_node_value[i]->node_value << endl;
		}
	}
	if (size_st > 0)
	{
		os << "#Source term or Neumann BC  (from " << m_file_name
		   << ".st file) " << endl;
		size_t n_st = 0;
		for (size_t i = 0; i < st_node_value.size(); i++)
			n_st += st_node_value[i].size();
		os << "#Total ST nodes  " << n_st << endl;
		os << "#Node index, x, y, z, name    value: " << endl;
		for (size_t i = 0; i < size_st; i++)
		{
			for (size_t ii = 0; ii < st_node_value[i].size(); ii++)
			{
				nindex = st_node_value[i][ii]->geo_node_number;
				//         anode = m_msh->nod_vector[nindex];
				//         os << nindex << "  " <<
				//         convertPrimaryVariableToString(
				//            st_node[i]->getProcessPrimaryVariable()) << " " <<
				//            std::setw(14)
				//            << anode->X() << " " << std::setw(14) <<
				//            anode->Y() << " "
				//            << std::setw(14) << anode->Z() << " " <<
				//            std::setw(14)
				//            << st_node_value[i]->node_value << endl;
				double const* const pnt(m_msh->nod_vector[nindex]->getData());
				os << nindex << "  "
				   << convertPrimaryVariableToString(
				          st_vector[i]->getProcessPrimaryVariable()) << " "
				   << std::setw(14) << pnt[0] << " " << std::setw(14) << pnt[1]
				   << " " << std::setw(14) << pnt[2] << " " << std::setw(14)
				   << st_node_value[i][ii]->node_value << endl;
			}
		}
	}
	os << "#STOP" << endl;
	os.close();
}

/*************************************************************************
   GeoSys-Function:
   Task: PI time step contorl
   Programming:
   08/2008 WW Implementation
   10/2008 WW Node value criteria (test)
   03/2009 WW Euclidean norm
 **************************************************************************/
void CRFProcess::PI_TimeStepSize()
{
	//----------------------------------------------------------------------
	//----------------------------------------------------------------------
	// Time step control
	double hmin;
	double hmax;
	double factor1;  // 1/hmin
	double factor2;  // 1/hmax
	double sfactor = 0.9;
	// WW double reject_factor;                          // BG

	double* u_n = _problem->GetBufferArray();

	double* eqs_x = NULL;
	if (m_num->nls_method == FiniteElement::NL_NEWTON)  // Newton-Raphson
	{
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
		eqs_x = eqs_new->GetGlobalSolution();
#else
#ifdef NEW_EQS
		eqs_x = eqs_new->getX();
#else
		eqs_x = eqs->x;
#endif
#endif
	}

	//
	//
	hmax = Tim->GetMaximumTSizeRestric();
	hmin = Tim->GetMinimumTSizeRestric();
	//
	if (hmin < DBL_MIN)
		factor1 = 5.0;
	else
		factor1 = 1.0 / hmin;
	if (hmax < DBL_MIN)
		factor2 = 0.166666666666666666667e+00;
	else
		factor2 = 1.0 / hmax;
	if (factor1 < 1.0e0) factor1 = 5.0;
	if (factor2 > 1.0e0) factor2 = 0.166666666666666666667e+00;
	//
	hmax = Tim->max_time_step;
	if (hmax < DBL_MIN) hmax = fabs(Tim->time_end - aktuelle_zeit);
	//
	// esitmate the error
	double hnew;
	double err, fac;
	double factorGus;
	double hacc = Tim->GetHacc();
	double erracc = Tim->GetErracc();

	err = 0.0;
	//
	int ii, nidx1;
	long g_nnodes, j, k, l, size_x;
	double x0, x1;
	double Rtol = Tim->GetRTol();
	double Atol = Tim->GetATol();
	double* u_k = _problem->GetBufferArray(true);

	size_x = 0;
	for (ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
		nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
		g_nnodes = m_msh->getNumNodesLocal();
#else
		g_nnodes = m_msh->GetNodesNumber(false);
#endif
		size_x += g_nnodes;

		if (m_num->nls_method == FiniteElement::NL_NEWTON)  // Newton-Raphson
		{
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
				k = j;
				l = pcs_number_of_primary_nvals * j + ii;
#else
				k = m_msh->Eqs2Global_NodeIndex[j];
				l = j + ii * g_nnodes;
#endif
				x0 = u_n[l];
				x1 = GetNodeValue(k, nidx1);
				err += pow((eqs_x[l]) / (Atol + Rtol * max(fabs(x0), fabs(x1))),
				           2);
			}
		}
		else
		{
			for (j = 0; j < g_nnodes; j++)
			{
#if defined(USE_PETSC)  // || defined(other parallel libs)//03.3012. WW
				k = j;
				l = pcs_number_of_primary_nvals * j + ii;
#else
				k = m_msh->Eqs2Global_NodeIndex[j];
				l = j + ii * g_nnodes;
#endif
				x0 = u_n[l];
				x1 = GetNodeValue(k, nidx1);
				err += pow(
				    (x1 - u_k[l]) / (Atol + Rtol * max(fabs(x0), fabs(x1))), 2);
			}
		}
	}

#if defined(USE_PETSC)  // || defined(other parallel libs)//04.3012. WW
	double err_l = err;
	MPI_Allreduce(&err_l, &err, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	long size_xloc = size_x;
	MPI_Allreduce(&size_xloc, &size_x, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
#endif
	err = sqrt(err / (double)size_x);

	//----------------------------------------------------------------------
	//
	// Set the reject factor for the first timestep BG
	// if (Tim->step_current == 1)
	//   Tim->reject_factor = 1;

	// compute hnew with the resriction: 0.2<=hnew/h<=6.0;
	fac = max(factor2, min(factor1, pow(err, 0.25) / sfactor));
	hnew = Tim->time_step_length / fac;  //*Tim->reject_factor;              //
	// BG, use the reject factor to lower
	// timestep after rejection BG

	// determine if the error is small enough
	if (err <= 1.0e0 && accepted)  // step is accept (unless Newton diverged!)
	{
		accept_steps++;
		if (Tim->GetPITimeStepCrtlType() ==
		    2)  // Mod. predictive step size controller (Gustafsson)
		{
			if (accept_steps > 1)
			{
				factorGus = (hacc / Tim->time_step_length) *
				            pow(err * err / erracc, 0.25) / sfactor;
				factorGus = max(factor2, min(factor1, factorGus));
				fac = max(fac, factorGus);
				hnew = Tim->time_step_length / fac;
			}
			hacc = Tim->time_step_length;
			erracc = max(1.0e-2, err);
			Tim->SetHacc(hacc);
			Tim->setErracc(erracc);
		}
		if (fabs(hnew) > hmax) hnew = hmax;
		if (!accepted) hnew = min(fabs(hnew), Tim->time_step_length);
		Tim->SetTimeStep(fabs(hnew));
		// store the used time steps for post-processing BG
		if (Tim->step_current == 1)  // BG
			Tim->time_step_vector.push_back(Tim->time_step_length);
		Tim->time_step_vector.push_back(hnew);
		// WW reject_factor = 1;
		// end of time storage BG
	}  // end if(err<=1.0e0)
	else
	{
		if (!accepted && err <= 1.0e0)
		{  // JT: Then error suggests success, but the iteration diverged.
			if (hnew / Tim->time_step_length >
			    0.99)  // Shock the system to escape the stagnation.
				hnew = Tim->time_step_length * 0.8;
		}
		//
		// WW Tim->reject_factor = 1;                     //BG; if the time step
		// is rejected the next timestep increase is reduced by the reject
		// factor (choose reject factor between 0.1 and 0.9); 1.0 means no
		// change
		reject_steps++;
		accepted = false;
		// WW hnew = hnew / Tim->reject_factor;           //BG
		Tim->SetTimeStep(hnew);
		Tim->time_step_vector.push_back(hnew);  // BG
		if (reject_steps > 100 && accept_steps == 0)
		{
			cout << "!!! More than 100 steps rejected and none of steps "
			        "accepted. Quit the simulation now" << endl;
			exit(1);
		}
		// Recover solutions
	}
}

#ifdef NEW_EQS  // WW
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  //For fluid momentum, WW
   Programming:
   01/2008 WW Implementation
 **************************************************************************/
void CRFProcess::EQSInitialize()
{
	eqs_new->Initialize();
}

#if defined(LIS)  // 11.03.2008 WW
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:  //For fluid momentum,
   Programming:
   02/2008 PCH Implementation
   03/2009 PCH option to tell if this is FLUID_MOMENTUM
 **************************************************************************/
void CRFProcess::EQSSolver(double* x)
{
	eqs_new->Solver(this->m_num);  // NW

	// OK411
	for (int i = 0; i < (int)m_msh->nod_vector.size(); ++i)
		x[i] = eqs_new->X(i);
}
#endif
#endif

#ifdef GEM_REACT
void CRFProcess::IncorporateSourceTerms_GEMS(void)
{
	// Initialization
	long it;       // iterator
	long N_Nodes;  // Number of Nodes
	int nDC;       // Number of mass transport components.
	int i = 0;     // index of the component

	// Get a vector pointing to the REACT_GEM class
	if (m_vec_GEM)
	{
		// Get needed informations.----------------------------
		// Number of Nodes
		N_Nodes = (long)m_msh->GetNodesNumber(false);
		// Number of DC
		nDC = m_vec_GEM->nDC;
		// Identify which PCS it is and its sequence in mcp.
		i = this->pcs_component_number;
		// ----------------------------------------------------

		// Loop over all the nodes-----------------------------
		for (it = 0; it < N_Nodes /*Number of Nodes*/; it++)
		{
// Adding the rate of concentration change to the right hand side of the
// equation.
#ifdef NEW_EQS  // 15.12.2008. WW
			eqs_new->getRHS()[it] -= m_vec_GEM->m_xDC_Chem_delta[it * nDC + i] /
			                  Tim->time_step_length;
#elif defined(USE_PETSC)
// eqs_new->getRHS()[it] -= m_vec_GEM->m_xDC_Chem_delta[it * nDC + i] /
// Tim->time_step_length;
#else
			eqs->b[it] -= m_vec_GEM->m_xDC_Chem_delta[it * nDC + i] /
			              Tim->time_step_length;
#endif
		}
		// ----------------------------------------------------
	}
}
#endif
/*************************************************************************
   ROCKFLOW - Function: CRFProcess::
   Task:
   Programming:
   10/2008 //WW/CB Implementation
   07/2010 TF substituted GEOGetPLYByName
   03/2010  WW Read binary file of precipitation
   06/2010  WW Output top surface flux and head to a raster file
 **************************************************************************/
void CRFProcess::UpdateTransientBC()
{
	//--------------- 24.03.2010. WW
	long i;
	CSourceTerm* precip;
	CSourceTerm* a_st;
	precip = NULL;

	for (i = 0; i < (long)st_vector.size(); i++)
	{
		a_st = st_vector[i];
	}
	// transient boundary condition
	if (bc_transient_index.size() == 0) return;

	bool valid = false;
	long end_i = 0;
	double t_fac = 0.;
	std::vector<double> node_value;

	for (size_t i = 0; i < bc_transient_index.size(); i++)
	{
		std::vector<double> interpolation_points;
		std::vector<double> interpolation_values;

		CBoundaryCondition* bc = bc_node[bc_transient_index[i]];
		long start_i = bc_transient_index[i];
		if (i == bc_transient_index.size() - 1)
			end_i = (long)bc_node.size();
		else
			end_i = bc_transient_index[i + 1];

		// fetch points (representing mesh nodes) along polyline for
		// interpolation
		std::vector<double> nodes_as_interpol_points;
		GEOLIB::Polyline const* ply(
		    static_cast<GEOLIB::Polyline const*>(bc->getGeoObj()));
		m_msh->getPointsForInterpolationAlongPolyline(ply,
		                                              nodes_as_interpol_points);

		valid = false;
		t_fac = 0.0;
		// Piecewise linear distributed.
		for (size_t i(0); i < bc->getDistribedBC().size(); i++)
		{
			for (size_t j = 0; j < ply->getNumberOfPoints(); j++)
				if (bc->getPointsWithDistribedBC()[i] ==
				    (int)ply->getPointID(j))
				{
					if (fabs(bc->getDistribedBC()[i]) < MKleinsteZahl)
						bc->getDistribedBC()[i] = 1.0e-20;
					interpolation_points.push_back(ply->getLength(j));
					interpolation_values.push_back(bc->getDistribedBC()[i]);

					CFunction* fct(FCTGet(bc->getPointsFCTNames()[i]));
					if (fct)
						t_fac = fct->GetValue(aktuelle_zeit, &valid);
					else
						std::cout << "Warning in CBoundaryConditionsGroup - no "
						             "FCT data"
						          << "\n";

					if (valid)
						interpolation_values[interpolation_values.size() - 1] *=
						    t_fac;

					break;
				}
		}
		std::vector<double> interpol_res;
		MathLib::PiecewiseLinearInterpolation(interpolation_points,
		                                      interpolation_values,
		                                      nodes_as_interpol_points,
		                                      interpol_res);

		for (long k = start_i; k < end_i; k++)
			bc_node_value[k]->node_value = interpol_res[k - start_i];
	}
}

/*************************************************************************
   ROCKFLOW - Function: CRFProcess::Add_GEMS_Water_ST
   Task: add one more record into the vector Water_ST_vec
   Programming:
   11/2008 //HS Implementation
 **************************************************************************/

void CRFProcess::Add_GEMS_Water_ST(long idx, double val)
{
	Water_ST_GEMS tmp_st;

	tmp_st.index_node = idx;
	tmp_st.water_st_value = val;

	Water_ST_vec.push_back(tmp_st);
}

void CRFProcess::Clean_Water_ST_vec()
{
	Water_ST_vec.clear();
}

void CRFProcess::configMaterialParameters()
{
	// Output material parameters
	// WW
	const size_t out_vector_size(out_vector.size());

	for (size_t i = 0; i < out_vector_size; i++)
	{
		COutput* out = out_vector[i];
		const size_t size(out->_nod_value_vector.size());
		for (size_t k = 0; k < size; k++)
			if (out->_nod_value_vector[k].find("PERMEABILITY_X1") !=
			    string::npos)
			{
				additioanl2ndvar_print = 1;
				break;
			}
		if (additioanl2ndvar_print == 1) break;
	}

	for (size_t i = 0; i < out_vector_size; i++)
	{
		COutput* out = out_vector[i];
		const size_t size(out->_nod_value_vector.size());
		for (size_t k = 0; k < size; k++)
		{
			if (out->_nod_value_vector[k].find("POROSITY") != string::npos)
			{
				if (additioanl2ndvar_print > 0)
					additioanl2ndvar_print = 2;
				else
					additioanl2ndvar_print = 3;
			}
			if (additioanl2ndvar_print > 1) break;
		}
		if (additioanl2ndvar_print > 1) break;
	}

	if (additioanl2ndvar_print > 0)  // WW
	{
		if (additioanl2ndvar_print < 3)
		{
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
			    "PERMEABILITY_X1";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] =
			    "1/m^2";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++;
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
			    "PERMEABILITY_Y1";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] =
			    "1/m^2";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++;
			if (max_dim == 2)  // 3D
			{
				pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
				    "PERMEABILITY_Z1";
				pcs_secondary_function_unit[pcs_number_of_secondary_nvals] =
				    "1/m^2";
				pcs_secondary_function_timelevel
				    [pcs_number_of_secondary_nvals] = 1;
				pcs_number_of_secondary_nvals++;
			}
		}
		if (additioanl2ndvar_print > 1)  // WW
		{
			pcs_secondary_function_name[pcs_number_of_secondary_nvals] =
			    "POROSITY";
			pcs_secondary_function_unit[pcs_number_of_secondary_nvals] = "-";
			pcs_secondary_function_timelevel[pcs_number_of_secondary_nvals] = 1;
			pcs_number_of_secondary_nvals++;
		}
	}
}

/*************************************************************************
   GeoSys - Function:
   06/2009 OK Implementation
 **************************************************************************/
bool PCSConfig()
{
	bool some_thing_done = false;
	for (int i = 0; i < (int)pcs_vector.size(); i++)  // OK
	{
		pcs_vector[i]->Config();
		some_thing_done = true;
	}
	return some_thing_done;
}
