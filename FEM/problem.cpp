/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "problem.h"

#include <cfloat>
#include <iostream>
#include <sstream>
#include <limits>

#ifdef USE_PETSC
#include <petsctime.h>
#endif

#include "makros.h"
#include "display.h"
#include "FileToolsRF.h"
#include "MemWatch.h"

#ifdef USE_PETSC
#include "PETSC/PETScLinearSolver.h"
#endif

#include "msh_lib.h"
#include "msh_node.h"

#include "ElementValue.h"
#include "fem_ele_std.h"
#include "files0.h"
#include "Output.h"
#include "rfmat_cp.h"
#include "rf_bc_new.h"
#include "rf_fluid_momentum.h"
#include "rf_kinreact.h"
#include "rf_node.h"
#include "rf_out_new.h"
#include "rf_pcs.h"
#include "rf_pcs_dm.h"
#include "rf_pcs_TH.h"
#include "rf_random_walk.h"
#include "rf_react.h"
#ifdef BRNS
#include "rf_REACT_BRNS.h"
#endif
#ifdef GEM_REACT
#include "rf_REACT_GEM.h"
#endif
#include "rf_st_new.h"
#include "rf_tim_new.h"
#include "tools.h"


namespace process
{
class CRFProcessDeformation;
}
using process::CRFProcessDeformation;

/**************************************************************************
   GeoSys - Function: Constructor
   Task:
   Programing:
   07/2008 WW Set it as an constructor of class problem based on the
            PreTimeloop
   Modification:
 ***************************************************************************/
Problem::Problem(char* filename)
    : dt0(0.),
      print_result(true),
      _geo_obj(new GEOLIB::GEOObjects),
      _geo_name(filename)
{
#if defined(USE_MPI) || defined(USE_PETSC)
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);
#endif

	if (filename != NULL)
	{
		// read data
		ReadData(filename, *_geo_obj, _geo_name);
#if defined(USE_MPI) || defined(USE_PETSC)
		MPI_Barrier(MPI_COMM_WORLD);
#endif
	}

	// set the link to Problem instance in CRFProcess objects
	for (size_t i = 0; i < pcs_vector.size(); i++)
		pcs_vector[i]->setProblemObjectPointer(this);

	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		hasAnyProcessDeactivatedSubdomains =
			(pcs_vector[i]->Deactivated_SubDomain.size() > 0);
		if (hasAnyProcessDeactivatedSubdomains) break;
	}
	//----------------------------------------------------------------------
	// Create ST
	// OK STCreateFromPNT();
	//----------------------------------------------------------------------
	GetHeterogeneousFields();  // OK/MB
	//----------------------------------------------------------------------
	// Test MSH-MMP //OK
	int g_max_mmp_groups = MeshLib::MSHSetMaxMMPGroups();

	if (g_max_mmp_groups > (int)mmp_vector.size())
	{
		std::cout << "Error: not enough MMP data";
		print_result = false;  // OK
		return;
	}
	//----------------------------------------------------------------------
	// Create PCS processes
	PCSCreate();
	if (!PCSCheck())  // OK4910 reactivated
	{
		print_result = false;  // OK
		return;
	}
	//
	// JT: Set to true to force node copy at end of loop
	force_post_node_copy = true;
	//
	// JT: Certain restrictions might be made if an external simulator is being
	// used
	external_coupling_exists = false;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		if (pcs_vector[i]->simulator.compare("GEOSYS") != 0)
			external_coupling_exists = true;
	}
#ifdef GEM_REACT
	external_coupling_exists = true;
#endif
#ifdef LIBPHREEQC
	external_coupling_exists = true;
#endif
#ifdef BRNS
	external_coupling_exists = true;
#endif
#ifdef CHEMAPP
	external_coupling_exists = true;
#endif
	//
	//......................................................................
	//#ifdef RESET_4410
	//  //if(pcs_vector[0]->pcs_type_name.compare("TWO_PHASE_FLOW")==0) //OK
	//  if(total_processes[3])  // 3: TWO_PHASE_FLOW. 12.12.2008. WW
	//    PCSCalcSecondaryVariables(); //OK
	//#endif
	//......................................................................
	// 09.07.2008 WW
	SetActiveProcesses();
	// OK if (!Check()) return; //OK
	//----------------------------------------------------------------------
	// REACTIONS
	// CB before the first time step
	// if(MASS_TRANSPORT_Process) //
	// if(MASS_TRANSPORT_Process&&NAPL_Dissolution) //CB Todo
	if (transport_processes.size() > 0)  // 12.12.2008. WW
	{
		// 1) set the id variable flow_pcs_type for Saturation and velocity
		// calculation
		//    in mass transport element matrices
		SetFlowProcessType();
		// 2) in case of Twophaseflow calculate NAPL- and the corresponding
		//    Water phase Saturation from the NAPL concentrations
		// if(pcs_vector[0]->pcs_type_name.compare("TWO_PHASE_FLOW")==0) //OK
		// WW CRFProcess *m_pcs = NULL;
		// WW if (m_pcs = PCSGet("TWO_PHASE_FLOW"))
		// CB: this fct will set the initial NAPL Density in case of
		// NAPL-Dissolution
		//     should this fct be executed in all cases?? --> CHP
		if (total_processes[3])  // 3: TWO_PHASE_FLOW. 12.12.2008. WW
			PCSCalcSecondaryVariables();
	}
//----------------------------------------------------------------------
// REACTIONS
// Initialization of REACT structure for rate exchange between MTM2 and
// Reactions

//--------------------------------------------------
// HB, for the GEM chemical reaction engine 05.2007
//--------------------------------------------------
#ifdef GEM_REACT
	m_vec_GEM = new REACT_GEM();
	GEMRead(FileName, m_vec_GEM);

	string path = "";  // to get the path of the file;
	path = FileName;   // first get full path and project name;
	int pos, npos;
	pos = 0;
	npos = (int)path.size();

// Get path
#ifdef _WIN32
	pos = (int)path.rfind("\\");  // HS keep this on windows
#else
	pos = (int)path.rfind("/");  // HS keep this on linux
#endif  // _WIN32
	if (pos < npos) path = path.substr(0, pos + 1);

	// now start initialization of GEMS
	if (m_vec_GEM->Init_Nodes(path) == 0)
	{
		if (m_vec_GEM->Init_RUN(path) == 0)
		{
			m_vec_GEM->initialized_flag = 1;
		}
		else  // something is wrong and we stop execution
		{
			cout << " GEMS: Error in Init_Nodes..check input " << endl;
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
			MPI_Finalize();  // make sure MPI exits
#endif

			exit(1);
		}
	}
	else  // something is wrong and we stop execution
	{
		cout << " GEMS: Error in Init_RUN..check input " << endl;
#if defined(USE_MPI_GEMS) || defined(USE_PETSC)
		MPI_Finalize();  // make sure MPI exits
#endif

		exit(1);
	}
#else  // GEM_REACT
	//---------------------------------------------------
	REACT* rc = NULL;  // SB
	//  rc->TestPHREEQC(); // Test if *.pqc file is present
	rc = rc->GetREACT();
	if (rc)  // OK
	{
		if (rc->flag_pqc)
		{
			if (cp_vec.size() > 0)
			{                       // OK
#ifdef REACTION_ELEMENT
				rc->CreateREACT();  // SB
				rc->InitREACT0();
				rc->ExecuteReactionsPHREEQC0();
				REACT_vec.clear();
				REACT_vec.push_back(rc);
#else

				rc->CreateREACT();  // SB
				rc->InitREACT();
// SB4501        rc->ExecuteReactions();

#ifdef LIBPHREEQC  // MDL: new functions with built-in phreeqc
				rc->ExecuteReactionsPHREEQCNewLib();
#else
				rc->ExecuteReactionsPHREEQCNew();
#endif  // LIBPHREEQC
				REACT_vec.clear();
				REACT_vec.push_back(rc);
#endif  // REACTION_ELEMENT
			}
		}
		//  delete rc;
	}
#endif  // GEM_REACT

#ifdef BRNS
	// Here to test BRNS; HB 02.11.2007
	// REACT_BRNS* pBRNS;
	// pBRNS = new REACT_BRNS();
	m_vec_BRNS = new REACT_BRNS();
	m_vec_BRNS->InitBRNS(this);
#endif

	//  delete rc;

	//----------------------------------------------------------------------
	// DDC
	size_t no_processes = pcs_vector.size();
	CRFProcess* m_pcs = NULL;
	// WW
	//----------------------------------------------------------------------
	PCSRestart();                        // SB
	if (transport_processes.size() > 0)  // WW. 12.12.2008
	{
		//----------------------------------------------------------------------
		KRConfig(*_geo_obj, _geo_name);
		//----------------------------------------------------------------------
		// Configure Data for Blobs (=>NAPL dissolution)
		KBlobConfig(*_geo_obj, _geo_name);
		KBlobCheck();
		// WW CreateClockTime();
	}
	OUTCheck();  // new SB
	//========================================================================
	// Controls for coupling. WW
	cpl_overall_max_iterations = 1;
	cpl_overall_min_iterations = 1;  // JT2012
	//========================================================================
	// WW
	char line[MAX_ZEILE];
	std::string line_string;
	std::ios::pos_type position;
	std::stringstream in_num;
	// File handling
	std::string num_file_name = FileName + NUM_FILE_EXTENSION;
	std::ifstream num_file(num_file_name.data(), std::ios::in);
	if (num_file.good())
	{
		num_file.seekg(0L, std::ios::beg);
		while (!num_file.eof())
		{
			num_file.getline(line, MAX_ZEILE);
			line_string = line;
			if (line_string.find("#STOP") != std::string::npos) break;
			if (line_string.find("$OVERALL_COUPLING") != std::string::npos)
			{
				in_num.str(GetLineFromFile1(&num_file));
				// JT// in_num >> max_coupling_iterations >> coupling_tolerance;
				// JT: coupling_tolerance is process dependent, cannot be here.
				// See m_num->cpl_tolerance (with $COUPLING_CONTROL)
				in_num >> cpl_overall_min_iterations >>
				    cpl_overall_max_iterations;
				break;
			}
		}
		num_file.close();
	}
	//========================================================================
	// For time stepping. WW
	CTimeDiscretization* m_tim = NULL;
	start_time = 1.0e+25;  // 1.e+8;  kg44 I need a really big time, as I have
	                       // starting times bigger than 3.e+13 (1 Million
	                       // years!!!)
	end_time = 0.;
	max_time_steps = 0;
	bool time_ctr = false;
	// Determine the start and end times from all available process related
	// data.
	for (size_t i = 0; i < time_vector.size(); i++)
	{
		m_tim = time_vector[i];
		m_tim->FillCriticalTime();
		if (m_tim->time_start < start_time) start_time = m_tim->time_start;
		if (m_tim->time_end > end_time) end_time = m_tim->time_end;
		if (max_time_steps < m_tim->time_step_vector.size())
			max_time_steps = m_tim->time_step_vector.size();
		if (m_tim->GetPITimeStepCrtlType() > 0) time_ctr = true;
		m_tim->last_active_time = start_time;  // NW
	}
	if (max_time_steps == 0)
		max_time_steps = std::numeric_limits<std::size_t>::max() -
		                 1;  // ULONG_MAX-1;  //kg44 increased the number to
	                         // maximum number (size_t)
	current_time = start_time;
	if (time_ctr)
	{
		// Modified on 14.02.2011. WW
		long maxi_eqs_dim = 0;
		for (size_t i = 0; i < no_processes; i++)
		{
			m_pcs = pcs_vector[i];
			if (m_pcs->size_unknowns > maxi_eqs_dim)
				maxi_eqs_dim = m_pcs->size_unknowns;
		}
		buffer_array = new double[maxi_eqs_dim];
	}
	else
		buffer_array = NULL;
	buffer_array1 = NULL;
	//========================================================================
	CRFProcessDeformation* dm_pcs = NULL;

	for (size_t i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];
		m_pcs->CalcSecondaryVariables(true);
		m_pcs->Extropolation_MatValue();
	}
	// Calculation of the initial stress and released load for excavation
	// simulation
	// 07.09.2007  WW
	// Excavation for defromation
	dm_pcs = (CRFProcessDeformation*)total_processes[12];
	if (dm_pcs) dm_pcs->CreateInitialState4Excavation();

#ifdef OGS_DELETE_EDGES_AFTER_INIT
	if (!fluid_mom_pcs)
	{
		for (size_t k = 0; k < fem_msh_vector.size(); k++)
		{
			fem_msh_vector[k]->FreeEdgeMemory();
		}
	}
#endif
}

/**************************************************************************
   GeoSys - Function: Desstructor
   Task:
   Programing:
   08/2008 WW Set it as an constructor of class problem based on the
            PreTimeloop

   Modification:
   12.2008  WW
 ***************************************************************************/
Problem::~Problem()
{
	if (_geo_obj) delete _geo_obj;
	delete[] active_processes;
	delete[] exe_flag;
	if (buffer_array) delete[] buffer_array;
	if (buffer_array1) delete[] buffer_array1;
	buffer_array = NULL;
	buffer_array1 = NULL;
	active_processes = NULL;
	exe_flag = NULL;
	//
	PCSDestroyAllProcesses();
	for (size_t i = 0; i < out_vector.size(); i++)
		delete out_vector[i];
	out_vector.clear();

#ifdef GEM_REACT
	delete m_vec_GEM;
#endif

#ifdef BRNS
	delete m_vec_BRNS;
#endif

	ScreenMessage("\nYour simulation is terminated normally\n");
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetActiveProcesses
   Task:
   total_processes:
    0: LIQUID_FLOW     | 1: GROUNDWATER_FLOW  | 2: RICHARDS_FLOW
    3: PS_GLOBAL   | 4: MULTI_PHASE_FLOW  | 5: COMPONENTAL_FLOW
    6: OVERLAND_FLOW   | 7: AIR_FLOW          | 8: HEAT_TRANSPORT
    9: FLUID_MOMENTUM  |10: RANDOM_WALK       |11: MASS_TRANSPORT
   12: DEFORMATION     |13: PTC_FLOW | 14: TH
   Return:
   Programming:
   07/2008 WW
   03/2009 PCH added PS_GLOBAL
   Modification:
   -------------------------------------------------------------------------*/
inline int Problem::AssignProcessIndex(CRFProcess* m_pcs, bool activefunc)
{
	if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
	{
		if (!activefunc) return 0;
		total_processes[0] = m_pcs;
		active_processes[0] = &Problem::LiquidFlow;
		return 0;
	}
	else if (m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
	{
		if (!activefunc) return 1;
		total_processes[1] = m_pcs;
		active_processes[1] = &Problem::GroundWaterFlow;
		return 1;
	}
	else if (m_pcs->getProcessType() == FiniteElement::RICHARDS_FLOW)
	{
		if (!activefunc) return 2;
		total_processes[2] = m_pcs;
		active_processes[2] = &Problem::RichardsFlow;
		return 2;
	}
	else if (m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW)
	{
		if (!activefunc) return 3;
		total_processes[3] = m_pcs;
		active_processes[3] = &Problem::TwoPhaseFlow;
		return 3;
	}
	else if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
	{
		if (!activefunc) return 4;
		total_processes[4] = m_pcs;
		active_processes[4] = &Problem::MultiPhaseFlow;
		return 4;
	}
	else if (m_pcs->getProcessType() == FiniteElement::AIR_FLOW)
	{
		if (!activefunc) return 7;
		total_processes[7] = m_pcs;
		active_processes[7] = &Problem::AirFlow;
		return 7;
	}
	else if (m_pcs->getProcessType() == FiniteElement::HEAT_TRANSPORT)
	{
		if (!activefunc) return 8;
		total_processes[8] = m_pcs;
		active_processes[8] = &Problem::HeatTransport;
		return 8;
	}
	else if (m_pcs->getProcessType() == FiniteElement::FLUID_MOMENTUM)
	{
		if (!activefunc) return 9;
		total_processes[9] = m_pcs;
		active_processes[9] = &Problem::FluidMomentum;
		return 9;
	}
#ifndef OGS_ONLY_TH
	else if (m_pcs->getProcessType() == FiniteElement::RANDOM_WALK)
	{
		if (!activefunc) return 10;
		total_processes[10] = m_pcs;
		active_processes[10] = &Problem::RandomWalker;
		return 10;
	}
#endif
	else if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		if (!activefunc) return 11;
		total_processes[11] = m_pcs;
		active_processes[11] = &Problem::MassTrasport;
		return 11;
	}
	else if (isDeformationProcess(m_pcs->getProcessType()))
	{
		if (!activefunc) return 12;
		total_processes[12] = m_pcs;
		active_processes[12] = &Problem::Deformation;
		return 12;
	}
	else if (m_pcs->getProcessType() == FiniteElement::PS_GLOBAL)
	{
		if (!activefunc) return 3;
		total_processes[3] = m_pcs;
		active_processes[3] = &Problem::PS_Global;
		return 3;
	}
	else if (m_pcs->getProcessType() == FiniteElement::TH_MONOLITHIC)
	{
		if (!activefunc) return 14;
		total_processes[14] = m_pcs;
		active_processes[14] = &Problem::TH_Monolithic;
		return 14;
	}
	std::cout << "Error: no process is specified. "
	          << "\n";
	return -1;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetActiveProcesses
   Task:
   total_processes:
    0: LIQUID_FLOW     | 1: GROUNDWATER_FLOW  | 2: RICHARDS_FLOW
    3: TWO_PHASE_FLOW  | 4: MULTI_PHASE_FLOW  | 5: COMPONENTAL_FLOW
    6: OVERLAND_FLOW   | 7: AIR_FLOW          | 8: HEAT_TRANSPORT
    9: FLUID_MOMENTUM  |10: RANDOM_WALK       |11: MASS_TRANSPORT
   12: DEFORMATION     |13: PS_GLOBAL         |
   Return:
   Programming:
   07/2008 WW
   03/2009 PCH add PS_GLOBAL
   Modification:
   --------------------------------------------------------------------*/
void Problem::SetActiveProcesses()
{
	int i;
	CRFProcess* m_pcs = NULL;
	const int max_processes = 15;  // PCH
	total_processes.resize(max_processes);
	active_processes = new ProblemMemFn[max_processes];
	coupled_process_index.resize(max_processes);
	exe_flag = new bool[max_processes];
	//
	for (i = 0; i < max_processes; i++)
	{
		total_processes[i] = NULL;
		active_processes[i] = NULL;
		coupled_process_index[i] = -1;
	}
	//
	for (i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		AssignProcessIndex(m_pcs);
	}
	//
	for (i = 0; i < max_processes; i++)
		if (total_processes[i])
		{
			// JT: Check for a coupled variable or process (variable is probably
			// necessary for multiple component transport situations.
			// First try for a variable, because this is more general for mass
			// transport
			m_pcs = PCSGet(total_processes[i]->m_num->cpl_variable, true);
			if (m_pcs)
			{
				m_pcs->pcs_is_cpl_underling = true;
				total_processes[i]->pcs_is_cpl_overlord = true;
				//
				coupled_process_index[i] = AssignProcessIndex(m_pcs, false);
				m_pcs->cpl_overlord = total_processes[i];
				total_processes[i]->cpl_underling = m_pcs;
			}
			else
			{  // Try for a process, because it may have been assigned this way
				m_pcs = PCSGet(total_processes[i]->m_num->cpl_process);
				if (m_pcs)
				{
					m_pcs->pcs_is_cpl_underling = true;
					total_processes[i]->pcs_is_cpl_overlord = true;
					//
					coupled_process_index[i] = AssignProcessIndex(m_pcs, false);
					m_pcs->cpl_overlord = total_processes[i];
					total_processes[i]->cpl_underling = m_pcs;
				}
			}
			active_process_index.push_back(i);
		}
	// Transport  porcesses
	for (size_t k = 0; k < pcs_vector.size(); k++)
	{
		m_pcs = pcs_vector[k];
		//		if (m_pcs->pcs_type_name.compare("MASS_TRANSPORT") == 0)
		// TF
		if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			transport_processes.push_back(m_pcs);
		//		if (m_pcs->pcs_type_name.compare("TWO_PHASE_FLOW") == 0)
		////09.01.2008. WW
		// TF
		if ((m_pcs->getProcessType() == FiniteElement::TWO_PHASE_FLOW) ||
		    (m_pcs->getProcessType() ==
		     FiniteElement::MULTI_PHASE_FLOW))  // BG 04/2011
			multiphase_processes.push_back(m_pcs);

		if ((m_pcs->getProcessType() == FiniteElement::GROUNDWATER_FLOW) ||
		    (m_pcs->getProcessType() ==
		     FiniteElement::LIQUID_FLOW))  // BG 04/2011
			singlephaseflow_process.push_back(m_pcs);
	}
}

/************************************************************************** <
   ROCKFLOW - Function: PCSCreate
   Task:
   Programing:
   02/2003 OK Implementation
   03/2003 OK H processes
   04/2003 OK C processes
   05/2003 OK T processes
   05/2003 OK TH processes
   08/2004 OK PCS2
   08/2004 WW The new creation of the deformation process
   10/2004 OK H gas processes
   01/2005 WW New element calculation
   01/2005 OK H unsatutated process
   02/2005 MB switch case in config()
   06/2005 OK MMP2PCSRelation
   07/2008 WW Capsulated into class Problem
   Modification:
   05/2010 TF formated source code
 ***************************************************************************/
void Problem::PCSCreate()
{
	ScreenMessage("---------------------------------------------\n");
	ScreenMessage("Create PCS processes\n");

	size_t no_processes = pcs_vector.size();
	// OK_MOD if(pcs_deformation>0) Init_Linear_Elements();
	for (size_t i = 0; i < no_processes; i++)
	{
		pcs_vector[i]->pcs_type_number = i;
		pcs_vector[i]->Config();  // OK
	}

#ifndef WIN32
	BaseLib::MemWatch mem_watch;
#endif

#ifdef NEW_EQS
	CreateEQS_LinearSolver();  // WW
#endif

	for (size_t i = 0; i < no_processes; i++)
	{
		ScreenMessage(".............................................\n");
		FiniteElement::ProcessType pcs_type(pcs_vector[i]->getProcessType());
		ScreenMessage(
		    "Create: %s\n",
		    FiniteElement::convertProcessTypeToString(pcs_type).c_str());
		//		if (!pcs_vector[i]->pcs_type_name.compare("MASS_TRANSPORT")) {
		// YS   // TF
		if (pcs_type != FiniteElement::MASS_TRANSPORT &&
		    pcs_type != FiniteElement::FLUID_MOMENTUM &&
		    pcs_type != FiniteElement::RANDOM_WALK)
		{
			ScreenMessage(" for %s pcs_component_number %d\n",
			              pcs_vector[i]->pcs_primary_function_name[0],
			              pcs_vector[i]->pcs_component_number);
		}
		ScreenMessage("-> TIM_TYPE: %s\n",
		              FiniteElement::convertTimTypeToString(
		                  pcs_vector[i]->tim_type).c_str());
		pcs_vector[i]->Create();
	}

#if defined(USE_PETSC)  // || defined(other solver libs)//03.3012. WW
	ScreenMessage("---------------------------------------------\n");
	ScreenMessage("Create linear solvers\n");
	CreateEQS_LinearSolver();
#endif

#ifndef WIN32
	ScreenMessage("\tcurrent memory: %d MB\n",
	              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif

	for (size_t i = 0; i < no_processes; i++)
		MMP2PCSRelation(pcs_vector[i]);

	for (size_t i = 0; i < no_processes; i++)  // WW

		pcs_vector[i]->ConfigureCouplingForLocalAssemblier();

	for (size_t i = 0; i < out_vector.size(); i++)
		// initialize process and mesh attributes of COutput objects
		out_vector[i]->init();
}

/*-------------------------------------------------------------------------
   ROCKFLOW - Function: PCSRestart
   Task: Insert process to list
   Programming:
   06/2003 OK Implementation
   07/2008 WW Capsulated into class Problem
   Modification:
   -------------------------------------------------------------------------*/
void Problem::PCSRestart()
{
	const size_t no_processes(pcs_vector.size());
	if (no_processes == 0) return;  // OK41

	int ok = 0;  // = ReadRFRRestartData(file_name_base);

	if (ok == 0)
	{
		ScreenMessage("RFR: no restart data\n");
		return;
	}

	// WW int nidx0; //, nidx1;
	for (size_t i = 0; i < no_processes; i++)
	{
		CRFProcess* m_pcs = pcs_vector[i];
		for (size_t j = 0; j < m_pcs->GetPrimaryVNumber(); j++)
		{
			// timelevel=0;
			// WW nidx0 = m_pcs->GetNodeValueIndex(m_pcs->GetPrimaryVName(j));
			// timelevel= 1;
			// WW nidx1 = nidx0 + 1;
			// OK411      CopyNodeVals(nidx1,nidx0);
		}
	}
}

/*-------------------------------------------------------------------------
   GeoSys - Function: SetTimeActiveProcesses
   Task: For independent time stepping. Set processes to active in time.
   Note: Minimum time step allowance in handled in CalcTimeStep()
   Return:
   Programming:
   03/2012 JT
   Modification:
   --------------------------------------------------------------------*/
void Problem::SetTimeActiveProcesses()
{
	size_t ii;
	double tval, next_time, lowest_next_active;
	CTimeDiscretization* m_tim = NULL;
	CTimeDiscretization* inactive_tim = NULL;
	next_time = current_time + dt;
	//
	lowest_next_active = DBL_MAX;
	for (ii = 0; ii < active_process_index.size(); ii++)
	{
		m_tim = total_processes[active_process_index[ii]]->Tim;
		m_tim->time_active = true;  // activate
		if (m_tim->time_independence && m_tim->next_active_time > next_time)
		{  // Process is then not active this time step
			m_tim->time_active = false;  // deactivate
			//
			// store the lowest next active time of inactive processes
			if (m_tim->next_active_time < lowest_next_active)
			{
				lowest_next_active = m_tim->next_active_time;
				inactive_tim = m_tim;
			}
		}
	}
	//
	// Check if we should shift time step slightly to hit a non-active process
	if (inactive_tim)
	{
		tval = next_time + dt / 1.0e3;  // a small dt increase is better than a
		                                // miniscule dt on the next step
		if (tval > lowest_next_active)
		{  // allow this slight increase
			inactive_tim->time_active = true;
			dt = lowest_next_active - current_time;
			next_time = current_time + dt;
		}
		else if (tval + dt > lowest_next_active)
		{  // Try to smooth to the target time from 2 time steps away
			dt = (lowest_next_active - current_time) / 2.0;
			next_time = current_time + dt;
		}
	}
	//
	// Set times for all active processes
	for (ii = 0; ii < active_process_index.size(); ii++)
	{
		m_tim = total_processes[active_process_index[ii]]->Tim;
		if (m_tim->time_active)
		{
			m_tim->time_step_length = next_time - m_tim->last_active_time;
			m_tim->last_active_time = next_time;
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   07/2008 WW Implementation
   01/2009 WW Update
   03/2012 JT Many changes. Allow independent time stepping.
**************************************************************************/
bool Problem::Euler_TimeDiscretize()
{
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
#endif
	long accepted_times = 0;
	long rejected_times = 0;
	double dt_rec;
	int i;
	bool force_output;
	last_dt_accepted =
	    false;  // JT: false first. Thus copy node values after first dt.
	//
	CTimeDiscretization* m_tim = NULL;
	aktueller_zeitschritt = 0;
	ScreenMessage("\n\n---------------------------------------------\n");
	ScreenMessage("Start time steps\n");

	// Output zero time initial values
	ScreenMessage("-> outputting initial values... \n");
	OUTData(current_time, aktueller_zeitschritt, true);

	// check if this is a steady state simulation
	bool isSteadySimulation = true;
	for (i = 0; i < (int)active_process_index.size(); i++)
	{
		if (total_processes[active_process_index[i]]->tim_type !=
		    FiniteElement::TIM_STEADY)
		{
			isSteadySimulation = false;
			break;
		}
	}

//
// ------------------------------------------
// PERFORM TRANSIENT SIMULATION
// ------------------------------------------
#ifdef USE_PETSC
	PetscLogDouble v1, v2;
#endif
	double previous_rejected_dt = .0;
	;
	while (end_time > current_time)
	{
#if defined(USE_MPI) || defined(USE_PETSC)
		MPI_Barrier(MPI_COMM_WORLD);
#if (PETSC_VERSION_NUMBER >= 3040)
		PetscTime(&v1);
#else
		PetscGetTime(&v1);
#endif
#endif

		// Get time step
		double dt_pre = dt;
		dt = dt_rec = DBL_MAX;
		for (i = 0; i < (int)active_process_index.size(); i++)
		{
			CRFProcess* pcs = total_processes[active_process_index[i]];
			m_tim = pcs->Tim;
			if (!m_tim->time_active) continue;  // JT
			m_tim->dt_pre = dt_pre;
			dt = MMin(dt, m_tim->CalcTimeStep(current_time));
			dt_rec =
			    MMin(dt_rec, m_tim->recommended_time_step);  // to know if we
			                                                 // have a critical
			                                                 // time alteration
			if (!m_tim->repeat && m_tim->isPIDControl())
			{
				pcs->e_pre2 = pcs->e_pre;
				pcs->e_pre = pcs->e_n;
			}
		}

		if (!last_dt_accepted && dt == previous_rejected_dt)
		{
			ScreenMessage(
			    "Stop this simulation. New time step size is same as the "
			    "rejected one.\n");
			break;
		}

		SetTimeActiveProcesses();  // JT2012: Activate or deactivate processes
                                   // with independent time stepping
		//
		// Update time settings
		aktueller_zeitschritt++;
		current_time += dt;
		aktuelle_zeit = current_time;
//
// Print messsage
#if defined(USE_MPI) || defined(USE_PETSC)
		if (myrank == 0)
		{
#endif
			std::cout << "\n\n#################################################"
			             "############\n";
			std::cout << "Time step: " << aktueller_zeitschritt
			          << "|  Time: " << current_time
			          << "|  Time step size: " << dt << "\n";
			if (dt_rec > dt)
			{
				std::cout << "This time step size was modified to match a "
				             "critical time!"
				          << "\n";
			}
#if defined(USE_MPI) || defined(USE_PETSC)
		}
#endif
		if (CouplingLoop())
		{
			// ---------------------------------
			// TIME STEP ACCEPTED
			// ---------------------------------
			last_dt_accepted = true;
			ScreenMessage("This step is accepted.\n");
			PostCouplingLoop();
#ifndef WIN32
			ScreenMessage("\tcurrent mem: %d MB\n",
			              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
			if (print_result)
			{
				ScreenMessage("-> output results\n");
				if (current_time < end_time)
					force_output = false;
				else  // JT: Make sure we printout on last time step
					force_output = true;
				//
				OUTData(current_time, aktueller_zeitschritt, force_output);
			}
			accepted_times++;
			for (i = 0; i < (int)active_process_index.size(); i++)
			{
				CRFProcess* pcs = total_processes[active_process_index[i]];
				m_tim = pcs->Tim;
				if (m_tim->time_active)
				{
					m_tim->accepted_step_count++;
					m_tim->dt_pre = dt;
				}
			}
#ifdef USE_PETSC
#if (PETSC_VERSION_NUMBER >= 3040)
			PetscTime(&v2);
#else
			PetscGetTime(&v2);
#endif
			ScreenMessage("\telapsed time for this time step: %g s\n", v2 - v1);
#endif
#ifndef WIN32
			ScreenMessage("\tcurrent mem: %d MB\n",
			              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
		}
		else if (isSteadySimulation)
		{
			ScreenMessage("This step is rejected.\n");
		}
		else
		{
			// ---------------------------------
			// TIME STEP FAILED
			// ---------------------------------
			last_dt_accepted = false;
			ScreenMessage(
			    "This step is rejected: Redo, with a new time step.\n");
			rejected_times++;
			current_time -= dt;
			aktuelle_zeit = current_time;
			aktueller_zeitschritt--;
			previous_rejected_dt = dt;
			//
			// decrement active dt, and increment count
			for (i = 0; i < (int)active_process_index.size(); i++)
			{
				m_tim = total_processes[active_process_index[i]]->Tim;
				if (!m_tim->time_active) continue;
				m_tim->rejected_step_count++;
				m_tim->last_active_time -= dt;
				m_tim->last_rejected_timestep = aktueller_zeitschritt + 1;
				m_tim->step_current--;
				m_tim->repeat = true;
				//
				// Copy nodal values in reverse
				if (isDeformationProcess(
				        total_processes[active_process_index[i]]
				            ->getProcessType()))
					continue;
				total_processes[active_process_index[i]]->CopyTimestepNODValues(
				    false);
				// JT: This wasn't done before. Is it needed? //
				// total_processes[active_process_index[i]]->CopyTimestepELEValues(false);
			}
			for (i = 0; i < (int)total_processes.size(); i++)
			{
				if (!active_processes[i] && total_processes[i] &&
				    total_processes[i]->tim_type == FiniteElement::TIM_STEADY)
				{
					m_tim = total_processes[i]->Tim;
					m_tim->step_current--;
					m_tim->repeat = true;
				}
			}
		}
		ScreenMessage(
		    "\n#############################################################"
		    "\n");
		if (aktueller_zeitschritt >= max_time_steps) break;

		//		// executing only one time step for profiling
		//		current_time = end_time;
	}
#if defined(USE_MPI) || defined(USE_PETSC)  // JT2012
	if (myrank == 0)
	{
#endif
		std::cout << "\n----------------------------------------------------\n";
		for (i = 0; i < (int)active_process_index.size(); i++)  // JT2012
		{
			m_tim = total_processes[active_process_index[i]]->Tim;
			std::cout << "\nFor process: "
			          << convertProcessTypeToString(
			                 total_processes[active_process_index[i]]
			                     ->getProcessType()) << "\n";
			if (m_tim->time_control_type == TimeControlType::INVALID)
			{
				std::cout << "No time control for this process."
				          << "\n";
			}
			else
			{
				std::cout << "Accepted time steps:                "
				          << m_tim->accepted_step_count << "\n";
				std::cout << "Rejected time steps:                "
				          << m_tim->rejected_step_count << "\n";
			}
			if (total_processes[active_process_index[i]]
			        ->m_num->nls_max_iterations > 1)
			{
				std::cout << "Number of non-converged iterations: "
				          << total_processes[active_process_index[i]]
				                 ->num_notsatisfied << "\n";
				std::cout
				    << "Number of stagnated iterations:     "
				    << total_processes[active_process_index[i]]->num_diverged
				    << "\n";
			}
		}
		std::cout << "\n----------------------------------------------------\n";
#if defined(USE_MPI) || defined(USE_PETSC)  // JT2012
	}
#endif

	return last_dt_accepted;
}

/*-----------------------------------------------------------------------
   GeoSys - Function: Coupling loop
   Task:
   Return: error
   Programming:
   07/2008 WW
   Modification:
   12.2008 WW Update
   03.2012 JT All new. Different strategy, generalized tolerance criteria
   -------------------------------------------------------------------------*/
bool Problem::CouplingLoop()
{
	int i, index, cpl_index;
	double max_outer_error, max_inner_error;  //, error;
	bool run_flag[20];
	int outer_index, inner_index, inner_max;  //, inner_min;
	//
	CRFProcess* a_pcs = NULL;
	CRFProcess* b_pcs = NULL;
	CTimeDiscretization* m_tim = NULL;
	//
	print_result = false;
	int acounter = 0;
	//
	for (i = 0; i < (int)pcs_vector.size(); i++)
		pcs_vector[i]->UpdateTransientBC();

	for (i = 0; i < (int)total_processes.size(); i++)
	{
		if (active_processes[i] && total_processes[i]->selected)
		{
			exe_flag[i] = true;
			m_tim = total_processes[i]->Tim;
			total_processes[i]->SetDefaultTimeStepAccepted();
			acounter++;
			m_tim->step_current++;
			// initilize
			total_processes[i]->iter_nlin_max = 0;
			total_processes[i]->iter_lin_max = 0;
		}
		else
		{  // 21.05.2010.  WW
			if (total_processes[i] &&
			    total_processes[i]->tim_type == FiniteElement::TIM_STEADY)
			{
				acounter++;
				m_tim = total_processes[i]->Tim;
				m_tim->step_current++;  // NW increment needed to get correct
				                        // time step length in
				                        // CTimeDiscretization::CalcTimeStep()
			}
			exe_flag[i] = false;
		}
	}
	int num_processes = (int)active_process_index.size();
	//
	// JT: All active processes must run on the overall loop. Strange this
	// wasn't the case before.
	for (i = 0; i < (int)total_processes.size(); i++)
	{
		run_flag[i] = exe_flag[i];
	}
	for (i = 0; i < num_processes; i++)
	{
		index = active_process_index[i];
		total_processes[index]->first_coupling_iteration = true;
	}
	//
	// To do
	// SB->WW I do not understand this condition, why switch off output?
	// WW Reason:
	/// Make output when all defined processes are activated.
	// JT->WW->SB:  I agree with SB. Just b/c one process is deactivated doesn't
	// mean we don't want output for the others.
	// if(acounter == num_processes)
	print_result = true;
	//
	bool accept = true;
	max_outer_error = 0.0;
	for (outer_index = 0; outer_index < cpl_overall_max_iterations;
	     outer_index++)
	{
		// JT: All active processes must run on the overall loop. Strange this
		// wasn't the case before.
		for (i = 0; i < num_processes; i++)
		{
			index = active_process_index[i];
			run_flag[index] = exe_flag[index];
		}
		for (i = 0; i < num_processes; i++)
		{
			index = active_process_index[i];
			m_tim = total_processes[index]->Tim;
			if (!m_tim->time_active) run_flag[index] = false;
		}

		max_outer_error = 0.0;  // NW reset error for each iteration
		for (i = 0; i < num_processes; i++)
		{
			index = active_process_index[i];
			if (!run_flag[index])
				continue;  // JT: may have been turned off after an inner loop!
			cpl_index = coupled_process_index[index];
			//
			// PERFORM AN INNER COUPLING
			// ---------------------------------------
			if (cpl_index >= 0 && run_flag[cpl_index])
			{
				a_pcs = total_processes[index];
				b_pcs = total_processes[cpl_index];
				//
				inner_max = a_pcs->m_num->cpl_max_iterations;
				//				inner_min = a_pcs->m_num->cpl_min_iterations; //
				// variable set but never used
				//
				a_pcs->iter_outer_cpl = outer_index + 1;
				b_pcs->iter_outer_cpl = outer_index + 1;
				//
				max_inner_error = 0.0;
				for (inner_index = 0;
				     inner_index < a_pcs->m_num->cpl_max_iterations;
				     inner_index++)
				{
					a_pcs->iter_inner_cpl = inner_index + 1;
					b_pcs->iter_inner_cpl = inner_index + 1;
					//
					// FIRST PROCESS
					loop_process_number = i;
					if (a_pcs->first_coupling_iteration) PreCouplingLoop(a_pcs);
					//					error = Call_Member_FN(this,
					// active_processes[index])();
					Call_Member_FN(this, active_processes[index])();
					if (!a_pcs->TimeStepAccept())
					{
						accept = false;
						break;
					}
					//
					// COUPLED PROCESS
					loop_process_number = i + 1;
					if (b_pcs->first_coupling_iteration) PreCouplingLoop(b_pcs);
					//					error = Call_Member_FN(this,
					// active_processes[cpl_index])();
					Call_Member_FN(this, active_processes[cpl_index])();
					if (!b_pcs->TimeStepAccept())
					{
						accept = false;
						break;
					}
					//
					// Check for break criteria
					max_inner_error = MMax(a_pcs->cpl_max_relative_error,
					                       b_pcs->cpl_max_relative_error);
					a_pcs->first_coupling_iteration =
					    false;  // No longer true (JT: these are important, and
					            // are also used elswhere).
					b_pcs->first_coupling_iteration = false;  // No longer true.
					//
					// Store the outer loop error
					if (inner_index == 0)
						max_outer_error =
						    MMax(max_outer_error, max_inner_error);
					//
					std::cout << "\n==========================================="
					             "===========\n";
					std::cout << "Inner coupling loop " << inner_index + 1
					          << "/" << inner_max << " complete."
					          << "\n";
					std::cout << "Max coupling error (relative to tolerance): "
					          << max_inner_error << "\n";
					std::cout << "============================================="
					             "=========\n";
					//
					// Coupling convergence criteria (use loop minimum from
					// a_pcs because this is where the coupled process was
					// called)
					if (max_inner_error <= 1.0 &&
					    inner_index + 2 >
					        a_pcs->m_num->cpl_min_iterations)  // JT: error is
						                                       // relative to
						                                       // the tolerance.
						break;
				}
				run_flag[cpl_index] = false;  // JT: CRUCIAL!!
			}
			else
			{
				// PERFORM AN OUTER COUPLING
				// ---------------------------------------
				a_pcs = total_processes[index];
				a_pcs->iter_outer_cpl = outer_index + 1;
				a_pcs->iter_inner_cpl = 0;
				//
				loop_process_number = i;
				if (a_pcs->first_coupling_iteration) PreCouplingLoop(a_pcs);
				//				error = Call_Member_FN(this,
				// active_processes[index])(); // TF: error set, but never used
				Call_Member_FN(this, active_processes[index])();
				a_pcs->first_coupling_iteration = false;  // No longer true.
				// Check for break criteria
				max_outer_error = MMax(max_outer_error, a_pcs->cpl_max_relative_error);
				if (outer_index>0)
					ScreenMessage("coupling error (relative to tolerance): %e\n", a_pcs->cpl_max_relative_error);
				if (!a_pcs->TimeStepAccept())
				{
					if (a_pcs->tim_type == FiniteElement::TIM_STEADY)
					{
#if 0
						ScreenMessage("***ERR: Cannot change time step size because the process is a steady state\n");
						ScreenMessage("***ERR: Terminate this simulation\n");
#ifdef USE_PETSC
						PetscFinalize();
#endif
						exit(1);
#endif
					}
					else
					{
						ScreenMessage(
						    "***WARN: The process rejected this time step.\n");
						accept = false;
						break;
					}
				}
			}
			if (!accept) break;
		}
		//
		if (cpl_overall_max_iterations > 1)
		{
			ScreenMessage("\n======================================================\n");
			if (outer_index == 0)
				ScreenMessage("Outer coupling loop %d/%d\n", outer_index + 1, cpl_overall_max_iterations);
			else
				ScreenMessage("Outer coupling loop %d/%d: err=%g\n", outer_index + 1, cpl_overall_max_iterations, max_outer_error);
			ScreenMessage("======================================================\n");
		}
		else
		{
			ScreenMessage("\n");
		}
		if (!accept)
		{
			ScreenMessage("\n");
			break;
		}
		// Coupling convergence criteria
		if (max_outer_error <= 1.0 &&
		    outer_index + 2 > cpl_overall_min_iterations)  // JT: error is
			                                               // relative to the
			                                               // tolerance.
			break;
	}
	if (max_outer_error > 1.0 && outer_index == cpl_overall_max_iterations &&
	    cpl_overall_max_iterations > 1)
		accept = false;
	//
	return accept;
}

/*-----------------------------------------------------------------------
   GeoSys - Function: pre Coupling loop
   Task: Process solution is beginning. Perform any pre-loop configurations
   Programming:
   03/2012 JT
   Modification:
-------------------------------------------------------------------------*/
void Problem::PreCouplingLoop(CRFProcess* m_pcs)
{
	if (!last_dt_accepted || force_post_node_copy)  // if last time step not
		                                            // accepted or values were
		                                            // already copied.
		return;
	//
	/*For mass transport this routine is only called once (for the overall
	  transport process)
	  and so we need to copy for all transport components*/
	if (m_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
	{
		CRFProcess* c_pcs = NULL;
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			c_pcs = pcs_vector[i];
			if (c_pcs->getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				c_pcs->CopyTimestepNODValues();
				c_pcs->CopyTimestepELEValues();
			}
		}
	}
	else
	{  // Otherwise, just copy this process
		m_pcs->CopyTimestepNODValues();
		m_pcs->CopyTimestepELEValues();
	}
}

/*-----------------------------------------------------------------------
   GeoSys - Function: post Coupling loop
   Task:
   Return: error
   Programming:
   08/2008 WW
   Modification:
   -------------------------------------------------------------------------*/
void Problem::PostCouplingLoop()
{
	CRFProcess* m_pcs = NULL;
	if (total_processes[12])
	{
		CRFProcessDeformation* dm_pcs =
		    (CRFProcessDeformation*)(total_processes[12]);
		if (pcs_vector.size() > 1)
			dm_pcs->ResetTimeStep();
		dm_pcs->Extropolation_GaussValue();
	}

	// CB new NAPL and Water Saturations after reactions for Two_Phase_Flow and
	// NAPL-Dissolution
	// WW if(MASS_TRANSPORT_Process) //
	// if(MASS_TRANSPORT_Process&&NAPL_Dissolution) //CB Todo
	// 12.2008. WW
	if (transport_processes.size() > 0 && total_processes[3])
		if (KNaplDissCheck())  // Check if NAPLdissolution is modeled
			CalcNewNAPLSat(total_processes[3]);

	/* CB 21/09 The next fct. was necessary in 4.08. Still needed here? I think
	   so
	   // for TWO_PHASE_FLOW the new time step results for secondary variables
	   // PRESSURE2 and SATURATION1 are not copied below in the function
	   // CopyTimestepNODValues(); but I can do it here:
	   if (m_pcs = PCSGet("TWO_PHASE_FLOW"))
	   CopyTimestepNODValuesSVTPhF();
	 */

	//  Update the results
	for (int i = 0; i < (int)pcs_vector.size(); i++)
	{
		m_pcs = pcs_vector[i];
		if (hasAnyProcessDeactivatedSubdomains)  // NW
			m_pcs->CheckMarkedElement();
		m_pcs->WriteSolution();  // WW
#ifdef GEM_REACT
		if (i == 0)  // for GEM_REACT we also need information on porosity
					 // (node porosity internally stored in Gems
					 // process)!....do it only once and it does not matter
					 // for which process ! ....we assume that the first pcs
					 // process is the flow process...if reload not defined
					 // for every process, restarting with gems will not
					 // work in any case

			if ((m_pcs->reload == 1 || m_pcs->reload == 3) &&
				!((aktueller_zeitschritt % m_pcs->nwrite_restart) > 0))
				m_vec_GEM->WriteReloadGem();
#endif

		m_pcs->Extropolation_MatValue();         // WW
		if (m_pcs->cal_integration_point_value)  // WW
			m_pcs->Extropolation_GaussValue();
		// JT: Now done in PreCouplingLoop() // m_pcs->CopyTimestepNODValues();
		// //MB
		if (force_post_node_copy)
		{  // JT: safety valve. Set this value to true (in Problem()) and values
			// will be copied here.
			m_pcs->CopyTimestepNODValues();
			m_pcs->CopyTimestepELEValues();
			for (size_t ie = 0; ie < ele_gp_value.size(); ie++)
			{
				ElementValue* ev = ele_gp_value[ie];
				ev->Velocity0 = ev->Velocity;
			}
		}
	}
	LOPCalcELEResultants();
}

const GEOLIB::GEOObjects* Problem::getGeoObj() const
{
	return _geo_obj;
}

const std::string& Problem::getGeoObjName() const
{
	return _geo_name;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: LiquidFlow
   Task: Similate liquid flow
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::LiquidFlow()
{
	double error = 0.;
	CRFProcess* m_pcs = total_processes[0];
	if (m_pcs->tim_type == FiniteElement::TIM_STEADY &&
	    aktueller_zeitschritt > 1)
		m_pcs->selected = false;
	if (!m_pcs->selected) return error;
	//  error = m_pcs->Execute();
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	PCSCalcSecondaryVariables();
	m_pcs->CalIntegrationPointValue();

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: RichardsFlow
   Task: Similate Richards flow
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::RichardsFlow()
{
	double error = 0.;
	CRFProcess* m_pcs = total_processes[2];
	if (!m_pcs->selected) return error;
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept())
	{
		// WW
		m_pcs->CalcSecondaryVariablesUnsaturatedFlow();
		CalcVelocities = true;
		conducted = true;  // WW
	}
	if (m_pcs->TimeStepAccept()) m_pcs->CalIntegrationPointValue();  // WW
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: TwoPhaseFlow
   Task: Similate twp-phase flow
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Update
   01.2008 WW Add phases
   -------------------------------------------------------------------------*/
inline double Problem::TwoPhaseFlow()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[3];
	if (!m_pcs->selected) return error;  // 12.12.2008 WW
	//
	// 08.01.2009. WW
	for (int i = 0; i < (int)multiphase_processes.size(); i++)
	{
		m_pcs = multiphase_processes[i];
		error = m_pcs->ExecuteNonLinear(loop_process_number);
		if (m_pcs->TimeStepAccept())
		{
			PCSCalcSecondaryVariables();
			m_pcs->CalIntegrationPointValue();
			// CB 12/09 (first time added on 010808) Velocity at
			// CenterOfGravity, required for NAPL dissolution
			if (i == 0)  // is 0 in all cases the correct index?
				m_pcs->CalcELEVelocities();
		}
	}
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: MultiPhaseFlow()
   Task: Similate multi-phase flow by p-p scheme
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Update
   -------------------------------------------------------------------------*/
inline double Problem::MultiPhaseFlow()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[4];
	if (!m_pcs->selected) return error;  // 12.12.2008 WW
	// initialize density and viscosity if the CO2 phase transition is used
	if (m_pcs->Phase_Transition_Model == 1)
	{
		//		if (m_pcs->Tim->step_current == 1)
		//		{
		//			std::cout << " The Viscosity is not calculated yet!!!" <<
		//"\n";
		//			m_pcs->CalculateFluidDensitiesAndViscositiesAtNodes(m_pcs);
		//		}
		//		else
		//			m_pcs->Phase_Transition_CO2(m_pcs, 1);
	}

	// m_pcs->CalculateFluidDensitiesAndViscositiesAtNodes(m_pcs);

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept()) m_pcs->CalIntegrationPointValue();

	// CO2-Phase_Transition BG, NB
	if (m_pcs->Phase_Transition_Model == 1)
	{
		// check if mfp-model for density and viscosity is 18
		if (m_pcs->Tim->step_current == 1)
		{
			CFluidProperties* FluidProp;

			FluidProp = MFPGet("LIQUID");
			if ((FluidProp->density_model != 18) ||
			    (FluidProp->viscosity_model != 18))
			{
				std::cout << "If the Phase_Transition_Model is used the "
				             "density model and the viscosity model should be "
				             "18!"
				          << "\n";
				std::cout << "The run is terminated now ..."
				          << "\n";
				// system("Pause");
				exit(0);
			}
			FluidProp = MFPGet("GAS");
			if ((FluidProp->density_model != 18) ||
			    (FluidProp->viscosity_model != 18))
			{
				std::cout << "If the Phase_Transition_Model is used the "
				             "density model and the viscosity model should be "
				             "18!"
				          << "\n";
				std::cout << "The run is terminated now ..."
				          << "\n";
				// system("Pause");
				exit(0);
			}
		}
		//		if (m_pcs->Phase_Transition_Model == 1)
		//		{
		//			m_pcs->CalculateFluidDensitiesAndViscositiesAtNodes(m_pcs);
		//			m_pcs->Phase_Transition_CO2(m_pcs, 0);
		//		}
	}

	if (m_pcs->tim_type == FiniteElement::TIM_STEADY) m_pcs->selected = false;


	if (m_pcs->OutputMassOfGasInModel == true)  // 05/2012 BG
		OutputMassOfGasInModel(m_pcs);

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: OutputMassOfGasInModel()
   Task: provides the sum of CO2 in the model domain in output files
   Return: nothing
   Programming:
   02/2011 BG
   -------------------------------------------------------------------------*/
void Problem::OutputMassOfGasInModel(CRFProcess* m_pcs)
{
	// Testoutput amount of co2 in model domain calculated at nodes
	CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CMediumProperties* m_mat_mp = NULL;
	std::ostringstream temp;
	double mass_Gas_gas, mass_Gas_water, mass_Gas;
	int index;
	double saturation_Gas;
	double saturation_water;
	double node_volume;
	double time;
	std::string tempstring;
	std::vector<std::string> vec_string;
	// int position;
	std::string path;
	double density_Gas;
	double porosity = 0.0;
	int variable_index;
	double concentration_Gas_water;
	int indexConcentration_Gas = -1;
	int ProcessIndex_GasInLiquid;
	CRFProcess* n_pcs = NULL;
	int group;
	std::string transport_process_name, Processname;
	double Molweight_Gas;
	double V_model;
	// bool Windows_System;
	int position;
	std::string filename;

	path = m_pcs->file_name_base;

#ifdef _WIN32
	{
		position = int(path.find_last_of("\\"));
		if (position < 0)
			filename = "Sum_Gas_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "\\Sum_Gas_nodes.csv";
		}
	}
#else
	{
		position = int(path.find_last_of("/"));
		if (position < 0)
			filename = "Sum_Gas_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "/Sum_Gas_nodes.csv";
		}
	}
#endif

	transport_process_name = "H2";
	Molweight_Gas = 2;
	Processname = "PS_GLOBAL";

	ProcessIndex_GasInLiquid = -1;
	if (ProcessIndex_GasInLiquid == -1)
		for (int i = 0; i < int(pcs_vector.size()); i++)
		{
			n_pcs = pcs_vector[i];
			// identify your process and store idx of pcs-vector
			if (n_pcs->nod_val_name_vector[0] == transport_process_name)
				ProcessIndex_GasInLiquid = i;
		}

	if (m_pcs->Tim->step_current == 1)
		// Header of the file
		vec_string.push_back(
		    "Time, massGas_gas, massGas_water, massGas, porosity, Modelvolume, "
		    "Pressure2, Density2");
#ifndef OGS_ONLY_TH
	else
	{
		//		CReadTextfiles_ECL TextFile;
		//		TextFile.Read_Text(filename);
		//
		//		for (int i = 0; i < TextFile.NumberOfRows; i++)
		//			vec_string.push_back(TextFile.Data[i]);
	}
#endif

	V_model = mass_Gas = mass_Gas_gas = mass_Gas_water = 0;
	// +1: new timelevel
	if (ProcessIndex_GasInLiquid > -1)
		indexConcentration_Gas =
		    pcs_vector[ProcessIndex_GasInLiquid]->GetNodeValueIndex(
		        pcs_vector[ProcessIndex_GasInLiquid]
		            ->pcs_primary_function_name[0]) +
		    1;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i];  // get element
		node_volume = 0;
		saturation_Gas = 0;
		if (mfp_vector[1]->density_model == 18)
		{
			variable_index = m_pcs->GetNodeValueIndex("DENSITY2");
			density_Gas = m_pcs->GetNodeValue(i, variable_index);
		}
		else
		{
			variable_index = m_pcs->GetNodeValueIndex("PRESSURE2");
			double p2 = m_pcs->GetNodeValue(i, variable_index);
			double dens_arg[3];
			dens_arg[0] = p2;
			density_Gas = mfp_vector[1]->Density(dens_arg);
			// cout << "Knoten: " << i << " Dichte: " << density_Gas << " P2: "
			// << p2 << endl;
		}

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			// CB Now provides also heterogeneous porosity, model 11
			porosity = m_mat_mp->Porosity(m_ele->GetIndex(), 1);
			// node_volume = node_volume + m_ele->GetVolume() / 8 * porosity;
			if (m_ele->GetElementType() == 2)
				node_volume =
				    node_volume +
				    m_ele->GetVolume() / m_ele->GetNodesNumber(true) * porosity;
			else
				node_volume = node_volume +
				              m_ele->GetVolume() /
				                  m_ele->GetNodesNumber(false) * porosity;
		}

		V_model = V_model + node_volume;

		//+1... new time level
		if (Processname == "PS_GLOBAL")
			index = m_pcs->GetNodeValueIndex("SATURATION1");
		else
			index = m_pcs->GetNodeValueIndex("SATURATION1") + 1;
		saturation_water = m_pcs->GetNodeValue(i, index);
		// if (saturation_water < 1)
		saturation_Gas = 1 - saturation_water;
		if (ProcessIndex_GasInLiquid > -1)
		{
			concentration_Gas_water =
			    pcs_vector[ProcessIndex_GasInLiquid]->GetNodeValue(
			        i, indexConcentration_Gas);
			mass_Gas_water = mass_Gas_water +
			                 node_volume * saturation_water *
			                     concentration_Gas_water * Molweight_Gas *
			                     0.001;
		}
		else
			mass_Gas_water = 0;

		mass_Gas_gas =
		    mass_Gas_gas + node_volume * saturation_Gas * density_Gas;
		// cout << " Node: " << i << " saturation: " << saturation_water << "
		// Density_CO2: " << density_CO2 << " node_volume: " << node_volume <<
		// endl;
	}
	mass_Gas = mass_Gas_gas + mass_Gas_water;
	// calculating time
	time = 0;
	for (int k = 0; k < m_pcs->Tim->step_current; k++)
		time += m_pcs->Tim->time_step_vector[k];
	temp.str("");
	temp.clear();
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_Gas_gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_Gas_water;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << mass_Gas;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << porosity;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << V_model;
	tempstring += temp.str();

	// Test output pressure and density
	variable_index = m_pcs->GetNodeValueIndex("PRESSURE2");
	double p2 = m_pcs->GetNodeValue(1, variable_index);
	double dens_arg[3];
	dens_arg[0] = p2;
	density_Gas = mfp_vector[1]->Density(dens_arg);
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << p2;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(12);
	temp << density_Gas;
	tempstring += temp.str();
	// Test output pressure and density

	vec_string.push_back(tempstring);

#ifndef OGS_ONLY_TH
//	//within the first timestep create file and write header
//	CWriteTextfiles_ECL TextFile;
//	TextFile.Write_Text(filename, vec_string);
#endif
}

/*-------------------------------------------------------------------------
GeoSys - Function: OutputMassOfComponentInModel()
Task: provides the sum of each transport component in the model domain in output
files
Return: nothing
Programming:
05/2011 BG
-------------------------------------------------------------------------*/
void Problem::OutputMassOfComponentInModel(std::vector<CRFProcess*> flow_pcs,
                                           CRFProcess* transport_pcs)
{
	// Testoutput mass of component in model domain calculated at nodes
	MeshLib::CFEMesh* m_msh = fem_msh_vector[0];  // SB: ToDo hart gesetzt
	MeshLib::CElem* m_ele = NULL;
	MeshLib::CNode* m_node = NULL;
	CMediumProperties* m_mat_mp = NULL;
	std::ostringstream temp;
	double ComponentMass;
	// //int index;
	double saturation_water;
	double node_volume;
	double time;
	std::string tempstring;
	std::vector<std::string> vec_string;
	// int position;
	std::string path;
	//	double density_water; // 2012-08 TF, variable ‘density_water’ set but
	// not used [-Wunused-but-set-variable]
	double porosity;
	//	int variable_index; // 2012-08 TF, variable set but not used
	double ComponentConcentration, TotalVolume, SourceTerm;
	int indexComponentConcentration;
	// CRFProcess *n_pcs = NULL;
	int group;
	std::string filename;
	// bool Windows_System;
	int position;

	path = flow_pcs[0]->file_name_base;
// position = int(path.find_last_of("\\"));
// if (position >= 0)
//	Windows_System = true;
// else
//	Windows_System = false;
// if (Windows_System == true)
#ifdef _WIN32
	{
		position = int(path.find_last_of("\\"));
		if (position < 0)
			filename =
			    "Mass_" + transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "\\Mass_" +
			           transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		}
	}
#else
	{
		position = int(path.find_last_of("/"));
		if (position < 0)
			filename =
			    "Mass_" + transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		else
		{
			std::string path_new;
			path_new = path.substr(0, position);
			filename = path_new + "/Mass_" +
			           transport_pcs->nod_val_name_vector[0] + "_nodes.csv";
		}
	}
#endif
	if (flow_pcs[0]->Tim->step_current == 1)
	{
		// Header of the file
		vec_string.push_back("Time, ComponentMass, TotalVolume");
	}
#ifndef OGS_ONLY_TH
	else
	{
		//		//read file and store data
		//		CReadTextfiles_ECL TextFile;
		//		TextFile.Read_Text(filename);
		//
		//		for (int i = 0; i < TextFile.NumberOfRows; i++)
		//		{
		//			vec_string.push_back(TextFile.Data[i]);
		//		}
	}
#endif
	ComponentMass = TotalVolume = SourceTerm = 0;

	for (long i = 0; i < (long)m_msh->nod_vector.size(); i++)
	{
		m_node = m_msh->nod_vector[i];  // get element
		node_volume = 0;
		if (mfp_vector[0]->density_model == 18)
		{
			//			variable_index =
			//flow_pcs[0]->GetNodeValueIndex("DENSITY1");
			//// // 2012-08 TF, variable set but not used
			//			density_water = flow_pcs[0]->GetNodeValue(i,
			// variable_index); // 2012-08 TF, variable ‘density_water’ set but
			// not used [-Wunused-but-set-variable]
		}
		//		else
		//			density_water = mfp_vector[0]->Density(); // 2012-08 TF,
		// variable ‘density_water’ set but not used [-Wunused-but-set-variable]

		for (int j = 0; j < int(m_node->getConnectedElementIDs().size()); j++)
		{
			m_ele = m_msh->ele_vector[m_node->getConnectedElementIDs()[j]];

			// get the phase volume of current element elem
			group = m_ele->GetPatchIndex();
			m_mat_mp = mmp_vector[group];
			porosity = m_mat_mp->Porosity(
			    m_ele->GetIndex(),
			    1);  // CB Now provides also heterogeneous porosity, model 11
			node_volume = node_volume +
			              porosity * m_ele->GetVolume() * m_ele->GetFluxArea() /
			                  m_ele->GetNodesNumber(false);
			// cout << m_ele->GetNodesNumber(false) << " " << endl;
		}

		TotalVolume += node_volume;
		if (flow_pcs[0]->getProcessType() == FiniteElement::GROUNDWATER_FLOW)
			saturation_water =
			    1;  // saturation does not exist for example in groundwater flow
		else
		{
			int index = flow_pcs[0]->GetNodeValueIndex("SATURATION1") +
			            1;  //+1... new time level
			saturation_water = flow_pcs[0]->GetNodeValue(i, index);
		}
		indexComponentConcentration =
		    transport_pcs->GetNodeValueIndex(
		        transport_pcs->pcs_primary_function_name[0]) +
		    1;  // +1: new timelevel
		ComponentConcentration =
		    transport_pcs->GetNodeValue(i, indexComponentConcentration);

		ComponentMass = ComponentMass +
		                node_volume * saturation_water * ComponentConcentration;
		// cout << " Node: " << i << " saturation: " << saturation_water << "
		// Density_water: " << density_water << " node_volume: " << node_volume
		// << " Concentration: " << ComponentConcentration << endl;
	}

	// calculating source term
	for (long i = 0; i < (long)transport_pcs->st_node_value.size(); i++)
	{
		for (size_t j = 0; j < transport_pcs->st_node_value[i].size(); j++)
			SourceTerm = transport_pcs->st_node_value[i][j]->node_value;
	}
	SourceTerm *= transport_pcs->Tim->this_stepsize;

	// calculating time
	time = 0;
	for (int k = 0; k < transport_pcs->Tim->step_current; k++)
	{
		time += transport_pcs->Tim->time_step_vector[k];
	}
	temp.str("");
	temp.clear();
	temp.precision(14);
	temp << time;
	tempstring = temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(14);
	temp << ComponentMass;
	tempstring += temp.str();
	tempstring += ", ";
	temp.str("");
	temp.clear();
	temp.precision(14);
	temp << TotalVolume;
	tempstring += temp.str();

	vec_string.push_back(tempstring);

#ifndef OGS_ONLY_TH
//	//within the first timestep create file and write header
//	CWriteTextfiles_ECL TextFile;
//	TextFile.Write_Text(filename, vec_string);
#endif
}

/*-------------------------------------------------------------------------
   GeoSys - Function: PS_Global()
   Task: Similate multi-phase flow by p-p scheme
   Return: error
   Programming:
   03/2009 PCH Implementation
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::PS_Global()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[3];
	if (!m_pcs->selected) return error;
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	if (m_pcs->TimeStepAccept()) m_pcs->CalIntegrationPointValue();

	if (m_pcs->OutputMassOfGasInModel == true)  // 05/2012 BG
		OutputMassOfGasInModel(m_pcs);

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: GroundWaterFlow()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Update
   05.2009 WW For surface-soil-ground coupled model
   -------------------------------------------------------------------------*/
inline double Problem::GroundWaterFlow()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[1];
	if (!m_pcs->selected) return error;  // 12.12.2008 WW
	error = m_pcs->ExecuteNonLinear(loop_process_number);
	//................................................................
	// Calculate secondary variables
	// NOD values
	conducted = true;  // WW
	// std::cout << "      Calculation of secondary NOD values" << "\n";
	if (m_pcs->TimeStepAccept())
	{
		PCSCalcSecondaryVariables();  // PCS member function
		// std::cout << "      Calculation of secondary GP values" << "\n";
		m_pcs->CalIntegrationPointValue();  // WW
		m_pcs->cal_integration_point_value =
		    true;  // WW Do not extropolate Gauss velocity

	}
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: AirFlow()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::AirFlow()
{
	double error = 1.e8;
	CRFProcess* m_pcs = total_processes[7];
	if (!m_pcs->selected) return error;  // 12.12.2008 WW

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	m_pcs->CalIntegrationPointValue();           // WW
	m_pcs->cal_integration_point_value = false;  // AKS
	m_pcs->CalcELEVelocities();                  // OK
	//
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: HeatTransport
   Task: Similate heat transport
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::HeatTransport()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[8];
	if (!m_pcs->selected) return error;  // 12.12.2008 WW

	error = m_pcs->ExecuteNonLinear(loop_process_number);
	// if(m_pcs->non_linear)
	//  error = m_pcs->ExecuteNonLinear();
	// else
	//  error = m_pcs->Execute();
	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: MassTrasport
   Task: Similate heat transport
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW Immigrtate the new functionalities  from loop_pcs.cpp
   -------------------------------------------------------------------------*/
inline double Problem::MassTrasport()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[11];
	//
	if (!m_pcs->selected) return error;  // 12.12.2008 WW

	for (int i = 0; i < (int)transport_processes.size(); i++)
	{
		m_pcs = transport_processes[i];  // 18.08.2008 WW
		// Component Mobile ?
		if (CPGetMobil(m_pcs->GetProcessComponentNumber()) > 0)
			error = m_pcs->ExecuteNonLinear(
			    loop_process_number);  // NW. ExecuteNonLinear() is called to
		                               // use the adaptive time step scheme

		int component = m_pcs->pcs_component_number;
		CompProperties* m_cp = cp_vec[component];

		if (m_cp->OutputMassOfComponentInModel == 1)
		{  // 05/2012 BG
			if (singlephaseflow_process.size() > 0)
				OutputMassOfComponentInModel(singlephaseflow_process, m_pcs);
			else
				OutputMassOfComponentInModel(multiphase_processes, m_pcs);
		}
	}

	// Calculate Chemical reactions, after convergence of flow and transport
	// Move inside iteration loop if couplingwith transport is implemented
	// SB:todo
	// SB:todo move into Execute Reactions	  if((aktueller_zeitschritt % 1) ==
	// 0)
	// REACT *rc = NULL; //OK
	// rc = REACT_vec[0]; //OK
	//				if(rc->flag_pqc) rc->ExecuteReactions();
	//				delete rc;
	if (KinReactData_vector.size() >
	    0)  // WW moved the following lines into this curly braces. 12.12.2008
		// SB4900    ClockTimeVec[0]->StopTime("Transport");
		// SB4900    ClockTimeVec[0]->StartTime();
		// Calculate Chemical reactions, after convergence of flow and transport
		// Move inside iteration loop if couplingwith transport is implemented
		// SB:todo
		// First calculate kinetic reactions
		KinReactData_vector[0]->ExecuteKinReact();
	// SB4900 ClockTimeVec[0]->StopTime("KinReactions");
	// SB4900 ClockTimeVec[0]->StartTime();
	if (REACT_vec.size() > 0)  // OK
	{
		if (REACT_vec[0]->flag_pqc)
		{
#ifdef REACTION_ELEMENT
			REACT_vec[0]->ExecuteReactionsPHREEQC0();
#else
// REACT_vec[0]->ExecuteReactions();

#ifdef LIBPHREEQC
			// MDL: built-in phreeqc
			REACT_vec[0]->ExecuteReactionsPHREEQCNewLib();
#else
			REACT_vec[0]->ExecuteReactionsPHREEQCNew();
#endif  // LIBPHREEQC
#endif  // REACTION_ELEMENT
		}
	}
#ifdef GEM_REACT
	else  // WW moved these pare of curly braces inside  ifdef GEM_REACT
	    if (m_vec_GEM->initialized_flag == 1)  // when it was initialized.
	{
		int m_time =
		    1;  // 0-previous time step results; 1-current time step results

		// Check if the Sequential Iterative Scheme needs to be intergrated
		if (m_pcs->m_num->cpl_max_iterations > 1)
			m_vec_GEM->flag_iterative_scheme =
			    1;  // set to standard iterative scheme;
		// Move current xDC to previous xDC
		m_vec_GEM->CopyCurXDCPre();
		// Get info from MT
		// m_vec_GEM->ConvPorosityNodeValue2Elem(); //
		// second arguments should be one if we work with concentrations
		m_vec_GEM->GetReactInfoFromMassTransport(m_time);
		// m_vec_GEM->ConcentrationToMass();
		m_vec_GEM->Run_MainLoop();  // Run GEM
		m_vec_GEM->CopyCurBPre();

		// m_vec_GEM->MassToConcentration();
		// Calculate the different of xDC
		m_vec_GEM->UpdateXDCChemDelta();
		// Set info in MT
		m_vec_GEM->SetReactInfoBackMassTransport(m_time);
		// m_vec_GEM->ConvPorosityNodeValue2Elem(); // update element porosity
		// and push back values
	}
#endif  // GEM_REACT

#ifdef BRNS
	if (m_vec_BRNS->init_flag == true)
		m_vec_BRNS->RUN(dt /*time value in seconds*/);
#endif

	// if(KinReactData_vector.size() > 0)  //12.12.2008 WW
	// SB4900    ClockTimeVec[0]->StopTime("EquiReact");

	return error;
}

/*-------------------------------------------------------------------------
   GeoSys - Function: FluidMomentum()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::FluidMomentum()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[9];
	//
	if (!m_pcs->selected) return error;  // 12.12.2008 WW

	CFluidMomentum* fm_pcs = (CFluidMomentum*)m_pcs;

	fm_pcs->Execute(loop_process_number);

	// Switch off rechard flow if
	if (m_pcs->num_type_name.compare("STEADY") == 0 &&
	    aktueller_zeitschritt > 1)
	{
		// Turn off FLUID_MOMENTUM
		m_pcs->selected = false;
		// Turn off RICHARDS_FLOW
		m_pcs = PCSGet("RICHARDS_FLOW");
		if (m_pcs) m_pcs->selected = false;
		// Turn off LIQUID_FLOW
		m_pcs = PCSGet("LIQUID_FLOW");
		if (m_pcs) m_pcs->selected = false;
		// Turn off GROUNDWATER_FLOW
		m_pcs = PCSGet("GROUNDWATER_FLOW");
		if (m_pcs) m_pcs->selected = false;
	}
	//
	// error = 0.0 // JT... in unsteady flow, setting error=0.0 corresponds to
	// error_cpl=0.0, and the coupling loop ceases before RWPT is performed
	//            // What is the correct way to handle this, rather than setting
	//            error=1.e8???
	return error;
}

#ifndef OGS_ONLY_TH
/*-------------------------------------------------------------------------
   GeoSys - Function: RandomWalker()
   Task:
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   12.2008 WW
   -------------------------------------------------------------------------*/
inline double Problem::RandomWalker()
{
	double error = 1.0e+8;
	//
	CRFProcess* m_pcs = total_processes[10];
	//
	if (!m_pcs || !m_pcs->selected) return error;  // 12.12.2008 WW
	//
	// CFEMesh* m_msh = NULL;

	if (m_pcs && m_pcs->selected)
	{
		lop_coupling_iterations = 1;

		RandomWalk* rw_pcs = (RandomWalk*)m_pcs;

// Do I need velocity fileds solved by the FEM?
#if 0
		if(m_pcs->tim_type_name.compare("PURERWPT") == 0)
		{
			rw_pcs->PURERWPT = 1;
			char* dateiname = NULL;
			int sizeOfWord = 100;
			dateiname = (char*)malloc(sizeOfWord * sizeof(char ));

			std::string filename = FileName;
			for(int i = 0; i <= (int)filename.size(); ++i)
				dateiname[i] = filename[i];

			rw_pcs->ReadInVelocityFieldOnNodes(dateiname);

			free(dateiname);
		}
#endif

		// Set the mode of the RWPT method
		if (m_pcs->num_type_name.compare("HETERO") == 0)
		{
			rw_pcs->RWPTMode = 1;  // Set it for heterogeneous media
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_ADVECTION") == 0)
		{
			rw_pcs->RWPTMode = 2;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_ADVECTION") == 0)
		{
			rw_pcs->RWPTMode = 3;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_DISPERSION") == 0)
		{
			rw_pcs->RWPTMode = 4;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_DISPERSION") == 0)
		{
			rw_pcs->RWPTMode = 5;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_FDM") == 0)
		{
			rw_pcs->RWPTMode = 1;  // Set it for heterogeneous media
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_ADVECTION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 2;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_ADVECTION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 3;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HOMO_DISPERSION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 4;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else if (m_pcs->num_type_name.compare("HETERO_DISPERSION_FDM") == 0)
		{
			rw_pcs->RWPTMode = 5;
			std::cout << "RWPT is on " << m_pcs->num_type_name << " mode."
			          << "\n";
		}
		else  // HOMO Advection + Dispersion
		{
			rw_pcs->RWPTMode = 0;
			std::cout << "RWPT is on HOMO_ADVECTION_DISPERSION mode."
			          << "\n";
		}

		if (m_pcs->num_type_name.find("FDM") != std::string::npos)
		{
			rw_pcs->PURERWPT = 2;
			if (rw_pcs->FDMIndexSwitch == 0)
			{
				rw_pcs->buildFDMIndex();
				// Switch off
				rw_pcs->FDMIndexSwitch = 1;
			}
		}

		if (rwpt_numsplits < 0)
			rwpt_numsplits = 10;  // JT 2010 set default value, unless specified
		                          // in .tim input file

		rw_pcs->AdvanceBySplitTime(dt, rwpt_numsplits);
		//	rw_pcs->TraceStreamline(); // JT, no longer needed
		print_result = true;
		rw_pcs->RandomWalkOutput(aktuelle_zeit, aktueller_zeitschritt);
	}

	return 0.0;
}
#endif

/*-------------------------------------------------------------------------
   GeoSys - Function: Deformation
   Task: Similate deformation
   Return: error
   Programming:
   07/2008 WW Extract from LOPTimeLoop_PCS();
   Modification:
   -------------------------------------------------------------------------*/
inline double Problem::Deformation()
{
	CRFProcessDeformation* dm_pcs = NULL;
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[12];
	//
	dm_pcs = (CRFProcessDeformation*)(m_pcs);
	error = dm_pcs->Execute(loop_process_number);
	// Error
	if (dm_pcs->type / 10 == 4)
	{
		m_pcs->cal_integration_point_value = true;
		dm_pcs->CalIntegrationPointValue();

		if (dm_pcs->type == 42)  // H2M. 07.2011. WW
			dm_pcs->CalcSecondaryVariablesUnsaturatedFlow();
	}
	return error;
}

inline double Problem::TH_Monolithic()
{
	double error = 1.0e+8;
	CRFProcess* m_pcs = total_processes[14];
	if (!m_pcs->selected) return error;
	CRFProcessTH* th_pcs = (CRFProcessTH*)m_pcs;
	error = th_pcs->Execute(loop_process_number);
	th_pcs->CalIntegrationPointValue();

	return error;
}


/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2003 SB Implementation
   01/2004 MX k_eff
   11/2004 OK fluid mass fluxes
   08/2008 WW Extract from LOPTimeLoop_PCS();
   last modification:
**************************************************************************/
void Problem::LOPCalcELEResultants()
{
	size_t no_processes = pcs_vector.size();
	CRFProcess* m_pcs = NULL;

	for (size_t p = 0; p < no_processes; p++)
	{
		m_pcs = pcs_vector[p];
		if (!m_pcs->selected)
			continue;
		std::string pcs_type_name(
		    convertProcessTypeToString(m_pcs->getProcessType()));
		switch (pcs_type_name[0])
		{
			default:
				break;
			case 'L':  // Liquid flow
				m_pcs->CalcELEVelocities();
				break;
			case 'G':  // Groundwater flow
				m_pcs->CalcELEVelocities();
				break;
			case 'A':  // Gas flow
				m_pcs->CalcELEVelocities();
				// m_pcs->CalcELEMassFluxes();			// BG
				break;
			case 'T':  // Two-phase flow
				break;
			case 'C':  // Componental flow
				break;
			case 'H':  // Heat transport
				break;
			case 'M':  // Mass transport
				break;
			case 'D':  // Deformation
				break;
			case 'R':  // Richards flow
				m_pcs->CalcELEVelocities();
				break;
			case 'F':  // Fluid Momentum
				break;
		}
	}
}

/**************************************************************************
   ROCKFLOW - Funktion: ASMCalcNodeWDepth

   Task:
   Berechnung und Speichern der Knotenfl?se
   Parameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E: long i: node index
   Result:
   - void -

   Programmaenderungen:
   11/2002   MB/OK  Implementation
   10/2004   MB     PCS
   12/2008   WW     Encapsulate to this class
**************************************************************************/
inline void Problem::ASMCalcNodeWDepth(CRFProcess* m_pcs)
{
	int nidx, nidy, nidz;
	// OK411 int timelevel = 1;
	double WDepth;

	nidx = m_pcs->GetNodeValueIndex("HEAD") + 1;
	nidy = m_pcs->GetNodeValueIndex("WDEPTH");
	nidz = m_pcs->GetNodeValueIndex("COUPLING");

	const size_t n_nodes(m_pcs->m_msh->nod_vector.size());
	for (size_t nn = 0; nn < n_nodes; nn++)
	{
		WDepth = m_pcs->GetNodeValue(nn, nidx) -
		         m_pcs->m_msh->nod_vector[nn]->getData()[2];
		// JOD only needed for GREEN_AMPT source term
		m_pcs->SetNodeValue(nn, nidz, m_pcs->GetNodeValue(nn, nidz + 1));
		if (WDepth < 0.0) WDepth = 0.0;
		m_pcs->SetNodeValue(nn, nidy, WDepth);
	}
}

/**************************************************************************/
/* ROCKFLOW - Funktion: PCSCalcSecondaryVariables
 */
/* Aufgabe:
   Berechung von secondary variables w?rend der Zeitschleife
   Abfrage je nach Prozess, der aktiv ist
 */
/* Programmaenderungen:
   08/2003   SB   Implementation
   01/2006   YD   add dual porosity
   01/2007 OK Two-phase flow
 */
/**************************************************************************/
void Problem::PCSCalcSecondaryVariables()
{
	int i, ptype;
	CRFProcess* m_pcs = NULL;
	bool NAPLdiss = false;
	NAPLdiss = KNaplDissCheck();
	/* go through all processes */
	int no_processes = (int)pcs_vector.size();
	for (i = 0; i < no_processes; i++)
	{
		m_pcs = pcs_vector[i];  // JOD
		if (m_pcs != NULL)
		{
			ptype = m_pcs->GetObjType();
			switch (ptype)
			{
				case 12: /* Multi-phase flow process */
					MMPCalcSecondaryVariablesNew(m_pcs, NAPLdiss);
					// MMPCalcSecondaryVariablesNew(m_pcs);
					break;
				default:
					break;
			}
		}  // If
	}      // while
}

/**************************************************************************
   FEMLib-Method:
   05/2009 OK Implementation
**************************************************************************/
bool Problem::Check()
{
	CRFProcess* m_pcs = NULL;
	for (int i = 0; i < (int)total_processes.size(); i++)
	{
		m_pcs = total_processes[i];
		if (!m_pcs->Check()) return false;
	}
	return true;
}

#ifdef BRNS

// BRNS-Coupling: For writing spatially resolved reaction rates at the final
// iteration,
// we need to get the timing information.

double Problem::getCurrentTime()
{
	return current_time;
}

double Problem::getEndTime()
{
	return end_time;
}
#endif  // BRNS
