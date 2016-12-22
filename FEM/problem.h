/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef problem_INC
#define problem_INC

#include <vector>

class CRFProcess;

#include "GEOObjects.h"

//---------------------------------------------------------------------
// Pointers to member functions
class Problem;
typedef double (Problem::*ProblemMemFn)(void);
#define Call_Member_FN(object, ptrToMember) ((object)->*(ptrToMember))
//---------------------------------------------------------------------

class Problem
{
public:
	explicit Problem(char* filename);
	~Problem();

	bool Euler_TimeDiscretize();
	void RosenBrock_TimeDiscretize() {}
	//
	void SetActiveProcesses();
	void SetTimeActiveProcesses();
	void PCSRestart();
	//
	bool CouplingLoop();
	void PostCouplingLoop();
	void PreCouplingLoop(CRFProcess* m_pcs = NULL);

	// Copy u_n for auto time stepping
	double* GetBufferArray(const bool is_x_k = false)
	{
		if (is_x_k)
			return buffer_array1;
		else
			return buffer_array;
	}
	int GetCPLMaxIterations() { return cpl_overall_max_iterations; }

	/**
	 * get the geometric objects stored in GEOLIB::GEOObjects
	 * @return a pointer to an instance of class GEOLIB::GEOObjects
	 */
	const GEOLIB::GEOObjects* getGeoObj() const;
	/**
	 * Get the name of the project. The name is used by GEOLIB::GEOObjects
	 * to access the geometric data.
	 * @return the name to acces geometric data
	 */
	const std::string& getGeoObjName() const;

private:
	// Time:
	double start_time;
	double end_time;
	double current_time;
	double* buffer_array;
	double* buffer_array1;
	int step_control_type;
	bool last_dt_accepted;      // JT
	bool force_post_node_copy;  // JT
	// Mixed time step WW
	double dt0;  // Save the original time step size

	// Controls
	bool external_coupling_exists;
	int cpl_overall_max_iterations;
	int cpl_overall_min_iterations;
	int loop_process_number;
	size_t max_time_steps;
	//
	int lop_coupling_iterations;
	bool CalcVelocities;
	bool conducted;

	// Print flag
	bool print_result;
	// Processes
	std::vector<CRFProcess*> total_processes;
	std::vector<CRFProcess*> transport_processes;
	std::vector<CRFProcess*> multiphase_processes;
	std::vector<CRFProcess*> singlephaseflow_process;
	ProblemMemFn* active_processes;
	std::vector<int> active_process_index;
	std::vector<int> coupled_process_index;
	//
	bool* exe_flag;
	int AssignProcessIndex(CRFProcess* m_pcs, bool activefunc = true);
	//
	void PCSCreate();
	// Perform processes:
	double LiquidFlow();
	double RichardsFlow();
	double TwoPhaseFlow();
	double MultiPhaseFlow();
	double PS_Global();
	double GroundWaterFlow();
	double AirFlow();
	double HeatTransport();
	double FluidMomentum();
	double RandomWalker();
	double MassTrasport();
	double Deformation();
	double TH_Monolithic();
	// Accessory
	void LOPCalcELEResultants();
	void PCSCalcSecondaryVariables();
	bool Check();
	void OutputMassOfComponentInModel(std::vector<CRFProcess*> flow_pcs,
	                                  CRFProcess* transport_pcs);  // BG
	void OutputMassOfGasInModel(CRFProcess* m_pcs);                // BG

	/**
	 * pointer to an instance of class GEOObjects,
	 * that manages geometric entities
	 */
	GEOLIB::GEOObjects* _geo_obj;
	/**
	 * project/file name for geometry file,
	 * used to access data in data manager GEOObjects
	 */
	std::string _geo_name;
};

#endif
