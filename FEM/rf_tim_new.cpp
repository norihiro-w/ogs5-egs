/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: TIM
   Task:
   Programing:
   08/2004 OK Implementation
   last modified:
**************************************************************************/
#include "rf_tim_new.h"

#include <cfloat>

#include "display.h"
#include "FileToolsRF.h"
#include "makros.h"

#include "fem_ele_std.h"
#include "mathlib.h"
#include "Output.h"
#include "rfmat_cp.h"
#include "rf_mmp_new.h"
#include "tools.h"

double aktuelle_zeit;
size_t aktueller_zeitschritt = 0;
double dt = 0.0;
int rwpt_numsplits = -1;
std::vector<CTimeDiscretization*> time_vector;

/**************************************************************************
   FEMLib-Method:
   Task: OBJ constructor
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CTimeDiscretization::CTimeDiscretization(void)
    : Write_tim_discrete(false), tim_discrete(NULL)  // YD
{
	step_current = 0;
	time_start = std::numeric_limits<double>::max();
	time_end = 0.0;
	time_type_name = "CONSTANT";                   // OK
	time_control_type = TimeControlType::INVALID;  // kg44//JT
	time_unit = "SECOND";
	max_time_step = 1.e10;  // YD
	min_time_step =
	    DBL_EPSILON;  // YD//JT Minimum allowed timestep, this process
	initial_time_step = DBL_EPSILON;
	repeat = false;                     // OK/YD
	step_current = 0;                   // WW
	this_stepsize = 0.;                 // WW
	dt_sum = .0;                        // WW
	relative_error = 1.e-4;             // 26.08.2008. WW
	absolute_error = 1.e-10;            // 26.08.2008. WW
	h_max = 6;                          // 27.08.2008. WW
	h_min = 0.2;                        // 27.08.2008. WW
	hacc = 0.;                          // 27.08.2008. WW
	erracc = 0.;                        // 27.08.2008. WW
	PI_tsize_ctrl_type = -1;            // 27.08.2008. WW
	minimum_dt_reached = false;         // JT
	time_active = true;                 // JT
	time_independence = false;          // JT
	dt_failure_reduction_factor = 1.0;  // JT
	accepted_step_count = 0;            // JT
	rejected_step_count = 0;            // JT
	last_active_time = 0.0;             // JT
	next_active_time = 0.0;             // JT
	dynamic_time_buffer = 0;            // JT
	for (size_t ii = 0; ii < DOF_NUMBER_MAX + 1; ii++)
	{
		dynamic_control_tolerance[ii] = -1.0;
	}
	last_rejected_timestep = 0;
	iteration_type = IterationType::LINEAR;
	last_dt_accepted = true;
	usePIDControlInSelfAdaptive = false;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ destructor
   Programing:
   08/2004 OK Implementation
**************************************************************************/
CTimeDiscretization::~CTimeDiscretization(void)
{
	if (tim_discrete)  // YD
	{
		tim_discrete->close();
		if (tim_discrete) delete tim_discrete;
		time_step_vector.clear();
		time_adapt_tim_vector.clear();
		time_adapt_coe_vector.clear();
	}
}


/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation
   11/2004 OK string streaming by SB for lines
   10/2005 YD Time Controls
   08/2008 WW General classic time step size control (PI control)
**************************************************************************/
std::ios::pos_type CTimeDiscretization::Read(std::ifstream* tim_file)
{
	std::string sub_line;
	std::string line_string;
	std::string delimiter(" ");
	bool new_keyword = false;
	std::string hash("#");
	std::ios::pos_type position;
	std::string sub_string;
	bool new_subkeyword = false;
	std::string dollar("$");
	int no_time_steps = 0;
	double time_step_length;
	std::ios::pos_type position_subkeyword;
	std::stringstream line;
	std::string line_complete;
	int iter_times;        // YD
	double multiply_coef;  // YD
	int i;
	CRFProcess* m_pcs = NULL;
	//    m_pcs = PCSGet("RICHARDS_FLOW");
	m_pcs = PCSGet("GROUNDWATER_FLOW");  // kg44 changed default

	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		if (new_subkeyword) tim_file->seekg(position, std::ios::beg);
		new_subkeyword = false;
		position = GetNextSubKeyword(tim_file, &line_string, &new_keyword);
		if (new_keyword) return position;
		/*
		    position = tim_file->tellg();
		    if(new_subkeyword)
		      tim_file->seekg(position_subkeyword,ios::beg);
		    new_subkeyword = false;
		    tim_file->getline(buffer,MAX_ZEILE);
		    line_string = buffer;
		   if(line_string.size()<1) // empty line
		      continue;
		    if(Keyword(line_string))
		      return position;
		 */
		//....................................................................

		// subkeyword found
		if (line_string.find("$PCS_TYPE") != std::string::npos)
		{
			line.str(GetLineFromFile1(tim_file));
			line >> pcs_type_name;
			line.clear();
			m_pcs = PCSGet(pcs_type_name);  // kg44 inserted to overwrite
			                                // default Richards_flow
			// this works only of pcs_type is read before adaption
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIME_START") != std::string::npos)
		{
			line.str(GetLineFromFile1(tim_file));
			line >> time_start;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIME_END") != std::string::npos)
		{
			line.str(GetLineFromFile1(tim_file));
			line >> time_end;
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIME_UNIT") != std::string::npos)
		{
			*tim_file >> time_unit >> std::ws;  // WW unit of time
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$INDEPENDENT") != std::string::npos)  // JT2012
		{
			line.str(GetLineFromFile1(tim_file));
			line >> time_independence;
			// =0: Process will adopt minimum time step of all processes.
			// =1: Process will have it's own time step (may not execute on
			// every time step). It will be independent.
			line.clear();
			continue;
		}
		//....................................................................
		/* //WW
		   if(line_string.find("$TIME_FIXED_POINTS")!=string::npos) { //
		   subkeyword found
		   int no_fixed_points;
		   double fixed_point;
		   line.str(GetLineFromFile1(tim_file));
		   line >> no_fixed_points;
		   line.clear();
		   for(i=0;i<no_fixed_points;i++) {
		    line.str(GetLineFromFile1(tim_file));
		    line >> fixed_point;
		   fixed_point_vector.push_back(fixed_point);
		   line.clear();
		   }
		   continue;
		   }
		 */
		//....................................................................
		// subkeyword found
		if (line_string.find("$TIME_STEPS") != std::string::npos)
			while ((!new_keyword) || (!new_subkeyword) || (!tim_file->eof()))
			{
				position = tim_file->tellg();
				line_string = GetLineFromFile1(tim_file);
				if (line_string.find("#") != std::string::npos) return position;
				if (line_string.find("$") != std::string::npos)
				{
					new_subkeyword = true;
					break;
				}
				line.str(line_string);
				line >> no_time_steps;
				line >> time_step_length;
				for (i = 0; i < no_time_steps; i++)
					time_step_vector.push_back(time_step_length);
				line.clear();
			}
		// subkeyword found
		if (line_string.find("$TIME_SPLITS") != std::string::npos)
		{
			line.str(GetLineFromFile1(tim_file));
			line >> rwpt_numsplits;
			line.clear();
			continue;
		}
		// 25.08.2008. WW
		if (line_string.find("$CRITICAL_TIME") != std::string::npos)
			while ((!new_keyword) || (!new_subkeyword) || (!tim_file->eof()))
			{
				position = tim_file->tellg();
				line_string = GetLineFromFile1(tim_file);
				if (line_string.find("#") != std::string::npos) return position;
				if (line_string.find("$") != std::string::npos)
				{
					new_subkeyword = true;
					break;
				}
				line.str(line_string);
				double crtime;
				line >> crtime;
				critical_time.push_back(crtime);
				line.clear();
			}
		// subkeyword found
		if (line_string.find("$TIME_CONTROL") != std::string::npos)
			while ((!new_keyword) || (!new_subkeyword) || (!tim_file->eof()))
			{
				position = tim_file->tellg();
				line_string = GetLineFromFile1(tim_file);

				if (line_string.find("#") != std::string::npos) return position;
				if (line_string.find("$") != std::string::npos)
				{
					new_subkeyword = true;
					break;
				}
				line.str(line_string);
				std::string time_control_name;
				line >> time_control_name;
				line.clear();
				time_control_type = convertTimeControlType(time_control_name);

				if (time_control_type ==
				    TimeControlType::PI_AUTO_STEP_SIZE)  // 26.08.2008. WW
				{
					line.str(GetLineFromFile1(tim_file));
					line >> PI_tsize_ctrl_type >> relative_error >>
					    absolute_error >> this_stepsize;
					// 13.03.2008. WW
					int real_type = (int)(PI_tsize_ctrl_type / 10);
					if (real_type < 10 && real_type > 0)  //
					{
						PI_tsize_ctrl_type = real_type;
						line >> h_min >> h_max >> max_time_step;
					}
					else
						max_time_step = 0.0;
					line.clear();
				}
				else if (time_control_type == TimeControlType::PID_CONTROL)
				{
					usePIDControl = true;
					line.str(GetLineFromFile1(tim_file));
					line >> pid_error >> initial_time_step >> min_time_step >>
					    max_time_step;
					line.clear();
				}
				else if (time_control_type ==
				         TimeControlType::DYNAMIC_VARIABLE)  // JT2012
				{
					// DYNAMIC TIME STEP SERIES
					line.str(GetLineFromFile1(tim_file));
					int num_tolerances = DOF_NUMBER_MAX;
					double include_third_variable = -1.0;
					double third_variable_tolerance = -1.0;
					//
					line >> dynamic_control_error_method;  // Corresponds to
					// FiniteElement::ErrorMethod.
					// Defines how tolerance
					// is applied.
					line >> time_step_length;           // initial_dt
					line >> min_time_step;              // min_dt
					line >> dynamic_failure_threshold;  // threshold to force a
					                                    // time failure
					                                    // (recommend 1.1-2.0.
					                                    // If tolerance is
					                                    // exceeeded on first
					                                    // iteration by this
					                                    // factor, dt will be
					                                    // decreased and a
					                                    // failure forced)
					line >> dynamic_control_tolerance
					            [DOF_NUMBER_MAX];  // max_increase_factor (dt
					                               // never allowed to increase
					                               // faster than this factor
					                               // (i.e. 1.5)
					//
					// tolerances[1:dof]: One tolerance for each degree of
					// freedom in process (or only one tolerance for certain
					// methods)
					switch (FiniteElement::convertErrorMethod(
					    dynamic_control_error_method))
					{
						case FiniteElement::ENORM:  // only 1 tolerance required
							num_tolerances = 1;
							line >> dynamic_control_tolerance[0];
							//
							// Next entry is OPTIONAL (may be left blank)!
							line >> include_third_variable;  // if >0.0, include
							                                 // a 3rd variable
							                                 // in the global
							                                 // norm (PwSnw
							                                 // model: Pc --
							                                 // PcPnw model: Snw
							                                 // -- PwPnw model:
							                                 // Snw)
							break;
						//
						case FiniteElement::ERNORM:  // only 1 tolerance
							                         // required
							num_tolerances = 1;
							line >> dynamic_control_tolerance[0];
							//
							// Next entry is OPTIONAL (may be left blank)!
							line >> include_third_variable;  // if >0.0, include
							                                 // a 3rd variable
							                                 // in the global
							                                 // norm (PwSnw
							                                 // model: Pc --
							                                 // PcPnw model: Snw
							                                 // -- PwPnw model:
							                                 // Snw)
							break;
						//
						case FiniteElement::EVNORM:  // 1 tolerance for each
							                         // primary variable (for
							                         // Deformation, only 1
							                         // tolerance required.
							                         // Applies to x,y,z)
							for (int i = 0; i < num_tolerances; i++)
								line >> dynamic_control_tolerance[i];
							//
							// Next entry is OPTIONAL (may be left blank)!
							line >> third_variable_tolerance;  // tolerance of a
							                                   // third variable
							                                   // (PwSnw model:
							                                   // Pc -- PcPnw
							                                   // model: Snw --
							                                   // PwPnw model:
							                                   // Snw)
							break;
						//
						case FiniteElement::LMAX:  // 1 tolerance for each
							                       // primary variable (for
							                       // Deformation, only 1
							                       // tolerance required.
							                       // Applies to x,y,z)
							for (int i = 0; i < num_tolerances; i++)
								line >> dynamic_control_tolerance[i];
							//
							// Next entry is OPTIONAL (may be left blank)!
							line >> third_variable_tolerance;  // tolerance of a
							                                   // third variable
							                                   // (PwSnw model:
							                                   // Pc -- PcPnw
							                                   // model: Snw --
							                                   // PwPnw model:
							                                   // Snw)
							break;
						//
						case FiniteElement::BNORM:
							ScreenMessage(
							    "ERROR in TIMRead. BNORM not configured for "
							    "time control.\n");
							ScreenMessage(
							    "We suggest ENORM as a valid companion for "
							    "NEWTON couplings.\n");
							exit(1);
							break;
						//
						default:
							ScreenMessage(
							    "ERROR in TIMRead. Invalid error method "
							    "selected for dynamic time control.\n");
							exit(1);
							break;
					}
					//
					if (num_tolerances > 1)
						dynamic_control_tolerance[DOF_NUMBER_MAX] =
						    third_variable_tolerance;
					else
						dynamic_control_tolerance[DOF_NUMBER_MAX] =
						    include_third_variable;
					line.clear();
				}
				else if (time_control_type ==
				         TimeControlType::DYNAMIC_COURANT)  // JT2012
				{
					std::cout << "Josh apologizes (especially to Marc), but "
					             "promises DYNAMIC_COURANT will be merged into "
					             "the next release."
					          << "\n";
					std::cout << "emoticon:: sad face"
					          << "\n";
					exit(1);
					//
					// DYNAMIC TIME STEP SERIES
					line.str(GetLineFromFile1(tim_file));
					//
					line >> time_step_length;  // initial_dt
					line >> min_time_step;     // min_dt
					line >> dynamic_control_tolerance
					            [DOF_NUMBER_MAX];  // max_increase_factor (dt
					                               // never allowed to increase
					                               // faster than this factor
					                               // (i.e. 1.5)
					line >>
					    dynamic_control_tolerance[0];  // desired courant number
					//
					// ADDITIONAL OPTIONS TO RESTRICT CALCULATION TO CERTAIN
					// ELEMENTS
					if (time_control_name.find("CONCENTRATION") !=
					    std::string::npos)
					{  // DYNAMIC_COURANT_CONCENTRATION
						line >> dynamic_control_tolerance[1];  // Concentration
						                                       // threshold
						                                       // (elements with
						                                       // concentration
						                                       // beneath this
						                                       // value are not
						                                       // included in
						                                       // Courant
						                                       // restriction)
					}
					else if (time_control_name.find("TEMPERATURE") !=
					         std::string::npos)
					{  // DYNAMIC_COURANT_TEMPERATURE
						line >> dynamic_control_tolerance[1];  // Temperature
						                                       // threshold
						                                       // (elements with
						                                       // temperature
						                                       // BENEATH this
						                                       // value are not
						                                       // included in
						                                       // Courant
						                                       // restriction)
						line >> dynamic_control_tolerance[2];  // Temperature
						                                       // threshold
						                                       // (elements with
						                                       // temperature
						                                       // ABOVE   this
						                                       // value are not
						                                       // included in
						                                       // Courant
						                                       // restriction)
					}
					//
					line.clear();
				}
				else if (time_control_type ==
				         TimeControlType::DYNAMIC_PRESSURE)  // JT2012
				{
					// DYNAMIC TIME STEP SERIES
					line.str(GetLineFromFile1(tim_file));
					//
					line >> time_step_length;  // initial_dt
					line >> min_time_step;     // min_dt
					line >> dynamic_control_tolerance
					            [DOF_NUMBER_MAX];  // max_increase_factor (dt
					                               // never allowed to increase
					                               // faster than this factor
					                               // (i.e. 1.5)
					line >> dynamic_control_tolerance[0];  // pressure tolerance
					                                       // (mean pressure)
					//
					line.clear();
				}
				else if (time_control_type ==
				         TimeControlType::STEP_SIZE_RESTRICTION)  // 26.08.2008.
				                                                  // WW
				{
					line.str(GetLineFromFile1(tim_file));
					line >> h_min >> h_max;
					line.clear();
				}
				else if (time_control_type == TimeControlType::NEUMANN)
				{
					line.clear();
				}
				else if (time_control_type ==
				         TimeControlType::ERROR_CONTROL_ADAPTIVE)
				{
					m_pcs->adaption = true;
					line.clear();
				}
				else if (time_control_type == TimeControlType::SELF_ADAPTIVE)
				{
					if (time_control_name.find("_LINEAR") != std::string::npos)
						iteration_type = IterationType::LINEAR;
					else if (time_control_name.find("_NONLINEAR") !=
					         std::string::npos)
						iteration_type = IterationType::NONLINEAR;
					else if (time_control_name.find("_COUPLING") !=
					         std::string::npos)
						iteration_type = IterationType::COUPLING;
					ScreenMessage(
					    "-> SELF_ADAPTIVE is chosen for time stepping based on "
					    "%s iteration numbers\n",
					    convertIterationTypeToString(iteration_type).c_str());

					while ((!new_keyword) || (!new_subkeyword) ||
					       (!tim_file->eof()))
					{
						position = tim_file->tellg();
						line_string = GetLineFromFile1(tim_file);
						if (line_string.find("#") != std::string::npos)
							return position;
						if (line_string.find("$") != std::string::npos)
						{
							new_subkeyword = true;
							break;
						}
						if (line_string.find("MAX_TIME_STEP") !=
						    std::string::npos)
						{
							*tim_file >> line_string;
							max_time_step = strtod(line_string.data(), NULL);
							line.clear();
							// kg44 should not break break;
						}
						else if (line_string.find("MIN_TIME_STEP") !=
						         std::string::npos)
						{
							*tim_file >> line_string;
							min_time_step = strtod(line_string.data(), NULL);
							line.clear();
							// kg44 should not break break;
						}
						else if (line_string.find("INITIAL_TIME_STEP") !=
						         std::string::npos)
						{
							*tim_file >> line_string;
							initial_time_step =
							    strtod(line_string.data(), NULL);
							line.clear();
						}
						else if (line_string.find("PID") != std::string::npos)
						{
							usePIDControlInSelfAdaptive = true;
							*tim_file >> pid_variable >> pid_error;
							line.clear();
						}
						/*  //WW
						   if(line_string.find("MINISH")!=string::npos){
						   *tim_file >> line_string;
						   minish = strtod(line_string.data(),NULL);
						   line.clear();
						   }
						 */
						else if (line_string.find("M") == std::string::npos)
						{
							line.str(line_string);
							line >> iter_times;
							line >> multiply_coef;
							time_adapt_tim_vector.push_back(iter_times);
							time_adapt_coe_vector.push_back(multiply_coef);
							line.clear();
						}
					}  // end of while loop adaptive
					   // end of if "SELF_ADAPTIVE"
				}
				else
				{
					ScreenMessage("ERROR: Unrecognized time control type.\n");
					exit(1);
				}
			}  // end of while
		// end of "TIME_CONTROL"
		//....................................................................
		/* //WW
		   if(line_string.find("$SUBSTEPS")!=string::npos) { // subkeyword found
		   JOD 4.7.10
		   *tim_file>>sub_steps>>ws;
		   continue;
		   }
		 */
		//....................................................................
	}  // end of while(!new_keyword)

	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   08/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
**************************************************************************/
bool TIMRead(std::string file_base_name)
{
	//----------------------------------------------------------------------
	// OK  TIMDelete();
	//----------------------------------------------------------------------
	CTimeDiscretization* m_tim = NULL;
	char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::ios::pos_type position;
	//========================================================================
	// File handling
	std::string tim_file_name = file_base_name + TIM_FILE_EXTENSION;
	std::ifstream tim_file(tim_file_name.data(), std::ios::in);
	if (!tim_file.good()) return false;
	tim_file.seekg(0L, std::ios::beg);
	//========================================================================
	// Keyword loop
	ScreenMessage("TIMRead\n");
	while (!tim_file.eof())
	{
		tim_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos) return true;
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#TIME_STEPPING") != std::string::npos)
		{
			m_tim = new CTimeDiscretization();
			position = m_tim->Read(&tim_file);
			m_tim->time_current = m_tim->time_start;
			//----------------------------------------------------------------------
			if (m_tim->Write_tim_discrete)  // YD Write out Time Steps &
			                                // Iterations
			{
				std::string m_file_name = file_base_name + "_TimeDiscrete.txt";
				m_tim->tim_discrete = new std::fstream(
				    m_file_name.c_str(), std::ios::trunc | std::ios::out);
				std::fstream* tim_dis = m_tim->tim_discrete;
				*tim_dis << " No. Time  Tim-Disc  Iter"
				         << "\n";
				if (!m_tim->tim_discrete->good())
					std::cout << "Warning : Time-Discrete files are not found"
					          << "\n";
			}
			//----------------------------------------------------------------------
			time_vector.push_back(m_tim);
			tim_file.seekg(position, std::ios::beg);
		}  // keyword found
	}      // eof
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task: master write function
   01/2005 OK Implementation
   06/2009 OK Write only existing data
**************************************************************************/
void TIMWrite(std::string base_file_name)
{
	//----------------------------------------------------------------------
	if ((int)time_vector.size() < 1) return;
	//----------------------------------------------------------------------
	CTimeDiscretization* m_tim = NULL;
	std::string sub_line;
	std::string line_string;
	//========================================================================
	// File handling
	std::string tim_file_name = base_file_name + TIM_FILE_EXTENSION;
	std::fstream tim_file(tim_file_name.data(),
	                      std::ios::trunc | std::ios::out);
	tim_file.setf(std::ios::scientific, std::ios::floatfield);
	tim_file.precision(12);
	if (!tim_file.good()) return;
	tim_file.seekg(0L, std::ios::beg);
	//========================================================================
	// OUT vector
	tim_file
	    << "GeoSys-TIM: ------------------------------------------------\n";
	for (int i = 0; i < (int)time_vector.size(); i++)
	{
		m_tim = time_vector[i];
		m_tim->Write(&tim_file);
	}
	tim_file << "#STOP";
	tim_file.close();
}

/**************************************************************************
   FEMLib-Method:
   01/2004 OK Implementation
   05/2009 OK $TIME_CONTROL
**************************************************************************/
void CTimeDiscretization::Write(std::fstream* tim_file)
{
	int i;
	//--------------------------------------------------------------------
	// KEYWORD
	*tim_file << "#TIME_STEPPING"
	          << "\n";
	//--------------------------------------------------------------------
	// PCS_TYPE
	*tim_file << " $PCS_TYPE"
	          << "\n";
	*tim_file << "  " << pcs_type_name << "\n";
	//--------------------------------------------------------------------
	*tim_file << " $TIME_START"
	          << "\n";
	*tim_file << "  " << time_start << "\n";
	//--------------------------------------------------------------------
	*tim_file << " $TIME_END"
	          << "\n";
	*tim_file << "  " << time_end << "\n";
	//--------------------------------------------------------------------
	if (time_control_type == TimeControlType::INVALID)
	{
		*tim_file << " $TIME_STEPS"
		          << "\n";
		for (i = 0; i < (int)time_step_vector.size(); i++)
			*tim_file << "  " << 1 << " " << time_step_vector[i] << "\n";
	}
	else
	{
		*tim_file << " $TIME_CONTROL"
		          << "\n";
		*tim_file << "  " << convertTimeControlTypeToString(time_control_type)
		          << "\n";
		//		if(time_control_name == "COURANT_MANIPULATE")
		//		{
		//			*tim_file << "  " << time_control_name << "\n";
		//			*tim_file << "   " << time_control_manipulate << "\n";
		//		}
		if (time_control_type == TimeControlType::PI_AUTO_STEP_SIZE)
		{
			*tim_file << "   " << PI_tsize_ctrl_type << " " << relative_error
			          << " " << absolute_error << " " << this_stepsize << "\n";
		}
		else if (time_control_type == TimeControlType::STEP_SIZE_RESTRICTION)
		{
			*tim_file << "   " << h_min << " " << h_max << "\n";
		}
		else if (time_control_type == TimeControlType::SELF_ADAPTIVE)
		{
			*tim_file << "  MAX_TIME_STEP " << max_time_step << "\n";
			*tim_file << "  MIM_TIME_STEP " << min_time_step << "\n";
		}
	}
	//--------------------------------------------------------------------
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2008 WW Force t+dt be indentical to the time for output or other special
time
        WW Auto time step size control
   03/2012 JT Modify critical time method, clean
**************************************************************************/
double CTimeDiscretization::CalcTimeStep(double current_time)
{
	int no_time_steps = (int)time_step_vector.size();
	//
	// TIME STEP VECTOR
	// -----------------------------------
	if (no_time_steps > 0)
	{
		time_step_length = time_step_vector[0];
		if (step_current < no_time_steps)
			time_step_length = time_step_vector[step_current];
	}
	//
	// TIME CONTROL METHODS
	// -----------------------------------
	if (time_control_type == TimeControlType::NEUMANN ||
	    time_control_type == TimeControlType::SELF_ADAPTIVE)
	{
		if (aktuelle_zeit < MKleinsteZahl && repeat == false)
		{
			time_step_length = FirstTimeStepEstimate();
		}
		else if (time_control_type == TimeControlType::NEUMANN)
		{
			time_step_length = NeumannTimeControl();
		}
		else if (time_control_type == TimeControlType::SELF_ADAPTIVE)
		{
			time_step_length = SelfAdaptiveTimeControl();
		}
	}
	else if (time_control_type == TimeControlType::ERROR_CONTROL_ADAPTIVE)
	{
		if (aktuelle_zeit < MKleinsteZahl)
		{
			time_step_length = AdaptiveFirstTimeStepEstimate();
		}
		else
		{
			time_step_length = ErrorControlAdaptiveTimeControl();
		}
	}
	else if (time_control_type == TimeControlType::PI_AUTO_STEP_SIZE)
	{
		time_step_length = this_stepsize;
	}
	else if (time_control_type == TimeControlType::DYNAMIC_COURANT ||
	         time_control_type == TimeControlType::DYNAMIC_PRESSURE ||
	         time_control_type == TimeControlType::DYNAMIC_VARIABLE)
	{  // JT2012: Soon to come.
		if (!last_dt_accepted)
		{
			time_step_length *= dt_failure_reduction_factor;
			dynamic_minimum_suggestion = time_step_length;
		}
		else if (accepted_step_count > 1)
		{   // initial time step is otherwise maintained for first 2 active time
			// steps
		}
	}
	else if (time_control_type == TimeControlType::PID_CONTROL)
	{
		if (!last_dt_accepted)
		{
			CRFProcess* m_pcs = PCSGet(pcs_type_name);
			double fac = std::min(1. / m_pcs->e_n, 0.8);
			double v = std::max(fac * dt_pre, min_time_step);
			time_step_length = v * v / dt_pre;
		}
		else if (aktuelle_zeit < MKleinsteZahl && !repeat)
			time_step_length = initial_time_step;
		//		else
		//			time_step_length = this_stepsize;
	}
	else if (no_time_steps == 0)
	{  // Processes without time control
		time_step_length =
		    DBL_MAX;  // Large, thus other processes will control the step
	}
	// Restrict by minimum
	if (time_step_length < min_time_step)
	{   // Default value of min_time_step is DBL_EPSILON, unless entered
		// otherwise in the .tim read
		time_step_length = min_time_step;
	}
	// JT: the recommended time step, before critical alteration (for dt
	// control)
	//     otherwise the critical time will govern time change, rather than
	//     primary variable rates.
	recommended_time_step = time_step_length;

	// WW. Critical time match (JT2012 modified)
	// ------------------------------------------------------
	const double tol = time_step_length * 1e-3;
	for (int i = 0; i < (int)critical_time.size(); i++)
	{
		if (current_time < critical_time[i])
		{
			double next = current_time + time_step_length;
			double tval = next + time_step_length / 1.0e3;  // JT2012. A tiny
			                                                // increase in dt is
			                                                // better than a
			                                                // miniscule dt on
			                                                // the next step
			if (tval > critical_time[i])
			{  // Critical time is hit
				bool modify = false;
				if (next - critical_time[i] > .0)
				{  // otherwise, match is already exact
					ScreenMessage(
					    "-> suggested time step exceeds a critical time: %g > "
					    "%g (tol=%g)\n",
					    next, critical_time[i], tol);
					modify = true;
				}
				else if (critical_time[i] - next < tol)
				{
					ScreenMessage(
					    "-> round suggested time step to a critical time: %g > "
					    "%g (tol=%g)\n",
					    next, critical_time[i], tol);
					modify = true;
				}
				if (modify)
				{
					double new_dt = critical_time[i] - current_time;
					ScreenMessage(
					    "-> modify the time step size to match a critical "
					    "time: %g -> %g\n",
					    time_step_length, new_dt);
					time_step_length = new_dt;
				}
				break;
			}
			else if (tval + time_step_length > critical_time[i])
			{   // We can hit the critical time in 2 time steps, smooth the
				// transition
				if (next + time_step_length != critical_time[i])
				{  // otherwise, match is already exact
					ScreenMessage(
					    "-> modify the time step size to hit the critical time "
					    "in 2 time steps: %g -> %g\n",
					    time_step_length,
					    (critical_time[i] - current_time) / 2.0);
					time_step_length = (critical_time[i] - current_time) / 2.0;
				}
				break;
			}
			break;
		}
	}
	//
	next_active_time = current_time + time_step_length;
	return time_step_length;
}

/**************************************************************************
   FEMLib-Method: Operator
   Task:
   Programing:
   08/2008 WW Implementation
**************************************************************************/
CTimeDiscretization::CTimeDiscretization(const CTimeDiscretization& a_tim,
                                         std::string pcsname)
{
	int i;
	safty_coe = a_tim.safty_coe;
	dt_sum = a_tim.dt_sum;
	this_stepsize = a_tim.this_stepsize;
	file_base_name = a_tim.file_base_name;
	time_start = a_tim.time_start;
	time_end = a_tim.time_end;
	time_current = a_tim.time_current;
	time_control_manipulate = a_tim.time_control_manipulate;
	step_current = a_tim.step_current;
	repeat = a_tim.repeat;
	pcs_type_name = pcsname;  // by argument
	time_type_name = a_tim.time_type_name;
	time_control_type = a_tim.time_control_type;
	time_unit = a_tim.time_unit;
	iter_times = a_tim.iter_times;
	multiply_coef = a_tim.multiply_coef;
	max_time_step = a_tim.max_time_step;
	min_time_step = a_tim.min_time_step;
	Write_tim_discrete = a_tim.Write_tim_discrete;
	tim_discrete = a_tim.tim_discrete;
	nonlinear_iteration_error = a_tim.nonlinear_iteration_error;
	//
	time_independence = a_tim.time_independence;
	minimum_dt_reached = a_tim.minimum_dt_reached;
	time_independence = false;
	time_active = true;
	last_active_time = 0.0;
	next_active_time = 0.0;  // JT
	//
	time_step_vector.clear();
	time_adapt_tim_vector.clear();
	time_adapt_coe_vector.clear();
	for (i = 0; i < (int)a_tim.time_step_vector.size(); i++)
		time_step_vector.push_back(a_tim.time_step_vector[i]);
	for (i = 0; i < (int)a_tim.time_adapt_tim_vector.size(); i++)
		time_adapt_tim_vector.push_back(a_tim.time_adapt_tim_vector[i]);
	for (i = 0; i < (int)a_tim.time_adapt_coe_vector.size(); i++)
		time_adapt_coe_vector.push_back(a_tim.time_adapt_coe_vector[i]);
	last_rejected_timestep = 0;
	last_dt_accepted = true;
	usePIDControlInSelfAdaptive = false;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2004 OK Implementation
   last modified:
**************************************************************************/
// kg44 const string gave trouble for me
CTimeDiscretization* TIMGet(const std::string& pcs_type_name)
{
	CTimeDiscretization* m_tim = NULL;
	int i;
	int no_times = (int)time_vector.size();
	for (i = 0; i < no_times; i++)
	{
		m_tim = time_vector[i];
		if (m_tim->pcs_type_name.compare(pcs_type_name) == 0)
			return time_vector[i];
	}
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void TIMDelete()
{
	long i;
	int no_tim = (int)time_vector.size();
	for (i = 0; i < no_tim; i++)
		delete time_vector[i];
	time_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void TIMDelete(std::string pcs_type_name)
{
	long i;
	CTimeDiscretization* m_tim = NULL;
	int no_tim = (int)time_vector.size();
	for (i = 0; i < no_tim; i++)
	{
		m_tim = TIMGet(pcs_type_name);
		if (!m_tim)  // OK
			continue;
		if (m_tim->pcs_type_name.compare(pcs_type_name) == 0)
		{
			delete time_vector[i];
			time_vector.erase(time_vector.begin() + i);
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Neumann estimation
   Programing:
   10/2005 YD Implementation
**************************************************************************/
double CTimeDiscretization::FirstTimeStepEstimate(void)
{
	CMediumProperties* m_mmp = NULL;
	CRFProcess* m_pcs = NULL;
	MeshLib::CElem* elem = NULL;
	int idxS;
	long group;
	double GP[3];
	static double Node_Sat[8];
	double buffer;
	// WW int no_time_steps;
	// WW  int no_processes =(int)pcs_vector.size();
	//	CFluidProperties* m_mfp = NULL; // 2012-08 TF not used
	//	m_mfp = MFPGet("LIQUID");             //WW
	//	double density_fluid = m_mfp->Density(); //WW // TF: set, but never used

	initial_time_step = std::max(initial_time_step, min_time_step);

	for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
	{
		m_pcs = pcs_vector[n_p];
		CFiniteElementStd* fem = m_pcs->GetAssembler();

		time_step_length = initial_time_step;  // take min time step as
		                                       // conservative best guess for
		                                       // testing
		//		switch (m_pcs->pcs_type_name[0]) {
		switch (m_pcs->getProcessType())  // TF
		{
			//		case 'G': // kg44 groudnwater flow ---if steady state, time
			// step should be greater zero...transient flow does not work with
			// adaptive stepping
			case FiniteElement::GROUNDWATER_FLOW:  // TF, if steady state, time
				                                   // step should be greater
				                                   // zero...transient flow does
				                                   // not work with adaptive
				                                   // stepping
				time_step_length = initial_time_step;  // take min time step as
				                                       // conservative best
				                                       // guess for testing
				break;
			//		case 'L': // kg44 liquid flow ---if steady state, time step
			// should be greater zero...transient flow does not work with
			// adaptive stepping
			case FiniteElement::LIQUID_FLOW:  // TF, if steady state, time step
				                              // should be greater
				                              // zero...transient flow does not
				                              // work with adaptive stepping
				time_step_length = initial_time_step;  // take min time step as
				                                       // conservative best
				                                       // guess for testing
				break;
			//		case 'M': // kg44 Mass transport ---if steady state, time
			//step
			// should be greater zero..
			case FiniteElement::MASS_TRANSPORT:  // TF, if steady state, time
				                                 // step should be greater
				                                 // zero..
				time_step_length = initial_time_step;  // take min time step as
				                                       // conservative best
				                                       // guess for testing
				break;
			//		case 'R': // Richards
			case FiniteElement::RICHARDS_FLOW:  // TF
			{
				idxS = m_pcs->GetNodeValueIndex("SATURATION1");
				// WW no_time_steps = 1000000000;           //OK (int)(1.0e10);
				time_step_length = 1.e10;
				size_t mmp_vector_size = mmp_vector.size();
				for (size_t m = 0; m < mmp_vector_size; m++)
					m_mmp = mmp_vector[m];
				const size_t size(m_pcs->m_msh->ele_vector.size());
				for (size_t i = 0; i < size; i++)
				{
					elem = m_pcs->m_msh->ele_vector[i];
					if (elem->GetMark())  // Element selected
					{
						// Activated Element
						group = elem->GetPatchIndex();
						m_mmp = mmp_vector[group];
						m_mmp->m_pcs = m_pcs;
						MshElemType::type EleType = elem->GetElementType();
						// Triangle
						if (EleType == MshElemType::TRIANGLE)
						{
							GP[0] = GP[1] = 0.1 / 0.3;
							GP[2] = 0.0;
						}
						else if (EleType == MshElemType::TETRAHEDRON)
							GP[0] = GP[1] = GP[2] = 0.25;
						else
							GP[0] = GP[1] = GP[2] = 0.0;
					}
					const int vertex_number(elem->GetVertexNumber());
					for (int j = 0; j < vertex_number; j++)
						Node_Sat[j] =
						    m_pcs->GetNodeValue(elem->GetNodeIndex(j), idxS);
					// JT: dSdP now returns actual sign (<0)
					buffer = -m_mmp->PressureSaturationDependency(
					    fem->interpolate(Node_Sat),
					    true);  // JT: now returns correct sign.
					buffer *= 0.5 * elem->GetVolume() * elem->GetVolume();
					buffer *= m_mmp->porosity_model_values[0] *
					          mfp_vector[0]->Viscosity() /
					          m_mmp->permeability_tensor[0];
					buffer /= m_pcs->time_unit_factor;
					time_step_length = MMin(time_step_length, buffer);
				}  // ele_vector

				if (time_step_length < MKleinsteZahl)
				{
					std::cout << "Warning : Time Control Step Wrong, dt = 0.0 "
					          << "\n";
					time_step_length = 1.e-6;
				}
				std::cout << "Neumann Time Step: " << time_step_length << "\n";
				time_step_length_neumann = 1.e10;
				time_step_length = MMin(time_step_length, max_time_step);
				if (Write_tim_discrete)
					*tim_discrete << aktueller_zeitschritt << "  "
					              << aktuelle_zeit << "   " << time_step_length
					              << "  " << m_pcs->iter_lin << "\n";
				break;
			}
			default:
				ScreenMessage(
				    "CTimeDiscretization::FirstTimeStepEstimate default "
				    "case\n");
				break;
		}
	}
	return time_step_length;
}

/**************************************************************************
   FEMLib-Method:
   Task: Control based on primary varable change
   Programing:
   02/2012 JT
**************************************************************************/
double CTimeDiscretization::DynamicVariableTimeControl(void)
{
	long node, gnodes;
	int ii, idx0, idx1, idx[DOF_NUMBER_MAX + 1];
	double error, tol, delta, edof, val, ndof, nerror,
	    dof_error[DOF_NUMBER_MAX + 1];
	double suggested_time_step_change, suggested_time_step;
	CRFProcess* m_pcs = PCSGet(pcs_type_name);
	int num_variables = m_pcs->GetDOF();
	gnodes = (long)m_pcs->m_msh->GetNodesNumber(false);
	error = nerror = 0.0;
	//
	// Indices of primaries
	for (ii = 0; ii < num_variables; ii++)
	{
		idx[ii] =
		    m_pcs->GetNodeValueIndex(m_pcs->pcs_primary_function_name[ii]);
	}
	// Is a third variable wished to be controlled (and is this allowed)
	if (dynamic_control_tolerance[DOF_NUMBER_MAX] > 0.0 &&
	    m_pcs->isPCSMultiFlow)
	{
		if (m_pcs->getProcessType() == FiniteElement::MULTI_PHASE_FLOW)
			idx[num_variables] = m_pcs->GetNodeValueIndex("SATURATION2");
		else
			idx[num_variables] = m_pcs->GetNodeValueIndex("PRESSURE2");
		//
		dynamic_control_tolerance[num_variables] =
		    dynamic_control_tolerance[DOF_NUMBER_MAX];  // shift to 3rd variable
		                                                // slot
		num_variables++;
	}
	//
	// Get error
	switch (FiniteElement::convertErrorMethod(dynamic_control_error_method))
	{
		default:
			return 0.0;
		//
		case FiniteElement::LMAX:  // "LMAX"
			for (ii = 0; ii < num_variables; ii++)
			{
				idx0 = idx[ii];   // old time
				idx1 = idx0 + 1;  // new time
				tol = dynamic_control_tolerance[ii];
				edof = 0.0;
				//
				for (node = 0; node < gnodes; node++)
				{
					delta = fabs(m_pcs->GetNodeValue(node, idx1) -
					             m_pcs->GetNodeValue(node, idx0)) /
					        tol;
					edof = MMax(edof, delta);
				}
				error = MMax(error, edof);
				dof_error[ii] = edof;
			}
			break;
		//
		case FiniteElement::ENORM:  // Error of global norm (single tolerance)
			tol = dynamic_control_tolerance[0];
			//
			for (ii = 0; ii < num_variables; ii++)
			{
				idx0 = idx[ii];   // old time
				idx1 = idx0 + 1;  // new time
				edof = 0.0;
				//
				for (node = 0; node < gnodes; node++)
				{
					delta = m_pcs->GetNodeValue(node, idx1) -
					        m_pcs->GetNodeValue(node, idx0);
					edof += delta * delta;
				}
				error += edof;
				dof_error[ii] = sqrt(edof) / tol;
			}
			error = sqrt(error) / tol;
			break;
		//
		case FiniteElement::EVNORM:  // Error of global norm (DOF specific
			                         // tolerance)
			for (ii = 0; ii < num_variables; ii++)
			{
				idx0 = idx[ii];   // old time
				idx1 = idx0 + 1;  // new time
				tol = dynamic_control_tolerance[0];
				edof = 0.0;
				//
				for (node = 0; node < gnodes; node++)
				{
					delta = m_pcs->GetNodeValue(node, idx1) -
					        m_pcs->GetNodeValue(node, idx0);
					edof += delta * delta;
				}
				dof_error[ii] = sqrt(edof) / tol;
				error = MMax(error, dof_error[ii]);
			}
			break;
		//
		case FiniteElement::ERNORM:  // Error of relative global norm (single
			                         // tolerance)
			tol = dynamic_control_tolerance[0];
			ndof = 0.0;
			//
			for (ii = 0; ii < num_variables; ii++)
			{
				idx0 = idx[ii];   // old time
				idx1 = idx0 + 1;  // new time
				edof = ndof = 0.0;

				//
				for (node = 0; node < gnodes; node++)
				{
					val = m_pcs->GetNodeValue(node, idx1);
					delta = val - m_pcs->GetNodeValue(node, idx0);
					edof += delta * delta;
					ndof += val * val;
				}
				error += edof;
				nerror += ndof;
				dof_error[ii] = (sqrt(edof) / (sqrt(ndof) + DBL_EPSILON)) / tol;
			}
			error = (sqrt(error) / (sqrt(nerror) + DBL_EPSILON)) / tol;
			break;
	}
	// Get suggested time step
	if (error < DBL_EPSILON) error = DBL_EPSILON;
	suggested_time_step_change = time_step_length * (1.0 / error - 1.0);
	suggested_time_step = suggested_time_step_change + time_step_length;
	//
	std::cout << "Dynamic Variable suggested time step:  "
	          << suggested_time_step << "\n";
	if (num_variables > 1)
	{
		for (ii = 0; ii < num_variables; ii++)
		{
			edof = time_step_length +
			       time_step_length * (1.0 / dof_error[ii] - 1.0);
			std::cout << "--> For DOF #:  " << ii
			          << "  suggested time step is: " << edof << "\n";
		}
	}
	//
	// Smooth the time step suggestion
	time_step_length = DynamicTimeSmoothing(suggested_time_step_change);
	return time_step_length;
}

/**************************************************************************
   FEMLib-Method:
   Task: Smooth a time step suggestion from any calculation method
   Programing:
   02/2012 JT
**************************************************************************/
double CTimeDiscretization::DynamicTimeSmoothing(
    double suggested_time_step_change)
{
	double ddt, val, time_step;
	int number_of_time_steps_to_smooth = 2;
	//
	ddt = fabs(suggested_time_step_change);
	val =
	    1.0 -
	    ddt / (time_step_length +
	           ddt);  // how different is the suggestion from the current value
	                  //
	if (suggested_time_step_change > 0.0)  // increase suggested in dt (take at
		                                   // least 1% but not more than 50%)
		val = MRange(0.01, val, 0.5);
	else  // decrease suggested in dt (take at least 20% but not more than 90%)
		val = MRange(0.2, val, 0.9);
	time_step = ddt * val +
	            time_step_length;  // the suggested time step (after weighting)

	// Buffer the time step change over "number_of_time_steps_to_smooth" time
	// steps (taking the worst case)
	dynamic_time_buffer++;
	if (dynamic_time_buffer == 1)  // Not yet tested 2 steps
		dynamic_minimum_suggestion = time_step;
	else
		dynamic_minimum_suggestion =
		    MMin(time_step, dynamic_minimum_suggestion);

	if (dynamic_time_buffer < number_of_time_steps_to_smooth)
	{
		time_step = recommended_time_step;  // recommended_time_step is the
		                                    // m_tim variable of last suggested
		                                    // time step before any critical
		                                    // changes
	}
	else
	{
		dynamic_time_buffer = 0;  // reset
		time_step = dynamic_minimum_suggestion;
		// do not allow to increase faster than this user input factor
		// note, minimum time step is checked in CalcTimeStep()
		time_step =
		    MMin(time_step, recommended_time_step *
		                        dynamic_control_tolerance[DOF_NUMBER_MAX]);
	}
	return time_step;
}

/**************************************************************************
   FEMLib-Method:
   Task: Nuemann Control
   Programing:
   10/2005 YD Implementation
**************************************************************************/
double CTimeDiscretization::NeumannTimeControl(void)
{
	CRFProcess* m_pcs = NULL;

	for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
	{
		m_pcs = pcs_vector[n_p];
		//		switch (m_pcs->pcs_type_name[0]) {
		//				case 'R': // Richards
		switch (m_pcs->getProcessType())  // TF
		{
			case FiniteElement::RICHARDS_FLOW:
				time_step_length = time_step_length_neumann;
				break;
			default:
				std::cout << "Fatal error: No valid PCS type"
				          << "\n";
				break;
		}
	}

	std::cout << "Neumann Time Step: " << time_step_length << "\n";
	time_step_length_neumann = 1.e10;
	if (Write_tim_discrete)
		*tim_discrete << aktueller_zeitschritt << "  " << aktuelle_zeit << "   "
		              << time_step_length << "  " << m_pcs->iter_lin << "\n";
	return time_step_length;
}

/**************************************************************************
   FEMLib-Method:
   Task: Self adaptive method
   Programing:
   10/2005 YD Implementation
   03/2008 HS KG Implementation for Groundwater flow and mass transport
   10/2010 KG updates
**************************************************************************/
double CTimeDiscretization::SelfAdaptiveTimeControl(void)
{
	ScreenMessage("-> calculate the next time step for %s\n",
	              pcs_type_name.c_str());


	// find pcs
	CRFProcess* m_pcs = NULL;
	const FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(pcs_type_name));
	for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
	{
		if (pcs_vector[n_p]->getProcessType() == pcs_type)
			m_pcs = pcs_vector[n_p];
	}
	if (!m_pcs)
	{
		ScreenMessage("-> ERROR in SelfAdaptiveTimeControl(): PCS not found\n");
		return 0.0;
	}

	// get iteration number
	int n_itr = 0;
	if (iteration_type == IterationType::LINEAR)
	{
		n_itr = m_pcs->iter_lin_max;
	}
	else if (iteration_type == IterationType::NONLINEAR)
	{
		n_itr = m_pcs->iter_nlin_max;
	}
	else if (iteration_type == IterationType::COUPLING)
	{
		n_itr = m_pcs->iter_outer_cpl;
	}

	// get the multiplier
	double multiplier = 1.0;
	if (!time_adapt_coe_vector.empty())
		multiplier = time_adapt_coe_vector.back();
	for (std::size_t i = 0; i < time_adapt_tim_vector.size(); i++)
	{
		if (n_itr <= time_adapt_tim_vector[i])
		{
			multiplier = time_adapt_coe_vector[i];
			break;
		}
	}
	//
	if (!m_pcs->accepted) multiplier = time_adapt_coe_vector.back();

	if (usePIDControlInSelfAdaptive && n_itr <= time_adapt_tim_vector.front())
	{
		if (aktueller_zeitschritt > 2)
		{
			double kP = 0.075, kI = 0.175, kD = 0.01;  // TODO
			multiplier = std::pow(m_pcs->e_pre / m_pcs->e_n, kP) *
			             std::pow(1 / m_pcs->e_n, kI) *
			             std::pow(std::pow(m_pcs->e_pre, 2) /
			                          (m_pcs->e_n * m_pcs->e_pre2),
			                      kD);
			ScreenMessage("-> PID: e_n-2=%g, e_n-1=%g, e_n=%g, fac=%g\n",
			              m_pcs->e_pre2, m_pcs->e_pre, m_pcs->e_n, multiplier);
		}
		else
			multiplier = 1.0;
		if (!m_pcs->accepted)
		{
			double fac = std::min(1. / m_pcs->e_n, 0.8);
			double v = std::max(fac * time_step_length, min_time_step);
			multiplier = v * v / time_step_length / time_step_length;
			ScreenMessage("-> PID: rejected. fac=%g\n", multiplier);
		}
	}

	// update the time step length
	time_step_length *= multiplier;
	time_step_length = std::min(time_step_length, max_time_step);
	time_step_length = std::max(time_step_length, min_time_step);

	ScreenMessage(
	    "  ->Self adaptive time stepping: nr. itr=%d, multiplier=%g, next "
	    "dt=%g\n",
	    n_itr, multiplier, time_step_length);

	if (Write_tim_discrete)
#if defined(USE_PETSC)
		if (myrank == 0)
#endif
			*tim_discrete << aktueller_zeitschritt << "  " << aktuelle_zeit
			              << "   " << time_step_length << "  " << n_itr << "\n";

	//	if (time_step_length<=min_time_step) {
	//		ScreenMessage("-> ***ERROR*** Next time step size %g is less than or
	// equal to the given minimum size %g. The simulation is aborted.\n",
	// time_step_length, min_time_step);
	//		exit(1);
	//	}
	return time_step_length;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:CMCD 03/2006
**************************************************************************/
double CTimeDiscretization::CheckCourant(void)
{
	long index;
	long group;
	double velocity[3] = {0., 0., 0.};
	double porosity, vg, advective_velocity, length, courant;
	CRFProcess* m_pcs = NULL;
	m_pcs = PCSGetFluxProcess();
	if (!m_pcs) return 0.0;
	int pcs_no = m_pcs->pcs_number;
	CMediumProperties* m_mmp = NULL;
	MeshLib::CElem* elem = NULL;
	ElementValue* gp_ele;
	long critical_element_no = -1;
	double recommended_time_step = 0.0;
	double stable_time_step = 0.0;
//  int edx;

#ifndef OGS_ONLY_TH
	for (index = 0; index < (long)m_pcs->m_msh->ele_vector.size(); index++)
	{
		elem = m_pcs->m_msh->ele_vector[index];
		length = elem->GetRepLength();
		group = elem->GetPatchIndex();
		m_mmp = mmp_vector[group];
		m_mmp->m_pcs = m_pcs;
		porosity = m_mmp->Porosity(m_mmp->Fem_Ele_Std);
		gp_ele = ele_gp_value[index];  // to get gp velocities
		gp_ele->getIPvalue_vec(pcs_no, velocity);
		vg = MBtrgVec(velocity, 3);
		advective_velocity = vg / porosity;
		// kg44 avoid zero velocity..otherwise stable_time_step is a problem
		if (advective_velocity < DBL_EPSILON) advective_velocity = DBL_EPSILON;
		courant = dt * advective_velocity / length;
		elem->SetCourant(courant);
		//    edx = m_pcs->GetElementValueIndex("COURANT"); //kg44 does this
		//    work?
		//    m_pcs->SetElementValue(index,edx,courant);    // kg44 seems not to
		//    work
		stable_time_step = (1. / courant) * dt;
		if (index == 0) recommended_time_step = stable_time_step;
		if (stable_time_step < recommended_time_step)
		{
			recommended_time_step = stable_time_step;
			critical_element_no = index;
		}
	}
#endif
	std::cout << "Courant time step control, critical element = "
	          << critical_element_no << " Recomended time step "
	          << recommended_time_step << "\n";
	return recommended_time_step;
}

/**************************************************************************
   FEMLib-Method:
   Task: Neumann estimation
   Programing:
   04/2006 YD Implementation
**************************************************************************/
double CTimeDiscretization::AdaptiveFirstTimeStepEstimate(void)
{
	CNumerics* m_num(num_vector[0]);
	MeshLib::CElem* elem = NULL;
	static double Node_p[8];
	double p_ini, buff = 0.0;
	// WW int no_time_steps;
	safty_coe = 5.0;
	p_ini = 1.0e-10;

	for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
	{
		CRFProcess* m_pcs = pcs_vector[n_p];
		CFiniteElementStd* fem = m_pcs->GetAssembler();

		//  switch(m_pcs->pcs_type_name[0]){
		switch (m_pcs->getProcessType())  // TF
		{
			//  case 'R': // Richards
			case FiniteElement::RICHARDS_FLOW:  // TF
			{
				int idxp = m_pcs->GetNodeValueIndex("PRESSURE1") + 1;
				// WW no_time_steps = 1000000000;           //OK (int)(1e10);
				time_step_length = 1.e10;
				for (size_t i = 0; i < m_pcs->m_msh->ele_vector.size(); i++)
				{
					elem = m_pcs->m_msh->ele_vector[i];
					for (int j = 0; j < elem->GetVertexNumber(); j++)
						Node_p[j] =
						    m_pcs->GetNodeValue(elem->GetNodeIndex(j), idxp);
					p_ini = MMax(fabs(fem->interpolate(Node_p)), p_ini);
				}
				buff = safty_coe * sqrt(m_num->nls_error_tolerance[0] / p_ini);
				buff /= m_pcs->time_unit_factor;
				time_step_length = MMin(time_step_length, buff);
				if (time_step_length < MKleinsteZahl)
				{
					std::cout << "Warning : Time Control Step Wrong, dt = 0.0 "
					          << "\n";
					time_step_length = 1.0e-8;
				}
				std::cout << "Error Control Time Step: " << time_step_length
				          << "\n";
				if (Write_tim_discrete)
					*tim_discrete << aktueller_zeitschritt << "  "
					              << aktuelle_zeit << "   " << time_step_length
					              << "  " << m_pcs->iter_lin << "\n";
				break;
			}
			//		case 'M': // kg44 mass transport
			case FiniteElement::MASS_TRANSPORT:    // TF
				time_step_length = min_time_step;  // take min time step as
				                                   // conservative best guess
				                                   // for testing
				break;
			//		case 'G': // kg44 groudnwater flow ---if steady state, time
			// step should be greater zeor...transient flow does not work with
			// adaptive stepping
			case FiniteElement::GROUNDWATER_FLOW:  // if steady state, time step
				                                   // should be greater zero ...
				                                   // transient flow does not
				                                   // work with adaptive
				                                   // stepping
				time_step_length = min_time_step;  // take min time step as
				                                   // conservative best guess
				                                   // for testing
				break;
			default:
				break;
		}
	}
	return time_step_length;
}

/**************************************************************************
   FEMLib-Method:
   Task: Error control adaptive method
   Programing:
   04/2006 YD Implementation
**************************************************************************/
double CTimeDiscretization::ErrorControlAdaptiveTimeControl(void)
{
	CRFProcess* m_pcs = NULL;
	double rmax = 5.0;
	double rmin = 0.5;
	double safty_coe = 0.8;

	for (size_t n_p = 0; n_p < pcs_vector.size(); n_p++)
	{
		m_pcs = pcs_vector[n_p];
		//		switch (m_pcs->pcs_type_name[0]) {
		switch (m_pcs->getProcessType())  // TF
		{
			default:
				std::cout << "Fatal error: No valid PCS type"
				          << "\n";
				break;
			//		case 'R': // Richards, accepted and refused time step
			case FiniteElement::RICHARDS_FLOW:  // accepted and refused time
				                                // step
				// nonlinear_iteration_error = m_pcs->nonlinear_iteration_error;
				if (repeat)
					time_step_length *= MMax(
					    safty_coe * sqrt(m_pcs->m_num->nls_error_tolerance[0] /
					                     nonlinear_iteration_error),
					    rmin);
				else
					time_step_length *= MMin(
					    safty_coe * sqrt(m_pcs->m_num->nls_error_tolerance[0] /
					                     nonlinear_iteration_error),
					    rmax);
				std::cout << "Error_Self_Adaptive Time Step: "
				          << time_step_length << "\n";
				if (Write_tim_discrete)
					*tim_discrete << aktueller_zeitschritt << "  "
					              << aktuelle_zeit << "   " << time_step_length
					              << "  " << m_pcs->iter_lin << "\n";
		}
	}
	return time_step_length;
}

/**************************************************************************
FEMLib-Method:
Task: Allow for time step retry if system changed rapidly and non-linearly
    :: Currently will only allow this if DynamicVariable time control is used
Programing:
02.2011 JT implementation
**************************************************************************/
bool CTimeDiscretization::isDynamicTimeFailureSuggested(CRFProcess* m_pcs)
{
	//
	if (time_control_type == TimeControlType::DYNAMIC_VARIABLE)
	{
		if (m_pcs->iter_nlin > 0 || this->minimum_dt_reached)
			return false;  // only checking on zeroth iteration (also do not
		                   // fail if already at the minimum allowed dt)
		//
		FiniteElement::ErrorMethod method(FiniteElement::convertErrorMethod(
		    this->dynamic_control_error_method));
		if (method == m_pcs->m_num->getNonLinearErrorMethod())
		{  // Same as NL method, can re-use the values
			for (int ii = 0; ii < m_pcs->pcs_num_dof_errors; ii++)
			{
				if (this->dynamic_failure_threshold >
				    m_pcs->pcs_absolute_error[ii] /
				        this->dynamic_control_tolerance[ii])
				{
					return true;
				}
			}
		}
		else if (method == m_pcs->m_num->getCouplingErrorMethod())
		{  // Alright, is different from NLS method. How about CPL method?
			for (int ii = 0; ii < m_pcs->cpl_num_dof_errors; ii++)
			{
				if (this->dynamic_failure_threshold >
				    m_pcs->cpl_absolute_error[ii] /
				        this->dynamic_control_tolerance[ii])
				{
					return true;
				}
			}
		}
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03.3012. WW
		else
		{   // Alright, this is annoying. Unfortunately we have to recalculate
			// the node errors.
			m_pcs->CalcIterationNODError(method, false, false);
			for (int ii = 0; ii < m_pcs->temporary_num_dof_errors; ii++)
			{
				if (this->dynamic_failure_threshold >
				    m_pcs->temporary_absolute_error[ii] /
				        this->dynamic_control_tolerance[ii])
				{
					return true;
				}
			}
		}
#endif
	}
	else if (time_control_type == TimeControlType::SELF_ADAPTIVE)
	{
		int n_itr = 0;
		if (iteration_type == IterationType::LINEAR)
		{
			n_itr = m_pcs->iter_lin_max;
		}
		else if (iteration_type == IterationType::NONLINEAR)
		{
			n_itr = m_pcs->iter_nlin_max;
		}
		if (usePIDControlInSelfAdaptive)
		{
			CRFProcess* m_pcs = PCSGet(pcs_type_name);
			if (pid_variable == "VELOCITY1")
			{
				m_pcs->e_n = m_pcs->CalcVelocityChanges();
			}
			else
			{
				int ii = m_pcs->GetNodeValueIndex(pid_variable);
				if (ii < 0)
				{
					static bool done = false;
					if (!done)
					{
						done = true;
						ScreenMessage("Given PID variable name not found: %s\n",
						              pid_variable.c_str());
					}
				}
				if (ii >= 0) m_pcs->e_n = m_pcs->CalcNodeValueChanges(ii);
			}
			m_pcs->e_n /= pid_error;
			ScreenMessage("-> PID: e_n=%g\n", m_pcs->e_n);
			if (aktueller_zeitschritt > 2 && m_pcs->e_n > 1.0) return true;
		}
		if (n_itr >= time_adapt_tim_vector.back()) return true;
	}
	else if (time_control_type == TimeControlType::PID_CONTROL)
	{
		return !PIDControl();
	}

	//
	return false;
}

/**************************************************************************
   FEMLib-Method:
   Task:  Used to force time steps matching the times requried by output or
       boundary
   Programing:
   08/2008 WW Implementation
**************************************************************************/
void CTimeDiscretization::FillCriticalTime()
{
	for (size_t i = 0; i < out_vector.size(); i++)
	{
		COutput* a_out = out_vector[i];
		for (size_t j = 0; j < a_out->getTimeVector().size(); j++)
		{
			bool done = false;
			for (size_t k = 0; k < critical_time.size(); k++)
				if (fabs(critical_time[k] - a_out->getTimeVector()[j]) <
				    DBL_MIN)
				{
					done = true;
					break;
				}
			if (!done) critical_time.push_back(a_out->getTimeVector()[j]);
		}
	}
	// Sort
	for (size_t i = 0; i < critical_time.size(); i++)
		for (size_t j = i; j < critical_time.size(); j++)
			if (critical_time[i] > critical_time[j])
				//				double val = critical_time[i];
				//				critical_time[i] = critical_time[j];
				//				critical_time[j] = val;
				std::swap(critical_time[i], critical_time[j]);
}

bool CTimeDiscretization::PIDControl()
{
	//-----------------------------------------------------
	if (aktueller_zeitschritt == 1)
		dt_pre = initial_time_step;
	else
		dt_pre = time_step_length;
	CRFProcess* m_pcs = PCSGet(pcs_type_name);
	int ii = 0;
	if (m_pcs->getProcessType() == FiniteElement::LIQUID_FLOW)
	{
		m_pcs->e_n = m_pcs->CalcVelocityChanges();
	}
	else
	{
		m_pcs->e_n = m_pcs->CalcNodeValueChanges(ii);
	}
	m_pcs->e_n /= pid_error;
	// std::cout << "-> e_n = " << m_pcs->e_n << "\n";

	//-----------------------------------------------------
	bool accept = true;
	double fac = 1.0;
	//	if (m_pcs->e_n > pid_error)
	//	{
	//		accept = false;
	//		fac *= pid_error / m_pcs->e_n;
	//	}
	//	else
	if (aktueller_zeitschritt > 2)
	{
		accept = true;
		double kP = 0.075, kI = 0.175, kD = 0.01;
		fac = std::pow(m_pcs->e_pre / m_pcs->e_n, kP) *
		      std::pow(1 / m_pcs->e_n, kI) *
		      std::pow(std::pow(m_pcs->e_pre, 2) / (m_pcs->e_n * m_pcs->e_pre2),
		               kD);
		// std::cout << "-> fac = " << fac << "\n";
	}

	//-----------------------------------------------------
	dt = dt_pre * fac;
	dt = std::max(time_step_length, min_time_step);
	dt = std::min(time_step_length, max_time_step);
	//	time_step_length = std::max(time_step_length, dt_pre*l);
	//	time_step_length = std::min(time_step_length, dt_pre*L);

	return accept;
}
