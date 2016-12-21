/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_num_new.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <string>

#include "makros.h"
#include "memory.h"
#include "display.h"
#include "FileToolsRF.h"
#include "StringTools.h"

#include "mathlib.h"

#include "rf_pcs.h"
#include "tools.h"

using namespace std;

vector<CNumerics*> num_vector;

/**************************************************************************
   FEMLib-Method:
   Task: constructor
   Programing:
   11/2004 OK Implementation
   10/2005 OK pcs_type_name
   07/2007 OK DEFORMATION
**************************************************************************/
CNumerics::CNumerics(string const& name)
{
	pcs_type_name = name;

	//
	// LS - Linear Solver
	ls_method = 2;
	ls_max_iterations = 1000;
	ls_error_method = 1;
	ls_error_tolerance = 1e-12;
	ls_theta = 1.0;
	ls_precond = 1;
	ls_storage_method = 2;
	ls_extra_arg = "";
#ifdef USE_PETSC
	petsc_split_fields = false;
	petsc_use_snes = false;
#endif
	//
	// NLS - Nonlinear Solver
	nls_method_name = "PICARD";
	nls_method = FiniteElement::NL_LINEAR;
	nls_error_method = 1;
	nls_max_iterations = 1;
	nls_relaxation = 0.0;
	for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
		nls_error_tolerance[i] = -1.0;  // JT2012: should not default this.
	                                    // Should always be entered by user!
	nls_jacobian_level = 0;             // full
	//
	// CPL - Coupled processes
	cpl_error_specified = false;
	cpl_master_process = false;
	cpl_process = "INVALID_PROCESS";  // JT2012: do not couple with any process,
	                                  // unless indicated
	cpl_variable = "NONE";
	cpl_variable_JOD = "FLUX";
	cpl_max_iterations = 1;
	cpl_min_iterations = 1;
	for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
		cpl_error_tolerance[i] = -1.0;  // JT2012: should not default this.
	                                    // Should always be entered by user!
	//
	// ELE
	ele_gauss_points = 3;
	ele_mass_lumping = 0;
	ele_upwind_method = 0;
	ele_upwinding = 0;
	ele_supg_method = 0;
	ele_supg_method_length = 0;
	ele_supg_method_diffusivity = 0;
	fct_method = -1;
	fct_prelimiter_type = 0;
	fct_const_alpha = -1.0;
	//----------------------------------------------------------------------
	// Deformation
	DynamicDamping = NULL;
	//
	_pcs_cpl_error_method = FiniteElement::LMAX;
	_pcs_nls_error_method = FiniteElement::LMAX;
	nls_plasticity_local_tolerance = 1e-12;
}

/**************************************************************************
   FEMLib-Method:
   Task: deconstructor
   Programing:
   11/2004 OK Implementation
**************************************************************************/
CNumerics::~CNumerics(void)
{
	delete[] DynamicDamping;
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
bool NUMRead(string file_base_name)
{
	//----------------------------------------------------------------------
	// OK  NUMDelete();
	//----------------------------------------------------------------------
	CNumerics* m_num = NULL;
	char line[MAX_ZEILE];
	bool overall_coupling_exists = false;  // JT
	string sub_line;
	string line_string;
	ios::pos_type position;
	//========================================================================
	// File handling
	string num_file_name = file_base_name + NUM_FILE_EXTENSION;
	ifstream num_file(num_file_name.data(), ios::in);
	if (!num_file.good()) return false;
	num_file.seekg(0L, ios::beg);
	//========================================================================
	// Keyword loop
	ScreenMessage("NUMRead\n");
	while (!num_file.eof())
	{
		num_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos) return true;
		//
		if (line_string.find("$OVERALL_COUPLING") != string::npos)
		{
			overall_coupling_exists = true;  // JT: for error checking
		}
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#NUMERICS") != string::npos)
		{
			m_num = new CNumerics("default");
			position = m_num->Read(&num_file);
			num_vector.push_back(m_num);
			num_file.seekg(position, ios::beg);
			m_num->NumConfigure(overall_coupling_exists);  // JT2012
		}                                                  // keyword found
	}                                                      // eof
	return true;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 WW Implementation
**************************************************************************/
bool CNumerics::CheckDynamic()
{
	if (DynamicDamping)
		return true;
	else
		return false;
}

/**************************************************************************
FEMLib-Method:
Task: After a read of each m_num, configure any defaults here
Programing:
05/2011 JT Implementation
**************************************************************************/
void CNumerics::NumConfigure(bool overall_coupling_exists)
{
	// Overall coupling check
	if (overall_coupling_exists && !cpl_error_specified)
	{
		std::cout << "ERROR in NUMRead. Overall coupling requested, but ";
		std::cout << this->pcs_type_name << " was not\n";
		std::cout << "supplied with coupling tolerance. See "
					 "$COUPLING_CONTROL keyword to enter this.\n";
		exit(1);
	}
	//
	// Check master processes
	if (cpl_master_process && !cpl_error_specified)
	{
		std::cout << "ERROR in NUMRead. Process coupling requested, but ";
		std::cout << this->pcs_type_name << " was not\n";
		std::cout << "supplied with coupling tolerance. See $COUPLING_CONTROL "
		             "keyword to enter this.\n";
		exit(1);
	}
	//
	// We are ok. Now check the tolerances.
	if (this->nls_method == FiniteElement::NL_LINEAR)
	{  // linear solution
		if (cpl_error_specified)
		{  // A coupling error was entered. Adopt this for error calculations.
			for (size_t i = 0; i < DOF_NUMBER_MAX; i++)
			{
				nls_error_tolerance[i] = cpl_error_tolerance[i];
			}
			setNonLinearErrorMethod(getCouplingErrorMethod());
		}
		else
		{  // We have no error tolerances for non-linear or coupled simulations.
			// Force some defaults.
			setNonLinearErrorMethod(FiniteElement::LMAX);
			setCouplingErrorMethod(FiniteElement::LMAX);
			nls_error_tolerance[0] = cpl_error_tolerance[0] = 1.0;
		}
	}
	// Default CPL error method to NLS method. Just so error is not checked
	// twice
	if (!cpl_error_specified)
	{
		setCouplingErrorMethod(getNonLinearErrorMethod());
	}
	//
	// Default all NLS tolerances to the previous DOF, if they were not entered.
	for (size_t i = 1; i < DOF_NUMBER_MAX; i++)
	{
		if (nls_error_tolerance[i] < 0.0)
			nls_error_tolerance[i] = nls_error_tolerance[i - 1];
	}
	//
	// Default all CPL tolerances to the previous DOF, if they were not entered.
	for (size_t i = 1; i < DOF_NUMBER_MAX; i++)
	{
		if (cpl_error_tolerance[i] < 0.0)
			cpl_error_tolerance[i] = cpl_error_tolerance[i - 1];
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: OBJ read function
   Programing:
   11/2004 OK Implementation
**************************************************************************/
ios::pos_type CNumerics::Read(ifstream* num_file)
{
	string line_string;
	std::string error_method_name;
	std::string coupling_target;
	bool new_keyword = false;
	bool new_subkeyword = false;
	ios::pos_type position;
	ios::pos_type position_subkeyword;
	std::stringstream line;
	//========================================================================
	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		if (new_subkeyword) num_file->seekg(position, ios::beg);
		new_subkeyword = false;
		position = GetNextSubKeyword(num_file, &line_string, &new_keyword);
		if (new_keyword) return position;
		//....................................................................
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> pcs_type_name;
			line.clear();
			continue;
		}
		//....................................................................
		// JT->WW: Local tolerance previously found in $NON_LINEAR_SOLVER for
		// NEWTON. Moved here for now.
		if (line_string.find("$PLASTICITY_TOLERANCE") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> nls_plasticity_local_tolerance;
		}
		//....................................................................
		// subkeyword found ($NON_LINEAR_ITERATION  -or- $NON_LINEAR_ITERATIONS)
		if (line_string.find("$NON_LINEAR_ITERATION") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> nls_method_name;
			line >> error_method_name;
			line >> nls_max_iterations;
			line >> nls_relaxation;
			//
			setNonLinearErrorMethod(
			    FiniteElement::convertErrorMethod(error_method_name));
			switch (getNonLinearErrorMethod())
			{
				case FiniteElement::ENORM:  // only 1 tolerance required
					line >> nls_error_tolerance[0];
					break;
				//
				case FiniteElement::ERNORM:  // only 1 tolerance required
					line >> nls_error_tolerance[0];
					break;
				//
				case FiniteElement::EVNORM:  // 1 tolerance for each primary
					                         // variable (for Deformation, only
					                         // 1 tolerance required. Applies to
					                         // x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> nls_error_tolerance[i];
					break;
				//
				case FiniteElement::LMAX:  // 1 tolerance for each primary
					                       // variable (for Deformation, only 1
					                       // tolerance required. Applies to
					                       // x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> nls_error_tolerance[i];
					break;
				//
				case FiniteElement::BNORM:  // only 1 tolerance required
					line >> nls_error_tolerance[0];
					break;
				//
				default:
					ScreenMessage(
					    "ERROR in NUMRead. Invalid non-linear iteration error "
					    "method selected.\n");
					exit(1);
					break;
			}

			if (nls_method_name.find("LINEAR") != string::npos)
				nls_method = FiniteElement::NL_LINEAR;
			else if (nls_method_name.find("PICARD") != string::npos)
				nls_method = FiniteElement::NL_PICARD;
			else if (nls_method_name.find("NEWTON") != string::npos)
				nls_method = FiniteElement::NL_NEWTON;
			//
			line.clear();
			continue;
		}
		else if (line_string.find("$NON_LINEAR_SOLVER") != string::npos)
		{
			ScreenMessage(" --\n Using old $NON_LINEAR_SOLVER keyword.\n");
			ScreenMessage(
			    " Eventually this will be obsolete. Consider switching to\n");
			ScreenMessage(
			    " $NON_LINEAR_ITERATIONS for better results and greater "
			    "flexibility.\n");
			ScreenMessage(" --\n");
			//
			// JT:	in >> method_name
			//		in >> tolerance
			//		if(NEWTON) in >> tolerance_local
			//		in >> max iter
			//		in >> relaxation
			//
			//
			line.str(GetLineFromFile1(num_file));
			line >> nls_method_name;
			//
			if (nls_method_name.find("PICARD") != string::npos)
				nls_method = FiniteElement::NL_PICARD;
			else if (nls_method_name.find("NEWTON") != string::npos)
				nls_method = FiniteElement::NL_NEWTON;
			//
			if (FiniteElement::isNewtonKind(nls_method))
			{
				line >> nls_error_tolerance[0];
				line >> nls_plasticity_local_tolerance;
				error_method_name =
				    "BNORM";  // JT: this is hardwired in old version
			}
			else
			{
				line >> nls_error_tolerance[0];
				error_method_name =
				    "LMAX";  // JT: this is hardwired in old version
			}
			setNonLinearErrorMethod(
			    FiniteElement::convertErrorMethod(error_method_name));
			//
			line >> nls_max_iterations;
			line >> nls_relaxation;
			line.clear();
			continue;
		}
		//....................................................................
		// JT->WW: Local tolerance previously found in $NON_LINEAR_SOLVER for
		// NEWTON. Moved here for now.
		if (line_string.find("$JACOBIAN_LEVEL") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> nls_jacobian_level;
			ScreenMessage("-> $JACOBIAN_LEVEL = %d\n", nls_jacobian_level);
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$LINEAR_SOLVER") != string::npos)
		{
			std::string str_buf = GetLineFromFile1(num_file);
			line.str(str_buf);
			if (str_buf.find("petsc") != string::npos)
			{
				line >> str_buf >> ls_solver_name >> ls_precond_name >>
				    ls_error_tolerance >> ls_max_iterations >> ls_theta;
			}
			else
			{
				line >> ls_method;
				line >> ls_error_method;
				line >> ls_error_tolerance;
				line >> ls_max_iterations;
				line >> ls_theta;
				line >> ls_precond;
				line >> ls_storage_method;
#ifdef USE_PETSC
				ls_solver_name = "bcgs";
				ls_precond_name = "bjacobi";
				if (ls_method == 4)
					ls_solver_name = "bcgs";
				else if (ls_method == 9)
					ls_solver_name = "gmres";
				if (ls_precond == 1)
					ls_precond_name = "bjacobi";
				else if (ls_precond == 2) // ILU
					ls_precond_name = "bjacobi";
				else if (ls_precond == 9) // ILUT
					ls_precond_name = "bjacobi";

				ScreenMessage("-> configure linear solver from Lis options (solver=%s, precon=%s)\n", ls_solver_name.data(), ls_precond_name.data());
#endif
			}
			line.clear();
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$PETSC_LINEAR_SOLVER") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ls_solver_name >> ls_precond_name >> ls_error_tolerance >> ls_max_iterations;
			line.clear();
			continue;
		}
		//....................................................................
		// JT subkeyword found
		if (line_string.find("$COUPLING_ITERATIONS") != string::npos)
		{
			ScreenMessage("$COUPLING_ITERATIONS keyword obsolete.\n");
			ScreenMessage(
			    "Use $COUPLING_CONTROL and $COUPLED_PROCESS for process "
			    "couplings.\n");
			exit(1);
		}
		//....................................................................
		// JT subkeyword found
		if (line_string.find("$COUPLING_CONTROL") !=
		    string::npos)  // JT: For this process, how is coupling error
		                   // handled?
		{
			// JT:	in >> error_method_name
			//		in >> tolerance[1:dof]
			//
			line.str(GetLineFromFile1(num_file));
			line >> error_method_name;
			//
			cpl_error_specified = true;
			setCouplingErrorMethod(
			    FiniteElement::convertErrorMethod(error_method_name));
			switch (getCouplingErrorMethod())
			{
				case FiniteElement::ENORM:  // only 1 tolerance required
					line >> cpl_error_tolerance[0];
					break;
				//
				case FiniteElement::ERNORM:  // only 1 tolerance required
					line >> cpl_error_tolerance[0];
					break;
				//
				case FiniteElement::EVNORM:  // 1 tolerance for each primary
					                         // variable (for Deformation, only
					                         // 1 tolerance required. Applies to
					                         // x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> cpl_error_tolerance[i];
					break;
				//
				case FiniteElement::LMAX:  // 1 tolerance for each primary
					                       // variable (for Deformation, only 1
					                       // tolerance required. Applies to
					                       // x,y,z)
					for (int i = 0; i < DOF_NUMBER_MAX; i++)
						line >> cpl_error_tolerance[i];
					break;
				//
				case FiniteElement::BNORM:
					ScreenMessage(
					    "ERROR in NUMRead. BNORM not configured for process "
					    "couplings.\n");
					ScreenMessage(
					    "We suggest ENORM as a valid companion for NEWTON "
					    "couplings.\n");
					exit(1);
					break;
				//
				default:
					ScreenMessage(
					    "ERROR in NUMRead. Invalid coupling error method "
					    "selected.\n");
					exit(1);
					break;
			}
			//
			line.clear();
			continue;
		}
		//....................................................................
		// JT subkeyword found
		if (line_string.find("$COUPLED_PROCESS") !=
		    string::npos)  // JT: Is this process coupled to another process in
		                   // an inner loop?
		{
			// in >> process name >> min iter >> max iter
			//
			line.str(GetLineFromFile1(num_file));
			line >> coupling_target;  // name of coupled process -OR- process
			                          // variable
			line >> cpl_min_iterations;
			line >> cpl_max_iterations;
			//
			cpl_master_process = true;
			//
			// Is coupling through a process name or a primary variable?
			if (FiniteElement::convertPrimaryVariable(coupling_target) !=
			    FiniteElement::INVALID_PV)
			{  // Then a valid process VARIABLE is entered. Use this.
				cpl_variable = coupling_target;
			}
			else if (PCSGet(coupling_target))
			{  // Then a valid process is entered
				cpl_process = coupling_target;
			}
			else
			{
				ScreenMessage(
				    "WARNING. $COUPLED_PROCESS keyword encountered, but a "
				    "valid process OR primary variable was not found.\n");
				cpl_master_process = false;
			}
			//
			line.clear();
			continue;
		}
		//....................................................................
		if (line_string.find("$EXTERNAL_SOLVER_OPTION") !=
		    string::npos)  // subkeyword found
		{
			ls_extra_arg = GetLineFromFile1(num_file);
			trim(ls_extra_arg);
			continue;
		}
		//....................................................................
		// subkeyword found
		if (line_string.find("$ELE_GAUSS_POINTS") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ele_gauss_points;  // probably element-type-wise
			line.clear();
			continue;
		}
		// subkeyword found
		if (line_string.find("$ELE_MASS_LUMPING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> ele_mass_lumping;
			line.clear();
			ScreenMessage("-> Mass Lumping method is selected.\n");
			continue;
		}
		// subkeyword found
		if (line_string.find("$ELE_UPWINDING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			// CB now read also upwinding method
			line >> ele_upwinding >> ele_upwind_method;
			line.clear();
			continue;
		}
		// subkeyword found
		if (line_string.find("$ELE_SUPG") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			// NW
			line >> ele_supg_method >> ele_supg_method_length >>
			    ele_supg_method_diffusivity;
			line.clear();
			ScreenMessage("-> SUPG method is selected.\n");
			continue;
		}
		// subkeyword found
		if (line_string.find("$DYNAMIC_DAMPING") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));  // WW
			DynamicDamping = new double[3];
			// Default
			DynamicDamping[0] = 0.515;
			DynamicDamping[1] = 0.51;
			DynamicDamping[2] = 0.51;
			line >> DynamicDamping[0] >> DynamicDamping[1] >> DynamicDamping[2];
			line.clear();
			continue;
		}
		// Flux corrected transport by Kuzmin (2009)
		// NW
		if (line_string.find("$FEM_FCT") != string::npos)
		{
			line.str(GetLineFromFile1(num_file));
			line >> fct_method;  // 1: linearized FCT
			line >>
			    fct_prelimiter_type;  // 0: just cancel, 1: minmod, 2: superbee
			line >>
			    fct_const_alpha;  //-1: off, [0.0,1.0] 0: Upwind, 1: Galerkin
			line.clear();
			ScreenMessage("-> FEM_FCT method is selected.");
			continue;
		}
#ifdef USE_PETSC
		if (line_string.find("$PETSC_SPLIT_FIELDS") != string::npos)
		{
			this->petsc_split_fields = true;
			ScreenMessage("-> Use field splits in PETSc\n");
			continue;
		}
#endif

		//....................................................................
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   Task: master write function
   Programing:
   11/2004 OK Implementation
   last modification:
**************************************************************************/
void NUMWrite(string base_file_name)
{
	CNumerics* m_num = NULL;
	string sub_line;
	string line_string;
	//========================================================================
	// File handling
	string num_file_name = base_file_name + NUM_FILE_EXTENSION;
	fstream num_file(num_file_name.data(), ios::trunc | ios::out);
	num_file.setf(ios::scientific, ios::floatfield);
	num_file.precision(12);
	if (!num_file.good()) return;
	num_file.seekg(0L, ios::beg);
	//========================================================================
	num_file << "GeoSys-NUM: Numerics "
	            "------------------------------------------------\n";
	//========================================================================
	// OUT vector
	int num_vector_size = (int)num_vector.size();
	int i;
	for (i = 0; i < num_vector_size; i++)
	{
		m_num = num_vector[i];
		m_num->Write(&num_file);
	}
	num_file << "#STOP";
	num_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task: write function
   Programing:
   11/2004 OK Implementation
   last modification:
**************************************************************************/
void CNumerics::Write(fstream* num_file)
{
	// KEYWORD
	*num_file << "#NUMERICS" << endl;
	//--------------------------------------------------------------------
	/*OK
	   *num_file << " $METHOD" << endl;
	   *num_file << method_name << endl;
	   if(method_name.find("LAGRANGE")!=string::npos){
	   *num_file << lag_quality << " " << lag_max_steps << " " << lag_local_eps
	   << " ";
	   *num_file << lag_time_weighting << " " << lag_min_weight << " ";
	   *num_file << lag_use_matrix << " " << lag_vel_method;
	   *num_file << endl;
	   }
	 */
	//--------------------------------------------------------------------
	*num_file << " $PCS_TYPE" << endl;
	*num_file << "  " << pcs_type_name << endl;
	//--------------------------------------------------------------------
	*num_file << " $NON_LINEAR_SOLVER" << endl;
	*num_file << "  " << nls_method_name;
	*num_file << " " << nls_error_tolerance;
	*num_file << " " << nls_max_iterations;
	*num_file << " " << nls_relaxation;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $LINEAR_SOLVER" << endl;
	*num_file << "  " << ls_method;
	*num_file << " " << ls_error_method;
	*num_file << " " << ls_error_tolerance;
	*num_file << " " << ls_max_iterations;
	*num_file << " " << ls_theta;
	*num_file << " " << ls_precond;
	*num_file << " " << ls_storage_method;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $ELE_GAUSS_POINTS" << endl;
	*num_file << "  " << ele_gauss_points;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $ELE_MASS_LUMPING" << endl;
	*num_file << "  " << ele_mass_lumping;
	*num_file << endl;
	//--------------------------------------------------------------------
	*num_file << " $ELE_UPWINDING" << endl;
	*num_file << "  " << ele_upwinding;
	*num_file << endl;
	//--------------------------------------------------------------------
}


int GetNumericsGaussPoints(int element_dimension)
{
	int m_gaussian_points = 3;
	int g_gaussian_points = 3;
	switch (element_dimension)
	{
		case 1:
			m_gaussian_points = 1;
			break;
		case 2:
			m_gaussian_points = g_gaussian_points;
			break;
		case 3:
			m_gaussian_points = g_gaussian_points;
			break;
		case 4:
			m_gaussian_points = g_gaussian_points;
			break;
	}
	return m_gaussian_points;
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void NUMDelete()
{
	long i;
	int no_num = (int)num_vector.size();
	for (i = 0; i < no_num; i++)
		delete num_vector[i];
	num_vector.clear();
}
