/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef rf_num_new_INC
#define rf_num_new_INC

#include <fstream>
#include <string>
#include <vector>

#include "makros.h"
#include "FEMEnums.h"

#define NUM_FILE_EXTENSION ".num"

class CNumerics
{
public:
	explicit CNumerics(std::string const& name);
	~CNumerics(void);

	std::ios::pos_type Read(std::ifstream*);
	void Write(std::fstream*);

	const char* getLinearSolverName() const { return ls_sover_name.c_str(); }
	const char* getPreconditionerName() const { return ls_precond_name.c_str(); }

	void setNonLinearErrorMethod(FiniteElement::ErrorMethod err_method)
	{
		_pcs_nls_error_method = err_method;
	}
	FiniteElement::ErrorMethod getNonLinearErrorMethod() const
	{
		return _pcs_nls_error_method;
	}

	void setCouplingErrorMethod(FiniteElement::ErrorMethod cpl_method)
	{
		_pcs_cpl_error_method = cpl_method;
	}
	FiniteElement::ErrorMethod getCouplingErrorMethod() const
	{
		return _pcs_cpl_error_method;
	}

	void NumConfigure(bool overall_coupling_exists);

private:
	FiniteElement::ErrorMethod _pcs_nls_error_method;
	FiniteElement::ErrorMethod _pcs_cpl_error_method;
	std::string ls_sover_name;
	std::string ls_precond_name;

public:

	std::string pcs_type_name;

	//
	// LS - Linear Solver
	int ls_method;
	int ls_max_iterations;
	int ls_error_method;
	double ls_error_tolerance;
	double time_theta;
	int ls_precond;
	int ls_storage_method;
	std::string ls_extra_arg;
#ifdef USE_PETSC
	bool petsc_split_fields;
	bool petsc_use_snes;
#endif


	// NLS - Non-linear Solver
	std::string nls_method_name;
	FiniteElement::NonlinearSolverType nls_method;
	int nls_error_method;
	int nls_max_iterations;
	double nls_relaxation;
	double nls_error_tolerance[DOF_NUMBER_MAX];
	double nls_plasticity_local_tolerance;
	int nls_jacobian_level;

	// CPL
	std::string cpl_variable;
	std::string cpl_process;
	std::string cpl_variable_JOD;
	int cpl_max_iterations;
	int cpl_min_iterations;
	double cpl_error_tolerance[DOF_NUMBER_MAX];
	bool cpl_error_specified;
	bool cpl_master_process;

	// Integration
	int ele_gauss_points;

	// Mass lumping
	int ele_mass_lumping;

	// Upwind
	int ele_upwind_method;
	double ele_upwinding;

	// SUPG
	int ele_supg_method;
	int ele_supg_method_length;
	int ele_supg_method_diffusivity;

	// FEM-FCT
	int fct_method;
	unsigned int fct_prelimiter_type;
	double fct_const_alpha;
};

extern std::vector<CNumerics*> num_vector;
extern bool NUMRead(std::string);
extern void NUMWrite(std::string);
extern void NUMDelete();

extern int GetNumericsGaussPoints(int element_dimension);
#endif
