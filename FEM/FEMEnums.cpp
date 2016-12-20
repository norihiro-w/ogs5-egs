/*
 * FEMEnums.cpp
 *
 *  Created on: Sep 2, 2010
 *      Author: TF
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FEMEnums.h"
#include <cstdlib>
#include <iostream>

namespace FiniteElement
{
ProcessType convertProcessType(const std::string& pcs_type_string)
{
	if (pcs_type_string.compare("LIQUID_FLOW") == 0) return LIQUID_FLOW;
	if (pcs_type_string.compare("TWO_PHASE_FLOW") == 0) return TWO_PHASE_FLOW;
	if (pcs_type_string.compare("RICHARDS_FLOW") == 0) return RICHARDS_FLOW;
	if (pcs_type_string.compare("GROUNDWATER_FLOW") == 0)
		return GROUNDWATER_FLOW;
	if (pcs_type_string.compare("HEAT_TRANSPORT") == 0) return HEAT_TRANSPORT;
	if (pcs_type_string.compare("DEFORMATION") == 0) return DEFORMATION;
	if (pcs_type_string.compare("DEFORMATION_FLOW") == 0)
		return DEFORMATION_FLOW;
	if (pcs_type_string.compare("DEFORMATION_DYNAMIC") == 0)
		return DEFORMATION_DYNAMIC;
	if (pcs_type_string.compare("MASS_TRANSPORT") == 0) return MASS_TRANSPORT;
	if (pcs_type_string.compare("MULTI_PHASE_FLOW") == 0)
		return MULTI_PHASE_FLOW;
	if (pcs_type_string.compare("DEFORMATION_H2") == 0) return DEFORMATION_H2;
	if (pcs_type_string.compare("AIR_FLOW") == 0) return AIR_FLOW;
	if (pcs_type_string.compare("FLUID_MOMENTUM") == 0) return FLUID_MOMENTUM;
	if (pcs_type_string.compare("RANDOM_WALK") == 0) return RANDOM_WALK;
	if (pcs_type_string.compare("FLUX") == 0) return FLUX;
	if (pcs_type_string.compare("PS_GLOBAL") == 0) return PS_GLOBAL;
	if (pcs_type_string.compare("TH_MONOLITHIC") == 0) return TH_MONOLITHIC;
	if (pcs_type_string.compare("NO_PCS") == 0) return NO_PCS;
	// else
	// std::cout << "WARNING in convertProcessType: process type #" <<
	// pcs_type_string <<
	//"# unknown" << "\n";
	return INVALID_PROCESS;
}

std::string convertProcessTypeToString(ProcessType pcs_type)
{
	if (pcs_type == LIQUID_FLOW) return "LIQUID_FLOW";
	if (pcs_type == TWO_PHASE_FLOW) return "TWO_PHASE_FLOW";
	if (pcs_type == RICHARDS_FLOW) return "RICHARDS_FLOW";
	if (pcs_type == GROUNDWATER_FLOW) return "GROUNDWATER_FLOW";
	if (pcs_type == HEAT_TRANSPORT) return "HEAT_TRANSPORT";
	if (pcs_type == DEFORMATION) return "DEFORMATION";
	if (pcs_type == DEFORMATION_FLOW) return "DEFORMATION_FLOW";
	if (pcs_type == DEFORMATION_DYNAMIC) return "DEFORMATION_DYNAMIC";
	if (pcs_type == MASS_TRANSPORT) return "MASS_TRANSPORT";
	if (pcs_type == MULTI_PHASE_FLOW) return "MULTI_PHASE_FLOW";
	if (pcs_type == DEFORMATION_H2) return "DEFORMATION_H2";
	if (pcs_type == AIR_FLOW) return "AIR_FLOW";
	if (pcs_type == FLUID_MOMENTUM) return "FLUID_MOMENTUM";
	if (pcs_type == RANDOM_WALK) return "RANDOM_WALK";
	if (pcs_type == FLUX) return "FLUX";
	if (pcs_type == PS_GLOBAL) return "PS_GLOBAL";
	if (pcs_type == TH_MONOLITHIC) return "TH_MONOLITHIC";
	if (pcs_type == NO_PCS) return "NO_PCS";
	return "INVALID_PROCESS";
}

bool isFlowProcess(ProcessType pcs_type)
{
	if (pcs_type == LIQUID_FLOW ||
	    pcs_type == RICHARDS_FLOW || pcs_type == GROUNDWATER_FLOW ||
	    pcs_type == PS_GLOBAL || pcs_type == MULTI_PHASE_FLOW ||
	    pcs_type == DEFORMATION_FLOW || pcs_type == DEFORMATION_H2 ||
		pcs_type == TWO_PHASE_FLOW ||
	    pcs_type == AIR_FLOW || pcs_type == TH_MONOLITHIC)
		return true;
	return false;
}

bool isMultiFlowProcess(ProcessType pcs_type)
{
	if (pcs_type == PS_GLOBAL || pcs_type == MULTI_PHASE_FLOW ||
	    pcs_type == TWO_PHASE_FLOW || pcs_type == DEFORMATION_H2)
		return true;
	return false;
}

bool isDeformationProcess(ProcessType pcs_type)
{
	if (pcs_type == DEFORMATION || pcs_type == DEFORMATION_H2 ||
	    pcs_type == DEFORMATION_FLOW || pcs_type == DEFORMATION_DYNAMIC)
		return true;
	return false;
}

const std::list<std::string> getAllProcessNames()
{
	size_t count(1);
	std::list<std::string> enum_names;

	while (count != PROCESS_END)
	{
		enum_names.push_back(
		    convertProcessTypeToString(static_cast<ProcessType>(count++)));
	}
	return enum_names;
}

PrimaryVariable convertPrimaryVariable(const std::string& pcs_pv_string)
{
	if (pcs_pv_string.compare("PRESSURE1") == 0) return PRESSURE;
	if (pcs_pv_string.compare("PRESSURE2") == 0) return PRESSURE2;
	if (pcs_pv_string.compare("PRESSURE_RATE1") == 0) return PRESSURE_RATE1;
	if (pcs_pv_string.compare("SATURATION1") == 0) return SATURATION;
	if (pcs_pv_string.compare("SATURATION2") == 0) return SATURATION2;
	if (pcs_pv_string.compare("TEMPERATURE1") == 0) return TEMPERATURE;
	if (pcs_pv_string.compare("DISPLACEMENT_X1") == 0) return DISPLACEMENT_X;
	if (pcs_pv_string.compare("DISPLACEMENT_Y1") == 0) return DISPLACEMENT_Y;
	if (pcs_pv_string.compare("DISPLACEMENT_Z1") == 0) return DISPLACEMENT_Z;
	if (pcs_pv_string.compare("CONCENTRATION1") == 0) return CONCENTRATION;
	if (pcs_pv_string.compare("HEAD") == 0) return HEAD;
	if (pcs_pv_string.compare("VELOCITY_DM_X") == 0) return VELOCITY_DM_X;
	if (pcs_pv_string.compare("VELOCITY_DM_Y") == 0) return VELOCITY_DM_Y;
	if (pcs_pv_string.compare("VELOCITY_DM_Z") == 0) return VELOCITY_DM_Z;
	if (pcs_pv_string.compare("VELOCITY1_X") == 0) return VELOCITY1_X;
	if (pcs_pv_string.compare("VELOCITY1_Y") == 0) return VELOCITY1_Y;
	if (pcs_pv_string.compare("VELOCITY1_Z") == 0) return VELOCITY1_Z;
	if (pcs_pv_string.compare("STRESS_XX") == 0) return STRESS_XX;
	if (pcs_pv_string.compare("STRESS_XY") == 0) return STRESS_XY;
	if (pcs_pv_string.compare("STRESS_XZ") == 0) return STRESS_XZ;
	if (pcs_pv_string.compare("STRESS_YY") == 0) return STRESS_YY;
	if (pcs_pv_string.compare("STRESS_YZ") == 0) return STRESS_YZ;
	if (pcs_pv_string.compare("STRESS_ZZ") == 0) return STRESS_ZZ;
	if (pcs_pv_string.compare("ACCELERATION_X1") == 0) return ACCELERATION_X1;
	if (pcs_pv_string.compare("ACCELERATION_Y1") == 0) return ACCELERATION_Y1;
	if (pcs_pv_string.compare("ACCELERATION_Z1") == 0) return ACCELERATION_Z1;
	if (pcs_pv_string.compare("EXCAVATION") == 0) return EXCAVATION;
	if (pcs_pv_string.compare("STRAIN_XX") == 0) return STRAIN_XX;
	if (pcs_pv_string.compare("STRAIN_XY") == 0) return STRAIN_XY;
	if (pcs_pv_string.compare("STRAIN_XZ") == 0) return STRAIN_XZ;
	if (pcs_pv_string.compare("STRAIN_YY") == 0) return STRAIN_YY;
	if (pcs_pv_string.compare("STRAIN_YZ") == 0) return STRAIN_YZ;
	if (pcs_pv_string.compare("STRAIN_ZZ") == 0) return STRAIN_ZZ;
	if (pcs_pv_string.compare("STRAIN_PLS") == 0) return STRAIN_PLS;
	// else
	//{
	// std::cout << "convertPrimaryVariable #" << pcs_pv_string << "# not found"
	// << "\n";
	// exit (1);
	//}
	return INVALID_PV;
}

std::string convertPrimaryVariableToString(PrimaryVariable pcs_pv)
{
	if (pcs_pv == PRESSURE) return "PRESSURE1";
	if (pcs_pv == PRESSURE2) return "PRESSURE2";
	if (pcs_pv == PRESSURE_RATE1) return "PRESSURE_RATE1";
	if (pcs_pv == SATURATION) return "SATURATION1";
	if (pcs_pv == SATURATION2) return "SATURATION2";
	if (pcs_pv == TEMPERATURE) return "TEMPERATURE1";
	if (pcs_pv == DISPLACEMENT_X) return "DISPLACEMENT_X1";
	if (pcs_pv == DISPLACEMENT_Y) return "DISPLACEMENT_Y1";
	if (pcs_pv == DISPLACEMENT_Z) return "DISPLACEMENT_Z1";
	if (pcs_pv == CONCENTRATION) return "CONCENTRATION1";
	if (pcs_pv == HEAD) return "HEAD";
	if (pcs_pv == VELOCITY_DM_X) return "VELOCITY_DM_X";
	if (pcs_pv == VELOCITY_DM_Y) return "VELOCITY_DM_Y";
	if (pcs_pv == VELOCITY_DM_Z) return "VELOCITY_DM_Z";
	if (pcs_pv == VELOCITY1_X) return "VELOCITY1_X";
	if (pcs_pv == VELOCITY1_Y) return "VELOCITY1_Y";
	if (pcs_pv == VELOCITY1_Z) return "VELOCITY1_Z";
	if (pcs_pv == STRESS_XX) return "STRESS_XX";
	if (pcs_pv == STRESS_XY) return "STRESS_XY";
	if (pcs_pv == STRESS_XZ) return "STRESS_XZ";
	if (pcs_pv == STRESS_YY) return "STRESS_YY";
	if (pcs_pv == STRESS_YZ) return "STRESS_YZ";
	if (pcs_pv == STRESS_ZZ) return "STRESS_ZZ";
	if (pcs_pv == STRAIN_XX) return "STRAIN_XX";
	if (pcs_pv == STRAIN_XY) return "STRAIN_XY";
	if (pcs_pv == STRAIN_XZ) return "STRAIN_XZ";
	if (pcs_pv == STRAIN_YY) return "STRAIN_YY";
	if (pcs_pv == STRAIN_YZ) return "STRAIN_YZ";
	if (pcs_pv == STRAIN_ZZ) return "STRAIN_ZZ";
	if (pcs_pv == STRAIN_PLS) return "STRAIN_PLS";
	if (pcs_pv == ACCELERATION_X1) return "ACCELERATION_X1";
	if (pcs_pv == ACCELERATION_Y1) return "ACCELERATION_Y1";
	if (pcs_pv == ACCELERATION_Z1) return "ACCELERATION_Z1";
	if (pcs_pv == EXCAVATION) return "EXCAVATION";
	return "INVALID_PRIMARY_VARIABLE";
}

const std::list<std::string> getAllPrimaryVariableNames()
{
	size_t count(1);
	std::list<std::string> enum_names;

	while (count != PV_END)
	{
		enum_names.push_back(convertPrimaryVariableToString(
		    static_cast<PrimaryVariable>(count++)));
	}
	return enum_names;
}

DistributionType convertDisType(const std::string& dis_type_string)
{
	if (dis_type_string.compare("CONSTANT") == 0) return CONSTANT;
	if (dis_type_string.compare("ANALYTICAL") == 0) return ANALYTICAL;
	if (dis_type_string.compare("AVERAGE") == 0) return AVERAGE;
	if (dis_type_string.compare("CONSTANT_GEO") == 0) return CONSTANT_GEO;
	if (dis_type_string.compare("GRADIENT") == 0) return GRADIENT;
	if (dis_type_string.compare("RESTART") == 0) return RESTART;
	if (dis_type_string.compare("LINEAR") == 0) return LINEAR;
	if (dis_type_string.compare("POINT") == 0) return POINT;
	if (dis_type_string.compare("CONSTANT_NEUMANN") == 0)
		return CONSTANT_NEUMANN;
	if (dis_type_string.compare("LINEAR_NEUMANN") == 0) return LINEAR_NEUMANN;
	if (dis_type_string.compare("NORMALDEPTH") == 0) return NORMALDEPTH;
	if (dis_type_string.compare("CRITICALDEPTH") == 0) return CRITICALDEPTH;
	if (dis_type_string.compare("GREEN_AMPT") == 0) return GREEN_AMPT;
	if (dis_type_string.compare("SYSTEM_DEPENDENT") == 0)
		return SYSTEM_DEPENDENT;
	if (dis_type_string.compare("PRECIPITATION") == 0) return PRECIPITATION;
	if (dis_type_string.compare("DIRECT") == 0) return DIRECT;
	if (dis_type_string.compare("FUNCTION") == 0)
		return FUNCTION;  // 24.08.2011. WW
	if (dis_type_string.compare("ELEMENT") == 0) return ELEMENT;
	if (dis_type_string.compare("INITIAL") == 0)
		return INITIAL;
	else
	{
		std::cout << "convertDisType #" << dis_type_string << "# not found"
		          << "\n";
		exit(1);
	}
	return INVALID_DIS_TYPE;
}

std::string convertDisTypeToString(DistributionType dis_type)
{
	if (dis_type == ANALYTICAL) return "ANALYTICAL";
	if (dis_type == AVERAGE) return "AVERAGE";
	if (dis_type == CONSTANT) return "CONSTANT";
	if (dis_type == CONSTANT_GEO) return "CONSTANT_GEO";
	if (dis_type == GRADIENT) return "GRADIENT";
	if (dis_type == RESTART) return "RESTART";
	if (dis_type == LINEAR) return "LINEAR";
	if (dis_type == POINT) return "POINT";
	if (dis_type == CONSTANT_NEUMANN) return "CONSTANT_NEUMANN";
	if (dis_type == LINEAR_NEUMANN) return "LINEAR_NEUMANN";
	if (dis_type == NORMALDEPTH) return "NORMALDEPTH";
	if (dis_type == CRITICALDEPTH) return "CRITICALDEPTH";
	if (dis_type == GREEN_AMPT) return "GREEN_AMPT";
	if (dis_type == SYSTEM_DEPENDENT) return "SYSTEM_DEPENDENT";
	if (dis_type == PRECIPITATION) return "PRECIPITATION";
	if (dis_type == DIRECT) return "DIRECT";
	if (dis_type == FUNCTION) return "FUNCTION";  // 24.08.2011. WW
	if (dis_type == INITIAL) return "INITIAL";
	if (dis_type == ELEMENT) return "ELEMENT";

	return "INVALID_DIS_TYPE";
}

const std::list<std::string> getAllDistributionNames()
{
	size_t count(1);
	std::list<std::string> enum_names;

	while (count != DIS_END)
	{
		enum_names.push_back(
		    convertDisTypeToString(static_cast<DistributionType>(count++)));
	}
	return enum_names;
}

ErrorMethod convertErrorMethod(const std::string& error_method_string)
{
	if (error_method_string.compare("LMAX") == 0) return LMAX;
	if (error_method_string.compare("ENORM") == 0) return ENORM;
	if (error_method_string.compare("EVNORM") == 0) return EVNORM;
	if (error_method_string.compare("ERNORM") == 0) return ERNORM;
	if (error_method_string.compare("BNORM") == 0)
		return BNORM;
	else
	{
		std::cout << "convertErrorMethod #" << error_method_string
		          << "# not found"
		          << "\n";
		exit(1);
	}
	return INVALID_ERROR_METHOD;
}

ComparisonOperatorType convertComparisonOperatorType(const std::string& str)
{
	if (str.compare("<") == 0) return LT;
	if (str.compare("<=") == 0) return LE;
	if (str.compare("=") == 0) return EQ;
	if (str.compare("!=") == 0) return NE;
	if (str.compare(">") == 0) return GT;
	if (str.compare(">=") == 0)
		return GE;
	else
	{
		std::cout << "convertComparisonOperatorType #" << str << "# not found"
		          << "\n";
		exit(1);
	}
	return INVALID_OPERATOR_TYPE;
}

SourceTermType convertSTType(const std::string& str)
{
	if (str.compare("SOURCE") == 0) return SOURCE;
	if (str.compare("NEUMANN") == 0) return NEUMANN;
	return INVALID_ST_TYPE;
}

std::string convertSTTypeToString(SourceTermType st_type)
{
	if (st_type == SOURCE) return "SOURCE";
	if (st_type == NEUMANN) return "NEUMANN";
	return "INVALID_ST_TYPE";
}

TimType convertTimType(const std::string& str)
{
	if (str.compare("STEADY") == 0) return TIM_STEADY;
	if (str.compare("TRANSIENT") == 0) return TIM_TRANSIENT;
	return INVALID_TIM_TYPE;
}

std::string convertTimTypeToString(TimType st_type)
{
	if (st_type == TIM_STEADY) return "STEADY";
	if (st_type == TIM_TRANSIENT) return "TRANSIENT";
	return "INVALID_TIM_TYPE";
}

}  // end namespace FiniteElement

std::string convertIterationTypeToString(IterationType::type itr_type)
{
	if (itr_type == IterationType::LINEAR)
		return "LINEAR";
	else if (itr_type == IterationType::NONLINEAR)
		return "NONLINEAR";
	else if (itr_type == IterationType::COUPLING)
		return "COUPLING";
	return "INVALID";
}

TimeControlType::type convertTimeControlType(const std::string& str)
{
	if (str == "PI_AUTO_STEP_SIZE")
		return TimeControlType::PI_AUTO_STEP_SIZE;
	else if (str == "DYNAMIC_VARIABLE")
		return TimeControlType::DYNAMIC_VARIABLE;
	else if (str.find("DYNAMIC_COURANT") != std::string::npos)
		return TimeControlType::DYNAMIC_COURANT;
	else if (str == "DYNAMIC_PRESSURE")
		return TimeControlType::DYNAMIC_PRESSURE;
	else if (str == "STEP_SIZE_RESTRICTION")
		return TimeControlType::STEP_SIZE_RESTRICTION;
	else if (str == "NEUMANN")
		return TimeControlType::NEUMANN;
	else if (str == "ERROR_CONTROL_ADAPTIVE")
		return TimeControlType::ERROR_CONTROL_ADAPTIVE;
	else if (str.find("SELF_ADAPTIVE") != std::string::npos)
		return TimeControlType::SELF_ADAPTIVE;
	else if (str.find("PID_CONTROL") != std::string::npos)
		return TimeControlType::PID_CONTROL;

	return TimeControlType::INVALID;
}

std::string convertTimeControlTypeToString(TimeControlType::type tc_type)
{
	if (tc_type == TimeControlType::PI_AUTO_STEP_SIZE)
		return "PI_AUTO_STEP_SIZE";
	else if (tc_type == TimeControlType::DYNAMIC_VARIABLE)
		return "DYNAMIC_VARIABLE";
	else if (tc_type == TimeControlType::DYNAMIC_COURANT)
		return "DYNAMIC_COURANT";
	else if (tc_type == TimeControlType::DYNAMIC_PRESSURE)
		return "DYNAMIC_PRESSURE";
	else if (tc_type == TimeControlType::STEP_SIZE_RESTRICTION)
		return "STEP_SIZE_RESTRICTION";
	else if (tc_type == TimeControlType::NEUMANN)
		return "NEUMANN";
	else if (tc_type == TimeControlType::ERROR_CONTROL_ADAPTIVE)
		return "ERROR_CONTROL_ADAPTIVE";
	else if (tc_type == TimeControlType::SELF_ADAPTIVE)
		return "SELF_ADAPTIVE";
	else if (tc_type == TimeControlType::PID_CONTROL)
		return "PID_CONTROL";
	return "INVALID";
}
