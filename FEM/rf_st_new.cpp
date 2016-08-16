/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
 FEMLib - Object: Source Terms ST
 Task:
 Programing:
 01/2004 OK Implementation
 last modified
 **************************************************************************/

#include "rf_st_new.h"

#include <algorithm>
#include <cfloat>
#include <iostream>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "makros.h"
#include "display.h"
#include "FileToolsRF.h"
#include "FileTools.h"
#include "MemWatch.h"
#include "quicksort.h"

#include "InterpolationAlgorithms/InverseDistanceInterpolation.h"
#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "PointWithID.h"
#include "GEOObjects.h"

#include "GeoIO.h"
#include "ProcessIO.h"
#include "readNonBlankLineFromInputStream.h"

#include "DistributionTools.h"
#include "fem_ele.h"
#include "mathlib.h"
#include "matrix_class.h"
#include "rfmat_cp.h"
#include "rf_mfp_new.h"
#include "rf_mmp_new.h"
#include "rf_node.h"
#include "rf_pcs.h"
#include "rf_tim_new.h"
#include "SourceTerm.h"
#include "tools.h"


using FiniteElement::CElement;
using MeshLib::CElem;
using MeshLib::CEdge;
using MeshLib::CNode;
using Math_Group::vec;

#ifndef GRAVITY_CONSTANT
#define GRAVITY_CONSTANT 9.81
#endif

std::vector<CSourceTerm*> st_vector;
std::list<CSourceTermGroup*> st_group_list;
std::vector<std::string> analytical_processes;
std::vector<std::string> analytical_processes_polylines;
std::vector<NODE_HISTORY*> node_history_vector;

/**************************************************************************
 FEMLib-Method:
 Task: ST constructor
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
CSourceTerm::CSourceTerm()
    : ProcessInfo(),
      GeoInfo(),
      _coupled(false),
      _sub_dom_idx(-1),
      dis_linear_f(NULL),
      GIS_shape_head(NULL)
// 07.06.2010, 03.2010. WW
{
	CurveIndex = -1;
	// KR critical_depth = false;
	//	COUPLING_SWITCH = false;
	geo_node_value = 0.0;
	nodes = NULL;        // OK
	analytical = false;  // CMCD
	//  display_mode = false; //OK
	this->TimeInterpolation = 0;  // BG
	is_transfer_bc = false;
	st_id = -1;
	has_constrain = false;
}

// KR: Conversion from GUI-ST-object to CSourceTerm
CSourceTerm::CSourceTerm(const SourceTerm* st)
    : ProcessInfo(st->getProcessType(), st->getProcessPrimaryVariable(), NULL),
      GeoInfo(st->getGeoType(), st->getGeoObj()),
      DistributionInfo(st->getProcessDistributionType())
{
	setProcess(PCSGet(this->getProcessType()));
	this->geo_name = st->getGeoName();
	const std::vector<size_t> dis_nodes = st->getDisNodes();
	const std::vector<double> dis_values = st->getDisValues();

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
	{
		this->geo_node_value = dis_values[0];
	}
	else if (this->getProcessDistributionType() ==
	         FiniteElement::LINEAR_NEUMANN)
	{
		for (size_t i = 0; i < dis_values.size(); i++)
		{
			this->PointsHaveDistribedBC.push_back(
			    static_cast<int>(dis_nodes[i]));
			this->DistribedBC.push_back(dis_values[i]);
		}
	}
	else if (this->getProcessDistributionType() == FiniteElement::DIRECT)
	{
		// variable "fname" needs to be set, this must be done from outside!
	}
	else
		std::cout << "Error in CBoundaryCondition() - DistributionType \""
		          << FiniteElement::convertDisTypeToString(
		                 this->getProcessDistributionType())
		          << "\" currently not supported."
		          << "\n";
}

/**************************************************************************
 FEMLib-Method:
 Task: BC deconstructor
 Programing:
 04/2004 OK Implementation
 **************************************************************************/
CSourceTerm::~CSourceTerm()
{
	DeleteHistoryNodeMemory();
	//    dis_file_name.clear();
	node_number_vector.clear();
	node_value_vector.clear();
	node_renumber_vector.clear();
	PointsHaveDistribedBC.clear();
	DistribedBC.clear();
	element_st_vector.clear();
	// WW----------22.02.2007-------------------
	// TF 06/2010
	size_t size(normal2surface.size());
	for (size_t i = 0; i < size; i++)
		delete normal2surface[i];
	size = pnt_parameter_vector.size();
	for (size_t i = 0; i < size; i++)
		delete pnt_parameter_vector[i];
	if (GIS_shape_head)  // 07.06.2010. WW
	{
		delete[] GIS_shape_head;
		GIS_shape_head = NULL;
	}
	// WW
	if (dis_linear_f) delete dis_linear_f;
	dis_linear_f = NULL;

	// WW---------------------------------------

	for (size_t i = 0; i < this->st_boundary_elements.size(); i++)
		delete this->st_boundary_elements[i];
}

const std::string& CSourceTerm::getGeoName() const
{
	return geo_name;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 11/2004 MB neues read Konzept
 02/2005 MB River condition
 03/2005 WW Node force released by excavation
 11/2005 CMCD Analytical source for matrix
 04/2006 OK CPL
 04/2006 OK MSH_TYPE
06/2010 TF modification of the signature, added geo_obj and unique_name
**************************************************************************/
std::ios::pos_type CSourceTerm::Read(std::ifstream* st_file,
                                     const GEOLIB::GEOObjects& geo_obj,
                                     const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string, sub_string;
	bool new_keyword = false;

	std::stringstream in;

	// JOD 4.10.01
	channel = 0, node_averaging = 0, no_surface_water_pressure = 0;
	std::ios::pos_type position;

	// read loop
	while (!new_keyword)
	{
		position = st_file->tellg();
		if (!GetLineFromFile(line, st_file)) break;
		line_string = line;
		if (line_string.find("#") != std::string::npos)
		{
			new_keyword = true;
			break;
		}
		remove_white_space(&line_string);  // OK

		/* search for keywords */
		// subkeyword found
		if (line_string.find("$PCS_TYPE") != std::string::npos)
		{
			FileIO::ProcessIO::readProcessInfo(*st_file, _pcs_type);
			//         in.str(GetLineFromFile1(st_file));
			//         std::string tmp;
			//         in >> tmp;
			//         setProcessType (convertProcessType (tmp));
			//         in.clear();
			continue;
		}

		// subkeyword found
		if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			//         in.str(GetLineFromFile1(st_file));
			std::string tmp;
			in >> tmp;

			if (this->getProcessType() == FiniteElement::MASS_TRANSPORT)
			{
				// HS set the pointer to MCP based on component name.
				// first do a check whether this name is existing and unique.
				if (cp_name_2_idx.count(tmp) == 1)
				{
					setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
					setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
				}
				else
				{
					DisplayErrorMsg(
					    "Error: In reading ST file, the input component names "
					    "are not found in MCP file!!!");
					exit(1);
				}
			}
			else
			{
				setProcess(PCSGet(this->getProcessType()));
				setProcessPrimaryVariable(
				    FiniteElement::convertPrimaryVariable(tmp));
			}
			in.clear();
			continue;
		}

		if (line_string.find("$COMP_NAME") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			//         in.str(GetLineFromFile1(st_file));
			std::string tmp;
			in >> tmp;
			// HS set the pointer to MCP based on component name.
			// first do a check whether this name is existing and unique.
			if (cp_name_2_idx.count(tmp) == 1)
			{
				setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
				setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
			}
			else
			{
				DisplayErrorMsg(
				    "Error: In reading ST file, the input component names are "
				    "not found in MCP file!!!");
				exit(1);
			}
			in.clear();
			continue;
		}

		if (line_string.find("$GEO_TYPE") != std::string::npos)
		{
			ReadGeoType(st_file, geo_obj, unique_name);
			continue;
		}

		// 05.09.2008 WW
		if (line_string.find("$DIS_TYPE") != std::string::npos)
		{
			// 10.04.2008. WW  if(line_string.compare("$DIS_TYPE")==0) {
			if (line_string.find("CONDITION") != std::string::npos)
			{
				_coupled = true;
				ReadDistributionType(st_file);
				in.str(readNonBlankLineFromInputStream(*st_file));
				in >> line_string >> pcs_type_name_cond;
				in.clear();
				in.str(readNonBlankLineFromInputStream(*st_file));  //
				in >> pcs_pv_name_cond;
				in.clear();
				//            in.str(GetLineFromFile1(st_file));
				in.str(readNonBlankLineFromInputStream(*st_file));
				in >> _coup_leakance >> rill_height;
				in.clear();
			}  // 05.09.2008 WW
			else
			{
				ReadDistributionType(st_file);
				continue;
			}
		}

		if (line_string.find("$NODE_AVERAGING") != std::string::npos)
		{
			in.clear();
			node_averaging = true;
			continue;
		}
		// JOD 4.10.01
		if (line_string.find("$NEGLECT_SURFACE_WATER_PRESSURE") !=
		    std::string::npos)
		{
			in.clear();
			no_surface_water_pressure = true;
			continue;
		}

		if (line_string.find("$CHANNEL") != std::string::npos)
		{
			in.clear();
			in.str(readNonBlankLineFromInputStream(*st_file));
			in >> channel_width;
			channel = 1;
			continue;
		}

		if (line_string.find("$TIM_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			in >> tim_type_name;
			if (tim_type_name.find("CURVE") != std::string::npos)
			{
				//				dis_type = 0;
				in >> CurveIndex;
			}
			in.clear();
			continue;
		}

		// defines if time dependent source terms are use as piecewise constant
		// or linear interpolated; BG 05/2011
		if (line_string.find("$TIME_INTERPOLATION") != std::string::npos)
		{
			in.str(GetLineFromFile1(st_file));
			in >> interpolation_method;
			if (interpolation_method.find("LINEAR") != std::string::npos)
			{
				this->TimeInterpolation = 0;
			}
			if (interpolation_method.find("PIECEWISE_CONSTANT") !=
			    std::string::npos)
			{
				this->TimeInterpolation = 1;
			}
			in.clear();
			continue;
		}

		if (line_string.find("$FCT_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			in >> fct_name;  // sub_line
			// WW
			if (fct_name.find("METHOD") != std::string::npos) in >> fct_method;
			in.clear();
		}

		if (line_string.find("$MSH_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			std::string sub_string;
			in >> sub_string;  // sub_line
			msh_type_name = "NODE";
			if (sub_string.find("NODE") != std::string::npos)
			{
				in >> msh_node_number;
				in.clear();
			}
			continue;
		}

		if (line_string.find("$TRANSFER_COEFFICIENT") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*st_file));
			std::string str;
			in >> str;
			if (str.find("CONSTANT") != std::string::npos)
			{
				transfer_h_values.resize(1);
				in >> transfer_h_values[0];
			}
			else if (str.find("SUBDOMAIN") != std::string::npos)
			{
				in.clear();
				in.str(readNonBlankLineFromInputStream(*st_file));
				unsigned n_doms = 0;
				in >> n_doms;
				transfer_h_values.resize(n_doms);
				transfer_h_matId.resize(n_doms);
				for (unsigned ii = 0; ii < n_doms; ii++)
				{
					in.clear();
					in.str(readNonBlankLineFromInputStream(*st_file));
					in >> transfer_h_matId[ii] >> transfer_h_values[ii];
				}
			}
			is_transfer_bc = true;
			std::cout << "-> ST is defined as a transfer condition"
			          << "\n";
			in.clear();
			continue;
		}

		//....................................................................
		if (line_string.find("$CONSTRAIN") != std::string::npos)  // NW
		{
			// var > value
			in.str(readNonBlankLineFromInputStream(*st_file));
			in >> constrain_var_name >> line_string >> constrain_value;
			constrain_operator =
			    FiniteElement::convertComparisonOperatorType(line_string);
			in.clear();
			has_constrain = true;
			std::cout << "-> $CONSTRAIN is given for ST"
			          << "\n";
		}
	}  // end !new_keyword
	return position;
}

/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 02/2009 WW  Add a functionality to directly assign source terms to element
 nodes.
 **************************************************************************/
void CSourceTerm::ReadDistributionType(std::ifstream* st_file)
{
	std::stringstream in;
	// 03.2010 WW
	std::string aline;
	std::stringstream ss;
	int abuff, nLBC = 0;
	double bbuff;

	std::string dis_type_name;
	in.str(GetLineFromFile1(st_file));
	in >> dis_type_name;

	std::size_t found = dis_type_name.find("_NEUMANN");
	if (found != std::string::npos)
	{
		this->st_type = FiniteElement::NEUMANN;
		dis_type_name.erase(dis_type_name.begin() + found, dis_type_name.end());
		ScreenMessage("-> This ST is recognized as Neumann BC\n");
	}
	else
	{
		this->st_type = FiniteElement::SOURCE;
	}

	this->setProcessDistributionType(
	    FiniteElement::convertDisType(dis_type_name));

	if (dis_type_name.compare(
	        convertDisTypeToString(this->getProcessDistributionType())) != 0)
	{
		std::cerr
		    << "Error in CSourceTerm::ReadDistributionType (): dist_type_name #"
		    << dis_type_name << "#, new: "
		    << convertDisTypeToString(this->getProcessDistributionType())
		    << "\n";
		exit(1);
	}

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT ||
	    this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN ||
	    this->getProcessDistributionType() == FiniteElement::CONSTANT_GEO)
	{
		in >> geo_node_value;
		in.clear();
	}

	//	if (dis_type_name.find("ANALYTICAL") != std::string::npos) {
	if (this->getProcessDistributionType() == FiniteElement::ANALYTICAL)
	{
		in >> analytical_material_group;  // Which material group is it being
		                                  // applied to
		in >> analytical_diffusion;       // D value
		in >> analytical_porosity;        // n value of matrix
		in >> analytical_tortousity;      // t value of matrix
		in >> analytical_linear_sorption_Kd;  // Linear sorption coefficient
		in >> analytical_matrix_density;      // Density of solid
		in >> number_of_terms;  // no timesteps to consider in solution
		in >> resolution;       // every nth term will be considered
		in >> factor;           // to convert temperature to energy
		analytical = true;
		analytical_processes.push_back(
		    convertPrimaryVariableToString(getProcessPrimaryVariable()));
		//		if (geo_type_name.compare("POLYLINE") == 0)
		if (this->getGeoType() == GEOLIB::POLYLINE)
			analytical_processes_polylines.push_back(geo_name);
		in.clear();
	}

	// If a linear function is given. 25.08.2011. WW
	if (getProcessDistributionType() == FiniteElement::FUNCTION)
	{
		in.clear();
		dis_linear_f = new LinearFunctionData(*st_file);
	}

	if (this->getProcessDistributionType() == FiniteElement::LINEAR ||
	    this->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
	{
		in >> nLBC;
		in.clear();
		for (int i = 0; i < nLBC; i++)
		{
			in.str(GetLineFromFile1(st_file));
			in >> abuff >> bbuff;
			in.clear();
			PointsHaveDistribedBC.push_back(abuff);
			DistribedBC.push_back(bbuff);
		}

		//      Read LINENODES AND VALUES......
		in.clear();
	}

	if (this->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
	{
		// KR critical_depth = true;
		in >> geo_node_value;
		in.clear();
		in.str(GetLineFromFile1(st_file));
		in >> rill_height;
		in.clear();
		//		dis_type = 6;
	}

	if (this->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
	{
		dis_type_name = "NORMALDEPTH";
		in >> geo_node_value;
		in.clear();
		in.str(GetLineFromFile1(st_file));
		in >> normaldepth_slope >> rill_height;
		in.clear();
	}

	if (this->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
	{
		dis_type_name = "GREEN_AMPT";
		in >> geo_node_value;
		in.clear();
		in.str(GetLineFromFile1(st_file));
		in >> sorptivity >> constant >> rainfall >> moistureDeficit;
		in.clear();
	}
	// Soure terms are assign to element nodes directly. 23.02.2009. WW
	if (dis_type_name.find("DIRECT") != std::string::npos)
	{
		dis_type_name = "DIRECT";
		in >> fname;
		fname = FilePath + fname;
		in.clear();
	}

	// Soure terms from precipitation are assign to element nodes
	// directly.03.2010. WW
	if (dis_type_name.find("PRECIPITATION") != std::string::npos)
	{
		dis_type_name = "PRECIPITATION";
		in >> fname;
		fname = FilePath + fname;
		std::ifstream ins(fname.c_str());
		if (!ins.good())
		{
			std::cout << "Could not find file " << fname << "\n";
			exit(0);
		}
		double timess;
		GIS_shape_head = new double[6];  // 07.06.2010. WW
		for (int i = 0; i < 6; i++)
		{
			getline(ins, aline);
			ss.str(aline);
			ss >> aline >> GIS_shape_head[i];
			ss.clear();
		}
		while (!ins.eof())
		{
			getline(ins, aline);
			if (aline.find("#STOP") != std::string::npos) break;
			ss.str(aline);
			ss >> timess >> aline;
			precip_times.push_back(timess);
			precip_files.push_back(aline);

			ss.clear();
		}
		in.clear();
	}

	if (this->getProcessDistributionType() == FiniteElement::GRADIENT)
	{
		dis_type_name = "GRADIENT";
		in >> gradient_ref_depth;
		in >> gradient_ref_depth_value;
		in >> gradient_ref_depth_gradient;
		in.clear();
	}
}

/**************************************************************************
 FEMLib-Method:
 Task: for CSourceTerm::Read
 Programing:
 11/2007 JOD Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
void CSourceTerm::ReadGeoType(std::ifstream* st_file,
                              const GEOLIB::GEOObjects& geo_obj,
                              const std::string& unique_name)
{
	FileIO::GeoIO::readGeoInfo(this, *st_file, geo_name, geo_obj, unique_name);

	if (getProcessPrimaryVariable() == FiniteElement::EXCAVATION)  // WW
	{
		std::stringstream strstr;
		strstr.str(GetLineFromFile1(st_file));
		// size_t tmp_geo_type;
		std::string sub_string;
		strstr >> sub_string >> _sub_dom_idx;
		strstr.clear();
	}
}

/**************************************************************************
 FEMLib-Method:
 Task: ST read function
 Programing:
 01/2004 OK Implementation
 06/2010 TF modification of the signature, added geo_obj and unique_name
 **************************************************************************/
bool STRead(const std::string& file_base_name,
            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string, st_file_name;
	std::ios::pos_type position;

	// File handling
	st_file_name = file_base_name + ST_FILE_EXTENSION;
	std::ifstream st_file(st_file_name.data(), std::ios::in);

	ScreenMessage("STRead ... \n");
	if (!st_file.good())
	{
		ScreenMessage("-> No source term found !\n");
		return false;
	}

	// Keyword loop
	while (!st_file.eof())
	{
		st_file.getline(line, MAX_ZEILE);
		line_string = line;
		// Code included to make dynamic memory for analytical solution
		if (line_string.find("#STOP") != std::string::npos)
		{
			size_t no_source_terms(st_vector.size());
			size_t no_an_sol = 0, number_of_terms = 0;
			// Need to find the memory size limits for the anal. solution.
			for (size_t i = 0; i < no_source_terms; i++)
			{
				if (st_vector[i]->isAnalytical())
				{
					no_an_sol++;
					number_of_terms = std::max(st_vector[i]->getNumberOfTerms(),
					                           number_of_terms);
				}
			}
			if (no_an_sol > 0)
			{
				for (size_t i = 0; i < no_source_terms; i++)
				{
					st_vector[i]->setNumberOfAnalyticalSolutions(no_an_sol);
					st_vector[i]->setMaxNumberOfTerms(number_of_terms);
				}
			}
			break;
		}
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#SOURCE_TERM") != std::string::npos)
		{
			CSourceTerm* st(new CSourceTerm());
			std::ios::pos_type pos(st_file.tellg());
			position = st->Read(&st_file, geo_obj, unique_name);
			if (pos != position)
			{
				st_vector.push_back(st);
			}
			else
			{
				std::cerr << "WARNING: in STRead: could not read source term"
				          << "\n";
				delete st;
			}
			st_file.seekg(position, std::ios::beg);
		}  // keyword found
	}      // eof

	ScreenMessage("-> done, read %d source terms\n", st_vector.size());

	return true;
}


/**************************************************************************
 FEMLib-Method:
 Task: ST to mesh nodes
 Programing:
 01/2004 OK Implementation
 **************************************************************************/
// void ST2NOD()
//{
//	CGLVolume *m_volume = NULL;
//	double st_value;
//
//	CGLPolyline *ply (NULL);
//	size_t points_along_polyline;
//
//	// Nodes
//	size_t no_st (st_vector.size());
//	for (size_t j = 0; j < no_st; j++) {
//		CSourceTerm *m_st = st_vector[j];
//		switch (m_st->getGeoType()) {
//		case GS_POINT:
//			m_st->node_number_vector.push_back(m_st->getGeoObjIdx());
//			break;
//		case GS_POLYLINE:
//			CGLPolyline *ply = GEOGetPLYByName(m_st->geo_prop_name);
//			points_along_polyline = ply->point_vector.size();
//			for (size_t i = 0; i < points_along_polyline; i++) {
//				m_st->node_number_vector.push_back(ply->point_vector[i]->id);
//			}
//			break;
//		case GS_SURFACE:
//			Surface *m_surface (GEOGetSFCByName(m_st->geo_prop_name));
//			long points_in_surface;
//			long* nodes (NULL);
//			nodes = GetPointsIn(m_surface, &points_in_surface);
//			//MB patch areas
//			for (size_t i = 0; i < static_cast<size_t>(points_in_surface); i++)
//{
//				m_st->node_number_vector.push_back(nodes[i]);
//			}
//			delete [] nodes;
//			break;
//		case GS_VOLUME:
//			m_volume = GEOGetVOL(m_st->geo_prop_name);
//			break;
//		default:
//			break;
//		} // switch
//	} // while
//	//========================================================================
//	size_t st_point_number;
//	size_t geo_point_number;
//	// Values
//	for (size_t j = 0; j < no_st; j++) {
//		CSourceTerm *m_st = st_vector[j];
//		switch (m_st->dis_type) {
//		case CONSTANT:
//			st_point_number = m_st->node_number_vector.size();
//			for (size_t i = 0; i < st_point_number; i++) {
//				m_st->node_value_vector.push_back(m_st->dis_prop[0]);
//			}
//			break;
//		case LINEAR: // for polylines
//			CGLPolyline *ply = GEOGetPLYByName(m_st->geo_prop_name);//CC
//
//			geo_point_number = ply->point_vector.size();
//			//if(!(st_point_number==geo_point_number)) Warning !
//			for (size_t i = 0; i < geo_point_number; i++) {
//				st_value = ply->point_vector[i]->propert; // i.e. node property
//				m_st->node_value_vector.push_back(st_value);
//			}
//			break;
//		} // switch
//	} // while
//}

/**************************************************************************
 FEMLib-Method: STWrite
 Task: master write function
 Programing:
 04/2004 OK Implementation
 last modification:
 05/2010 TF
 **************************************************************************/
void STWrite(std::string base_file_name)
{
	// File handling
	std::string st_file_name = base_file_name + ST_FILE_EXTENSION;
	std::fstream st_file(st_file_name.data(), std::ios::trunc | std::ios::out);
	st_file.setf(std::ios::scientific, std::ios::floatfield);
	st_file.precision(12);
	if (!st_file.good()) return;

	st_file << "GeoSys-ST: Source Terms "
	           "------------------------------------------------\n";

	size_t no_st(st_vector.size());
	for (size_t j = 0; j < no_st; j++)
	{
		st_vector[j]->Write(&st_file);
	}
	st_file << "#STOP";
	st_file.close();
}

/**************************************************************************
 FEMLib-Method:
 Task: write function
 Programing:
 02/2004 OK Implementation
 04/2005 OK PRIMARY_VARIABLE
 06/2005 OK RIVER
 last modification:
 **************************************************************************/
void CSourceTerm::Write(std::fstream* st_file)
{
	// KEYWORD
	*st_file << "#SOURCE_TERM"
	         << "\n";
	//--------------------------------------------------------------------
	// NAME+NUMBER
	*st_file << " $PCS_TYPE"
	         << "\n";
	*st_file << "  ";
	//	*st_file << pcs_type_name << "\n";
	*st_file << convertProcessTypeToString(getProcessType()) << "\n";
	*st_file << " $PRIMARY_VARIABLE"
	         << "\n";
	*st_file << "  ";
	*st_file << convertPrimaryVariableToString(getProcessPrimaryVariable())
	         << "\n";
	//--------------------------------------------------------------------
	// GEO_TYPE
	if (this->getProcessDistributionType() != FiniteElement::DIRECT)
	{
		*st_file << " $GEO_TYPE"
		         << "\n";
		*st_file << "  ";
		*st_file << getGeoTypeAsString() << " " << geo_name << "\n";
	}
	//--------------------------------------------------------------------
	// TIM_TYPE
	if (tim_type_name.size() > 0)  // OK
	{
		*st_file << " $TIM_TYPE"
		         << "\n";
		*st_file << "  ";
		*st_file << tim_type_name << "\n";
	}
	//--------------------------------------------------------------------
	// DIS_TYPE
	*st_file << " $DIS_TYPE"
	         << "\n";
	*st_file << "  ";
	*st_file << convertDisTypeToString(this->getProcessDistributionType());
	switch (this->getProcessDistributionType())
	{
		case FiniteElement::CONSTANT:
			*st_file << " " << geo_node_value;
			*st_file << "\n";
			break;
		case FiniteElement::CONSTANT_NEUMANN:
			*st_file << " " << geo_node_value;
			*st_file << "\n";
			break;
		case FiniteElement::LINEAR:
			*st_file << " " << (int)PointsHaveDistribedBC.size() << "\n";
			for (long i = 0; i < (long)PointsHaveDistribedBC.size(); i++)
			{
				*st_file << "  " << PointsHaveDistribedBC[i] << " ";
				*st_file << "  " << DistribedBC[i] << "\n";
			}
			break;
		case FiniteElement::LINEAR_NEUMANN:
			*st_file << " " << PointsHaveDistribedBC.size() << "\n";
			for (size_t i = 0; i < PointsHaveDistribedBC.size(); i++)
			{
				*st_file << "  " << PointsHaveDistribedBC[i] << " ";
				*st_file << "  " << DistribedBC[i] << "\n";
			}
			break;
		case FiniteElement::DIRECT:
			*st_file << " " << this->fname << "\n";
			break;
		default:
			std::cerr << "this distributition type is not handled in "
			             "CSourceTerm::Write"
			          << "\n";
			break;
	}
	//--------------------------------------------------------------------
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 04/2004 OK Implementation
 11/2007 JOD Reaktivation
 last modification:
 **************************************************************************/
// void CSourceTerm::SetDISType()
//{
////	if (this->getProcessDistributionType() == CONSTANT)
////		dis_type = 1;
////////	if (dis_type_name.compare("CONSTANT_GEO") == 0)
////////		dis_type = 12; //SB flux is distributed along polyline. To do.
/// 4.10.06
////	if (this->getProcessDistributionType() == LINEAR)
////		dis_type = 2;
////	if (this->getProcessDistributionType() == CONSTANT_NEUMANN)
////		dis_type = 3;
////	if (this->getProcessDistributionType() == LINEAR_NEUMANN)
////		dis_type = 4;
////	if (this->getProcessDistributionType() == RIVER) {
////		dis_type = 5;
////	}
//
////	if (this->getProcessDistributionType() == CRITICALDEPTH)
////		dis_type = 6;
////	if (this->getProcessDistributionType() == SYSTEM_DEPENDENT)
////		dis_type = 7; //YD
////	if (this->getProcessDistributionType() == NORMALDEPTH)
////		dis_type = 8; //JOD MB
////	if (this->getProcessDistributionType() == ANALYTICAL)
////		dis_type = 9;//CMCD 02 2006
//////	if (dis_type_name.compare("PHILIP") == 0)
////////	if (this->getProcessDistributionType() == PHILIP)
//////		dis_type = 10; // JOD
////	if (this->getProcessDistributionType() == GREEN_AMPT)
////		dis_type = 11; // JOD
//	// if(dis_type_name.compare("CONSTANT")==0) dis_type = 0;
//	// if(dis_type_name.compare("LINEAR")  ==0) dis_type = 1;
//}

void setDistributionData(CSourceTerm* st, DistributionData& distData)
{
	distData.dis_type = st->getProcessDistributionType();
	distData.geo_obj = st->getGeoObj();
	distData.geo_name = st->getGeoName();
	distData.geo_type = st->getGeoType();
	switch (st->getProcessDistributionType())
	{
		case FiniteElement::CONSTANT:
		case FiniteElement::CONSTANT_NEUMANN:
			distData.dis_parameters.push_back(st->getGeoNodeValue());
			break;
		case FiniteElement::LINEAR:
		case FiniteElement::LINEAR_NEUMANN:
			distData._DistribedBC = st->getDistribedST();
			distData._PointsHaveDistribedBC = st->getPointsWithDistribedST();
			break;
		case FiniteElement::GRADIENT:
			distData.dis_parameters.push_back(st->gradient_ref_depth);
			distData.dis_parameters.push_back(st->gradient_ref_depth_value);
			distData.dis_parameters.push_back(st->gradient_ref_depth_gradient);
			break;
		case FiniteElement::FUNCTION:
			distData.linear_f = st->dis_linear_f;
			break;
		default:
			break;
	}
	distData.mesh_type_name = st->getMeshTypeName();
	if (st->getMeshTypeName() == "NODE")
	{
		distData.mesh_node_id = st->getMeshNodeNumber();
	}
}

/**************************************************************************
 FEMLib-Method:
 Task: set ST group member
 Programing:
 02/2004 OK Implementation
 09/2004 WW Face integration of Neumann boundary condition for all element type
 09/2004 WW Interpolation for piece-wise linear distributed source term or BC
 03/2005 OK LINE sources
 02/2005 MB River condition, CriticalDepth
 08/2006 WW Re-implementing edge,face and domain integration versatile for all
element types
 04/2006 OK MSH types
02/2009 WW Direct assign node source terms
**************************************************************************/
void CSourceTermGroup::Set(CRFProcess* m_pcs, const int ShiftInNodeVector,
                           std::string this_pv_name)
{
	if (!this_pv_name.empty()) pcs_pv_name = this_pv_name;
	const FiniteElement::ProcessType pcs_type =
	    FiniteElement::convertProcessType(pcs_type_name);
	const FiniteElement::PrimaryVariable pcs_pv_type =
	    FiniteElement::convertPrimaryVariable(pcs_pv_name);

	m_msh = FEMGet(pcs_type_name);
	if (!m_msh)
	{
		std::cout << "Warning in CSourceTermGroup::Set - no MSH data"
		          << "\n";
		return;
	}

	/// In case of P_U coupling monolithic scheme
	bool isQuadratic = false;
	if (m_pcs->type == 41)  // WW Mono
	{
		if (pcs_pv_name.find("DISPLACEMENT") != std::string::npos)
			isQuadratic = true;
	}
	else if (m_pcs->type == 4)
	{
		isQuadratic = true;
	}
	m_pcs->m_msh->SwitchOnQuadraticNodes(isQuadratic);

	//====================================================================
	const size_t no_st = st_vector.size();
	m_pcs->st_node_value.resize(no_st);
	for (size_t i = 0; i < no_st; i++)
	{
		CSourceTerm* st(st_vector[i]);

		if (st->getProcessType() != pcs_type) continue;
		if (st->getProcessPrimaryVariable() != pcs_pv_type &&
		    st->getProcessPrimaryVariable() != FiniteElement::CONCENTRATION)
			continue;
		if (st->getProcessDistributionType() == FiniteElement::PRECIPITATION)
			continue;

		if (st->getProcessType() == FiniteElement::MASS_TRANSPORT)
			if (cp_vec[cp_name_2_idx[convertPrimaryVariableToString(
			               st->getProcessPrimaryVariable())]]->getProcess() !=
			    m_pcs)
				continue;

		ScreenMessage("* %s on %s\n", FiniteElement::convertPrimaryVariableToString(st->getProcessPrimaryVariable()).data(), st->getGeoName().data());

		//-- 23.02.3009. WW
		if (st->getProcessDistributionType() == FiniteElement::DIRECT)
		{   // NB For climate ST, the source terms (recharge in this case) will
			// also be assigned directly to surface nodes
			st->DirectAssign(ShiftInNodeVector);
			continue;
		}

		if (st->is_transfer_bc)
		{
			std::vector<double> tmp_K(st->transfer_h_values);
			st->transfer_h_values.resize(mmp_vector.size());
			if (st->transfer_h_matId.size() > 0)
			{
				for (unsigned j = 0; j < st->transfer_h_matId.size(); j++)
					st->transfer_h_values[st->transfer_h_matId[j]] = tmp_K[j];
			}
			else
			{
				for (unsigned j = 0; j < st->transfer_h_values.size(); j++)
					st->transfer_h_values[j] = tmp_K[0];
			}
		}

		//------------------------------------------------------------------
		std::vector<long> nodes_vector;
		std::vector<double> node_value;

		if (st->isCoupled()) m_msh_cond = FEMGet(st->pcs_type_name_cond);

		DistributionData distData;
		setDistributionData(st, distData);
		//------------------------------------------------------------------
		// Detect mesh nodes for this ST
		//------------------------------------------------------------------
		nodes_vector.clear();
		getNodesOnDistribution(distData, *m_msh, nodes_vector);
		if (nodes_vector.empty())
		{
#ifndef USE_PETSC
			ScreenMessage("-> ***ERROR* No nodes found on %s %s\n", st->getGeoName().data());
#endif
			continue;
		}
		//------------------------------------------------------------------
		// Calculate ST values
		//------------------------------------------------------------------
		node_value.resize(nodes_vector.size());
		setDistribution(distData, *m_msh, nodes_vector, node_value);
		//------------------------------------------------------------------
		// Calculate integrated values
		//------------------------------------------------------------------
		if (st->getSTType() == FiniteElement::NEUMANN)
		{
			if (st->getGeoType() == GEOLIB::POINT && st->is_transfer_bc)
			{
				node_value[0] *= st->transfer_h_values[0];
			}
			else if (m_msh->GetMaxElementDim() == 2 &&
			         st->getGeoType() == GEOLIB::POLYLINE)
			{
				st->EdgeIntegration(m_msh, nodes_vector, node_value);
			}
			else if (m_msh->GetMaxElementDim() == 3 &&
			         st->getGeoType() == GEOLIB::SURFACE)
			{
				st->FaceIntegration(m_msh, nodes_vector, node_value);
			}
			else if (st->getGeoType() == GEOLIB::GEODOMAIN ||
			         (m_msh->GetMaxElementDim() == 1 &&
			          st->getGeoType() == GEOLIB::POLYLINE) ||
			         (m_msh->GetMaxElementDim() == 2 &&
			          st->getGeoType() == GEOLIB::SURFACE))
			{
				st->DomainIntegration(m_msh, nodes_vector, node_value);
			}
		}

		//------------------------------------------------------------------
		// create ST node
		//------------------------------------------------------------------
		std::vector<long> nodes_cond;
		ScreenMessaged("-> Create ST nodes%s\n", "");
		st->SetNodeValues(nodes_vector, nodes_cond, node_value,
		                  ShiftInNodeVector);
#ifndef WIN32
		BaseLib::MemWatch mem_watch;
		ScreenMessaged("\tcurrent mem: %d MB\n",
		               mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif

		//------------------------------------------------------------------
		// FCT
		//------------------------------------------------------------------
		if (st->fct_name.size() > 0) fct_name = st->fct_name;

		//------------------------------------------------------------------
		// SYSTEM_DEPENDENT
		//------------------------------------------------------------------
		if (st->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
		{
			m_pcs->compute_domain_face_normal = true;  // WW
			if (st->getGeoType() == GEOLIB::POINT)
			{
				for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
				{
					CElem* elem = m_msh->ele_vector[i];
					if (!elem->GetMark()) continue;
					int nn = elem->GetNodesNumber(m_msh->getOrder());
					for (int j = 0; j < nn; j++)
					{
						CNode* cnode = elem->GetNode(j);  // WW
						if (cnode->GetIndex() == (size_t)st->geo_node_number)
							st->element_st_vector.push_back(i);
					}
				}
			}
			else if (st->getGeoType() == GEOLIB::POLYLINE)
			{
				long no_face = (long)m_msh->face_vector.size();
				for (long i = 0; i < no_face; i++)
				{
					int node_on_line = 0;
					int no_vertex = m_msh->face_vector[i]->GetVertexNumber();
					for (long jj = 0; jj < no_vertex; jj++)
					{
						for (size_t kk = 0; kk < nodes_vector.size(); kk++)
						{
							if (nodes_vector[kk] ==
							    (m_msh->face_vector[i]->GetNodeIndex(jj)))
								node_on_line++;
						}  // end nodes
					}      // end vertices
					if (node_on_line == 2)
						st->element_st_vector.push_back(
						    m_msh->face_vector[i]->GetOwner()->GetIndex());
				}  // end faces
			}
		}

	}  // end st loop
}

/**************************************************************************
 ROCKFLOW - Funktion: FaceIntegration
 Task: Translate distributed Neumann boundary condition /source term on edges
 found on a polyline to nodes value for all kinds of element
 Programming:
 07/2005 WW Re-Implementation
 12/2005 WW Axismmetry
 **************************************************************************/
// void CSourceTerm::EdgeIntegration(CFEMesh* msh, vector<long>&nodes_on_ply,
//		vector<double>&node_value_vector) {
//	long i, j, k, l;
//	long this_number_of_nodes;
//	int elemsCnode;
//	int nedges, ii;
//	vec<CNode*> e_nodes(3);
//	vec<CEdge*> e_edges(12);
//
//	double Jac = 0.0;
//	double Weight = 0.0;
//	double eta = 0.0;
//	double v1, v2, radius = 0.0;
//	double Shfct[3];
//	bool Const = false;
//
//	if (dis_type_name.find("CONSTANT") != std::string::npos)
//		Const = true;
//
//	//CFEMesh* msh = m_pcs->m_msh;
//	//CFEMesh* msh;  // JOD
//	//msh = FEMGet(pcs_type_name);
//	CElem* elem = NULL;
//	CEdge* edge = NULL;
//	CNode* node = NULL;
//
//	int nSize = (long) msh->nod_vector.size();
//	this_number_of_nodes = (long) nodes_on_ply.size();
//	vector<long> G2L(nSize);
//	vector<double> NVal(this_number_of_nodes);
//
//	// Unmakr edges.
//	for (i = 0; i < (long) msh->edge_vector.size(); i++)
//		msh->edge_vector[i]->SetMark(false);
//	for (i = 0; i < nSize; i++) {
//		msh->nod_vector[i]->SetMark(false);
//		G2L[i] = -1;
//	}
//
//	// Search edges on polyline
//	for (i = 0; i < this_number_of_nodes; i++) {
//		NVal[i] = 0.0;
//		k = nodes_on_ply[i];
//		G2L[k] = i;
//		node = msh->nod_vector[k];
//		elemsCnode = (int) node->connected_elements.size();
//		for (j = 0; j < elemsCnode; j++) {
//			l = msh->nod_vector[k]->connected_elements[j];
//			elem = msh->ele_vector[l];
//			nedges = elem->GetEdgesNumber();
//			elem->GetEdges(e_edges);
//			for (ii = 0; ii < nedges; ii++) {
//				edge = e_edges[ii];
//				if (edge->GetMark())
//					continue;
//				edge->GetNodes(e_nodes);
//				// Edge A
//				if (*node == *e_nodes[0])
//					e_nodes[0]->SetMark(true);
//				// Edge B
//				if (*node == *e_nodes[1])
//					e_nodes[1]->SetMark(true);
//				if (msh->getOrder()) // Quadratic
//				{
//					if (*node == *e_nodes[2])
//						e_nodes[2]->SetMark(true);
//				}
//				if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()) {
//					if (msh->getOrder()) {
//						if (e_nodes[2]->GetMark())
//							edge->SetMark(true);
//					} else
//						edge->SetMark(true);
//				}
//
//			}// e_edges
//		}
//	}
//
//	for (i = 0; i < (long) msh->edge_vector.size(); i++) {
//		edge = msh->edge_vector[i];
//		if (!edge->GetMark())
//			continue;
//		edge->GetNodes(e_nodes);
//		if (msh->getOrder()) // Quad
//		{
//			if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()
//					&& e_nodes[2]->GetMark()) {
//				Jac = 0.5 * edge->Length();
//				v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
//				v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
//				if (Const && (!msh->isAxisymmetry())) {
//					NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1 / 3.0;
//					NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1 / 3.0;
//					NVal[G2L[e_nodes[2]->GetIndex()]] += 4.0 * Jac * v1 / 3.0;
//
//				} else {
//					for (k = 0; k < 3; k++) // Three nodes
//					{
//						// Numerical integration
//						for (l = 0; l < 3; l++) // Gauss points
//						{
//							Weight = Jac * MXPGaussFkt(3, l);
//							eta = MXPGaussPkt(3, l);
//							ShapeFunctionLineHQ(Shfct, &eta);
//							//Axisymmetical problem
//							if (msh->isAxisymmetry()) {
//								radius = 0.0;
//								for (ii = 0; ii < 3; ii++)
//									radius += Shfct[ii] * e_nodes[ii]->X();
//								Weight *= radius; //2.0*pai*radius;
//							}
//							NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
//									+ eta * (v2 - v1)) * Shfct[k] * Weight;
//						}
//					}
//				}
//			}
//		} else // Linear
//		{
//			if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()) {
//				Jac = 0.5 * edge->Length();
//				v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
//				v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
//				if (!msh->isAxisymmetry()) {
//					if (Const) {
//						NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1;
//						NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1;
//					} else {
//						NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * (2.0 * v1
//								+ v2) / 3.0;
//						NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * (v1 + 2.0
//								* v2) / 3.0;
//					}
//				} else // Axisymmetry
//				{
//
//					for (k = 0; k < 2; k++) // Three nodes
//					{
//						// Numerical integration
//						for (l = 0; l < 3; l++) // Gauss points
//						{
//							Weight = Jac * MXPGaussFkt(3, l);
//							eta = MXPGaussPkt(3, l);
//							ShapeFunctionLine(Shfct, &eta);
//							//Axisymmetical problem
//							if (msh->isAxisymmetry()) {
//								radius = 0.0;
//								for (ii = 0; ii < 2; ii++)
//									radius += Shfct[ii] * e_nodes[ii]->X();
//								Weight *= radius; //2.0*pai*radius;
//							}
//							NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
//									+ eta * (v2 - v1)) * Shfct[k] * Weight;
//						}
//					}
//				}// End of is (!axi)
//			}
//		}
//	}
//	for (i = 0; i < this_number_of_nodes; i++)
//		node_value_vector[i] = NVal[i];
//	for (i = 0; i < (long) msh->edge_vector.size(); i++)
//		msh->edge_vector[i]->SetMark(true);
//	for (i = 0; i < nSize; i++)
//		msh->nod_vector[i]->SetMark(true);
//	NVal.clear();
//	G2L.clear();
//	e_nodes.resize(0);
//	e_edges.resize(0);
//}

void CSourceTerm::EdgeIntegration(CFEMesh* msh,
                                  const std::vector<long>& nodes_on_ply,
                                  std::vector<double>& node_value_vector,
                                  std::vector<bool>* active_elements)
{
	long i, j, k, l;
	long this_number_of_nodes;
	int elemsCnode;
	int nedges, ii;
	vec<CNode*> e_nodes(3);
	vec<CEdge*> e_edges(12);

#if 0
   double Jac = 0.0;
   double Weight = 0.0;
   double eta = 0.0;
   double v1, v2, radius = 0.0;
   double Shfct[3];

   bool Const = false;
   if (this->getProcessDistributionType() == FiniteElement::CONSTANT || this->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
      Const = true;
#endif

	int Axisymm = 1;                         // ani-axisymmetry
	if (msh->isAxisymmetry()) Axisymm = -1;  // Axisymmetry is true
	CElement* fem = new CElement(Axisymm * msh->GetCoordinateFlag());

	// CFEMesh* msh = m_pcs->m_msh;
	// CFEMesh* msh;  // JOD
	// msh = FEMGet(pcs_type_name);
	CElem* elem = NULL;
	CEdge* edge = NULL;
	CNode* node = NULL;

	int nSize = (long)msh->nod_vector.size();
	this_number_of_nodes = (long)nodes_on_ply.size();
	std::vector<long> G2L(nSize);
	std::vector<double> NVal(this_number_of_nodes);

	// Unmakr edges.
	for (i = 0; i < (long)msh->edge_vector.size(); i++)
		msh->edge_vector[i]->SetMark(false);
	for (i = 0; i < nSize; i++)
	{
		msh->nod_vector[i]->SetMark(false);
		G2L[i] = -1;
	}

	// Search edges on polyline
	for (i = 0; i < this_number_of_nodes; i++)
	{
		NVal[i] = 0.0;
		k = nodes_on_ply[i];
		G2L[k] = i;
		node = msh->nod_vector[k];
		elemsCnode = (int)node->getConnectedElementIDs().size();
		for (j = 0; j < elemsCnode; j++)
		{
			l = msh->nod_vector[k]->getConnectedElementIDs()[j];
			elem = msh->ele_vector[l];
			nedges = elem->GetEdgesNumber();
			elem->GetEdges(e_edges);
			for (ii = 0; ii < nedges; ii++)
			{
				edge = e_edges[ii];
				if (edge->GetMark()) continue;
				edge->SetPatchIndex(elem->GetPatchIndex());  // NW
				edge->GetNodes(e_nodes);
				// Edge A
				if (*node == *e_nodes[0]) e_nodes[0]->SetMark(true);
				// Edge B
				if (*node == *e_nodes[1]) e_nodes[1]->SetMark(true);
				if (msh->getOrder())  // Quadratic
				{
					if (*node == *e_nodes[2]) e_nodes[2]->SetMark(true);
				}
				if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark())
				{
					if (msh->getOrder())
					{
						if (e_nodes[2]->GetMark()) edge->SetMark(true);
					}
					else
						edge->SetMark(true);
				}

			}  // e_edges
		}
	}

	for (i = 0; i < (long)msh->edge_vector.size(); i++)
	{
		edge = msh->edge_vector[i];
		if (!edge->GetMark()) continue;
		if (active_elements && !(*active_elements)[i]) continue;
		edge->GetNodes(e_nodes);
		if (msh->getOrder())
		{
			if (!(e_nodes[0]->GetMark() && e_nodes[1]->GetMark() &&
			      e_nodes[2]->GetMark()))
				continue;
		}
		else
		{
			if (!(e_nodes[0]->GetMark() && e_nodes[1]->GetMark())) continue;
		}
		CElem* edge_ele = new CElem(1);
		edge_ele->setElementProperties(MshElemType::LINE, true);
		const size_t nen = msh->getOrder() ? 3 : 2;
		std::vector<CNode*> edge_nodes(nen);
		for (size_t j = 0; j < nen; j++)
			edge_nodes[j] = edge->GetNode(j);
		edge_ele->setPatchIndex(edge->GetPatchIndex());
		edge_ele->setNodes(edge_nodes);
		edge_ele->SetOrder(msh->getOrder());
		edge_ele->ComputeVolume();
		if (!active_elements) st_boundary_elements.push_back(edge_ele);

		fem->setOrder(msh->getOrder() + 1);
		fem->ConfigElement(edge_ele, true);
		double nodesFVal[3] = {};
		nodesFVal[0] = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
		nodesFVal[1] = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
		if (msh->getOrder())
			nodesFVal[2] = node_value_vector[G2L[e_nodes[2]->GetIndex()]];
		fem->FaceIntegration(nodesFVal);
		if (this->is_transfer_bc)
		{
			for (k = 0; k < (signed)nen; k++)
				nodesFVal[k] *=
				    this->transfer_h_values[edge_ele->GetPatchIndex()];
		}
		for (k = 0; k < (signed)nen; k++)
		{
			NVal[G2L[edge->GetNode(k)->GetIndex()]] += nodesFVal[k];
		}

#if 0
      if (msh->getOrder())                        // Quad
      {
         if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark()
            && e_nodes[2]->GetMark())
         {
            Jac = 0.5 * edge->getLength();
            v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
            v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
            if (Const && (!msh->isAxisymmetry()))
            {
               NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1 / 3.0;
               NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1 / 3.0;
               NVal[G2L[e_nodes[2]->GetIndex()]] += 4.0 * Jac * v1 / 3.0;

            }
            else
            {
               for (k = 0; k < 3; k++)            // Three nodes
               {
                  // Numerical integration
                  for (l = 0; l < 3; l++)         // Gauss points
                  {
                     Weight = Jac * MXPGaussFkt(3, l);
                     eta = MXPGaussPkt(3, l);
                     ShapeFunctionLineHQ(Shfct, &eta);
                     //Axisymmetical problem
                     if (msh->isAxisymmetry())
                     {
                        radius = 0.0;
                        for (ii = 0; ii < 3; ii++)
                           radius += Shfct[ii] * e_nodes[ii]->getData()[0];
                        Weight *= radius;         //2.0*pai*radius;
                     }
                     NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
                        + eta * (v2 - v1)) * Shfct[k] * Weight;
                  }
               }
            }
         }
      } else                                      // Linear
      {
         if (e_nodes[0]->GetMark() && e_nodes[1]->GetMark())
         {
            Jac = 0.5 * edge->getLength();
            v1 = node_value_vector[G2L[e_nodes[0]->GetIndex()]];
            v2 = node_value_vector[G2L[e_nodes[1]->GetIndex()]];
            if (!msh->isAxisymmetry())
            {
               if (Const)
               {
                  NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * v1;
                  NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * v1;
               }
               else
               {
                  NVal[G2L[e_nodes[0]->GetIndex()]] += Jac * (2.0 * v1
                     + v2) / 3.0;
                  NVal[G2L[e_nodes[1]->GetIndex()]] += Jac * (v1 + 2.0
                     * v2) / 3.0;
               }
            } else                                // Axisymmetry
            {

               for (k = 0; k < 2; k++)            // Three nodes
               {
                  // Numerical integration
                  for (l = 0; l < 3; l++)         // Gauss points
                  {
                     Weight = Jac * MXPGaussFkt(3, l);
                     eta = MXPGaussPkt(3, l);
                     ShapeFunctionLine(Shfct, &eta);
                     //Axisymmetical problem
                     if (msh->isAxisymmetry())
                     {
                        radius = 0.0;
                        for (ii = 0; ii < 2; ii++)
                           radius += Shfct[ii] * e_nodes[ii]->getData()[0];
                        Weight *= radius;         //2.0*pai*radius;
                     }
                     NVal[G2L[e_nodes[k]->GetIndex()]] += 0.5 * (v1 + v2
                        + eta * (v2 - v1)) * Shfct[k] * Weight;
                  }
               }
            }                                     // End of is (!axi)
         }
      }
#endif
	}
	for (i = 0; i < this_number_of_nodes; i++)
		node_value_vector[i] = NVal[i];
	for (i = 0; i < (long)msh->edge_vector.size(); i++)
		msh->edge_vector[i]->SetMark(true);
	for (i = 0; i < nSize; i++)
		msh->nod_vector[i]->SetMark(true);
	NVal.clear();
	G2L.clear();
	e_nodes.resize(0);
	e_edges.resize(0);

	delete fem;
}

/**************************************************************************
 ROCKFLOW - Funktion: FaceIntegration
 Task: Translate distributed Neumann boundary condition /source term on faces
 found on a surface into nodes value for all kinds of element
 Programming:
 08/2005 WW Re-Implementation
 11/2005 WW/OK Layer optimization
 01/2010 NW improvement of efficiency to search faces
 **************************************************************************/

void CSourceTerm::FaceIntegration(CFEMesh* msh,
                                  std::vector<long> const& nodes_on_sfc,
                                  std::vector<double>& node_value_vector,
                                  std::vector<bool>* active_elements)
{
	if (!msh)
	{
		std::cout << "Warning in CSourceTerm::FaceIntegration: no MSH data, "
		             "function doesn't function";
		return;
	}

	long i, j, k, l;

	const long n_sfc_nodes = (long)nodes_on_sfc.size();
	const long n_msh_nodes = (long)msh->nod_vector.size();


	//----------------------------------------------------------------------
	// search element faces on the surface
	//----------------------------------------------------------------------
	#pragma omp parallel for
	for (i = 0; i < n_msh_nodes; i++)
	{
		msh->nod_vector[i]->SetMark(false);
	}

	// init
	const long n_ele = (long)msh->ele_vector.size();
	#pragma omp parallel for
	for (i = 0; i < n_ele; i++)
	{
		msh->ele_vector[i]->selected = 0;  // TODO can use a new variable
	}

#define UNIQUE_OMP
#ifdef UNIQUE_OMP
// unique set of node id on the surface
#ifdef _OPENMP
	std::cout << "[CSourceTerm::FaceIntegration] creating unique nodes ... "
	          << std::flush;
	double begin_unique = omp_get_wtime();
#endif
	std::vector<long> unique_nodes_on_sfc(nodes_on_sfc);
	std::sort(unique_nodes_on_sfc.begin(), unique_nodes_on_sfc.end());
	unique_nodes_on_sfc.erase(
	    std::unique(unique_nodes_on_sfc.begin(), unique_nodes_on_sfc.end()),
	    unique_nodes_on_sfc.end());
#ifdef _OPENMP
	std::cout << "done, took " << omp_get_wtime() - begin_unique << " s"
	          << "\n";
#endif

	// filtering elements: elements should have nodes on the surface
	// Notice: node-elements relation has to be constructed beforehand
	#pragma omp parallel for default(none) private(k, j, l) shared(msh, nodes_on_sfc)
	for (i = 0; i < n_sfc_nodes; i++)
	{
		k = nodes_on_sfc[i];
		CNode* nod = msh->nod_vector[k];
		const long n_k_conn_ele = nod->getConnectedElementIDs().size();
		for (j = 0; j < n_k_conn_ele; j++)
		{
			l = nod->getConnectedElementIDs()[j];
			CElem* e = msh->ele_vector[l];
			#pragma omp atomic
			e->selected++;  // remember how many nodes of an element are on the
			                // surface
		}
	}

	long n_selected_ele = 0;
	#pragma omp parallel for reduction(+ : n_selected_ele)
	for (i = 0; i < n_ele; i++)
	{
		if (msh->ele_vector[i]->selected > 0) n_selected_ele++;
	}

	std::vector<long> vec_possible_elements;
	vec_possible_elements.reserve(n_selected_ele);
	for (i = 0; i < n_ele; i++)
	{
		if (msh->ele_vector[i]->selected > 0)
			vec_possible_elements.push_back(i);
	}
#else
	std::set<long> set_nodes_on_sfc;  // unique set of node id on the surface
	for (i = 0; i < (long)nodes_on_sfc.size(); i++)
	{
		set_nodes_on_sfc.insert(nodes_on_sfc[i]);
	}

	// filtering elements: elements should have nodes on the surface
	// Notice: node-elements relation has to be constructed beforehand
	std::vector<long> vec_possible_elements;
	for (i = 0; i < this_number_of_nodes; i++)
	{
		k = nodes_on_sfc[i];
		for (j = 0;
		     j < (long)msh->nod_vector[k]->getConnectedElementIDs().size();
		     j++)
		{
			l = msh->nod_vector[k]->getConnectedElementIDs()[j];
			if (msh->ele_vector[l]->selected == 0)
				vec_possible_elements.push_back(l);
			msh->ele_vector[l]->selected +=
			    1;  // remember how many nodes of an element are on the surface
		}
	}
#endif

	//----------------------------------------------------------------------
	// face integration
	//----------------------------------------------------------------------
//#define ST_OMP
#ifdef _OPENMP
	std::cout << "[CSourceTerm::FaceIntegration] face integration ... "
	          << std::flush;
	//double begin_int = omp_get_wtime();
#endif
	std::vector<long> mshNodeId2sfcNodeId(n_msh_nodes, -1); // msh node id -> sfc node index
	#pragma omp parallel for private(k)
	for (i = 0; i < n_sfc_nodes; i++)
	{
		k = nodes_on_sfc[i];
		mshNodeId2sfcNodeId[k] = i;
	}
	std::vector<double> sfc_node_values(n_sfc_nodes, .0);

	const int Axisymm_sign = msh->isAxisymmetry() ? -1 : 1;
	CElement* fem = new CElement(Axisymm_sign * msh->GetCoordinateFlag());
	CElem face(1);
	face.SetFace();
	// search elements & face integration
	const long n_vec_possible_elements = vec_possible_elements.size();
	int face_nodes_local_index[8];
	double face_node_values[8];
//#define ST_OMP
#ifdef ST_OMP
#pragma omp parallel for default(none)                                         \
    shared(active_elements, NVal, G2L, msh, vec_possible_elements,             \
           unique_nodes_on_sfc,                                                \
           node_value_vector) private(elem, nodesFVal, e_nei, nodesFace, j, k, \
                                      nfaces, e_node) firstprivate(fem, face)
#endif
	for (i = 0; i < n_vec_possible_elements; i++)
	{
		CElem* elem = msh->ele_vector[vec_possible_elements[i]];
		if (!elem->GetMark())
			continue;
		if (elem->GetDimension() < 3)
			continue;
		if (active_elements && !(*active_elements)[i])
			continue;

		elem->SetOrder(msh->getOrder());

		const int nfaces = elem->GetFacesNumber();
		for (j = 0; j < nfaces; j++)
		{
			const int n_face_nodes = elem->GetElementFaceNodes(j, face_nodes_local_index);
			// 1st check
			if (elem->selected < n_face_nodes) continue;
			// 2nd check: if all nodes of the face are on the surface
			int count = 0;
			for (k = 0; k < n_face_nodes; k++)
			{
				CNode* face_node = elem->GetNode(face_nodes_local_index[k]);
#ifdef UNIQUE_OMP
				if (std::binary_search(unique_nodes_on_sfc.begin(),
				                       unique_nodes_on_sfc.end(),
				                       face_node->GetIndex()))
//            if (std::find(unique_nodes_on_sfc.begin(),
//            unique_nodes_on_sfc.end(), e_node->GetIndex()) !=
//            unique_nodes_on_sfc.end())
#else
				if (set_nodes_on_sfc.count(face_node->GetIndex()) > 0)
#endif
				{
					count++;
				}
			}
			if (count != n_face_nodes)
				continue;

			// face integration
			for (k = 0; k < n_face_nodes; k++)
			{
				CNode* face_node = elem->GetNode(face_nodes_local_index[k]);
				face_node_values[k] = node_value_vector[mshNodeId2sfcNodeId[face_node->GetIndex()]];
			}
			double fac = 1.0;
			// Not a surface face
			CElem* e_neighbor = elem->GetNeighbor(j);
			if (elem->GetDimension() == e_neighbor->GetDimension())
				fac = 0.5;

			CElem* face = new CElem(1);
			face->SetOrder(msh->getOrder());
			face->SetFace(elem, j);
			face->ComputeVolume();
			if (active_elements == NULL)
				st_boundary_elements.push_back(face);
			fem->setOrder(msh->getOrder() ? 2 : 1);
			fem->ConfigElement(face, true);
			fem->FaceIntegration(face_node_values);

			if (this->is_transfer_bc)
			{
				for (k = 0; k < n_face_nodes; k++)
					face_node_values[k] *= this->transfer_h_values[face->GetPatchIndex()];
			}

			for (k = 0; k < n_face_nodes; k++)
			{
				CNode* face_node = elem->GetNode(face_nodes_local_index[k]);
				long id = mshNodeId2sfcNodeId[face_node->GetIndex()];
				#pragma omp atomic
				sfc_node_values[id] += fac * face_node_values[k];
			}
		}
	}

	#pragma omp for
	for (i = 0; i < n_sfc_nodes; i++)
		node_value_vector[i] = sfc_node_values[i];
	#pragma omp for
	for (i = 0; i < n_msh_nodes; i++)
		msh->nod_vector[i]->SetMark(true);

	delete fem;
	// delete face;
}

/**************************************************************************
 ROCKFLOW - Funktion: DomainIntegration
 Task:  Translate distributed source term within elements into nodes value
 for all kinds of element
 Programming:
 08/2005 WW Re-Implementation
 09/2010 TF re structured some things
 **************************************************************************/
void CSourceTerm::DomainIntegration(
    CFEMesh* msh, const std::vector<long>& nodes_in_dom,
    std::vector<double>& node_value_vector) const
{
	double nodesFVal[8];

	int Axisymm = 1;                         // ani-axisymmetry
	if (msh->isAxisymmetry()) Axisymm = -1;  // Axisymmetry is true
	CElement* fem = new CElement(Axisymm * msh->GetCoordinateFlag());
	vec<CNode*> e_nodes(20);

	const size_t this_number_of_nodes(nodes_in_dom.size());
	const size_t nSize(msh->nod_vector.size());
	std::vector<long> G2L(nSize);
	std::vector<double> NVal(this_number_of_nodes);

	for (size_t i = 0; i < nSize; i++)
	{
		msh->nod_vector[i]->SetMark(false);
		G2L[i] = -1;
	}

	for (size_t i = 0; i < this_number_of_nodes; i++)
	{
		NVal[i] = 0.0;
		G2L[nodes_in_dom[i]] = i;
	}

	size_t count = 0;
	for (size_t i = 0; i < msh->ele_vector.size(); i++)
	{
		CElem* elem(msh->ele_vector[i]);
		if (!elem->GetMark()) continue;
		elem->GetNodes(e_nodes);
		size_t nn = elem->GetNodesNumber(msh->getOrder());
		if (nodes_in_dom.size() != msh->nod_vector.size())
		{
			count = 0;
			for (size_t j = 0; j < nn; j++)
			{
				for (size_t k = 0; k < this_number_of_nodes; k++)
				{
					if (*e_nodes[j] == *msh->nod_vector[nodes_in_dom[k]])
					{
						count++;
						break;
					}
				}
			}
			if (count != nn) continue;
		}
		for (size_t j = 0; j < nn; j++)
			nodesFVal[j] = node_value_vector[G2L[e_nodes[j]->GetIndex()]];
		fem->ConfigElement(elem, true);
		fem->setOrder(msh->getOrder() ? 2 : 1);
		fem->DomainIntegration(nodesFVal);
		for (size_t j = 0; j < nn; j++)
			NVal[G2L[e_nodes[j]->GetIndex()]] += nodesFVal[j];
	}

	for (size_t i = 0; i < this_number_of_nodes; i++)
		node_value_vector[i] = NVal[i];
	for (size_t i = 0; i < nSize; i++)
		msh->nod_vector[i]->SetMark(true);

	NVal.clear();
	G2L.clear();
	e_nodes.resize(0);
	delete fem;
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STDelete()
{
	long i;
	int no_st = (int)st_vector.size();
	for (i = 0; i < no_st; i++)
	{
		delete st_vector[i];
	}
	st_vector.clear();
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 01/2005 OK Implementation
 last modified:
 **************************************************************************/
void STGroupsDelete()
{
	CSourceTermGroup* m_st_group = NULL;
	std::list<CSourceTermGroup*>::const_iterator p = st_group_list.begin();
	while (p != st_group_list.end())
	{
		m_st_group = *p;
		delete m_st_group;
		// st_group_list.remove(*p);
		++p;
	}
	st_group_list.clear();
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 03/2005 OK Implementation
 last modification:
 05/2010 TF restructured a little bit
 **************************************************************************/
CSourceTerm* STGet(std::string geo_name)
{
	size_t st_vector_size(st_vector.size());
	for (size_t i = 0; i < st_vector_size; i++)
	{
		if (st_vector[i]->getGeoName().compare(geo_name) == 0)
			return st_vector[i];
	}
	return NULL;
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 04/2005 OK Implementation
 last modification:
 **************************************************************************/
CSourceTermGroup* STGetGroup(std::string pcs_type_name, std::string pcs_pv_name)
{
	CSourceTermGroup* m_st_group = NULL;
	std::list<CSourceTermGroup*>::const_iterator p_st_group =
	    st_group_list.begin();
	while (p_st_group != st_group_list.end())
	{
		m_st_group = *p_st_group;
		if (m_st_group->pcs_type_name.compare(pcs_type_name) == 0 &&
		    (m_st_group->pcs_pv_name.compare(pcs_pv_name) == 0))
			return m_st_group;
		++p_st_group;
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
void STGroupDelete(std::string pcs_type_name, std::string pcs_pv_name)
{
	CSourceTermGroup* m_st_group = NULL;
	std::list<CSourceTermGroup*>::const_iterator p = st_group_list.begin();
	while (p != st_group_list.end())
	{
		m_st_group = *p;
		if ((m_st_group->pcs_type_name.compare(pcs_type_name) == 0) &&
		    (m_st_group->pcs_pv_name.compare(pcs_pv_name) == 0))
		{
			delete m_st_group;
			st_group_list.remove(m_st_group);
			return;
		}
		++p;
	}
}

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 modification:
 05/2010
 **************************************************************************/
#if 0
void CSourceTermGroup::SetPNT(CRFProcess* pcs, CSourceTerm* st,
const int ShiftInNodeVector)
{
   std::cout << "-> Set ST: Point " << st->geo_name << " with " << convertDisTypeToString(st->getProcessDistributionType()) << "\n";
   CNodeValue *nod_val (new CNodeValue());
   nod_val->_st = st;

   // TF removed some checks - check validity of data while reading data

   nod_val->msh_node_number = m_msh->GetNODOnPNT(static_cast<const GEOLIB::Point*>(st->getGeoObj())) + ShiftInNodeVector;

   nod_val->CurveIndex = st->CurveIndex;
                                                  //WW
   nod_val->geo_node_number = nod_val->msh_node_number - ShiftInNodeVector;
   nod_val->node_value = st->geo_node_value;
   nod_val->tim_type_name = st->tim_type_name;

   std::cout << "-> node id " << nod_val->geo_node_number << " found for this ST\n";

   if (st->getProcessDistributionType() == FiniteElement::CRITICALDEPTH)
   {
      //	if (st->dis_type_name.compare("CRITICALDEPTH") == 0) {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      nod_val->node_area = 1.0;
   }

   if (st->getProcessDistributionType() == FiniteElement::NORMALDEPTH)
   {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      nod_val->node_area = 1.0;
   }

   //	if (st->dis_type_name.compare("PHILIP") == 0) { // JOD
   //		nod_val->node_distype = 10;
   //		nod_val->node_area = 1.0;
   //	}

   if (st->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      nod_val->node_area = 1.0;
   }

   if (st->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT)
   {
      nod_val->setProcessDistributionType (st->getProcessDistributionType());
      pcs->compute_domain_face_normal = true;     //WW
      CElem* elem = NULL;
      CNode* cnode = NULL;                        //WW
      for (size_t i = 0; i < m_msh->ele_vector.size(); i++)
      {
         elem = m_msh->ele_vector[i];
         if (!elem->GetMark())
            continue;
         int nn = elem->GetNodesNumber(m_msh->getOrder());
         for (long j = 0; j < nn; j++)
         {
            cnode = elem->GetNode(j);             //WW
            if (cnode->GetIndex() == (size_t)st->geo_node_number)
               st->element_st_vector.push_back(i);
         }
      }
   }
   if (st->is_transfer_bc) {
      nod_val->node_value *= st->transfer_h_values;
   }
   st->msh_node_number = nod_val->msh_node_number;
   st->geo_node_number = nod_val->geo_node_number;
   //WW        group_vector.push_back(m_node_value);
   //WW        st_group_vector.push_back(st); //OK
   pcs->st_node_value.push_back(nod_val);         //WW
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
#if 0
void CSourceTermGroup::SetLIN(CRFProcess* m_pcs, CSourceTerm* m_st,
const int ShiftInNodeVector)
{
   (void)m_pcs;
   (void)m_st;
   (void)ShiftInNodeVector;
   /*OK411
    long number_of_nodes;
    vector<long>lin_nod_vector;
    vector<double>lin_nod_val_vector;
    CGLLine* m_lin = NULL;
    CGLPolyline* m_ply = NULL;
    long *nodes = NULL;
    m_lin = m_lin->GEOGetLine(m_st->geo_id);

    if(m_lin){
    double* coordinates;
   m_ply = new CGLPolyline;
   m_ply->point_vector.push_back(m_lin->m_point1);
   m_ply->point_vector.push_back(m_lin->m_point2);
   nodes = MSHGetNodesClose(&number_of_nodes, m_ply);//CC
   lin_nod_val_vector.resize(number_of_nodes);
   for(long i = 0; i < number_of_nodes; i++){
   lin_nod_val_vector[i] =  m_st->geo_node_value / number_of_nodes;
   coordinates = new double[3];
   coordinates[0] = GetNodeX(nodes[i]);
   coordinates[1] = GetNodeY(nodes[i]);
   coordinates[2] = GetNodeZ(nodes[i]);
   m_lin->nodes_coor_vector.push_back(coordinates);
   }
   //InterpolationAlongPolyline(m_polyline,node_value_vector);
   for(long i=0; i < number_of_nodes; i++){
   CNodeValue* m_nod_val = NULL;
   m_nod_val = new CNodeValue();
   m_nod_val->msh_node_number = -1;
   m_nod_val->msh_node_number = nodes[i]+ShiftInNodeVector;
   m_nod_val->geo_node_number = nodes[i];
   m_nod_val->node_value = lin_nod_val_vector[i];
   m_nod_val->CurveIndex = m_st->CurveIndex;
   //WW        group_vector.push_back(m_node_value);
   //WW        st_group_vector.push_back(m_st); //OK
   m_pcs->st_node_value.push_back(m_nod_val);  //WW
   m_pcs->st_node.push_back(m_st); //WW
   }
   lin_nod_val_vector.clear();
   m_ply->point_vector.clear();
   delete m_ply;
   }
   else
   cout << "Warning - CSourceTermGroup::Set: LIN not found" << endl;
   */
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing
 07/2005 OK Implementation based on CSourceTermGroup::Set
 modified
 07/2010 TF substituted GEOGetPLYByName
 **************************************************************************/
#if 0
void CSourceTermGroup::SetPLY(CSourceTerm* st, int ShiftInNodeVector)
{
	CGLPolyline* old_ply (GEOGetPLYByName(st->geo_name));
	if (old_ply) {
		std::cout << "-> Set ST: Polyline " << st->geo_name << " with " << convertDisTypeToString(st->getProcessDistributionType()) << "\n";

		std::vector<long> ply_nod_vector;
		std::vector<long> ply_nod_vector_cond;
		std::vector<double> ply_nod_val_vector;

		double min_edge_length (m_msh->getMinEdgeLength());
		double new_epsilon = std::min(min_edge_length, old_ply->epsilon);
		m_msh->setMinEdgeLength (new_epsilon);
		m_msh->GetNODOnPLY(static_cast<const GEOLIB::Polyline*>(st->getGeoObj()), ply_nod_vector);
		m_msh->setMinEdgeLength (min_edge_length);

		if (st->isCoupled()) {
			SetPolylineNodeVectorConditional(st, ply_nod_vector, ply_nod_vector_cond);
		}

		SetPolylineNodeValueVector(st, ply_nod_vector, ply_nod_vector_cond, ply_nod_val_vector);

		st->SetNodeValues(ply_nod_vector, ply_nod_vector_cond, ply_nod_val_vector, ShiftInNodeVector);
	} // end polyline
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
#if 0
void CSourceTermGroup::SetDMN(CSourceTerm *m_st, const int ShiftInNodeVector)
{
   std::vector<long> dmn_nod_vector;
   std::vector<double> dmn_nod_val_vector;
   std::vector<long> dmn_nod_vector_cond;

   if (m_st->getProcessDistributionType() == FiniteElement::ANALYTICAL) {
      GEOGetNodesInMaterialDomain(m_msh, m_st->analytical_material_group,
         dmn_nod_vector, false);
      size_t number_of_nodes (dmn_nod_vector.size());
      dmn_nod_val_vector.resize(number_of_nodes);
      for (size_t i = 0; i < number_of_nodes; i++)
         dmn_nod_val_vector[i] = 0;
   } else {
      dmn_nod_vector.resize(m_msh->nod_vector.size());
      for (size_t i=0; i<m_msh->nod_vector.size(); i++)
         dmn_nod_vector[i] = i;
      dmn_nod_val_vector.resize(dmn_nod_vector.size());
      if (m_st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN)
      {
          for (size_t i=0; i<dmn_nod_vector.size(); i++)
        	  dmn_nod_val_vector[i] = m_st->geo_node_value;
      }
   }


   if (m_st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN || m_st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      m_st->DomainIntegration(m_msh, dmn_nod_vector, dmn_nod_val_vector);
   }

   m_st->SetNodeValues(dmn_nod_vector, dmn_nod_vector_cond,
      dmn_nod_val_vector, ShiftInNodeVector);

}
#endif


/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 09/2010 WW  For the case that nodes are directly given
 **************************************************************************/
#if 0
void CSourceTermGroup::SetSFC(CSourceTerm* m_st, const int ShiftInNodeVector)
{
   std::vector<long> sfc_nod_vector;
   std::vector<long> sfc_nod_vector_cond;
   std::vector<double> sfc_nod_val_vector;
   Surface* m_sfc = NULL;

   m_sfc = GEOGetSFCByName(m_st->geo_name);       //CC

   if (m_sfc)
   {
      std::cout << "-> Set ST: Surface " << m_st->geo_name << " with " << convertDisTypeToString(m_st->getProcessDistributionType()) << "\n";

// NW use TF's node searching function instead of SetSurfaceNodeVector
//      SetSurfaceNodeVector(m_sfc, sfc_nod_vector);
      GEOLIB::Surface const* sfc(static_cast<const GEOLIB::Surface*> (m_st->getGeoObj()));
      std::vector<size_t> msh_nod_vec;
      m_msh->GetNODOnSFC(sfc, msh_nod_vec);
      for (size_t k(0); k < msh_nod_vec.size(); k++) {
         sfc_nod_vector.push_back (msh_nod_vec[k]);
      }
      std::cout << "-> " << sfc_nod_vector.size() << " nodes are found for this ST" << "\n";
      if (m_st->isCoupled())
         m_st->SetSurfaceNodeVectorConditional(sfc_nod_vector,
            sfc_nod_vector_cond);
      //		m_st->SetDISType();
      SetSurfaceNodeValueVector(m_st, m_sfc, sfc_nod_vector,
         sfc_nod_val_vector);
      m_st->SetNodeValues(sfc_nod_vector, sfc_nod_vector_cond,
         sfc_nod_val_vector, ShiftInNodeVector);

   }                                              // end surface
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task:
 Programing:
 11/2007 JOD Implementation
 **************************************************************************/
void CSourceTerm::SetNOD()
{
	std::vector<long> nod_vector;
	std::vector<long> nod_vector_cond;
	std::vector<double> nod_val_vector;
	int ShiftInNodeVector;

	nod_vector.push_back(msh_node_number);
	nod_vector_cond.push_back(msh_node_number);
	nod_val_vector.push_back(geo_node_value);

	/*nod_vector[0] = msh_node_number;
	 nod_vector_cond[0] = msh_node_number;
	 nod_val_vector[0] =geo_node_value;*/
	ShiftInNodeVector = 0;

	SetNodeValues(nod_vector, nod_vector_cond, nod_val_vector,
	              ShiftInNodeVector);
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetSurfaceNodeVector(Surface* m_sfc,
                                            std::vector<long>& sfc_nod_vector)
{
	m_msh->GetNODOnSFC(m_sfc, sfc_nod_vector);
}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::SetPolylineNodeVectorConditional(
    CSourceTerm* st, std::vector<long>& ply_nod_vector,
    std::vector<long>& ply_nod_vector_cond)
{
	size_t assembled_mesh_node, number_of_nodes;

	if (st->node_averaging)
	{
		if (m_msh_cond)
		{
			if (pcs_type_name == "RICHARDS_FLOW")
			{
				//				m_msh_cond->GetNODOnPLY(m_ply,
				// ply_nod_vector_cond);
				m_msh_cond->GetNODOnPLY(
				    static_cast<const GEOLIB::Polyline*>(st->getGeoObj()),
				    ply_nod_vector_cond);
				number_of_nodes = ply_nod_vector_cond.size();
				assembled_mesh_node = ply_nod_vector[0];
				ply_nod_vector.resize(number_of_nodes);
				for (size_t i = 0; i < number_of_nodes; i++)
					ply_nod_vector[i] = assembled_mesh_node;
			}  // end richards
			else if (pcs_type_name == "GROUNDWATER_FLOW")  // JOD 4.10.01
			{
				number_of_nodes = ply_nod_vector.size();
				//				m_msh_cond->GetNODOnPLY(m_ply,
				// ply_nod_vector_cond);
				m_msh_cond->GetNODOnPLY(
				    static_cast<const GEOLIB::Polyline*>(st->getGeoObj()),
				    ply_nod_vector_cond);
				assembled_mesh_node = ply_nod_vector_cond[0];
				ply_nod_vector_cond.resize(number_of_nodes);
				for (size_t i = 0; i < number_of_nodes; i++)
					ply_nod_vector_cond[i] = assembled_mesh_node;
			}  // end overland, groundwater
			else
				std::cout << "Warning in "
				             "CSourceTermGroup::"
				             "SetPolylineNodeVectorConditional - no area "
				             "assembly for this process"
				          << "\n";
		}  // end mesh_cond
		else
			std::cout
			    << "Warning in CSourceTermGroup::SetPLY - no MSH_COND data"
			    << "\n";
	}  // end area_assembly
	else
	{
		number_of_nodes = ply_nod_vector.size();
		ply_nod_vector_cond.resize(number_of_nodes);
		st->SetNOD2MSHNOD(ply_nod_vector, ply_nod_vector_cond);
	}  // end !area_assembly
}


/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::InterpolatePolylineNodeValueVector(
    CGLPolyline* m_ply, std::vector<double>& Distribed,
    std::vector<double>& ply_nod_vector)
{
	for (long k = 0; k < (long)DistribedBC.size(); k++)
	{
		for (long l = 0; l < (long)m_ply->point_vector.size(); l++)
		{
			if (PointsHaveDistribedBC[k] == m_ply->point_vector[l]->id)
			{
				if (fabs(DistribedBC[k]) < MKleinsteZahl)
					DistribedBC[k] = 1.0e-20;
				m_ply->point_vector[l]->setPropert(Distribed[k]);
				break;
			}
		}
	}

	InterpolationAlongPolyline(m_ply, ply_nod_vector);
}

void CSourceTerm::InterpolatePolylineNodeValueVector(
    std::vector<double> const& nodes_as_interpol_points,
    std::vector<double>& node_values) const
{
	std::vector<double> interpolation_points;
	std::vector<double> interpolation_values;

	GEOLIB::Polyline const* ply(
	    dynamic_cast<GEOLIB::Polyline const*>(this->getGeoObj()));

	for (size_t i(0); i < DistribedBC.size(); i++)
	{
		for (size_t j = 0; j < ply->getNumberOfPoints(); j++)
		{
			if ((size_t)(PointsHaveDistribedBC[i]) == ply->getPointID(j))
			{
				interpolation_points.push_back(ply->getLength(j));
				if (fabs(DistribedBC[i]) < MKleinsteZahl)
					interpolation_values.push_back(1.0e-20);
				else
					interpolation_values.push_back(DistribedBC[i]);
				break;
			}
		}
	}

	MathLib::PiecewiseLinearInterpolation(
	    interpolation_points, interpolation_values, nodes_as_interpol_points,
	    node_values);
}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
#if 0
void CSourceTermGroup::SetPolylineNodeValueVector(CSourceTerm* st, CGLPolyline * old_ply,
		const std::vector<long>& ply_nod_vector,
		std::vector<long>& ply_nod_vector_cond,
		std::vector<double>& ply_nod_val_vector)
{
	long number_of_nodes = (long) ply_nod_vector.size();
	ply_nod_val_vector.resize(number_of_nodes);

	if (st->getProcessDistributionType() == FiniteElement::LINEAR
			|| st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN) {
		st->InterpolatePolylineNodeValueVector(old_ply, st->DistribedBC, ply_nod_val_vector);
	} else if (st->getProcessDistributionType() == FiniteElement::SYSTEM_DEPENDENT) {
		CRFProcess* m_pcs = NULL;
		m_pcs = PCSGet(pcs_type_name);
		m_pcs->compute_domain_face_normal = true; //WW
		long no_face = (long) m_msh->face_vector.size();
		for (long i = 0; i < no_face; i++) {
			int node_on_line = 0;
			int no_vertex = m_msh->face_vector[i]->GetVertexNumber();
			for (long jj = 0; jj < no_vertex; jj++) {
				for (long kk = 0; kk < number_of_nodes; kk++) {
					if (ply_nod_vector[kk]
							== m_msh->face_vector[i]->GetNodeIndex(jj)) node_on_line++;
				} // end nodes
			} // end vertices
			if (node_on_line == 2) st->element_st_vector.push_back(
					m_msh->face_vector[i]->GetOwner()->GetIndex());
		} // end faces
	} // end system dependent
	else //WW
	{
		for (long i = 0; i < number_of_nodes; i++) {
			ply_nod_val_vector[i] = st->geo_node_value;
			//			if (st->dis_type == 12)
			if (st->getProcessDistributionType() == FiniteElement::CONSTANT_GEO)
				ply_nod_val_vector[i] = st->geo_node_value / (double) number_of_nodes; // distribute flow to nodes along polyline. To do.. 4.10.06
		}
	}
	/*|| st->getProcessDistributionType() == PHILIP */
	if (st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN
			|| st->getProcessDistributionType()
					== FiniteElement::LINEAR_NEUMANN
			|| st->getProcessDistributionType() == FiniteElement::GREEN_AMPT) {
		if (m_msh->GetMaxElementDim() == 1) // 1D  //WW MB
			st->DomainIntegration(m_msh, ply_nod_vector, ply_nod_val_vector);
		else st->EdgeIntegration(m_msh, ply_nod_vector, ply_nod_val_vector);
	}

	if ( st->getProcessDistributionType() == FiniteElement::CRITICALDEPTH
			|| st->getProcessDistributionType() == FiniteElement::NORMALDEPTH
			|| st->getProcessDistributionType() == FiniteElement::ANALYTICAL) {
		st->node_value_vectorArea.resize(number_of_nodes);
		for (long i = 0; i < number_of_nodes; i++)
			st->node_value_vectorArea[i] = 1.0; //Element width !
		st->EdgeIntegration(m_msh, ply_nod_vector, st->node_value_vectorArea);
	}

	if (st->isCoupled() && st->node_averaging)
		AreaAssembly(st, ply_nod_vector_cond, ply_nod_val_vector);
}
#endif

#if 0
// 09/2010 TF
void CSourceTermGroup::SetPolylineNodeValueVector(CSourceTerm* st,
		std::vector<long> const & ply_nod_vector,
		std::vector<long>& ply_nod_vector_cond,
		std::vector<double>& ply_nod_val_vector)
{
	size_t number_of_nodes(ply_nod_vector.size());
	ply_nod_val_vector.resize(number_of_nodes);

	FiniteElement::DistributionType distype(st->getProcessDistributionType());

	// linear
	if (distype == FiniteElement::LINEAR || distype == FiniteElement::LINEAR_NEUMANN) {
		// fetch data for the linear interpolation
		GEOLIB::Polyline const* polyline (dynamic_cast<GEOLIB::Polyline const*>(st->getGeoObj()));
		if (polyline) {
			std::vector<double> nodes_as_interpol_points;
			m_msh->getPointsForInterpolationAlongPolyline (polyline, nodes_as_interpol_points);
			st->InterpolatePolylineNodeValueVector(nodes_as_interpol_points, ply_nod_val_vector);
		}
	} else if (distype == FiniteElement::SYSTEM_DEPENDENT) { //System Dependented YD
		CRFProcess* m_pcs(PCSGet(pcs_type_name));
		m_pcs->compute_domain_face_normal = true; //WW
		long no_face = (long) m_msh->face_vector.size();
		for (long i = 0; i < no_face; i++) {
			int node_on_line = 0;
			int no_vertex = m_msh->face_vector[i]->GetVertexNumber();
			for (long jj = 0; jj < no_vertex; jj++) {
				for (size_t kk = 0; kk < number_of_nodes; kk++) {
					if (ply_nod_vector[kk] == (m_msh->face_vector[i]->GetNodeIndex(jj)))
						node_on_line++;
				} // end nodes
			} // end vertices
			if (node_on_line == 2) st->element_st_vector.push_back(
					m_msh->face_vector[i]->GetOwner()->GetIndex());
		} // end faces
	} // end system dependent
   else if (distype == FiniteElement::FUNCTION) // 25.08.2011. WW
   {
      for (size_t i = 0; i < number_of_nodes; i++)
      {
         double const*const pnt (m_msh->nod_vector[ply_nod_vector[i]]->getData());
         ply_nod_val_vector[i] = st->dis_linear_f->getValue(pnt[0], pnt[1], pnt[2]);
      }
   }
   else if (distype == FiniteElement::GRADIENT) //NW
	{
		for (size_t i(0); i < number_of_nodes; i++) {
			double z = m_msh->nod_vector[ply_nod_vector[i]]->getData()[2];
			ply_nod_val_vector[i] = st->gradient_ref_depth_gradient * (st->gradient_ref_depth - z) + st->gradient_ref_depth_value;
		}
	}
	else //WW
	{
		for (size_t i = 0; i < number_of_nodes; i++) {
			ply_nod_val_vector[i] = st->geo_node_value;
			if (st->getProcessDistributionType() == FiniteElement::CONSTANT_GEO)
				ply_nod_val_vector[i] = st->geo_node_value / (double) number_of_nodes;
		}
	}
	if (st->is_transfer_bc) {
		for (size_t i = 0; i < number_of_nodes; i++) {
			ply_nod_val_vector[i] *= st->transfer_h_values;
		}
	}
	if (distype == FiniteElement::CONSTANT_NEUMANN
			|| distype == FiniteElement::LINEAR_NEUMANN
			|| distype == FiniteElement::GREEN_AMPT
			|| st->is_transfer_bc)
	{
		if (m_msh->GetMaxElementDim() == 1) // 1D  //WW MB
			st->DomainIntegration(m_msh, ply_nod_vector,
					ply_nod_val_vector);
		else st->EdgeIntegration(m_msh, ply_nod_vector, ply_nod_val_vector);
	}

	if (distype == FiniteElement::CRITICALDEPTH
			|| distype == FiniteElement::NORMALDEPTH
			|| distype == FiniteElement::ANALYTICAL) {
		st->node_value_vectorArea.resize(number_of_nodes);
		for (size_t i = 0; i < number_of_nodes; i++)
			st->node_value_vectorArea[i] = 1.0; //Element width !
		st->EdgeIntegration(m_msh, ply_nod_vector,
				st->node_value_vectorArea);
	}

	if (st->isCoupled() && st->node_averaging)
		AreaAssembly(st, ply_nod_vector_cond, ply_nod_val_vector);
}
#endif

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTermGroup::AreaAssembly(
    CSourceTerm* st,
    const std::vector<long>& ply_nod_vector_cond,
    std::vector<double>& ply_nod_val_vector)
{
	if (pcs_type_name == "RICHARDS_FLOW")
	{
		if (m_msh_cond->GetMaxElementDim() == 1)  // 1D  //WW MB
			st->DomainIntegration(m_msh_cond, ply_nod_vector_cond,
			                      ply_nod_val_vector);
		else
			st->EdgeIntegration(m_msh_cond, ply_nod_vector_cond,
			                    ply_nod_val_vector);
		double sum_node_value = 0;
		for (size_t i = 0; i < ply_nod_val_vector.size(); i++)
			sum_node_value += ply_nod_val_vector[i];
		for (size_t i = 0; i < ply_nod_val_vector.size(); i++)
			ply_nod_val_vector[i] /= sum_node_value;
	}
}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
#if 0
void CSourceTermGroup::SetSurfaceNodeValueVector(CSourceTerm* st,
		Surface* m_sfc, std::vector<long>&sfc_nod_vector,
		std::vector<double>&sfc_nod_val_vector)
{
   // CRFProcess* m_pcs = NULL;
   // m_pcs = PCSGet(pcs_type_name);
   long number_of_nodes = (long) sfc_nod_vector.size();
   sfc_nod_val_vector.resize(number_of_nodes);

   double node_value = st->geo_node_value;
   if (st->is_transfer_bc)
	   node_value *= st->transfer_h_values;
   for (long i = 0; i < number_of_nodes; i++)
      sfc_nod_val_vector[i] = node_value;
   // KR & TF - case not used
   //	if (m_st->dis_type == 12) //To do. 4.10.06
   //		for (long i = 0; i < number_of_nodes; i++)
   //			sfc_nod_val_vector[i] = m_st->geo_node_value
   //					/ (double) number_of_nodes;

   //	if (st->dis_type == 2 || st->dis_type == 4) { // Piecewise linear distributed, polygon-wise WW
   if (st->getProcessDistributionType() == FiniteElement::LINEAR || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN)
   {
      CGLPolyline* m_ply = NULL;
      std::vector<CGLPolyline*>::iterator p =
         m_sfc->polyline_of_surface_vector.begin();
      p = m_sfc->polyline_of_surface_vector.begin();
      while (p != m_sfc->polyline_of_surface_vector.end())
      {
         m_ply = *p;
         for (long k = 0; k < (long) st->DistribedBC.size(); k++)
         {
            for (long l = 0; l < (long) m_ply->point_vector.size(); l++)
            {
               if (st->PointsHaveDistribedBC[k]
                  == m_ply->point_vector[l]->id)
               {
                  if (fabs(st->DistribedBC[k]) < MKleinsteZahl)
                     st->DistribedBC[k] = 1.0e-20;
                  m_ply->point_vector[l]->setPropert (st->DistribedBC[k]);
                  break;
               }                                  // end l
            }                                     // end k
         }                                        // end polyline
         // InterpolationAlongPolyline(m_polyline, node_value_vector);
         p++;
      }                                           // end while
   }                                              // end linear

   // neumann, Green-Ampt, Philip
   //	if (st->dis_type == 3 || st->dis_type == 4 || st->dis_type == 10
   //				|| st->dis_type == 11) {
                                                  /*|| st->getProcessDistributionType() == PHILIP */
   if (st->getProcessDistributionType() == FiniteElement::CONSTANT_NEUMANN || st->getProcessDistributionType() == FiniteElement::LINEAR_NEUMANN
      || st->getProcessDistributionType() == FiniteElement::GREEN_AMPT)
   {
      if (m_msh->GetMaxElementDim() == 2)         // For all meshes with 1-D or 2-D elements
         st->DomainIntegration(m_msh, sfc_nod_vector, sfc_nod_val_vector);
      else if (m_msh->GetMaxElementDim() == 3)    // For all meshes with 3-D elements
         st->FaceIntegration(m_msh, sfc_nod_vector, sfc_nod_val_vector);
   }                                              // end neumann
  else if (st->getProcessDistributionType() == FiniteElement::FUNCTION) // 25.08.2011. WW
   {
      for (size_t j = 0; j < sfc_nod_vector.size(); j++)
      {
         double const*const pnt (m_msh->nod_vector[sfc_nod_vector[j]]->getData());
         sfc_nod_val_vector[j] = st->dis_linear_f->getValue(pnt[0], pnt[1], pnt[2]);
      }
   }

}
#endif

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::SetSurfaceNodeVectorConditional(
    std::vector<long>& sfc_nod_vector, std::vector<long>& sfc_nod_vector_cond)
{
	long number_of_nodes;
	number_of_nodes = (long)sfc_nod_vector.size();

	sfc_nod_vector_cond.resize(number_of_nodes);
	SetNOD2MSHNOD(sfc_nod_vector, sfc_nod_vector_cond);
}

/**************************************************************************
 MSHLib-Method:
 Task:
 Programing:
 11/2007 JOD
 last modification:
 **************************************************************************/
void CSourceTerm::SetNodeValues(const std::vector<long>& nodes,
                                const std::vector<long>& nodes_cond,
                                const std::vector<double>& node_values,
                                int ShiftInNodeVector)
{
	CNodeValue* m_nod_val = NULL;
	size_t number_of_nodes(nodes.size());

	size_t st_id = getID();
	_pcs->st_node_value[st_id].clear();

	for (size_t i = 0; i < number_of_nodes; i++)
	{
		m_nod_val = new CNodeValue();
		m_nod_val->_st = this;
		m_nod_val->msh_node_number = nodes[i] + ShiftInNodeVector;
		m_nod_val->geo_node_number = nodes[i];
		m_nod_val->setProcessDistributionType(getProcessDistributionType());
		m_nod_val->node_value = node_values[i];
		m_nod_val->CurveIndex = CurveIndex;
		if (_coupled)  // JOD 4.7.10
		{
			m_nod_val->msh_node_number_conditional = nodes_cond[i];
			// JOD 4.10.01
			if (getProcessType() == FiniteElement::GROUNDWATER_FLOW &&
			    node_averaging)
			{
				double weights = 0;
				for (size_t j = 0; j < number_of_nodes; j++)
				{
					m_nod_val->msh_node_numbers_averaging.push_back(nodes[j]);
					m_nod_val->msh_node_weights_averaging.push_back(
					    node_values[j]);
					weights += node_values[j];
				}
				for (size_t j = 0; j < number_of_nodes; j++)
					m_nod_val->msh_node_weights_averaging[j] /= weights;
			}
		}
		// WW        group_vector.push_back(m_node_value);
		// WW        st_group_vector.push_back(m_st); //OK
		//		if (getProcessDistributionType() == RIVER) {
		//			m_nod_val->node_value = node_value_vectorArea[i];
		//			m_nod_val->node_parameterA = node_value_vectorA[i];
		//			m_nod_val->node_parameterB = node_value_vectorB[i];
		//			m_nod_val->node_parameterC = node_value_vectorC[i];
		//			m_nod_val->node_parameterD = node_value_vectorD[i];
		//			m_nod_val->node_parameterE = node_value_vectorE[i];
		//		}
		//		if (dis_type == 6 || dis_type == 8 || dis_type == 9) // critical
		// depth, normal depth, analytical
		if (getProcessDistributionType() == FiniteElement::CRITICALDEPTH ||
		    getProcessDistributionType() == FiniteElement::NORMALDEPTH ||
		    getProcessDistributionType() == FiniteElement::ANALYTICAL)
		{
			m_nod_val->node_value = node_value_vectorArea[i];
			// CMCD bugfix on 4.9.06
			m_nod_val->node_area = node_value_vectorArea[i];
		}
		_pcs->st_node_value[st_id].push_back(m_nod_val);  // WW
	}                                                     // end nodes
}

void CSourceTerm::SetNOD2MSHNOD(std::vector<long>& nodes,
                                std::vector<long>& conditional_nodes)
{
	CFEMesh* m_msh_cond(FEMGet(pcs_type_name_cond));
	CFEMesh* m_msh_this(FEMGet(convertProcessTypeToString(getProcessType())));

	for (size_t i = 0; i < nodes.size(); i++)
	{
		const GEOLIB::Point pnt(m_msh_this->nod_vector[nodes[i]]->getData());

		conditional_nodes[i] = m_msh_cond->GetNODOnPNT(&pnt);
	}
}

void CSourceTerm::SetNOD2MSHNOD(const std::vector<size_t>& nodes,
                                std::vector<size_t>& conditional_nodes) const
{
	CFEMesh* m_msh_cond(FEMGet(pcs_type_name_cond));
	CFEMesh* m_msh_this(FEMGet(convertProcessTypeToString(getProcessType())));

	for (size_t i = 0; i < nodes.size(); i++)
	{
		const GEOLIB::Point pnt(m_msh_this->nod_vector[nodes[i]]->getData());

		conditional_nodes[i] = m_msh_cond->GetNODOnPNT(&pnt);
	}
}

int CSourceTerm::getID()
{
	if (st_id < 0)
	{
		for (st_id = 0; st_id < (int)st_vector.size(); st_id++)
			if (st_vector[st_id] == this) break;
	}
	return st_id;
}

/**************************************************************************
 GeoSys source term function:
 02/2009 WW Implementation
 **************************************************************************/
void CSourceTerm::DirectAssign(long ShiftInNodeVector)
{
	CRFProcess* m_pcs = PCSGet(convertProcessTypeToString(getProcessType()));

	{
		std::string line_string;
		std::string st_file_name;
		std::stringstream in;
		long n_index;
		double n_val;

		//========================================================================
		// File handling
		std::ifstream d_file(fname.c_str(), std::ios::in);
		// if (!st_file.good()) return;

		if (!d_file.good())
		{
			std::cout << "! Error in direct node source terms: Could not find "
			             "file:!\n" << fname << "\n";
			abort();
		}
		// Rewind the file
		d_file.clear();
		d_file.seekg(0L, std::ios::beg);
		//========================================================================
		while (!d_file.eof())
		{
			line_string = GetLineFromFile1(&d_file);
			if (line_string.find("#STOP") != std::string::npos) break;

			in.str(line_string);
			in >> n_index >> n_val;
			in.clear();
			//
			CNodeValue* m_nod_val(new CNodeValue());
			m_nod_val->msh_node_number = n_index + ShiftInNodeVector;
			m_nod_val->geo_node_number = n_index;
			m_nod_val->setProcessDistributionType(getProcessDistributionType());
			m_nod_val->node_value = n_val;
			m_nod_val->CurveIndex = CurveIndex;
			m_pcs->st_node_value[getID()].push_back(m_nod_val);
			//
		}  // eof
	}
}


/**************************************************************************
 FEMLib-Method:
 Task: Analytical diffusion in matrix. Replaces matrix. See paper to be issued.
 Programing:
 11/2005 CMCD Implementation
 04/2006 Moved from CSourceTermGroup and changed the arguments
 last modification:
 04/2006 CMCD Updated
 **************************************************************************/
// , CSourceTerm *m_st)
#if 0
double CSourceTerm::GetAnalyticalSolution(long location)
{
   int idx, n;
   int size, process_no;
   long i;
   long step, no_terms_included;
   double value, source, gradient, ref_value = 0.0;
   double timevalue;
   double fac = 1.0;
   double temp_time, temp_value;
   double pi = 3.1415926;
   double D = this->analytical_diffusion;
   double ne = this->analytical_porosity;
   double tort = this->analytical_tortousity;
   double Kd = this->analytical_linear_sorption_Kd;
   double rho = this->analytical_matrix_density;
   double Dtrans = (D * ne) / ((ne + Kd * rho) * tort);
   double Dsteady = D * ne / tort;
   double t0, tn, tnn, val1, val2, node_area;
   double tvol, vol, flux_area, tflux_area;
   double mass_solute_present, mass_to_remove;
   //WW  bool out = false;
   //WW  int dimension = this->analytical_material_group;
   std::string process;
   CRFProcess* m_pcs = NULL;
   m_pcs = PCSGet(convertProcessTypeToString(this->getProcessType()), convertPrimaryVariableToString(this->getProcessPrimaryVariable()));
   CFEMesh* m_msh = m_pcs->m_msh;                 //WW
   CElem* Ele = NULL;
   long node_number = location;                   //WW m_pcs->st_node_value[location]->msh_node_number;
   CNode* Node = m_msh->nod_vector[node_number];
   double area = m_pcs->st_node_value[location]->node_area;
   std::vector<double> time_history;
   std::vector<double> value_history;
   //Initialise
   time_history.clear();
   value_history.clear();
   t0 = tn = tnn = source = gradient = val1 = val2 = node_area = flux_area
      = 0.0;
   idx = m_pcs->GetNodeValueIndex(convertPrimaryVariableToString(this->getProcessPrimaryVariable()));
   value = m_pcs->GetNodeValue(node_number, idx);
   if (value < MKleinsteZahl)
      value = 0.0;
   timevalue = aktuelle_zeit;
   step = aktueller_zeitschritt;
   if (step < 10)
      step = 10;
   size = (int) analytical_processes.size();
   process_no = 0;
   //Domain or Polyline
   for (i = 0; i < size; i++)
   {
      if (analytical_processes[i] == convertPrimaryVariableToString(this->getProcessPrimaryVariable()))
      {
         if (this->getGeoType () == GEOLIB::POLYLINE)
         {
            if (this->getGeoName().compare(
               analytical_processes_polylines[i]) == 0)
               process_no = i;
         }
         //			if (this->geo_type_name.compare("DOMAIN") == 0)
         if (this->getGeoType () == GEOLIB::GEODOMAIN)
            process_no = i;
      }
   }
   //Identify process
   if (process_no == 1)
   {
      process_no = process_no;
   }
   process_no *= 2;                               //first column time, second column value, hence two columns per process;

   //If time step require new calculation of source term then start
   if ((aktueller_zeitschritt - 1) % this->resolution == 0) {
      //Save data in a vector attached to the nodes
      this->SetNodePastValue(node_number, process_no, 0, timevalue);
      this->SetNodePastValue(node_number, process_no + 1, 0, value);

      //Recall historical data
      ref_value = this->GetNodePastValueReference(node_number, (process_no
         / 2));
      if ((size_t)step > this->number_of_terms)
         no_terms_included = this->number_of_terms;
      else
         no_terms_included = step;
      for (i = 0; i < no_terms_included; i++)
      {
         temp_time = this->GetNodePastValue(node_number, process_no, i);
         temp_value
            = (this->GetNodePastValue(node_number, process_no + 1, i)
            - ref_value);
         time_history.push_back(temp_time);
         value_history.push_back(temp_value);
      }

      //Calculate individual terms and sum for gradient
      for (i = no_terms_included - 1; i > 0; i--)
      {
         t0 = time_history[0];
         if (i == no_terms_included - 1)
            tn = (t0 - time_history[i]) + (time_history[i - 1]
               - time_history[i]);
         else
            tn = t0 - time_history[i + 1];
         tnn = t0 - time_history[i];
         val1 = 1 / (sqrt(pi * Dtrans * tn));
         val2 = 1 / (sqrt(pi * Dtrans * tnn));
         gradient += ((val2 - val1) * value_history[i]);
      }
      tn = t0 - time_history[1];
      tnn = 0;
      val1 = 1 / (sqrt(pi * Dtrans * tn));
      gradient -= (val1 * value_history[0]);

      //Area calculations
      mass_solute_present = 1.e99;                //Initially set value very high

      //Area for lines, triangles and quads in domain.
      //  if (area < 0) {//set in set source terms function, domain area = -1 to start with
      if (area < DBL_MIN) { // HS 04.2008
			tvol = 0.0;
			tflux_area = 0.0;
			for (i = 0; i < (int) Node->getConnectedElementIDs().size(); i++) {
				Ele = m_msh->ele_vector[Node->getConnectedElementIDs()[i]];
				vol = Ele->GetVolume(); //Assuming 1m thickness
				flux_area = Ele->GetFluxArea(); //Real thickness for a fracture
				n = Ele->GetVertexNumber();
				tvol += (vol / n);
				tflux_area += (flux_area / n);
			}
			node_area = tvol * 2.; //* 2 because the diffusion is in two direction perpendicular to the fracture
			mass_solute_present = tflux_area * tvol * value;
		}
		//Area for polylines
		else
		 node_area = area;

      //factor for conversion to energy for temperature if necessary
      fac = this->factor;
      source = gradient * Dsteady * fac * node_area;
      mass_to_remove = fabs(source) * dt;
      if (mass_to_remove > mass_solute_present)
      {
         source *= (mass_solute_present / mass_to_remove);
      }
      this->SetNodeLastValue(node_number, (process_no / 2), source);

   }                                              // If new source term calculation not required
   else
      source = this->GetNodeLastValue(node_number, (process_no / 2));
   return source;
}
#endif

/**************************************************************************
 FEMLib-Method:
 Task: master write function
 Programing:
 11/2005 CMCD Implementation, functions to access and write the time history
 data
 of nodes
 last modification:
 **************************************************************************/
void CSourceTerm::SetNodePastValue(long n, int idx, int pos, double value)
{
	bool endstepreached = false;
	pos = pos;  // WW

	// Check whether this is the first call
	CRFProcess* m_pcs = NULL;
	m_pcs = PCSGet(convertProcessTypeToString(getProcessType()),
	               convertPrimaryVariableToString(getProcessPrimaryVariable()));
	if (!m_pcs)  // OK
	{
		std::cout << "Warning in SetNodePastValue - no PCS data"
		          << "\n";
		return;
	}

	size_t size1(m_pcs->nod_val_vector.size());

	// Create memory for the values
	if (node_history_vector.empty())
	{
		// WW     int number_of_terms = max_no_terms;
		for (size_t k = 0; k < size1; k++)
		{
			NODE_HISTORY* nh = new NODE_HISTORY;
			CreateHistoryNodeMemory(nh);
			node_history_vector.push_back(nh);
		}
		for (size_t k = 0; k < size1; k++)
		{
			for (size_t j = 0; j < _no_an_sol; j++)
				node_history_vector[k]->value_reference.push_back(-1.0);
		}
	}

	// Store the first set of values as reference values
	int flipflop = idx % 2;
	if (aktueller_zeitschritt == 1)
	{
		// if (size2 == idx)
		if (flipflop == 1)
			node_history_vector[n]->value_reference[(idx - 1) / 2] = value;
	}
	//  size2 = (int) node_history_vector[n]->value_reference.size();
	size_t steps = 10;
	if (aktueller_zeitschritt > steps) steps = aktueller_zeitschritt;
	if (_max_no_terms >= steps)
		steps = aktueller_zeitschritt;
	else
	{
		steps = _max_no_terms;
		endstepreached = true;
	}
	// Enter the value and push the other values back
	if (!endstepreached)
	{
		for (size_t k = steps - 1; k > 0; k--)
			node_history_vector[n]->value_store[idx][k] =
			    node_history_vector[n]->value_store[idx][k - 1];
		node_history_vector[n]->value_store[idx][0] = value;
	}
	double cutvalue = 0.0;
	double nextvalue = 0.0;
	long no_steps_past_cutof = 0;
	if (endstepreached)
	{
		no_steps_past_cutof = aktueller_zeitschritt - _max_no_terms;
		cutvalue = node_history_vector[n]->value_store[idx][steps - 1];
		nextvalue = node_history_vector[n]->value_store[idx][steps - 2];
		node_history_vector[n]->value_store[idx][steps - 1] =
		    (cutvalue * (double)(no_steps_past_cutof) + nextvalue) /
		    ((double)no_steps_past_cutof + 1);
		for (size_t k = steps - 2; k > 0; k--)
			node_history_vector[n]->value_store[idx][k] =
			    node_history_vector[n]->value_store[idx][k - 1];
		node_history_vector[n]->value_store[idx][0] = value;
	}
}

void CSourceTerm::SetNodeLastValue(long n, int idx, double value)
{
	size_t size(node_history_vector[n]->last_source_value.size());
	if (size == 0)
	{
		for (size_t i = 0; i < _no_an_sol; i++)
			node_history_vector[n]->last_source_value.push_back(0);
	}
	node_history_vector[n]->last_source_value[idx] = value;
}

double CSourceTerm::GetNodeLastValue(long n, int idx)
{
	double value = 0.0;
	value = node_history_vector[n]->last_source_value[idx];
	return value;
}

double CSourceTerm::GetNodePastValue(long n, int idx, int pos)
{
	double value;
	value = node_history_vector[n]->value_store[idx][pos];
	return value;
}

double CSourceTerm::GetNodePastValueReference(long n, int idx)
{
	double value;
	value = node_history_vector[n]->value_reference[idx];
	return value;
}

void CSourceTerm::CreateHistoryNodeMemory(NODE_HISTORY* nh)
{
	size_t s_col = _no_an_sol * 2;
	size_t s_row = number_of_terms;

	nh->value_store = new double* [s_col];
	for (size_t i = 0; i < s_col; i++)
		nh->value_store[i] = new double[s_row];

	for (size_t k = 0; k < s_col; k++)
	{
		for (size_t j = 0; j < s_row; j++)
			nh->value_store[k][j] = 0.0;
	}
}

void CSourceTerm::DeleteHistoryNodeMemory()
{
	size_t size(node_history_vector.size());
	size_t s_row = _no_an_sol * 2;

	if (size > 0)
	{
		for (size_t j = 0; j < size; j++)
		{
			for (size_t i = 0; i < s_row; i++)
				delete node_history_vector[j]->value_store[i];
			delete node_history_vector[j]->value_store;
		}
		node_history_vector.clear();
	}
}
