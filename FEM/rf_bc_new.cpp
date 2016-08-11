/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_bc_new.h"

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>

#include "display.h"
#include "FileToolsRF.h"
#include "makros.h"
#include "memory.h"
#include "readNonBlankLineFromInputStream.h"

#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

#include "GEOObjects.h"

#include "BoundaryCondition.h"
#include "BoundaryConditionIO.h"
#include "DistributionTools.h"
#include "GeoIO.h"
#include "mathlib.h"
#include "ProcessIO.h"
#include "rfmat_cp.h"
#include "rf_pcs.h"
#include "tools.h"


CBoundaryConditionNode::CBoundaryConditionNode()
{
	conditional = false;
	_bc = NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task: destructor
   Programing:
   08/2011 WW Implementation
**************************************************************************/
void CBoundaryConditionNode::Read(std::istream& is)
{
	is >> geo_node_number;
	is >> msh_node_number;
	is >> CurveIndex;
	is >> node_value;
	is >> std::ws;
}
/**************************************************************************
   FEMLib-Method:
   Task: destructor
   Programing:
   08/2011 WW Implementation
**************************************************************************/
void CBoundaryConditionNode::Write(std::ostream& os) const
{
	std::string deli = "  ";
	os << geo_node_number << deli;
	os << msh_node_number << deli;
	os << CurveIndex << deli;
	os << node_value << deli;
	os << "\n";
}

//==========================================================================
std::list<CBoundaryCondition*> bc_list;
std::vector<std::string> bc_db_head;
std::list<CBoundaryConditionsGroup*> bc_group_list;
std::vector<CBoundaryCondition*> bc_db_vector;

/**************************************************************************
   FEMLib-Method:
   Task: BC constructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::CBoundaryCondition()
    : GeoInfo(), geo_name(""), _curve_index(-1), dis_linear_f(NULL)
{
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
	// FCT
	conditional = false;
	time_dep_interpol = false;
	epsilon = 1e-9;         // NW
	time_contr_curve = -1;  // WX
	bcExcav = -1;           // WX
	MatGr = -1;             // WX
	is_MatGr_set = false;   // NW
	TimeInterpolation = 0;
	has_constrain = false;
	constrain_value = .0;
	constrain_operator = FiniteElement::INVALID_OPERATOR_TYPE;
	constrain_var_id = -1;
	ele_interpo_method = "AVERAGE";  // NW
}

// KR: Conversion from GUI-BC-object to CBoundaryCondition
CBoundaryCondition::CBoundaryCondition(const BoundaryCondition* bc)
    : ProcessInfo(bc->getProcessType(), bc->getProcessPrimaryVariable(), NULL),
      GeoInfo(bc->getGeoType(), bc->getGeoObj()),
      DistributionInfo(bc->getProcessDistributionType())
{
	setProcess(PCSGet(this->getProcessType()));
	this->geo_name = bc->getGeoName();
	const std::vector<size_t> dis_nodes = bc->getDisNodes();
	const std::vector<double> dis_values = bc->getDisValues();

	if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
	{
		this->geo_node_value = dis_values[0];
	}
	else if (this->getProcessDistributionType() == FiniteElement::LINEAR)
	{
		for (size_t i = 0; i < dis_values.size(); i++)
		{
			this->_PointsHaveDistribedBC.push_back(
			    static_cast<int>(dis_nodes[i]));
			this->_DistribedBC.push_back(dis_values[i]);
		}
	}
	else
		std::cout << "Error in CBoundaryCondition() - DistributionType \""
		          << FiniteElement::convertDisTypeToString(
		                 this->getProcessDistributionType())
		          << "\" currently not supported."
		          << "\n";

	is_MatGr_set = false;  // NW
}

/**************************************************************************
   FEMLib-Method:
   Task: BC deconstructor
   Programing:
   01/2004 OK Implementation
**************************************************************************/
CBoundaryCondition::~CBoundaryCondition()
{
	// DIS
	node_number_vector.clear();
	geo_node_number = -1;
	geo_node_value = 0.0;

	// WW
	if (dis_linear_f) delete dis_linear_f;
	dis_linear_f = NULL;
}

const std::string& CBoundaryCondition::getGeoName() const
{
	return geo_name;
}

/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   09/2004 OK POINTS method
   11/2004 MX stream string
**************************************************************************/
std::ios::pos_type CBoundaryCondition::Read(std::ifstream* bc_file,
                                            const GEOLIB::GEOObjects& geo_obj,
                                            const std::string& unique_fname,
                                            bool& valid)
{
	std::string line_string;
	bool new_keyword = false;
	std::ios::pos_type position;

	std::string sub_string, strbuff;
	int ibuff;     // pos,
	double dbuff;  // WW
	std::stringstream in;

	// Schleife ueber alle Phasen bzw. Komponenten
	while (!new_keyword)
	{
		position = bc_file->tellg();
		line_string = readNonBlankLineFromInputStream(*bc_file);
		if (line_string.size() < 1) break;
		if (line_string.find("#") != std::string::npos)
		{
			new_keyword = true;
			break;
		}

		if (line_string.find("$PCS_TYPE") != std::string::npos)
			if (!FileIO::ProcessIO::readProcessInfo(*bc_file, _pcs_type))
				valid = false;

		if (line_string.find("$PRIMARY_VARIABLE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			std::string tmp;
			in >> tmp;  // _pcs_pv_name;
			if (this->_pcs_type == FiniteElement::MASS_TRANSPORT)
			{
				// HS set the pointer to MCP based on component name.
				// a check whether this name is existing and unique.
				if (cp_name_2_idx.count(tmp) == 1)
				{
					setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
					setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
				}
				else
				{
					DisplayErrorMsg(
					    "Error: In reading BC file, the input component names "
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
		}

		// HS, this is new. later on we should stick to COMP_NAME,
		// PRIMARY_VARIABLE support will be removed.
		if (line_string.find("$COMP_NAME") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			std::string tmp;
			in >> tmp;  // _pcs_pv_name;
			if (this->_pcs_type == FiniteElement::MASS_TRANSPORT)
			{
				// HS set the pointer to MCP based on component name.
				// check whether this name is existing and unique.
				if (cp_name_2_idx.count(tmp) == 1)
				{
					setProcess(cp_vec[cp_name_2_idx[tmp]]->getProcess());
					setProcessPrimaryVariable(FiniteElement::CONCENTRATION);
				}
				else
				{
					DisplayErrorMsg(
					    "Error: In reading BC file, the input component names "
					    "are not found in MCP file!!!");
					exit(1);
				}
			}
			in.clear();
		}

		if (line_string.find("$GEO_TYPE") != std::string::npos)
			if (!FileIO::GeoIO::readGeoInfo(this, *bc_file, geo_name, geo_obj,
			                                unique_fname))
				valid = false;

		// PCH
		if (line_string.find("$DIS_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> line_string;  // sub_line
			_periodic = false;  // JOD

			// Source terms are assign to element nodes directly. 23.02.2009. WW
			if (line_string.find("DIRECT") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::DIRECT);
				in >> fname;
				fname = FilePath + fname;
				in.clear();
			}

			if (line_string.find("CONSTANT") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::CONSTANT);
				in >> geo_node_value;  // sub_line
				in.clear();
			}
			// If a linear function is given. 25.08.2011. WW
			if (line_string.find("FUNCTION") != std::string::npos)
			{
				setProcessDistributionType(FiniteElement::FUNCTION);
				in.clear();
				dis_linear_f = new LinearFunctionData(*bc_file);
			}
			if (line_string.find("LINEAR") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::LINEAR);
				// Distribuded. WW
				size_t nLBC;
				in >> nLBC;  // sub_line
				in.clear();

				for (size_t i = 0; i < nLBC; i++)
				{
					in.str(readNonBlankLineFromInputStream(*bc_file));
					in >> ibuff >> dbuff >> strbuff;
					in.clear();

					//           *bc_file>>ibuff>>dbuff;
					_PointsHaveDistribedBC.push_back(ibuff);
					_DistribedBC.push_back(dbuff);
					if (strbuff.size() > 0)
					{
						_PointsFCTNames.push_back(strbuff);
						time_dep_interpol = true;
					}
				}
				//        bc_file->ignore(MAX_ZEILE,'\n');
			}

			if (line_string.find("GRADIENT") !=
			    std::string::npos)  // 6/2012  JOD
			{
				this->setProcessDistributionType(FiniteElement::GRADIENT);
				in >> gradient_ref_depth;
				in >> gradient_ref_depth_value;
				in >> gradient_ref_depth_gradient;
				in.clear();
			}
			if (line_string.find("INITIAL") != std::string::npos)  // NW
			{
				this->setProcessDistributionType(FiniteElement::INITIAL);
			}
			if (line_string.find("ELEMENT") != std::string::npos)  // NW
			{
				this->setProcessDistributionType(FiniteElement::ELEMENT);
				in >> fname >> ele_interpo_method;
				fname = FilePath + fname;
				in.clear();
			}
		}

		// Time dependent function
		//..Time dependent curve ............................................
		if (line_string.find("$TIM_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> line_string;

			if (line_string.find("CURVE") != std::string::npos)
			{
				//				tim_type_name = "CURVE";
				this->setProcessDistributionType(FiniteElement::CONSTANT);
				in >> _curve_index;
				in.clear();

				//        pos1=pos2+1;
				//        sub_string = get_sub_string(buffer,"  ",pos1,&pos2);
				//		_curve_index = atoi(sub_string.c_str());
			}
			continue;
		}

		if (line_string.find("$FCT_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> fct_name;  // sub_line
			in.clear();
		}

		if (line_string.find("$MSH_TYPE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> sub_string;  // sub_line
			_msh_type_name = "NODE";
			if (sub_string.find("NODE") != std::string::npos)
			{
				in >> _msh_node_number;
				in.clear();
			}
		}

		if (line_string.find("$DIS_TYPE_CONDITION") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(
			    *bc_file));  // CONSTANT -21500.0
			in >> line_string;
			if (line_string.find("CONSTANT") != std::string::npos)
			{
				this->setProcessDistributionType(FiniteElement::CONSTANT);
				in >> geo_node_value;
				in.clear();
			}
			in.str(readNonBlankLineFromInputStream(
			    *bc_file));                // 0.0 IF HEAD > 0.04
			std::string pcs_pv_name_cond;  // 07/2010 TF temp string
			in >> node_value_cond >> line_string >> pcs_pv_name_cond >>
			    line_string >> condition;
			in.clear();
			in.str(readNonBlankLineFromInputStream(
			    *bc_file));  // PCS OVERLAND_FLOW
			std::string pcs_type_name_cond;
			in >> line_string >> pcs_type_name_cond;
			in.clear();
			conditional = true;
		}

		//....................................................................
		if (line_string.find("$CONSTRAIN") != std::string::npos)  // NW
		{
			// var > value
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> constrain_var_name >> line_string >> constrain_value;
			constrain_operator =
			    FiniteElement::convertComparisonOperatorType(line_string);
			in.clear();
			has_constrain = true;
			std::cout << "-> $CONSTRAIN is given for BC"
			          << "\n";
		}

		//....................................................................
		if (line_string.find("$EPSILON") != std::string::npos)  // NW
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> epsilon;
			in.clear();
		}
		//....................................................................
		// aktive state of the bc is time controlled  WX
		if (line_string.find("$TIME_CONTROLLED_ACTIVE") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> time_contr_curve;
			in.clear();
		}
		//....................................................................
		// bc for excated boundaries WX
		if (line_string.find("$EXCAVATION") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> bcExcav >> MatGr;
			in.clear();
		}
		//....................................................................
		// assignment of BC on mesh nodes connected to certain material elements
		// NW
		if (line_string.find("$MAT_ID") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
			in >> MatGr;
			is_MatGr_set = true;
			in.clear();
			continue;
		}
		if (line_string.find("$TIME_INTERPOLATION") != std::string::npos)
		{
			in.str(readNonBlankLineFromInputStream(*bc_file));
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
	}
	return position;
}

///**************************************************************************
// FEMLib-Method: CBoundaryCondition::Write
// 02/2004 OK Implementation
// 07/2007 OK LINEAR
// 10/2008 OK NOD
// 06/2009 OK MSH_TYPE off
// **************************************************************************/
// void CBoundaryCondition::Write(std::fstream* rfd_file) const
//{
//   //KEYWORD
//   *rfd_file << "#BOUNDARY_CONDITION" << "\n";
//   //--------------------------------------------------------------------
//   //NAME+NUMBER
//   *rfd_file << " $PCS_TYPE" << "\n";
//   *rfd_file << "  " << convertProcessTypeToString(getProcessType()) << "\n";
//   *rfd_file << " $PRIMARY_VARIABLE" << "\n";
//   *rfd_file << "  " <<
//   convertPrimaryVariableToString(this->getProcessPrimaryVariable()) << "\n";
//   //--------------------------------------------------------------------
//   //GEO_TYPE
//   *rfd_file << " $GEO_TYPE" << "\n";
//   *rfd_file << "  ";
//   *rfd_file << getGeoTypeAsString() << " " << geo_name << "\n";
//
//   //--------------------------------------------------------------------
//   /*OK4910
//    //MSH_TYPE
//    if(msh_node_number>0){
//    *rfd_file << " $MSH_TYPE" << endl;
//    *rfd_file << "  ";
//    *rfd_file << "NODE" << " " << msh_node_number << endl;
//    }
//    */
//   //--------------------------------------------------------------------
//   //DIS_TYPE
//   *rfd_file << " $DIS_TYPE" << "\n";
//   *rfd_file << "  ";
//   *rfd_file << convertDisTypeToString(this->getProcessDistributionType());
//   //switch (dis_type_name[0]) {
//   //case 'C': // Constant
//   if (this->getProcessDistributionType() == FiniteElement::CONSTANT)
//   {
//      *rfd_file << " " << geo_node_value;
//      *rfd_file << "\n";
//      //break;
//   }
//   //case 'L': // Linear
//   else if (this->getProcessDistributionType() == FiniteElement::LINEAR)
//   {
//      *rfd_file << " " << _PointsHaveDistribedBC.size() << "\n";
//      for (size_t i = 0; i < _PointsHaveDistribedBC.size(); i++)
//      {
//         *rfd_file << "  " << _PointsHaveDistribedBC[i] << " ";
//         *rfd_file << "  " << _DistribedBC[i] << "\n";
//      }
//      //break;
//   }
//
//   //FCT
//   if (fct_name.length() > 0)                     //OK4108
//   {
//      *rfd_file << " $FCT_TYPE" << "\n";
//      *rfd_file << "  ";
//      *rfd_file << fct_name << "\n";
//   }
//}

/**************************************************************************
   FEMLib-Method: CBoundaryCondition::Write
   Task: write function
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
void CBoundaryCondition::WriteTecplot(std::fstream* tec_file) const
{
	long i;
	CGLPolyline* m_polyline1 = NULL;
	CGLPolyline* m_polyline2 = NULL;
	// list<CGLPolyline*>::const_iterator p;
	std::vector<CGLPolyline*>::iterator p;
	Surface* m_surface = NULL;
	long no_points = 0;
	std::vector<CTriangle*> triangle_vector;

	*tec_file << "VARIABLES = X,Y,Z,V1"
	          << "\n";

	if (getGeoType() == GEOLIB::SURFACE)
	{
		m_surface = GEOGetSFCByName(geo_name);  // CC
		if (m_surface) switch (m_surface->type)
			{
				case 2:
					p = m_surface->polyline_of_surface_vector.begin();
					while (p != m_surface->polyline_of_surface_vector.end())
					{
						m_polyline1 = *p;
						++p;
						m_polyline2 = *p;
						break;
					}
					no_points = (long)m_polyline1->point_vector.size();
					/*
					   for(i=0;i<no_points-1;i++) {
					   m_triangle = new CTriangle;
					   m_triangle->x[0] = m_polyline1->point_vector[i]->x;
					   m_triangle->y[0] = m_polyline1->point_vector[i]->y;
					   m_triangle->z[0] = m_polyline1->point_vector[i]->z;
					   m_triangle->x[1] = m_polyline1->point_vector[i+1]->x;
					   m_triangle->y[1] = m_polyline1->point_vector[i+1]->y;
					   m_triangle->z[1] = m_polyline1->point_vector[i+1]->z;
					   m_triangle->x[2] = m_polyline2->point_vector[i+1]->x;
					   m_triangle->y[2] = m_polyline2->point_vector[i+1]->y;
					   m_triangle->z[2] = m_polyline2->point_vector[i+1]->z;
					   triangle_vector.push_back(m_triangle);
					   m_triangle = new CTriangle;
					   m_triangle->x[0] = m_polyline2->point_vector[i]->x;
					   m_triangle->y[0] = m_polyline2->point_vector[i]->y;
					   m_triangle->z[0] = m_polyline2->point_vector[i]->z;
					   m_triangle->x[1] = m_polyline2->point_vector[i+1]->x;
					   m_triangle->y[1] = m_polyline2->point_vector[i+1]->y;
					   m_triangle->z[1] = m_polyline2->point_vector[i+1]->z;
					   m_triangle->x[2] = m_polyline1->point_vector[i+1]->x;
					   m_triangle->y[2] = m_polyline1->point_vector[i+1]->y;
					   m_triangle->z[2] = m_polyline1->point_vector[i+1]->z;
					   triangle_vector.push_back(m_triangle);
					   }
					 */
					break;
			}
	}

	long no_nodes = 2 * no_points;
	// long no_elements = triangle_vector.size();
	long no_elements = 2 * (no_points - 1);
	// Write
	*tec_file << "ZONE T = " << geo_name << ", "
	          << "N = " << no_nodes << ", "
	          << "E = " << no_elements << ", "
	          << "F = FEPOINT"
	          << ", "
	          << "ET = TRIANGLE"
	          << "\n";
	if (m_polyline1)
		for (i = 0; i < no_points; i++)
			*tec_file << m_polyline1->point_vector[i]->x << " "
			          << m_polyline1->point_vector[i]->y << " "
			          << m_polyline1->point_vector[i]->z << " "
			          << geo_node_value << "\n";

	if (m_polyline2)
		for (i = 0; i < no_points; i++)
			*tec_file << m_polyline2->point_vector[i]->x << " "
			          << m_polyline2->point_vector[i]->y << " "
			          << m_polyline2->point_vector[i]->z << " "
			          << geo_node_value << "\n";

	for (i = 0; i < no_points - 1; i++)
		*tec_file << i + 1 << " " << i + 1 + 1 << " " << no_points + i + 1
		          << "\n";
	for (i = 0; i < no_points - 1; i++)
		*tec_file << no_points + i + 1 << " " << no_points + i + 1 + 1 << " "
		          << i + 1 + 1 << "\n";
}

/**************************************************************************
   FEMLib-Method:
   Task: BC read function
   Programing:
   01/2004 OK Implementation
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   05/2010 TF changes due to new GEOLIB integration, some improvements
**************************************************************************/
bool BCRead(std::string const& file_base_name,
            const GEOLIB::GEOObjects& geo_obj, const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string, bc_file_name;

	// File handling
	bc_file_name = file_base_name + BC_FILE_EXTENSION;

	std::ifstream bc_file(bc_file_name.data(), std::ios::in);
	if (!bc_file.good())
	{
		std::cout << "! Error in BCRead: No boundary conditions !"
		          << "\n";
		return false;
	}

	// Keyword loop
	ScreenMessage("BCRead ... \n");
	while (!bc_file.eof())
	{
		bc_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos)
		{
			ScreenMessage("-> done, read %d boundary conditions\n",
			              bc_list.size());
			return true;
		}
		if (line_string.find("#BOUNDARY_CONDITION") != std::string::npos)
		{
			CBoundaryCondition* bc(new CBoundaryCondition());
			bool valid(true);
			std::ios::pos_type position =
			    bc->Read(&bc_file, geo_obj, unique_name, valid);
			if (valid)
				bc_list.push_back(bc);
			else
				delete bc;
			bc_file.seekg(position, std::ios::beg);
		}  // keyword found
	}      // eof
	return true;
}

/**************************************************************************
   FEMLib-Method: BCWrite
   Task: master write function
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/

void BCWrite(std::string const& base_file_name)
{
	std::string sub_line;
	std::string line_string;

	// File handling
	std::string bc_file_name(base_file_name + BC_FILE_EXTENSION);
	std::fstream bc_file(bc_file_name.data(), std::ios::trunc | std::ios::out);
	bc_file.setf(std::ios::scientific, std::ios::floatfield);
	bc_file.precision(12);
	// OK string tec_file_name = base_file_name + ".tec";
	// OK fstream tec_file (tec_file_name.data(),ios::trunc|ios::out);
	// OK tec_file.setf(ios::scientific,ios::floatfield);
	// OK tec_file.precision(12);
	if (!bc_file.good()) return;
	bc_file.seekg(0L, std::ios::beg);  // rewind?
	bc_file << "GeoSys-BC: Boundary Conditions "
	           "------------------------------------------------\n";
	//========================================================================
	// BC list
	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
	while (p_bc != bc_list.end())
	{
		FileIO::BoundaryConditionIO::write(bc_file, *(*p_bc));
		++p_bc;
	}
	bc_file << "#STOP";
	bc_file.close();
	// OK tec_file.close();
}

/**************************************************************************
   FEMLib-Method:
   01/2004 OK Implementation
   07/2007 OK V2, global function
**************************************************************************/
// CBoundaryCondition* BCGet(const std::string &pcs_name, const std::string
// &geo_type_name,
//		const std::string &geo_name)
//{
//	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
//	while (p_bc != bc_list.end()) {
//		if (((*p_bc)->pcs_type_name.compare(pcs_name) == 0)
//				&& ((*p_bc)->geo_type_name.compare(geo_type_name) == 0)
//				&& ((*p_bc)->getGeoName().compare(geo_name) == 0))
//			return *p_bc;
//		++p_bc;
//	}
//	return NULL;
//}

/**************************************************************************
   GeoSys source term function:
   02/2009 WW Implementation
**************************************************************************/
inline void CBoundaryCondition::DirectAssign(long ShiftInNodeVector)
{
	std::string line_string;
	std::stringstream in;
	long n_index;
	double n_val;
	CRFProcess* m_pcs = NULL;
	CBoundaryConditionNode* m_node_value = NULL;

	m_pcs = PCSGet(convertProcessTypeToString(this->getProcessType()));

	//========================================================================
	// File handling
	std::ifstream d_file(fname.c_str(), std::ios::in);
	// if (!st_file.good()) return;

	if (!d_file.good())
	{
		std::cout
		    << "! Error in direct node source terms: Could not find file:!\n"
		    << fname << "\n";
		abort();
	}
	// Rewind the file
	d_file.clear();
	d_file.seekg(0L, std::ios::beg);
	//========================================================================
	while (!d_file.eof())
	{
		line_string = readNonBlankLineFromInputStream(d_file);
		if (line_string.find("#STOP") != std::string::npos) break;

		in.str(line_string);
		in >> n_index >> n_val;
		in.clear();
		//
		m_node_value = new CBoundaryConditionNode;
		m_node_value->conditional = false;
		m_node_value->msh_node_number = n_index + ShiftInNodeVector;
		m_node_value->geo_node_number = n_index;
		m_node_value->node_value = n_val;
		m_node_value->CurveIndex = _curve_index;
		m_pcs->bc_node.push_back(this);
		m_pcs->bc_node_value.push_back(m_node_value);
	}  // eof
}

/**************************************************************************
   GeoSys BC function:
   03/2009 WW Implementation
**************************************************************************/
inline void CBoundaryCondition::PatchAssign(long ShiftInNodeVector)
{
	std::string line_string;
	std::stringstream in;
	long n_index;
	std::vector<long> sfc_nodes;
	CBoundaryConditionNode* m_node_value = NULL;

	CRFProcess* pcs(PCSGet(convertProcessTypeToString(this->getProcessType())));
	Surface* surface(GEOGetSFCByName(geo_name));

	// File handling
	std::ifstream d_file(fname.c_str(), std::ios::in);

	if (!d_file.good())
	{
		std::cout
		    << "! Error in direct node source terms: Could not find file:!\n"
		    << fname << "\n";
		abort();
	}
	// Rewind the file
	d_file.clear();
	d_file.seekg(0L, std::ios::beg);

	while (!d_file.eof())
	{
		line_string = readNonBlankLineFromInputStream(d_file);
		if (line_string.find("#STOP") != std::string::npos) break;

		in.str(line_string);
		in >> n_index;
		in.clear();
		sfc_nodes.push_back(n_index);
	}

	if (surface) pcs->m_msh->GetNODOnSFC_PLY_XY(surface, sfc_nodes, true);

	for (size_t i = 0; i < sfc_nodes.size(); i++)
	{
		m_node_value = new CBoundaryConditionNode;
		m_node_value->conditional = false;
		n_index = sfc_nodes[i];
		m_node_value->msh_node_number = n_index + ShiftInNodeVector;
		m_node_value->geo_node_number = n_index;
		m_node_value->node_value = geo_node_value;
		m_node_value->CurveIndex = _curve_index;
		pcs->bc_node.push_back(this);
		pcs->bc_node_value.push_back(m_node_value);
	}  // eof
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2012 NW Implementation
**************************************************************************/
void CBoundaryCondition::SetByElementValues(long ShiftInNodeVector)
{
	// File handling
	std::ifstream d_file(fname.c_str(), std::ios::in);
	if (!d_file.is_open())
	{
		std::cout << "! Error in CBoundaryCondition::SetByElementValues(): "
		             "Could not find file " << fname << "\n";
		abort();
	}

	CRFProcess* pcs = this->getProcess();
	// MeshLib::CFEMesh* msh = pcs->m_msh;
	// read element values
	std::vector<long> bc_ele_ids;
	std::map<long, double> map_eleId_val;
	std::string line_string;
	std::stringstream in;
	size_t ele_id;
	double val;
	while (!d_file.eof())
	{
		line_string = GetLineFromFile1(&d_file);
		if (line_string.find("#STOP") != std::string::npos) break;

		in.str(line_string);
		in >> ele_id >> val;
		map_eleId_val[ele_id] = val;
		bc_ele_ids.push_back(ele_id);
		in.clear();
	}

	// get a list of nodes connecting elements
	std::vector<long> vec_nodes;
	pcs->m_msh->GetNODOnELE(bc_ele_ids, vec_nodes);

	// set BC
	CRFProcess* m_pcs =
	    PCSGet(convertProcessTypeToString(this->getProcessType()));
	const size_t n_nodes = vec_nodes.size();
	CBoundaryConditionNode* m_node_value = NULL;
	EleToNodeInterpolationMethod::type iterpo_type =
	    EleToNodeInterpolationMethod::VOLUME_WEIGHTED;
	if (ele_interpo_method.find("SHAPE") != std::string::npos)
	{
		iterpo_type = EleToNodeInterpolationMethod::GAUSS_EXTRAPOLATION;
	}
	for (size_t i = 0; i < n_nodes; i++)
	{
		long nod_id = vec_nodes[i];
		double v = getNodalValueFromElementValue(*pcs, map_eleId_val,
		                                         iterpo_type, nod_id);

		m_node_value = new CBoundaryConditionNode;
		m_node_value->conditional = false;
		m_node_value->msh_node_number = nod_id + ShiftInNodeVector;
		m_node_value->geo_node_number = nod_id;
		m_node_value->node_value = v;
		m_node_value->CurveIndex = _curve_index;
		m_pcs->bc_node.push_back(this);
		m_pcs->bc_node_value.push_back(m_node_value);
	}
}

CBoundaryConditionsGroup::CBoundaryConditionsGroup(void)
{
	msh_node_number_subst = -1;  //
	time_dep_bc = -1;
}

CBoundaryConditionsGroup::~CBoundaryConditionsGroup(void)
{
	/*
	   int group_vector_length = group_vector.size();
	   int i;
	   for(i=0;i<group_vector_length;i++)
	   group_vector.pop_back();
	 */
	//  group_vector.clear();
}

void setDistributionData(CBoundaryCondition* bc, DistributionData& distData)
{
	distData.dis_type = bc->getProcessDistributionType();
	distData.geo_obj = bc->getGeoObj();
	distData.geo_name = bc->getGeoName();
	distData.geo_type = bc->getGeoType();
	switch (bc->getProcessDistributionType())
	{
		case FiniteElement::CONSTANT:
			distData.dis_parameters.push_back(bc->getGeoNodeValue());
			break;
		case FiniteElement::LINEAR:
			distData._DistribedBC = bc->getDistribedBC();
			distData._PointsHaveDistribedBC = bc->getPointsWithDistribedBC();
			break;
		case FiniteElement::GRADIENT:
			distData.dis_parameters.push_back(bc->gradient_ref_depth);
			distData.dis_parameters.push_back(bc->gradient_ref_depth_value);
			distData.dis_parameters.push_back(bc->gradient_ref_depth_gradient);
			break;
		case FiniteElement::FUNCTION:
			distData.linear_f = bc->dis_linear_f;
			break;
		default:
			break;
	}
	distData.mesh_type_name = bc->getMeshTypeName();
	if (bc->getMeshTypeName() != "NODE")
	{
		distData.mesh_node_id = bc->getMeshNodeNumber();
	}
	if (bc->getExcav() > 0 || bc->isMatGrSet())
	{
		distData.mat_id = bc->getExcavMatGr();
	}
}

/**************************************************************************
   FEMLib-Method: CBoundaryCondition::Set
   Task: set boundary conditions
   Programing:
   02/2004 OK Implementation
   09/2004 WW Interpolation of piecewise linear BC
   02/2005 OK MSH types
   03/2005 OK MultiMSH, PNT
   08/2005 WW Changes due to the new geometry finite element.
   12/2005 OK FCT
   04/2006 WW New storage
   09/2006 WW Move linear interpolation to new MSH structure
   12/2007 WW Linear distributed BC in a surface
   10/2008 WW/CB SetTransientBCtoNodes
   last modification:
**************************************************************************/
void CBoundaryConditionsGroup::Set(CRFProcess* pcs, int ShiftInNodeVector,
                                   const std::string& this_pv_name)
{
	long* nodes = NULL;
	std::vector<long> nodes_vector;
	std::vector<double> node_value;
	CGLPolyline* m_polyline = NULL;
	CBoundaryConditionNode* m_node_value = NULL;
	group_name = _pcs_type_name;
	bool quadratic = false;

	if (!this_pv_name.empty()) _pcs_pv_name = this_pv_name;
	CFEMesh* m_msh = pcs->m_msh;
	// Tests //OK

	if (!m_msh)
		std::cout << "Warning in CBoundaryConditionsGroup::Set - no MSH data"
		          << "\n";
	// return;
	if (m_msh)  // WW
	{
		/// In case of P_U coupling monolithic scheme
		if (pcs->type == 41)  // WW Mono
		{
			// Deform
			if (_pcs_pv_name.find("DISPLACEMENT") != std::string::npos ||
			    _pcs_pv_name.find("VELOCITY_DM") != std::string::npos)
				quadratic = true;
			else
				quadratic = false;
		}
		else if (pcs->type == 4)
			quadratic = true;
		else
			quadratic = false;
		pcs->m_msh->SwitchOnQuadraticNodes(quadratic);
	}

	FiniteElement::PrimaryVariable primary_variable(
	    FiniteElement::convertPrimaryVariable(_pcs_pv_name));
	int pv_id = pcs->GetNodeValueIndex(_pcs_pv_name);
	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();

	clock_t start_time(clock());

	while (p_bc != bc_list.end())
	{
		CBoundaryCondition* bc(*p_bc);
		if (bc->time_dep_interpol)  // WW/CB
		{
			++p_bc;
			continue;
		}

		if (!(bc->getProcess() == pcs &&
		      bc->getProcessPrimaryVariable() == primary_variable))
		{
			++p_bc;
			continue;
		}

		//-- 23.02.3009. WW
		if (bc->getProcessDistributionType() == FiniteElement::DIRECT)
		{
			bc->DirectAssign(ShiftInNodeVector);
			++p_bc;
			continue;
		}
		//------------------------------------------------------------------
		if (bc->getProcessDistributionType() == FiniteElement::ELEMENT)  // NW
		{
			bc->SetByElementValues(ShiftInNodeVector);
			++p_bc;
			continue;
		}

		DistributionData distData;
		setDistributionData(bc, distData);
		//------------------------------------------------------------------
		// Detect mesh nodes for this BC
		//------------------------------------------------------------------
		nodes_vector.clear();
		getNodesOnDistribution(distData, *m_msh, nodes_vector);
		if (!nodes_vector.empty())
			ScreenMessage2d("-> %d nodes are found for this BC\n",
			                nodes_vector.size());
		//------------------------------------------------------------------
		// Calculate BC values
		//------------------------------------------------------------------
		node_value.resize(nodes_vector.size());
		setDistribution(distData, *m_msh, nodes_vector, node_value);
		if (bc->getProcessDistributionType() == FiniteElement::INITIAL)
		{
			for (size_t i = 0; i < nodes_vector.size(); i++)
				node_value[i] = pcs->GetNodeValue(nodes_vector[i], pv_id);
		}
		//------------------------------------------------------------------
		// create BC node
		//------------------------------------------------------------------
		const size_t nodes_vector_size(nodes_vector.size());
		for (size_t i(0); i < nodes_vector_size; i++)
		{
			m_node_value = new CBoundaryConditionNode();
			m_node_value->_bc = bc;
			m_node_value->msh_node_number = nodes_vector[i] + ShiftInNodeVector;
			m_node_value->geo_node_number = nodes_vector[i];
			m_node_value->node_value = node_value[i];
			m_node_value->CurveIndex = bc->getCurveIndex();
			m_node_value->pcs_pv_name = _pcs_pv_name;
			pcs->bc_node.push_back(bc);
			pcs->bc_node_value.push_back(m_node_value);
		}
		//------------------------------------------------------------------
		// FCT types //OK
		//------------------------------------------------------------------
		if (bc->fct_name.size() > 0)
		{
			for (size_t i = 0; i < pcs->bc_node_value.size(); i++)
			{
				pcs->bc_node_value[i]->fct_name = bc->fct_name;
				pcs->bc_node_value[i]->msh_node_number_subst =
				    msh_node_number_subst;
			}
		}

		++p_bc;
	}  // list

	clock_t end_time(clock());
	ScreenMessage2d("\t[BC] set BC took %g s\n",
	                (end_time - start_time) / (double)(CLOCKS_PER_SEC));

	start_time = clock();
	// SetTransientBCtoNodes  10/2008 WW/CB Implementation
	p_bc = bc_list.begin();
	while (p_bc != bc_list.end())
	{
		CBoundaryCondition* bc(*p_bc);
		if (!bc->time_dep_interpol)  // WW/CB
		{
			++p_bc;
			continue;
		}
		if (bc->getProcess() != pcs) continue;

		//................................................................
		if (bc->getGeoType() == GEOLIB::POLYLINE)
		{
			// CC
			m_polyline = GEOGetPLYByName(bc->geo_name);

			if (m_polyline)
			{
				// WW
				if (bc->getProcessDistributionType() == FiniteElement::LINEAR)
				{
					// TF
					double msh_min_edge_length = m_msh->getMinEdgeLength();
					m_msh->setMinEdgeLength(m_polyline->epsilon);
					std::vector<size_t> my_nodes_vector;
					GEOLIB::Polyline const* ply(
					    static_cast<GEOLIB::Polyline const*>(bc->getGeoObj()));
					m_msh->GetNODOnPLY(ply, my_nodes_vector);
					m_msh->setMinEdgeLength(msh_min_edge_length);

					nodes_vector.clear();
					for (size_t k(0); k < my_nodes_vector.size(); k++)
						nodes_vector.push_back(my_nodes_vector[k]);

					pcs->bc_transient_index.push_back(
					    (long)pcs->bc_node.size());
					for (size_t i = 0; i < nodes_vector.size(); i++)
					{
						m_node_value = new CBoundaryConditionNode();
						m_node_value->_bc = bc;
						m_node_value->msh_node_number = -1;
						m_node_value->msh_node_number =
						    nodes_vector[i] + ShiftInNodeVector;
						m_node_value->geo_node_number = nodes_vector[i];
						m_node_value->node_value = 0.0;
						// YD/WW
						m_node_value->pcs_pv_name = _pcs_pv_name;
						m_node_value->CurveIndex = bc->getCurveIndex();
						// WW
						pcs->bc_node.push_back(bc);
						// WW
						pcs->bc_node_value.push_back(m_node_value);
					}
					node_value.clear();
				}
				//................................................................
				// delete(values);
				Free(nodes);
			}  // if(m_ply)
		}
		//------------------------------------------------------------------
		// PCS
		++p_bc;
	}  // list
	   /* // Make the following as comment by WW
	      // Test
	      long no_bc = (long)pcs->bc_node_value.size();
	      if(no_bc<1)
	      cout << "Warning: no boundary conditions specified for " << pcs_type_name
	      << endl;
	    */
	end_time = clock();
	ScreenMessage2d("\t[BC] set transient BC took %g s.\n",
	                (end_time - start_time) / (double)(CLOCKS_PER_SEC));
}

/**************************************************************************
   FEMLib-Method: CBoundaryCondition::Get
   Task: set boundary conditions
   Programing:
   02/2004 OK Implementation
   last modification:
**************************************************************************/
CBoundaryConditionsGroup* CBoundaryConditionsGroup::Get(
    const std::string& pcs_name)
{
	CBoundaryConditionsGroup* m_bc_group = NULL;
	std::list<CBoundaryConditionsGroup*>::const_iterator p_bc_group =
	    bc_group_list.begin();
	while (p_bc_group != bc_group_list.end())
	{
		m_bc_group = *p_bc_group;
		if (m_bc_group->group_name.compare(pcs_name) == 0) return m_bc_group;
		++p_bc_group;
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
CBoundaryConditionsGroup* BCGetGroup(const std::string& pcs_type_name,
                                     const std::string& pcs_pv_name)
{
	std::list<CBoundaryConditionsGroup*>::const_iterator it =
	    bc_group_list.begin();
	while (it != bc_group_list.end())
	{
		if (((*it)->getProcessTypeName().compare(pcs_type_name) == 0) &&
		    ((*it)->getProcessPrimaryVariableName().compare(pcs_pv_name) == 0))
			return *it;
		++it;
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
void BCDelete()
{
	CBoundaryCondition* m_bc = NULL;
	std::list<CBoundaryCondition*>::const_iterator p = bc_list.begin();
	while (p != bc_list.end())
	{
		// bc_list.remove(*p);
		m_bc = *p;
		delete m_bc;
		++p;
	}
	bc_list.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void BCGroupDelete()
{
	CBoundaryConditionsGroup* m_bc_group = NULL;
	std::list<CBoundaryConditionsGroup*>::const_iterator p =
	    bc_group_list.begin();
	while (p != bc_group_list.end())
	{
		m_bc_group = *p;
		delete m_bc_group;
		// bc_group_list.remove(*p);
		++p;
	}
	bc_group_list.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void BCGroupDelete(const std::string& pcs_type_name,
                   const std::string& pcs_pv_name)
{
	std::list<CBoundaryConditionsGroup*>::iterator p = bc_group_list.begin();
	while (p != bc_group_list.end())
	{
		if (((*p)->getProcessTypeName().compare(pcs_type_name) == 0) &&
		    ((*p)->getProcessPrimaryVariableName().compare(pcs_pv_name) == 0))
		{
			delete *p;
			bc_group_list.erase(p);
			return;
		}
		++p;
	}
}

/**************************************************************************
   FEMLib-Method:
   07/2007 OK Implementation
**************************************************************************/
CBoundaryCondition* BCGet(const std::string& pcs_type_name)
{
	CBoundaryCondition* m_bc = NULL;

	std::list<CBoundaryCondition*>::const_iterator p_bc = bc_list.begin();
	while (p_bc != bc_list.end())
	{
		m_bc = *p_bc;
		if (m_bc->getProcessType() ==
		    FiniteElement::convertProcessType(pcs_type_name))
			return m_bc;
		++p_bc;
	}
	return NULL;
}

/**************************************************************************
   ROCKFLOW - Funktion:
   Programming:
   11/2007 WW Implementation
**************************************************************************/
void CBoundaryCondition::SurfaceInterpolation(
    CRFProcess* m_pcs,
    std::vector<long>& nodes_on_sfc,
    std::vector<double>& node_value_vector)
{
	long i, j, k, l;

	//----------------------------------------------------------------------
	// Interpolation of polygon values to nodes_on_sfc
	int nPointsPly = 0;
	double Area1, Area2;
	// NW. Default tolerance is 1e-9 but it can be changed in a BC file.
	double Tol = this->epsilon;
	bool Passed;
	double gC[3], p1[3], p2[3], vn[3], unit[3], NTri[3];
	//
	CGLPolyline* m_polyline = NULL;
	Surface* m_surface = NULL;
	m_surface = GEOGetSFCByName(geo_name);  // CC

	// list<CGLPolyline*>::const_iterator p =
	// m_surface->polyline_of_surface_list.begin();
	std::vector<CGLPolyline*>::iterator p =
	    m_surface->polyline_of_surface_vector.begin();

	for (j = 0; j < (long)nodes_on_sfc.size(); j++)
	{
		//      pn[0] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->X();
		//      pn[1] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->Y();
		//      pn[2] = m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->Z();
		double const* const pn(
		    m_pcs->m_msh->nod_vector[nodes_on_sfc[j]]->getData());
		node_value_vector[j] = 0.0;
		Passed = false;
		// nodes close to first polyline
		p = m_surface->polyline_of_surface_vector.begin();
		while (p != m_surface->polyline_of_surface_vector.end())
		{
			m_polyline = *p;
			// Gravity center of this polygon
			for (i = 0; i < 3; i++)
				gC[i] = 0.0;
			vn[2] = 0.0;
			nPointsPly = (int)m_polyline->point_vector.size();
			for (i = 0; i < nPointsPly; i++)
			{
				gC[0] += m_polyline->point_vector[i]->x;
				gC[1] += m_polyline->point_vector[i]->y;
				gC[2] += m_polyline->point_vector[i]->z;
				vn[2] += m_polyline->point_vector[i]->getPropert();
			}
			for (i = 0; i < 3; i++)
				gC[i] /= (double)nPointsPly;
			// BC value at center is an average of all point values of polygon
			vn[2] /= (double)nPointsPly;
			// Area of this polygon by the gravity center
			for (i = 0; i < nPointsPly; i++)
			{
				p1[0] = m_polyline->point_vector[i]->x;
				p1[1] = m_polyline->point_vector[i]->y;
				p1[2] = m_polyline->point_vector[i]->z;
				k = i + 1;
				if (i == nPointsPly - 1) k = 0;
				p2[0] = m_polyline->point_vector[k]->x;
				p2[1] = m_polyline->point_vector[k]->y;
				p2[2] = m_polyline->point_vector[k]->z;
				vn[0] = m_polyline->point_vector[i]->getPropert();
				vn[1] = m_polyline->point_vector[k]->getPropert();

				Area1 = fabs(ComputeDetTri(p1, gC, p2));

				Area2 = 0.0;
				// Check if pn is in the triangle by points (p1, gC, p2)
				Area2 = fabs(ComputeDetTri(p2, gC, pn));
				unit[0] = fabs(ComputeDetTri(gC, p1, pn));
				unit[1] = fabs(ComputeDetTri(p1, p2, pn));
				Area2 += unit[0] + unit[1];
				if (fabs(Area1 - Area2) < Tol)
				{
					// Interpolation within a triangle (p1,p2,gC)
					// Shape function
					for (l = 0; l < 2; l++)
						unit[l] /= Area1;
					ShapeFunctionTri(NTri, unit);
					for (l = 0; l < 3; l++)
						node_value_vector[j] += vn[l] * NTri[l];
					Passed = true;
					break;
				}
			}
			//
			p++;
			if (Passed) break;
		}  // while
	}      // j
}
