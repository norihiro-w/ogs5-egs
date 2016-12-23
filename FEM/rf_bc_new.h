/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   FEMLib - Object: Boundary Conditions
   Task: class implementation
   Programing:
   02/2004 OK Implementation
   last modified
**************************************************************************/
#ifndef rf_bc_new_INC
#define rf_bc_new_INC

//#include <list>
//#include <fstream>
//#include <string>
//#include <vector>

namespace FileIO
{
class BoundaryConditionIO;
}

// new GEOLIB
#include "DistributionInfo.h"  // TF
#include "GEOObjects.h"
#include "GeoInfo.h"             // TF
#include "LinearFunctionData.h"  // TF
#include "ProcessInfo.h"         // KR

// GEOLib
//#include "geo_ply.h"
// MSHLib
//#include "msh_lib.h"
// PCSLib
//#include "rf_pcs.h"
namespace MeshLib
{
class CFEMesh;
}

class BoundaryCondition;

class CBoundaryCondition : public ProcessInfo,
                           public GeoInfo,
                           public DistributionInfo
{
public:
	//	CBoundaryCondition(ProcessInfo const& process_info,
	//			GeoInfo const& geo_info,
	//			DistributionInfo const& distribution_info,
	//			);
	friend class CBoundaryConditionsGroup;
	friend class FileIO::BoundaryConditionIO;
	CBoundaryCondition();
	CBoundaryCondition(const BoundaryCondition* bc);

	~CBoundaryCondition();
	//      void Write(std::fstream*) const;
	void WriteTecplot(std::fstream*) const;

	/**
	 * reads a boundary condition from stream
	 * @param in input file stream for reading
	 * @param geo_obj pointer to the geometric object manager
	 * @param unique_fname the project name
	 * @param valid after return the variable valid contains the status of the
	 * object,
	 * valid is false if there occured an error while reading the data, else
	 * true
	 * @return the position in the stream after the boundary condition
	 */
	// TF
	std::ios::pos_type Read(std::ifstream* in,
	                        const GEOLIB::GEOObjects& geo_obj,
	                        const std::string& unique_fname, bool& valid);

	/**
	 * ToDo remove after transition to new GEOLIB - REMOVE CANDIDATE
	 * getGeoName returns a string used as id for geometric entity
	 * @return the value of attribute geo_name in case of
	 * geo_type_name == POLYLINE or geo_type_name = SURFACE
	 * If geo_type_name == POINT the id of the point is returned.
	 */
	const std::string& getGeoName() const;  // TF 05/2010

	int getCurveIndex() const  // TF 05/2010
	{
		return _curve_index;
	}

	bool isPeriodic() const  // TF 07/2010
	{
		return _periodic;
	}
	double getPeriodeTimeLength() const  // TF 07/2010
	{
		return _periode_time_length;
	}
	double getPeriodePhaseShift() const  // TF 07/2010
	{
		return _periode_phase_shift;
	}

	const std::vector<int>& getPointsWithDistribedBC() const
	{
		return _PointsHaveDistribedBC;
	}
	const std::vector<double>& getDistribedBC() const { return _DistribedBC; }
	std::vector<double>& getDistribedBC() { return _DistribedBC; }
	double getGeoNodeValue() const { return geo_node_value; }
	// KR

	const std::vector<std::string>& getPointsFCTNames() const
	{
		return _PointsFCTNames;
	}

	size_t getMeshNodeNumber() const { return _msh_node_number; }
	const std::string& getMeshTypeName() const { return _msh_type_name; }

	int getMaterialID() const
	{
		return MatGr;
	}
	int getTimeContrCurve()
	{
		return time_contr_curve;
	}

	bool isMatGrSet() const { return is_MatGr_set; }

private:
	std::vector<std::string> _PointsFCTNames;
	std::vector<int> _PointsHaveDistribedBC;
	std::vector<double> _DistribedBC;

	// GEO
	/**
	 * the id of the geometric object as string REMOVE CANDIDATE
	 */
	std::string geo_name;  // TF 05/2010
	std::string geo_type_name;

	std::string fname;               // 27.02.2009. WW
	int _curve_index;                // Time function index
	std::string ele_interpo_method;  // NW

	// DIS
	std::vector<long> node_number_vector;
	std::vector<double> node_value_vector;
	long geo_node_number;
	double geo_node_value;

	double _periode_phase_shift;  // JOD
	double _periode_time_length;  // JOD
	bool _periodic;               // JOD

public:
	double gradient_ref_depth;  // 6/2012 JOD
	double gradient_ref_depth_value;
	double gradient_ref_depth_gradient;
	std::string interpolation_method;
	int TimeInterpolation;

private:
	double node_value_cond;  // OK
	double condition;        // OK
	bool conditional;

	double epsilon;  // NW. temporally set here for surface interpolation
	bool time_dep_interpol;

public:
	bool has_constrain;
	std::string constrain_var_name;
	int constrain_var_id;
	double constrain_value;
	FiniteElement::ComparisonOperatorType constrain_operator;

private:
	// FCT
	std::string fct_name;

public:
	LinearFunctionData* dis_linear_f;  // 24.8.2011. WW
private:
	// WW
	void SurfaceInterpolation(CRFProcess* m_pcs,
	                          std::vector<long>& nodes_on_sfc,
	                          std::vector<double>& node_value_vector);
	inline void DirectAssign(long ShiftInNodeVector);
	// 19.03.2009. WW
	inline void PatchAssign(long ShiftInNodeVector);

	void SetByElementValues(long);  // NW
	// MSH
	long _msh_node_number;
	std::string _msh_type_name;  // OK4105

	int MatGr;
	// aktive state is controlled by time curve WX:01.2011
	int time_contr_curve;

	// NW
	bool is_MatGr_set;
};

class CBoundaryConditionNode  // OK raus
{
public:
	long geo_node_number;
	long msh_node_number;
	long msh_node_number_subst;  // WW

	double node_value;
	int CurveIndex;           // Time dependent function index
	std::string pcs_pv_name;  // YD/WW
	//
	std::string fct_name;  // WW
	// FCT
	int conditional;  // OK

	CBoundaryCondition* _bc;

	CBoundaryConditionNode();

	// 25.08.2011. WW
	void Read(std::istream& is);
	void Write(std::ostream& os) const;
};

class CBoundaryConditionsGroup
{
public:
	CBoundaryConditionsGroup(void);
	~CBoundaryConditionsGroup(void);

	void Set(CRFProcess* pcs, int ShiftInNodeVector,
	         const std::string& this_pv_name = "");
	CBoundaryConditionsGroup* Get(const std::string&);

	const std::string& getProcessTypeName() const { return _pcs_type_name; }
	void setProcessTypeName(const std::string& pcs_type_name)
	{
		_pcs_type_name = pcs_type_name;
	}
	const std::string& getProcessPrimaryVariableName() const
	{
		return _pcs_pv_name;
	}
	void setProcessPrimaryVariableName(const std::string& pcs_pv_name)
	{
		if (_pcs_type_name.find("MASS_TRANSPORT") == std::string::npos)
			_pcs_pv_name = pcs_pv_name;
		else
			_pcs_pv_name = "CONCENTRATION1";
	}
	long msh_node_number_subst;  // WW
	std::string fct_name;        // OK

	MeshLib::CFEMesh* m_msh;  // OK
	// WW std::vector<CBoundaryCondition*>bc_group_vector; //OK
	// WW double GetConditionalNODValue(int,CBoundaryCondition*); //OK
	int time_dep_bc;

private:
	std::string group_name;
	std::string _pcs_type_name;  // OK
	std::string _pcs_pv_name;    // OK
};

//========================================================================
#define BC_FILE_EXTENSION ".bc"
extern std::list<CBoundaryConditionsGroup*> bc_group_list;
extern CBoundaryConditionsGroup* BCGetGroup(const std::string& pcs_type_name,
                                            const std::string& pcs_pv_name);
extern std::list<CBoundaryCondition*> bc_list;

/**
 * read boundary conditions from file
 * @param file_base_name the base name of the file (without extension)
 * @param geo_obj the geometric object managing geometric entities
 * @param unique_name the (unique) name of the project
 * @return false, if the file can not opened, else true
 */
bool BCRead(std::string const& file_base_name,
            const GEOLIB::GEOObjects& geo_obj,
            const std::string& unique_name);

extern void BCWrite(std::string const&);
extern void BCDelete();
extern void BCGroupDelete(const std::string& pcs_type_name,
                          const std::string& pcs_pv_name);
extern void BCGroupDelete(void);
// OK
extern CBoundaryCondition* BCGet(const std::string&, const std::string&,
                                 const std::string&);
extern CBoundaryCondition* BCGet(std::string);  // OK

// ToDo
extern void ScalingDirichletBoundaryConditions(const double factor);
#endif
