/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef msh_mesh_INC
#define msh_mesh_INC

#include <algorithm>
#include <string>

#include "GEOObjects.h"
#include "geo_lib.h"
#include "geo_sfc.h"
#include "geo_vol.h"
#include "Grid.h"
#include "Point.h"
#include "Polyline.h"
#include "Surface.h"

#include "MeshNodesAlongPolyline.h"
#include "MSHEnums.h"
#include "msh_elem.h"
#include "OGSMeshIO.h"


namespace MeshLib
{
#ifndef NON_GEO
/*!
   Class to handle topologic relationship among grids
   Designed by WW
   First version:    13.05.2009
 */
class GridsTopo
{
private:
	long* local_indices;  // of border nodes in local node array of a grid
	// double *comm_data;   // data for communication
	std::string neighbor_name;
	long bnodes;
	friend class CFEMesh;

public:
	GridsTopo(std::istream& in, std::string sec_name);
	std::string getNeighbor_Name() const { return neighbor_name; }
	long* getBorderNodeIndicies() const { return local_indices; }
	long getBorderNodeNumber() const { return bnodes; }

	// void Write(std::ostream &os=cout);
	~GridsTopo();
};
#endif  //#ifndef NON_GEO

#if defined(USE_PETSC)
typedef struct
{
	int global_id = -1;
	int dom_id = -1;
	int eqs_id = -1;
	int eqs_id_Q = -1;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
} MeshNodes;

struct MeshHeader
{
	int n_dom_nodes_Q;
	int n_dom_nodes_L;
	int n_inner_elements;
	int n_ghost_elements;
	int n_inner_nodes_L;
	int n_inner_nodes_Q;
	int n_global_nodes_L;
	int n_global_nodes_Q;
	int n_global_elements;
	int n_element_integers;
	int n_ghost_element_integers;

	void set(int *array)
	{
		int i = 0;
		n_dom_nodes_Q = array[i++];
		n_dom_nodes_L = array[i++];
		n_inner_elements = array[i++];
		n_ghost_elements = array[i++];
		n_inner_nodes_L = array[i++];
		n_inner_nodes_Q = array[i++];
		n_global_nodes_L = array[i++];
		n_global_nodes_Q = array[i++];
		n_global_elements = array[i++];
		n_element_integers = array[i++];
		n_ghost_element_integers = array[i++];
	}
};

#endif

//------------------------------------------------------------------------
// Class definition
class CFEMesh
{
	friend class FileIO::OGSMeshIO;
public:
	/// Constructor using geometric information.
	CFEMesh(GEOLIB::GEOObjects* geo_obj = NULL,
	        std::string* unique_name = NULL);

	/// Copy-Constructor.
	/// Note that this is not a real copy-constructor. It copies only nodes and
	/// elements and calls ConstructGrid() afterwards.
	CFEMesh(CFEMesh const& mesh);

	/// Destructor
	~CFEMesh();

#ifndef NON_GEO
	GEOLIB::GEOObjects* getGEOObjects() const { return _geo_obj; }

	std::string* getProjectName() const { return _geo_name; }
#endif  //#ifndef NON_GEO

	/**
	 * sets the value for element type
	 * @param ele_type
	 */
	void setElementType(MshElemType::type ele_type);

	/**
	 * set the number of mesh layer
	 * @param n_msh_layer
	 */
	void setNumberOfMeshLayers(size_t n_msh_layer);

	/**
	 * returns the number of mesh layers
	 * @return the number of mesh layers
	 */
	size_t getNumberOfMeshLayers() const;

	/**
	 *
	 * @return
	 */
	bool hasCrossSection() const;

	size_t getNumberOfLines() const;
	size_t getNumberOfQuads() const;
	size_t getNumberOfHexs() const;
	size_t getNumberOfTris() const;
	size_t getNumberOfTets() const;
	size_t getNumberOfPrisms() const;
	size_t getNumberOfPyramids() const;

	CNode* const* getNodes() const
	{
		return &nod_vector[0];
	}

	/**
	 * do not use this method REMOVE CANDIDATE
	 * @param val
	 */
	void setMinEdgeLength(double val)
	{
		_min_edge_length = val;
	}

	double getMinEdgeLength() const;

	/**
	 * Access method to the search length for geometric search algorithms
	 * computed
	 * by method \c computeSearchLength().
	 * @return the search length
	 */
	double getSearchLength() const;

	/**
	 * @brief Compute the search length for geometric search algorithms.
	 *
	 * Let \f$\mu\f$ the mean value of all edge length and \f$s\f$ the
	 * standard deviation. The search length \f$\ell\f$ is computed by the
	 *formula
	 * \f$\ell = \mu - c \cdot s.\f$
	 * @param c (input) scaling constant, default value = 2
	 */
	void computeSearchLength(double c = 2);

	bool Read(std::ifstream*);

	void ConstructGrid();
	void GenerateHighOrderNodes();

#if defined(USE_PETSC)
	void ConfigHighOrderElements();

	/*!
	 Fill data for subdomain mesh
	 @param header  : mesh header
	 @param s_nodes : mesh nodes
	 */
	void setSubdomainNodes(MeshHeader const& header, const MeshNodes* s_nodes);
	/*!
	 Fill data for subdomain mesh
	 @param header    : mesh header
	 @param elem_info : element information
	 @param inside    : indicator for elements that are inside the subdomain
	 */
	void setSubdomainElements(MeshHeader const& header, const int* elem_info,
	                          const bool inside);

	int calMaximumConnectedLocalNodes(bool quadratic, std::vector<int> &d_nnz);
	int calMaximumConnectedGhostNodes(bool quadratic, std::vector<int> &o_nnz);
	/// Get number of nodes of the entire mesh
	int getNumNodesGlobal() const { return glb_NodesNumber_Linear; }
	/// Get number of nodes of the entire mesh of quadratic elements
	int getNumNodesGlobal_Q() const { return glb_NodesNumber_Quadratic; }
	/// Get number of nodes of the subdomain mesh
	int getNumNodesLocal() const { return loc_NodesNumber_Linear; }
	/// Get number of nodes of the subdomain mesh of quadratic elements
	int getNumNodesLocal_Q() const { return loc_NodesNumber_Quadratic; }
	/// is the node owned by this domain
	bool isNodeLocal(long node_id) const
	{
		// node numbering is assumed in the following order
		// 1. local linear nodes
		// 2. ghost linear nodes
		// 3. local quad nodes
		// 4. ghost quad nodes
		if (isHigherOrderNode(node_id))
		{
			const size_t loc_nnodes_quad_only =
			    loc_NodesNumber_Quadratic - loc_NodesNumber_Linear;
			return static_cast<size_t>(node_id) <
			       (NodesNumber_Linear + loc_nnodes_quad_only);
		}
		else
		{
			return node_id < loc_NodesNumber_Linear;
		}
		return true;
	}
#else
	/// is the node owned by this domain
	bool isNodeLocal(long) const { return true; }
#endif
	/// is the node used for higher-order elements
	bool isHigherOrderNode(long node_id) const
	{
		return static_cast<size_t>(node_id) >= NodesNumber_Linear;
	}

	bool hasHigherOrderNodes() const { return NodesNumber_Linear != NodesNumber_Quadratic; }

	int GetMaxElementDim() const { return max_ele_dim; }
	void SwitchOnQuadraticNodes(bool quad) { useQuadratic = quad; }
	bool getOrder() const { return useQuadratic; }
	bool isAxisymmetry() const { return _axisymmetry; }
	void isAxisymmetry(bool f) { _axisymmetry = f; }

	size_t GetNodesNumber(const bool quadr) const
	{
		if (quadr)
			return NodesNumber_Quadratic;
		else
			return NodesNumber_Linear;
	}

	size_t NodesInUsage() const
	{
		if (useQuadratic)
			return NodesNumber_Quadratic;
		else
			return NodesNumber_Linear;
	}

	/// Free the memory occupied by edges
	void FreeEdgeMemory();  // 09.2012. WW

	/**
	 * \brief depreciated method - uses old surface class
	 */
	void GetNODOnPLY_XY(CGLPolyline* m_ply, std::vector<long>& msh_nod_vector);
	//	/**
	//	 * \brief depreciated method - uses old surface class
	//	 */
	void GetNODOnSFC_Vertical(Surface* m_sfc,
	                          std::vector<long>& msh_nod_vector);

	/**
	 * \brief depreciated method
	 */
	void CreateLineELEFromPLY(CGLPolyline*);

	// GEO-SFC
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC(Surface*, std::vector<long>&);
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_PLY(Surface const*, std::vector<long>&) const;
	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_PLY_XY(
	    Surface* m_sfc, std::vector<long>& msh_nod_vector,
	    bool givenNodesOnSurface = false);  // givenNodeOnSurface by WW
	                                        /**
	                                         * \brief depreciated method
	                                         */
	// 02.2009/OK
	void GetNODOnSFC_PLY_Z(Surface*, std::vector<long>&);

	/**
	 * \brief depreciated method
	 */
	void GetNODOnSFC_TIN(Surface*, std::vector<long>&);
	/**
	 * \brief deprecated method
	 */
	void GetNodesOnCylindricalSurface(Surface* m_sfc,
	                                  std::vector<long>& NodesS);
	/**
	 * \brief depreciated method
	 */
	void CreateQuadELEFromSFC(Surface*);

	/**
	 * GetNODOnPNT searchs the nearest node to the geometric point
	 * */
	long GetNODOnPNT(const GEOLIB::Point* const pnt) const;
	/**
	 * GetNearestELEOnPNT searchs the nearest element (gravity center)
	 * to the geometric point
	 * */
	long GetNearestELEOnPNT(const GEOLIB::Point* const pnt) const;

	/**
	 * GetNODOnPLY search the nearest nodes along the Polyline object
	 * @param ply constant pointer to a constant Polyline object
	 * @param msh_nod_vector the mesh node indices are saved in this vector
	 * */
	void GetNODOnPLY(const GEOLIB::Polyline* const ply,
	                 std::vector<size_t>& msh_nod_vector);

	/**
	 *
	 * @param ply
	 * @return
	 */
	const MeshNodesAlongPolyline& GetMeshNodesAlongPolyline(
	    const GEOLIB::Polyline* const ply);

	/**
	 *
	 * @param ply
	 * @param points
	 */
	void getPointsForInterpolationAlongPolyline(
	    const GEOLIB::Polyline* const ply, std::vector<double>& points);

	/**
	 * GetNODOnPLY search the nearest nodes to the Polyline
	 * */
	void GetNODOnPLY(const GEOLIB::Polyline* const ply,
	                 std::vector<long>& msh_nod_vector);

	/**
	 * \brief gives the indices of CElement elements, which have an edge
	 * in common with the polyline.
	 */
	void GetELEOnPLY(const GEOLIB::Polyline*, std::vector<size_t>&,
	                 bool With1DElements);  // 11/2011 BG

	/**
	 * \brief gives the indices of nodes, which are contained in the surface
	 */
	void GetNODOnSFC(const GEOLIB::Surface* sfc,
	                 std::vector<size_t>& msh_nod_vector) const;


	/**
	 * \brief Find all nodes connecting the given elements
	 */
	void GetNODOnELE(const std::vector<long>& vec_ele,
	                 std::vector<long>& vec_nod) const;

	int GetCoordinateFlag() const { return coordinate_system; }

	void FillTransformMatrix();

	/**
	 * returns the vector storing pointers to all nodes (class CNode) of the
	 * mesh
	 * @return
	 */
	const std::vector<MeshLib::CNode*>& getNodeVector() const
	{
		return nod_vector;
	}

	std::vector<MeshLib::CNode*> nod_vector;
	std::vector<MeshLib::CEdge*> edge_vector;
	std::vector<MeshLib::CElem*> face_vector;
	std::vector<double*> face_normal;
	std::vector<MeshLib::CElem*> ele_vector;

	const std::vector<MeshLib::CElem*>& getElementVector() const
	{
		return ele_vector;
	}

	// To record eqs_index->global node index
	std::vector<long> Eqs2Global_NodeIndex;
	std::vector<long> Eqs2Global_NodeIndex_Q;

	void ConnectedNodes(bool quadratic);
	void ConnectedElements2Node(bool quadratic = false);
	void FaceNormal();

	std::vector<std::string> mat_names_vector;

	std::string pcs_name;
	std::string geo_name;
	std::string geo_type_name;

	size_t max_mmp_groups;
	int msh_max_dim;


	/**
	 * This is a getter method to access the private attribute _mesh_grid
	 * that is an instance of class Grid.
	 * @return
	 */
	GEOLIB::Grid<MeshLib::CNode> const* getGrid() const;

private:
	/**
	 * reference to object of class GEOObject, that manages the geometry data
	 */
	GEOLIB::GEOObjects* _geo_obj;
	/**
	 * identifier for geometry
	 */
	std::string* _geo_name;
	std::vector<MeshLib::MeshNodesAlongPolyline> _mesh_nodes_along_polylines;

	MshElemType::type _ele_type;
	size_t _n_msh_layer;
	bool _cross_section;
	size_t _msh_n_lines;
	size_t _msh_n_quads;
	size_t _msh_n_hexs;
	size_t _msh_n_tris;
	size_t _msh_n_tets;
	size_t _msh_n_prisms;
	size_t _msh_n_pyras;

	/**
	 * method initializes the minimal edge length that is used in search
	 * algorithms
	 */
	void computeMinEdgeLength();
	double _min_edge_length;

	/**
	 * length for geometrical search algorithms, calculated in method
	 * \c computeSearchLength() and returned by method \c getSearchLength()
	 */
	double _search_length;

	size_t NodesNumber_Linear;
	size_t NodesNumber_Quadratic;
#if defined(USE_PETSC)
	int glb_ElementsNumber;
	int glb_NodesNumber_Linear;
	int glb_NodesNumber_Quadratic;
	int loc_NodesNumber_Linear;  // index of shadow nodes starts from this
	                             // number
	int loc_NodesNumber_Quadratic;
#endif
	bool useQuadratic;
	bool _axisymmetry;
	bool top_surface_checked;

	// Coordinate indicator
	// 1:  X component only
	// 12: Y component only
	// 13: Z component only
	// 2:  X, Y component
	// 23:  X, Z component
	// 3:  X, Y, Z component
	int coordinate_system;
	bool has_multi_dim_ele;
	int max_ele_dim;
	int map_counter;


	/// Store border nodes among different grids.
public:
	std::vector<GridsTopo*> grid_neighbors;
private:
	void CreateLineElementsFromMarkedEdges(
	    CFEMesh* m_msh_ply,
	    std::vector<long>& ele_vector_at_ply);
public:
	void constructMeshGrid();

private:
	GEOLIB::Grid<MeshLib::CNode>* _mesh_grid;

#ifdef USE_PETSC
public:
	std::size_t getLocalNodeID(std::size_t global_id) const
	{
		CompareGlobalNodeID comp(global_id);
		std::vector<std::pair<std::size_t, std::size_t> >::const_iterator itr =
		    std::find_if(_global_local_nodeids.begin(),
		                 _global_local_nodeids.end(), comp);
		if (itr != _global_local_nodeids.end()) return itr->second;
		return -1;
	}

private:
	std::vector<std::pair<std::size_t, std::size_t> > _global_local_nodeids;
	std::vector<std::size_t > _vec_node_dom_ids;

	class CompareGlobalNodeID
	{
	public:
		CompareGlobalNodeID(std::size_t gid) : _global_id(gid) {}
		bool operator()(std::pair<std::size_t, std::size_t> const& element)
		{
			return element.first == _global_id;
		}

	private:
		std::size_t _global_id;
	};
#endif
};

}  // namespace MeshLib

#endif
