/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**
 * \file msh_mesh.h
 */

/**************************************************************************
   MSHLib - Object:
   Task:
   Programing:
   08/2005 WW/OK Encapsulation from rf_ele_msh
   last modified
**************************************************************************/
#ifndef msh_mesh_INC
#define msh_mesh_INC

#include <algorithm>
#include <string>

#include "geo_lib.h"
#include "geo_sfc.h"
#include "geo_vol.h"
#include "GEOObjects.h"
#include "Point.h"
#include "Polyline.h"
#include "Surface.h"
#include "Grid.h"

#include "MSHEnums.h"
#include "MeshNodesAlongPolyline.h"

#include "OGSMeshIO.h"

#include "msh_elem.h"

class RandomWalk;
class CFluidMomentum;
class Problem;
class CRFProcess;

#ifdef NEW_EQS
namespace Math_Group
{
class SparseTable;
}
using Math_Group::SparseTable;
#endif

//------------------------------------------------------------------------
namespace MeshLib
{

#if defined(USE_PETSC)
typedef struct
{
	int index;
	double x;
	double y;
	double z;
} MeshNodes;

#endif

//------------------------------------------------------------------------
// Class definition
class CFEMesh
{
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
	size_t getNumberOfMeshLayers() const;  // TF

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
	double getMinEdgeLength() const;
	CNode* const* getNodes() const  // WW 05.2012.
	{
		return &nod_vector[0];
	}
	/**
	 * do not use this method REMOVE CANDIDATE
	 * @param val
	 */
	void setMinEdgeLength(double val)  // TF temporary
	{
		_min_edge_length = val;
	}

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

	friend class FileIO::OGSMeshIO;
	std::ios::pos_type GMSReadTIN(std::ifstream*);
	//
	void ConstructGrid();
	void GenerateHighOrderNodes();
/// For parallel computing. 03.2012. WW
#if defined(USE_PETSC)  // || defined(other parallel solver libs)
	void ConfigHighOrderElements();

	/*!
	 Fill data for subdomain mesh
	 @param header  : mesh header
	 @param s_nodes : mesh nodes
	 */
	void setSubdomainNodes(int* header, const MeshNodes* s_nodes);
	/*!
	 Fill data for subdomain mesh
	 @param header    : mesh header
	 @param elem_info : element information
	 @param inside    : indicator for elements that are inside the subdomain
	 */
	void setSubdomainElements(int* header, const int* elem_info,
	                          const bool inside);
	int calMaximumConnectedNodes();
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

	//
	//         void RenumberNodesForGlobalAssembly();
	// For number of nodes
	int GetMaxElementDim() const { return max_ele_dim; }
	void SwitchOnQuadraticNodes(bool quad) { useQuadratic = quad; }
	bool getOrder() const { return useQuadratic; }
	bool isAxisymmetry() const { return _axisymmetry; }
	// Get number of nodes
	// CMCD int to long
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

	void InitialNodesNumber()  // WW
	{
		NodesNumber_Quadratic = NodesNumber_Linear = nod_vector.size();
	}

	void setNumberOfNodesFromNodesVectorSize()
	{
		NodesNumber_Linear = nod_vector.size();
	}
	/// Free the memory occupied by edges
	void FreeEdgeMemory();  // 09.2012. WW

#ifndef NON_GEO  // WW
	             /**
	              * @{
	              * */

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
	/**
	 * \brief depreciated method
	 */
	void CreateLayerPolylines(CGLPolyline*);  // OK

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

	/** @} */  // close doxygen group
#endif         // WW #ifndef NON_GEO
	           //#ifndef NON_GEO                             //  WW
	           //         /**
	           //          * Store border nodes among different grids.
	           //          */
	           //         std::vector<GridsTopo*> grid_neighbors;
	           //         friend class ::Problem;
	           //#endif

	/**
	 * \brief Find all nodes connecting the given elements
	 */
	void GetNODOnELE(const std::vector<long>& vec_ele,
	                 std::vector<long>& vec_nod) const;
//....................................................................
#ifdef ObsoleteGUI  // WW 03.2012
	// QUAD->HEX
	void CreateHexELEFromQuad(int, double);
	// QUAD->LINE
	void CreateLineELEFromQuad(int, double, int);
	void SetActiveElements(std::vector<long>&);
	// MB
	void SetMSHPart(std::vector<long>&, long);
	bool NodeExists(size_t node);
	// LINE->LINE
	void AppendLineELE();

	// TRI->PRIS
	void CreatePriELEFromTri(int, double);

#endif

	// Coordinate system
	int GetCoordinateFlag() const { return coordinate_system; }
	void FillTransformMatrix();

	void RenumberNodesForGlobalAssembly();
	/**
	 * returns the vector storing pointers to all nodes (class CNode) of the
	 * mesh
	 * @return
	 */
	const std::vector<MeshLib::CNode*>& getNodeVector() const
	{
		return nod_vector;
	}

	// All nodes - should be private!!!
	std::vector<MeshLib::CNode*> nod_vector;

	// All edges
	std::vector<MeshLib::CEdge*> edge_vector;
	// All surface faces
	std::vector<MeshLib::CElem*> face_vector;
	// All surface normal
	std::vector<double*> face_normal;  // YD

	const std::vector<MeshLib::CElem*>& getElementVector() const
	{
		return ele_vector;
	}
	/**
	 * all elements stored in this vector
	 * */
	std::vector<MeshLib::CElem*> ele_vector;
	// Nodes in usage
	// To record eqs_index->global node index
	std::vector<long> Eqs2Global_NodeIndex;

	// OK
	void PrismRefine(int Layer, int subdivision);

	void ConnectedNodes(bool quadratic) const;
	// WW
	void ConnectedElements2Node(bool quadratic = false);
	// OK
	std::vector<std::string> mat_names_vector;
	void DefineMobileNodes(CRFProcess*);  // OK
	void FaceNormal();                    // YD
	void SetNODPatchAreas();              // OK4310
	void SetNetworkIntersectionNodes();   // OK4319->PCH
#ifndef NON_PROCESS                       // 05.03.2010 WW
#ifdef NEW_EQS                            // 1.11.2007 WW
	// Compute the graph of the sparse matrix related to this mesh. 1.11.2007 WW
	void CreateSparseTable();
	// Get the sparse graph   1.11.2007 WW
	SparseTable* GetSparseTable(bool quad = false) const
	{
		if (!quad)
			return sparse_graph;
		else
			return sparse_graph_H;
	}
#endif
#endif

	std::string pcs_name;
	std::string geo_name;       // MB
	std::string geo_type_name;  // OK10_4310

	size_t max_mmp_groups;  // OKCC
	int msh_max_dim;

	RandomWalk* PT;  // PCH

	CFluidMomentum* fm_pcs;  // by PCH

	/// Import MODFlow grid. 10.2009 WW
	void ImportMODFlowGrid(std::string const& fname);
	/// Convert raster cells into grid. 12.2009 WW
	void ConvertShapeCells(std::string const& fname);
#ifdef USE_HydSysMshGen
	// Be activated if it is still needed.
	// WW
	/// Generate Column-surface grid system for the modeling of
	/// surafce-subsuface coupled processes
	// 15.05.2009. WW // removed useless const TF
	void HydroSysMeshGenerator(std::string fname, int nlayers, double thickness,
	                           int mapping);

#endif
	void MarkInterface_mHM_Hydro_3D();  // 07.06.2010. WW
	void mHM2NeumannBC();
	/// Comptute int {f} a dA on top surface.
	void TopSurfaceIntegration();

	//#ifndef NDEBUG
	/**
	 * This is a getter method to access the private attribute _mesh_grid
	 * that is an instance of class Grid.
	 * @return
	 */
	GEOLIB::Grid<MeshLib::CNode> const* getGrid() const;
	//#endif

private:
	// private attributes
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
	size_t _n_msh_layer;  // OK
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
	double _min_edge_length;  // TK

	/**
	 * length for geometrical search algorithms, calculated in method
	 * \c computeSearchLength() and returned by method \c getSearchLength()
	 */
	double _search_length;

	// Process friends
	friend class ::CRFProcess;

	size_t NodesNumber_Linear;
	size_t NodesNumber_Quadratic;
/// For parallel computing. 03.2012. WW
#if defined(USE_PETSC)  // || defined(using other parallel scheme)
	// int n_sub_elements;
	int glb_NodesNumber_Linear;
	int glb_NodesNumber_Quadratic;
	int loc_NodesNumber_Linear;  // index of shadow nodes starts from this
	                             // number
	int loc_NodesNumber_Quadratic;
#endif
	bool useQuadratic;
	bool _axisymmetry;
	bool top_surface_checked;  // 07.06.2010.  WW

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
	int map_counter;  // 21.01.2009 WW

	bool mapping_check;  // 23.01.2009 WW
	/// Import shape file. 16.03.2026. WW
	size_t ncols, nrows;
	/// (x_0, y_0): coordinate of the left down corner
	double x0, y0, csize, ndata_v;
	std::vector<double> zz;  // Elevation
	inline void ReadShapeFile(std::string const& fname);
	// 03.2010. WW
	inline void Precipitation2NeumannBC(std::string const& fname,
	                                    std::string const& ofname,
	                                    double ratio = 0.8);
#ifndef NON_GEO  //  WW
	friend class ::Problem;
#endif  // #ifndef NON_GEO
        //
        // Sparse graph of this mesh. 1.11.2007 WW
#ifdef NEW_EQS
	SparseTable* sparse_graph;
	SparseTable* sparse_graph_H;  // For high order interpolation
#endif

#ifndef NON_GEO  // WW
	void CreateLineElementsFromMarkedEdges(
	    CFEMesh* m_msh_ply,
	    std::vector<long>& ele_vector_at_ply);  // NW
#endif                                          // #ifndef NON_GEO //WW
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
