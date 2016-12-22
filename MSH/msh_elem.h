/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef msh_elem_INC
#define msh_elem_INC

#include <iostream>
#include <string>

#include "MSHEnums.h"
#include "msh_edge.h"

namespace Math_Group
{
class Matrix;
}

namespace MeshLib
{

class CElem : public CCore
{
public:
	CElem();
	explicit CElem(size_t Index);
	// For Faces: Face, local face index
	CElem(size_t Index, CElem* onwer, int Face);
	CElem(size_t Index, CElem* m_ele_parent);

	CElem(CElem const& elem);

	/**
	 * constructor for a triangle element
	 * @param t have to be MshElemType::TRIANGLE
	 * @param node0 first node of triangle
	 * @param node1 second node of triangle
	 * @param node2 third node of triangle
	 * @param mat number of material
	 * @return a mesh element object
	 */
	CElem(MshElemType::type t, size_t node0, size_t node1, size_t node2,
	      int mat);

	/**
	 * constructor for a quad element
	 * @param t have to be MshElemType::QUAD
	 * @param node0 first node of quad
	 * @param node1 second node of quad
	 * @param node2 third node of quad
	 * @param node3 fourth node of quad
	 * @param mat number of material
	 * @return a mesh element object
	 */
	CElem(MshElemType::type t, size_t node0, size_t node1, size_t node2,
	      size_t node3, int mat);

	~CElem();

	//------------------------------------------------------------------

	/**
	 * \brief Sets the default properties for the given element type.
	 * \param t The element type of this element.
	 * \param isFace Signals if the element is initialised as a face of an owner
	 * element, some properties of the element might change in this case.
	 */
	void setElementProperties(MshElemType::type t, bool isFace = false);

	/**
	 * Method (re)sets the patch index. Patch index is used to assign a
	 * material to the element.
	 * @param pidx an identifier for material
	 */
	void setPatchIndex(size_t pidx) { patch_index = pidx; }

	/// Depricated method kept for backward compatibility. Use
	/// setElementProperties(MshElemType::type t) instead.
	void Config(MshElemType::type t) { setElementProperties(t); }

	// Geometry
	int GetDimension() const { return ele_dim; }
	double const* GetGravityCenter() const { return gravity_center; }
	double* ComputeGravityCenter();
	size_t GetPatchIndex() const { return patch_index; }  // MatGroup
	void SetPatchIndex(int value) { patch_index = value; }
	void ComputeVolume();
	void SetFluxArea(double fluxarea)
	{
		area = fluxarea;
	}
	double GetFluxArea()
	{
		return area;
	}

	double calcVolume() const;

	double GetVolume() const { return volume; }
	void SetVolume(double Vol) { volume = Vol; }
#ifndef OGS_ONLY_TH
	// This will be activated after m_tim->CheckCourant() is ready to work
	// kg44 21042010 activated
	void SetCourant(double Cour)
	{
		courant = Cour;
	}
	double GetCourant()
	{
		return courant;
	}
	void SetNeumann(double Neum)
	{
		neumann = Neum;
	}
	double GetNeumann()
	{
		return neumann;
	}
	double GetRepLength()
	{
		return representative_length;
	}
#endif
	//------------------------------------------------------------------
	// ID
	MshElemType::type GetElementType() const { return geo_type; }
	void SetElementType(MshElemType::type type) { geo_type = type; }

	void MarkingAll(bool makop);
	std::string GetName() const;
	//------------------------------------------------------------------
	// Nodes
	void GetNodeIndeces(Math_Group::vec<long>& node_index) const
	{
		for (int i = 0; i < (int)nodes_index.Size(); i++)
			node_index[i] = nodes_index[i];
	}

	void getNodeIndices(std::vector<size_t>& node_indices)
	{
		for (size_t i = 0; i < nodes_index.Size(); i++)
			node_indices.push_back(nodes_index[i]);
	}

	Math_Group::vec<long>& getNodeIndices() { return nodes_index; }
	Math_Group::vec<long> const& getNodeIndices() const { return nodes_index; }
	/**
	 * const access to the vector nodes_index
	 * @return a const reference to the vector
	 */
	const Math_Group::vec<long>& GetNodeIndeces() const { return nodes_index; }

	long GetNodeIndex(int index) const { return nodes_index[index]; }
	void SetNodeIndex(int index, long g_index) { nodes_index[index] = g_index; }

	void GetNodes(Math_Group::vec<CNode*>& ele_nodes)
	{
		for (size_t i = 0; i < nodes.Size(); i++)
			ele_nodes[i] = nodes[i];
	}

	void GetNodes(std::vector<CNode*>& nodesVec)
	{
		for (size_t i = 0; i < nodes.Size(); i++)
			nodesVec.push_back(nodes[i]);
	}

	CNode* GetNode(int index) { return nodes[index]; }

	CNode const* GetNode(int index) const { return nodes[index]; }

	void SetNodes(Math_Group::vec<CNode*>& ele_nodes, bool ReSize = false);

	void setNodes(std::vector<CNode*> const& ele_nodes);

	int GetNodesNumber_H() const { return nnodesHQ; }
	size_t GetNodesNumber(bool quad) const
	{
		if (quad)
			return (size_t)nnodesHQ;
		else
			return (size_t)nnodes;
	}
	int GetVertexNumber() const { return nnodes; }
	void SetNodesNumber(int ivalue) { nnodes = ivalue; }
	CElem* GetOwner() { return owner; }

	// Initialize topological properties
	void InitializeMembers();
	//------------------------------------------------------------------
	// Edges
	void GetEdges(Math_Group::vec<CEdge*>& ele_edges)
	{
		for (size_t i = 0; i < nedges; i++)
			ele_edges[i] = edges[i];
	}
	CEdge* GetEdge(int index) { return edges[index]; }
	void SetEdges(Math_Group::vec<CEdge*>& ele_edges)
	{
		for (size_t i = 0; i < nedges; i++)
			edges[i] = ele_edges[i];
	}
	void SetEdgesOrientation(Math_Group::vec<int>& ori_edg)
	{
		for (size_t i = 0; i < nedges; i++)
			edges_orientation[i] = ori_edg[i];
	}
	void FreeEdgeMemory()
	{
		edges.resize(0);
		edges_orientation.resize(0);
	}
	void GetLocalIndicesOfEdgeNodes(const int Edge, int* EdgeNodes);
	size_t GetEdgesNumber() const { return nedges; }
	//------------------------------------------------------------------
	// Faces
	size_t GetFacesNumber() const { return nfaces; }

	void SetFace();
	void SetFace(CElem* onwer, const int Face);
	int GetSurfaceFacesNumber() const { return no_faces_on_surface; }
	int GetLocalFaceIndex() const { return face_index; }
	int GetElementFaceNodes(int Face, int* FacesNode);
	//------------------------------------------------------------------

	// Neighbors
	void SetNeighbors(Math_Group::vec<CElem*>& ele_neighbors)
	{
		for (size_t i = 0; i < nfaces; i++)
			neighbors[i] = ele_neighbors[i];
	}
	void SetNeighbor(const int LocalIndex, CElem* ele_neighbor)
	{
		neighbors[LocalIndex] = ele_neighbor;
	}
	void GetNeighbors(Math_Group::vec<CElem*>& ele_neighbors)
	{
		for (size_t i = 0; i < nfaces; i++)
			ele_neighbors[i] = neighbors[i];
	}
	CElem* GetNeighbor(int index) { return neighbors[index]; }

	//------------------------------------------------------------------
	// Coordinates transform
	void FillTransformMatrix();
	void FillTransformMatrix(int noneed);
	double getTransformTensor(int idx);
	Math_Group::Matrix const* getTransformTensor() const {return transform_tensor; }
#ifndef OGS_ONLY_TH
	void AllocateMeomoryforAngle()
	{
		if (!angle) angle = new double[3];
	}
	double GetAngle(int i) const { return angle[i]; }
	void SetAngle(int i, double value) { angle[i] = value; }
#endif
	//------------------------------------------------------------------
	// I/O
	void Read(std::istream& is = std::cin, int fileType = 0);

	void WriteIndex(std::ostream& os = std::cout) const;
	void WriteIndex_TEC(std::ostream& os = std::cout) const;
	void WriteAll(std::ostream& os = std::cout) const;
	void WriteNeighbors(std::ostream& os = std::cout) const;

	Math_Group::Vector mat_vector;

	int selected;

	void FaceNormal(const int index0, const int index1, double*);
	double* normal_vector;
	void SetNormalVector();

#ifndef OGS_ONLY_TH
	// For stability calculations
	double representative_length;
	double courant;
	double neumann;

	int GetExcavState() { return excavated; }
	void SetExcavState(const int ExcavState)
	{
		excavated = ExcavState;
	}
#endif
#if defined(USE_PETSC)
	bool isOverlapped() const { return g_index != NULL; }
	int* getGhostNodeIndices() { return g_index; }
	long GetGlobalIndex() const { return global_index; }
#else
	long GetGlobalIndex() const { return this->GetIndex(); }
#endif

private:
	MshElemType::type geo_type;
	CElem* owner;
	int ele_dim;

	int nnodes;
	int nnodesHQ;
	Math_Group::vec<CNode*> nodes;
	Math_Group::vec<long> nodes_index;
#if defined(USE_PETSC)
	int* g_index;
	long global_index = -1;
#endif

	size_t nedges;
	Math_Group::vec<CEdge*> edges;
	Math_Group::vec<int> edges_orientation;

	size_t nfaces;
	int no_faces_on_surface;
	int face_index;  // Local face index for the instance for face
	double volume;
	double gravity_center[3];
	size_t patch_index;
	double area;
	Math_Group::Matrix* transform_tensor;
	Math_Group::vec<CElem*> neighbors;
#ifndef OGS_ONLY_TH
	double* angle;
	int excavated;
#endif

	int GetElementFaces1D(int* FaceNode);
	int GetElementFacesTri(int Face, int* FaceNode);
	int GetElementFacesQuad(int Face, int* FaceNode);
	int GetElementFacesHex(int Face, int* FaceNode);
	int GetElementFacesTet(int Face, int* FaceNode);
	int GetElementFacesPri(int Face, int* FaceNode);
	int GetElementFacesPyramid(int Face, int* FaceNode);

	//-- Friends
	friend class CFEMesh;
};

}  // namespace MeshLib

#endif
