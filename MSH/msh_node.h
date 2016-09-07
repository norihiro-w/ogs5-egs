/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef msh_node_INC
#define msh_node_INC

#include <cassert>
#include <string>
#include <vector>

#include "msh_core.h"

namespace MeshLib
{
// Class definition
class CNode : public CCore
{
public:
#ifndef OGS_ONLY_TH
	// The vector to store the representative element index.
	// This can be used to extract the norm of the plane that the element lies
	// on.
	// Establishing this vector is done in the Fluid Momentum
	// since this is bounded by velocity.
	std::vector<long> connected_planes;  // PCH

	bool crossroad = false;
#endif

	/** constructor */
	CNode(size_t Index);
	CNode(double x, double y, double z);
	CNode(double const* const coords);
	CNode(size_t Index, double x, double y, double z = 0.0);
	CNode(size_t Index, double const* coordinates);
	CNode(size_t Index, const CNode* parent);  // NW
	~CNode() {}

	// Operator
	void operator=(const CNode& n);
	bool operator==(const CNode& n);

	// Get functions
	/** \brief const access operator
	 *  The access to the coordinates is like the access to a field.
	 */
	double const& operator[](size_t idx) const
	{
		assert(idx <= 2);
		return coordinate[idx];
	}
	/** \brief access operator (see book Effektiv C++ programmieren - subsection
	 * 1.3.2 ).
	 * \sa const T& operator[] (size_t idx) const
	 */
	double& operator[](size_t idx)
	{
		return const_cast<double&>(static_cast<const CNode&>(*this)[idx]);
	}

	/** write mesh node coordinates into stream (used from operator<<)
	 * \param os a standard output stream
	 */
	virtual void write(std::ostream& os) const
	{
		os << coordinate[0] << " " << coordinate[1] << " " << coordinate[2]
		   << std::flush;
	}

	//	 double X() const { return coordinate[0]; };
	//	 double Y() const { return coordinate[1]; };
	//	 double Z() const { return coordinate[2]; };

	//	 void Coordinates(double *xyz) const
	//	 {
	//		for (size_t i = 0; i < 3; i++)
	//		   xyz[i] = coordinate[i];
	//	 }

	// 04/2010 TF
	/**
	 * returns the coordinates of the mesh node
	 * @return the coordinates of this mesh node
	 */
	inline double const* getData() const { return coordinate; }

	// Set functions
	void SetX(double argX) { coordinate[0] = argX; }

	void SetY(double argY) { coordinate[1] = argY; }

	void SetZ(double argZ) { coordinate[2] = argZ; }

	void SetCoordinates(const double* argCoord);

	int GetEquationIndex(bool quadratic = false) const { return !quadratic ? eqs_index : eqs_index_quadratic; }

	void SetEquationIndex(long eqIndex, bool quadratic = false)
	{
		if (!quadratic)
			eqs_index = eqIndex;
		else
			eqs_index_quadratic = eqIndex;
	}

	int GetGlobalIndex() const
	{
#ifdef USE_PETSC
		return global_index;
#else
		return GetIndex();
#endif
	}

	void SetGlobalIndex(long index)
	{
		global_index = index;
	}

	// Output
	void Write(std::ostream& os = std::cout) const;

	std::vector<size_t> const& getConnectedElementIDs() const
	{
		return _connected_elements;
	}

	std::vector<size_t>& getConnectedElementIDs()
	{
		return _connected_elements;
	}

	std::vector<size_t> const& getConnectedNodes() const
	{
		return _connected_nodes;
	}

	std::vector<size_t>& getConnectedNodes() { return _connected_nodes; }

	size_t getNumConnectedNodes() const { return _connected_nodes.size(); }

private:
	double coordinate[3];
	long global_index = -1;
	long eqs_index = -1;
	long eqs_index_quadratic = -1;
	std::vector<size_t> _connected_nodes;
	std::vector<size_t> _connected_elements;
};

std::ostream& operator<<(std::ostream& os, MeshLib::CNode const& node);

}  // namespace MeshLib

#endif
