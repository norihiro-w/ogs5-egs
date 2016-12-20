/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ElementMatrix.h"

namespace FiniteElement
{

/**************************************************************************
   ElementMatrix::AllocateMemory

   Arguments:
    const int EleIndex:  Element index,
    int type          : type
                         used to get element type.
                         type = 0, for Possion type
                         type = 1, for Possion equation with deformation
coupling
                         type = 2, for Navier equation
                         type = 3, for Navier equation with pressure coupling
   type = 4, Monlithic scheme of u-p coupling
   type = 5, Mass Transport
   default = 0.

   Programmaenderungen:
   01/2005     WW

**************************************************************************/
void ElementMatrix::AllocateMemory(MeshLib::CElem* ele, int type)
{
	int nnodes, nnodesHQ, dim, size;
	size = 0;
	// The following two lines will be updated when new FEMGEO is ready
	nnodes = ele->GetVertexNumber();
	nnodesHQ = ele->GetNodesNumber_H();
	dim = ele->GetDimension();
	switch (type)
	{
		case 0:  // H || T Process
			Mass = new Matrix(nnodes, nnodes);
			//        Laplace = new SymMatrix(nnodes);
			Laplace = new Matrix(nnodes, nnodes);
			RHS = new Vec(nnodes);
			break;
		case 1:  // HM Partioned scheme, Flow
			Mass = new Matrix(nnodes, nnodes);
			//        Laplace = new SymMatrix(nnodes);
			Laplace = new Matrix(nnodes, nnodes);
			RHS = new Vec(nnodes);
			CouplingB = new Matrix(nnodes, dim * nnodesHQ);
			break;
		case 2:  // M_Process only
			size = dim * nnodesHQ;
			Stiffness = new Matrix(size, size);
			RHS = new Vec(size);
			break;
		case 3:  // MH Partioned scheme, M_Process
			size = dim * nnodesHQ;
			Stiffness = new Matrix(size, size);
			RHS = new Vec(size);
			CouplingA = new Matrix(dim * nnodesHQ, nnodes);
			break;
		case 4:  // HM monothlic scheme
			Mass = new Matrix(nnodes, nnodes);
			//        Laplace = new SymMatrix(nnodes);
			Laplace = new Matrix(nnodes, nnodes);
			size = dim * nnodesHQ;
			Stiffness = new Matrix(size, size);
			RHS = new Vec(size + nnodes);
			CouplingA = new Matrix(dim * nnodesHQ, nnodes);
			CouplingB = new Matrix(nnodes, dim * nnodesHQ);
			break;
		case 5:  // Mass Transport process
			Mass = new Matrix(nnodes, nnodes);
			Laplace = new Matrix(nnodes, nnodes);
			Advection = new Matrix(nnodes, nnodes);
			Storage = new Matrix(nnodes, nnodes);
			Content = new Matrix(nnodes, nnodes);
			RHS = new Vec(nnodes);
			break;
	}
}

/**************************************************************************
   ElementMatrix::ElementMatrix

   Arguments:
      const int EleIndex:  Element index,
                           used to get element type.
   Programmaenderungen:
   01/2005     WW

**************************************************************************/
ElementMatrix::~ElementMatrix()
{
	if (Mass) delete Mass;
	if (Laplace) delete Laplace;
	if (Advection) delete Advection;
	if (Storage) delete Storage;
	if (Content) delete Content;
	if (RHS) delete RHS;
	if (CouplingA) delete CouplingA;
	if (CouplingB) delete CouplingB;
	Mass = NULL;
	Laplace = NULL;
	Advection = NULL;
	Storage = NULL;
	Content = NULL;
	RHS = NULL;
	CouplingA = NULL;
	CouplingB = NULL;
}
}  // end namespace FiniteElement
