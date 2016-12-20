/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ELEMENT_MATRIX_INC
#define ELEMENT_MATRIX_INC

#include "msh_elem.h"

namespace Math_Group
{
class SymMatrix;
class Matrix;
typedef Matrix Vector;
}

using namespace Math_Group;

namespace FiniteElement
{

/*------------------------------------------------------------------
   Element matrices:
   All local matrices are stored for the purpose of reducing
   compatation  time when steady state of the problem is reached
   12.01.2005. WW
   ------------------------------------------------------------------*/
class ElementMatrix
{
public:
	ElementMatrix()
	    : Mass(NULL),
	      Mass_notsym(NULL),
	      Laplace(NULL),
	      Advection(NULL),
	      Storage(NULL),
	      Content(NULL),
	      CouplingA(NULL),
	      CouplingB(NULL),
	      Stiffness(NULL),
		  RHS(NULL)
	{
	}
	~ElementMatrix();

	// Allocate memory for strain coupling matrix
	void AllocateMemory(MeshLib::CElem* ele, int type = 0);

	// Set members
	void SetMass(Matrix* mass) { Mass = mass; }
	void SetMass_notsym(Matrix* mass) { Mass_notsym = mass; }
	void SetLaplace(Matrix* laplace) { Laplace = laplace; }
	void SetStiffness(Matrix* x) { Stiffness = x; }
	void SetAdvection(Matrix* x) { Advection = x; }
	void SetStorage(Matrix* x) { Storage = x; }
	void SetContent(Matrix* x) { Content = x; }
	void SetCouplingMatrixA(Matrix* cplM) { CouplingA = cplM; }
	void SetCouplingMatrixB(Matrix* cplM) { CouplingB = cplM; }
	void SetRHS(Vector* rhs) { RHS = rhs; }
	// Get members
	Matrix* GetMass() { return Mass; }
	Matrix* GetMass_notsym() { return Mass_notsym; }
	Matrix* GetLaplace() { return Laplace; }
	Matrix* GetStiffness() { return Stiffness; }
	Matrix* GetAdvection()
	{
		return Advection;
	}
	Matrix* GetStorage()
	{
		return Storage;
	}
	Matrix* GetContent()
	{
		return Content;
	}
	Matrix* GetCouplingMatrixA() { return CouplingA; }
	Matrix* GetCouplingMatrixB() { return CouplingB; }
	Vector* GetRHS() { return RHS; }

private:
	Matrix* Mass;
	Matrix* Mass_notsym;
	Matrix* Laplace;
	Matrix* Advection;
	Matrix* Storage;
	Matrix* Content;
	Matrix* CouplingA;  // Pressure coupling for M_Process
	Matrix* CouplingB;  // Strain coupling gor H_Process
	Matrix* Stiffness;
	Vector* RHS;
};
}  // end namespace

#endif
