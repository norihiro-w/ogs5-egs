/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef fem_INC
#define fem_INC

#include "prototyp.h"
#include "MSHEnums.h"
#include "matrix_class.h"

namespace MeshLib
{
class CElem;
class CNode;
class CEdge;
}
class CRFProcess;
namespace process
{
class CRFProcessDeformation;
}

namespace FiniteElement
{
struct ExtrapolationMethod
{
	enum type
	{
		EXTRAPO_LINEAR,
		EXTRAPO_NEAREST,
		EXTRAPO_AVERAGE
	};
};

class CElement
{
public:
	CElement(int CoordFlag, const int order = 1);
	virtual ~CElement();
	//
	void ConfigElement(MeshLib::CElem* MElement, bool FaceIntegration = false);
	void setOrder(const int order);
	// Set Gauss point
	void SetGaussPoint(const int gp, int& gp_r, int& gp_s, int& gp_t);
	// Get Gauss integration information
	double GetGaussData(int gp, int& gp_r, int& gp_s, int& gp_t);

	// Compute values of shape function at integral point unit
	void ComputeShapefct(const int order);
	// Compute the Jacobian matrix. Return its determinate
	double computeJacobian(const int order);

	// Compute values of the derivatives of shape function at integral point
	void ComputeGradShapefct(const int order);
	// Compute the real coordinates from known unit coordinates
	void RealCoordinates(double* realXYZ);
	// Compute the unit coordinates from known unit coordinates
	void UnitCoordinates(double* realXYZ);
	// For axisymmetrical problems
	void CalculateRadius();
	//
	void setUnitCoordinates(double* u)
	{
		for (int i = 0; i < 3; i++)
			unit[i] = u[i];
	}

	// Finite element matrices and vectors
	// Compute the local finite element matrices
	void LocalAssembly(const long, const int) {}
	// Set the number of Gauss points
	// 26.03.2007 WW
	void SetGaussPointNumber(const int nGuassP) { nGauss = nGuassP; }
	// Get values;
	int GetNumGaussPoints() const { return nGaussPoints; }
	int GetNumGaussSamples() const { return nGauss; }
	int Dim() const { return ele_dim; }
	double Getdshapefct(int in) { return dshapefct[in]; }

	// Integrate Neumman type BC
	void FaceIntegration(double* NodeVal);
	void DomainIntegration(double* NodeVal);
	void CalcFaceMass(double* mass_matrix);

	// Coupling
	//
	bool isTemperatureCoupling() const { return T_Flag; }
	bool isFluidPressureCoupling() const { return F_Flag; }
	int isDeformationCoupling() const { return D_Flag; }
	int isConcentrationCoupling() const { return C_Flag; }

	// Interpolate Gauss values
	double interpolate(double* nodalVal, const int order = 1) const;
	double interpolate(const int idx, CRFProcess* m_pcs, const int order = 1);
	// double elemnt_average (const int idx, const int order =1);
	double elemnt_average(const int idx, CRFProcess* m_pcs,
	                      const int order = 1);

	void SetCenterGP();
	int GetGPindex() const { return gp; }
	int GetElementIndex() const { return Index; }
	MeshLib::CElem* GetMeshElement() const
	{
		return MeshElement;
	}

	// For extropolating gauss value to node
	int GetLocalIndex(const int gp_r, const int gp_s, int gp_t);

	void SetRWPT(const int idx)
	{
		PT_Flag = idx;
	}

protected:
	MeshLib::CElem* MeshElement;

	friend class ::CRFProcess;
	friend class process::CRFProcessDeformation;

	// Coordinate indicator
	// 10:  X component only
	// 11: Y component only
	// 12: Z component only
	// 20:  X, Y component
	// 22:  X, Z component
	// 32:  X, Y, Z component
	int coordinate_system;
	bool axisymmetry;
	// Order of shape functions
	// Displacement, 2. Others, 1. Default, 1
	int Order;
	size_t ele_dim;          // Dimension of element
	size_t dim;              // Dimension of real dimension
	int nGaussPoints;        // Number of Gauss points
	int nGauss;              // Number of sample points for Gauss integration
	int gp;                  // Gauss point index.
	mutable double unit[4];  // Local coordintes
	double* Jacobian;        // Jacobian matrix
	double* invJacobian;     // Inverse of Jacobian matrix.
	double* shapefct;        // Results of linear shape function at Gauss points
	double* shapefctHQ;  // Results of quadratic shape function at Gauss points
	// Results of derivatives of linear shape function at Gauss points
	double* dshapefct;
	// Results of derivatives of quadratic shape function at Gauss points
	double* dshapefctHQ;
	//
	double x1buff[3], x2buff[3], x3buff[3], x4buff[3];
	// Pointer to the linear interpolation function
	VoidFuncDXCDX ShapeFunction;
	// Pointer to the quadratic interpolation function
	VoidFuncDXCDX ShapeFunctionHQ;
	// Pointer to the gradient of linear interpolation function
	VoidFuncDXCDX GradShapeFunction;
	// Pointer to the gradient of Quadratic interpolation function
	VoidFuncDXCDX GradShapeFunctionHQ;
	// Coupling
	int NodeShift[5];
	// Displacement column indeces in the node value table
	int Idx_dm0[3];
	int Idx_dm1[3];

	int idx_c0, idx_c1;

	// Coupling flag
	bool T_Flag;   // Temperature
	bool C_Flag;   // Concentration
	bool F_Flag;   // Fluid
	int D_Flag;    // Deformation
	int PT_Flag;   // Particle Tracking Random Walk
	bool RD_Flag;  // Dual Richards
	// For extropolation
	double Xi_p;
	void SetExtropoGaussPoints(const int i);  // 25.2.2007 WW
	double CalcAverageGaussPointValues(double* GpValues);
	double CalcXi_p();

	// Buffer
	int Index;
	int nNodes;
	int nnodes;
	int nnodesHQ;
	double time_unit_factor;
	double Radius;  // For axisymmetrical problems
	long nodes[20];
	long eqs_number[20];
	double dShapefct[27];  // Auxullary
	double X[20];
	double Y[20];
	double Z[20];
	double node_val[20];
	double dbuff[20];

#if defined(USE_PETSC)
	int* idxm;          //> global indices of local matrix rows
	int* idxn;          //> global indices of local matrix columns
#endif
	ExtrapolationMethod::type extrapo_method;
	ExtrapolationMethod::type GetExtrapoMethod() { return extrapo_method; }

private:
	void ConfigNumerics(MshElemType::type elem_type);
};

}  // end namespace

#endif
