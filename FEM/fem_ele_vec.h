/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef fem_dm_INC
#define fem_dm_INC

#include "fem_ele.h"
#include "matrix_class.h"

namespace SolidProp
{
class CSolidProperties;
}

class CRFProcess;
class CRFProcessDeformation;
class CFluidProperties;
class CMediumProperties;

namespace MeshLib
{
class CElem;
}
using namespace SolidProp;
using namespace Math_Group;
using namespace MeshLib;

namespace FiniteElement
{

class ElementValue_DM;

class CFiniteElementVec : public CElement
{
public:
	CFiniteElementVec(CRFProcessDeformation* dm_pcs,
	                  const int C_Sys_Flad, const int order = 2);
	virtual ~CFiniteElementVec();

	void SetMaterial();

	void AssembleLinear();
	void AssembleResidual();
	void AssembleJacobian();

	void UpdateStressStrain();


private:
	void Init();
	void SetMemory();

	bool GlobalAssembly();

	void ComputeStrain();

private:
	CRFProcessDeformation* pcs;
	CRFProcess* h_pcs;
	CRFProcess* t_pcs;
	//
	int ns;  // Number of stresses components
	// Flow coupling
	int Flow_Type;

	// Primary value indeces
	// Column index in the node value table
	int idx_P, idx_P0, idx_P1, idx_P1_0, idx_P2;
	int idx_T0, idx_T1;
	int idx_pls;
	// Displacement column indeces in the node value table
	int* Idx_Stress;
	int* Idx_Strain;

	// B matrix
	Matrix* B_matrix;
	Matrix* B_matrix_T;
	std::vector<Matrix*> vec_B_matrix;
	std::vector<Matrix*> vec_B_matrix_T;

	//------ Material -------
	CSolidProperties* smat;
	CFluidProperties* m_mfp;
	CMediumProperties* m_mmp;

	double CalDensity();

	// Elastic constitutive matrix
	Matrix* De;
	// Consistent tangential matrix
	Matrix* ConsistDep;

	// Local matricies and vectors
	Matrix* AuxMatrix;
	Matrix* AuxMatrix2;
	Matrix* Stiffness;
	Matrix* PressureC;
	Vector* RHS;

	// Global RHS.
	double* b_rhs;

	//  Stresses:
	//  s11, s22, s33, s12, s13, s23
	double* dstress;
	//  Straines:
	//  s11, s22, s33, s12, s13, s23
	double* dstrain;
	double* strain_ne;
	double* stress_ne;
	double* stress0;
	// Results, displacements
	//  u_x1, u_x2, u_x3, ..., u_xn,
	//  u_y1, u_y2, u_y3, ..., u_yn,
	//  u_z1, u_z2, u_z3, ..., u_zn
	double* Disp;

	// Temperatures of nodes
	double* dT, Tem;
	double* T1;

	// Element value
	ElementValue_DM* eleV_DM;

	// Compute principle stresses
	double ComputePrincipleStresses(const double* Stresses, double* pr_stress);

	// Form B matric
	void setB_Matrix(const int LocalIndex);
	// Form the tanspose of B matric
	void setTransB_Matrix(const int LocalIndex);
	//
	void ComputeMatrix_RHS(const double fkt, const Matrix* p_D);

	// Temporarily used variables
	double* Sxx, *Syy, *Szz, *Sxy, *Sxz, *Syz, *pstr;

	/// Extropolation
	bool RecordGuassStrain(const int gp, const int gp_r, const int gp_s,
	                       int gp_t);
	// Effictive strain
	double CalcStrain_v();
	void ExtropolateGaussStrain();
	void ExtropolateGaussStress();
	double CalcStress_eff();

	// Compute the local finite element matrices
	void LocalAssembly_Linear();

	// Assembly local stiffness matrix
	void GlobalAssembly_Stiffness();
	void GlobalAssembly_PressureCoupling(Matrix* pCMatrix, double fct,
	                                     const int phase = 0);
	void GlobalAssembly_RHS();
#ifdef USE_PETSC
	void add2GlobalMatrixII();
#endif
	void assembleGlobalVector();
	void assembleGlobalMatrix();

	friend class ::CRFProcessDeformation;

	// Auxillarary vector
	double* AuxNodal0;
	double* AuxNodal;
	double* AuxNodal1;
	double* AuxNodal2;

};
}  // end namespace

#endif
