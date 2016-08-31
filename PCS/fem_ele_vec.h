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

#include <valarray>

#include "fem_ele.h"
#include "matrix_class.h"

class CRFProcess;
class CRFProcessDeformation;
class CFluidProperties;
class CMediumProperties;
class CSolidProperties;

namespace MeshLib
{
class CElem;
}
using namespace Math_Group;
using namespace MeshLib;

namespace FiniteElement
{

class ElementValue_DM;

class CFiniteElementVec : public CElement
{
	friend class ::CRFProcessDeformation;

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

	double CalDensity();

	// Compute principle stresses
	double ComputePrincipleStresses(const double* Stresses, double* pr_stress);

	// Form B matric
	void setB_Matrix(const int LocalIndex);
	// Form the tanspose of B matric
	void setTransB_Matrix(const int LocalIndex);
	//
	void ComputeMatrix_RHS(const double fkt, const Matrix* p_D);

	// Temporarily used variables
	double* Sxx = nullptr, *Syy = nullptr, *Szz = nullptr, *Sxy = nullptr, *Sxz = nullptr, *Syz = nullptr, *pstr = nullptr;

	// Extropolation
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


private:
	CRFProcessDeformation* pcs = nullptr;
	CRFProcess* h_pcs = nullptr;
	CRFProcess* t_pcs = nullptr;

	int ns;  // Number of stresses components

	// Primary value indeces
	// Column index in the node value table
	int idx_P1 = -1;
	int idx_T0 = -1, idx_T1 = -1;
	int idx_pls = -1;
	// Displacement column indeces in the node value table
	int* Idx_Stress = nullptr;
	int* Idx_Strain = nullptr;

	// B matrix
	Matrix* B_matrix = nullptr;
	Matrix* B_matrix_T = nullptr;
	std::vector<Matrix*> vec_B_matrix;
	std::vector<Matrix*> vec_B_matrix_T;

	//------ Material -------
	CSolidProperties* m_msp = nullptr;
	CFluidProperties* m_mfp = nullptr;
	CMediumProperties* m_mmp = nullptr;


	// Elastic constitutive matrix
	Matrix* De = nullptr;
	// Consistent tangential matrix
	Matrix* ConsistDep = nullptr;

	// Local matricies and vectors
	Matrix* AuxMatrix = nullptr;
	Matrix* AuxMatrix2 = nullptr;
	Matrix* Stiffness = nullptr;
	Matrix* PressureC = nullptr;
	Vector* RHS = nullptr;

	//  Stresses:
	//  s11, s22, s33, s12, s13, s23
	std::valarray<double> dstress;
	//  Straines:
	//  s11, s22, s33, s12, s13, s23
	std::valarray<double> dstrain;
	std::valarray<double> strain_ne;
	std::valarray<double> stress_ne;
	std::valarray<double> stress0;
	// Results, displacements
	//  u_x1, u_x2, u_x3, ..., u_xn,
	//  u_y1, u_y2, u_y3, ..., u_yn,
	//  u_z1, u_z2, u_z3, ..., u_zn
	double* Disp = nullptr;

	// Temperatures of nodes
	double* dT = nullptr;
	double* T1 = nullptr;

	// Element value
	ElementValue_DM* eleV_DM = nullptr;

	// Auxillarary vector
	std::valarray<double> AuxNodal;
	std::valarray<double> AuxNodal1;

};
}  // end namespace

#endif
