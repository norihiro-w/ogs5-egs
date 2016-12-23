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

#include "matrix_class.h"

#include "ElementValueDM.h"
#include "fem_ele.h"

namespace SolidProp
{
class CSolidProperties;
}

class CRFProcess;
class CFluidProperties;
class CMediumProperties;

namespace process
{
class CRFProcessDeformation;
}
namespace MeshLib
{
class CElem;
}
using namespace SolidProp;
using namespace Math_Group;

namespace FiniteElement
{

// Derived element for deformation caculation
class CFiniteElementVec : public CElement
{
public:
	CFiniteElementVec(process::CRFProcessDeformation* dm_pcs,
	                  const int C_Sys_Flad, const int order = 2);
	~CFiniteElementVec();

private:
	void SetMemory();
	void LocalAssembly(const int update);
	bool GlobalAssembly();
	void ComputeStrain();
	void SetMaterial();
	double* GetStrain() const { return dstrain; }

private:
	process::CRFProcessDeformation* pcs;
	::CRFProcess* h_pcs;
	::CRFProcess* t_pcs;
	//
	int ns;  // Number of stresses components
	// Flow coupling
	int Flow_Type;

	// Primary value indeces
	// Column index in the node value table
	int idx_P, idx_P0, idx_P1, idx_P1_0, idx_P2;
	int idx_T0, idx_T1;
	int idx_S0, idx_S, idx_Snw;
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
	Matrix* AuxMatrix2;  // NW
	Matrix* Stiffness;
	Matrix* PressureC;
	Matrix* PressureC_S;     // Function of S
	Matrix* PressureC_S_dp;  // Function of S and ds_dp
	Vector* RHS;
	// Global RHS. 08.2010. WW
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
	double* Temp, Tem;
	double* T1;
	double S_Water;

	// Element value
	ElementValue_DM* eleV_DM;

	// principle stresses
	double* pr_stress;
	// Compute principle stresses
	double ComputePrincipleStresses(const double* Stresses);

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
	void LocalAssembly_continuum(const int update);

	// Assembly local stiffness matrix
	void GlobalAssembly_Stiffness();
	void GlobalAssembly_PressureCoupling(Matrix* pCMatrix, double fct,
	                                     const int phase = 0);
	void GlobalAssembly_RHS();
#ifdef USE_PETSC
	void add2GlobalMatrixII();
#endif

	friend class process::CRFProcessDeformation;

	// Auxillarary vector
	double* AuxNodal0;
	double* AuxNodal;
	double* AuxNodal_S0;
	double* AuxNodal_S;
	double* AuxNodal1;
	double* AuxNodal2;
};
}  // end namespace

extern std::vector<FiniteElement::ElementValue_DM*> ele_value_dm;
#endif
