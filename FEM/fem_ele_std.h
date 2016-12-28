/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
   Class element declaration
   class for finite element.
   Designed and programmed by WW/OK, 12/2004
   modified by TF, 10/2010
*/

#ifndef fem_std_INC
#define fem_std_INC

//#include "FEMEnums.h"

#include "fem_ele.h"
#include "matrix_class.h"

// Problems
#include "rf_mfp_new.h"
//#include "rf_msp_new.h"
#include "rf_out_new.h"  //OK

//-----------------------------------------------------
// Process type
// L: Liquid flow
// U: Unconfined flow
// G: Groundwater flow
// T: Two-phase flow
// C: Componental flow
// H: heat transport
// M: Mass transport
// R: Richards flow
// F: Fluid momentum
// A: Gas flow
enum EnumProcessType
{
	L,
	U,
	G,
	T,
	C,
	H,
	M,
	R,
	F,
	A,
	V,
	P,
	TH,
	EnumProcessType_INVALID
};
//-----------------------------------------------------

namespace process
{
class CRFProcessDeformation;
}
namespace SolidProp
{
class CSolidProperties;
}
class CMediumProperties;
class CFluidProperties;
class CRFProcess;

using namespace Math_Group;

namespace FiniteElement
{

class CFiniteElementStd : public CElement
{
public:
	CFiniteElementStd(CRFProcess* Pcs,
	                  const int C_Sys_Flad,
	                  const int order = 1,
	                  const EnumProcessType pcsType = EnumProcessType_INVALID);
	~CFiniteElementStd();

	// Set material data
	void SetMaterial(const int phase = 0);
	// Set memory for local matrices
	void SetMemory();
	// Set variable  YD
	void SetVariable();
	//
	void Config();  // NW made it public for outputting MMP

	// Set coupling information
	void ConfigureCoupling(CRFProcess* pcs, const int* Shift);

	// Element claculation
	// 1. Mass matrix
	void CalcMass();
	void CalcMass2();
	void CalcMassPSGLOBAL();  // PCH
	// 2. Lumped mass matrix
	void CalcLumpedMass();
	void CalcLumpedMass2();
	void CalcLumpedMassPSGLOBAL();  // PCH
	// 3. Laplace matrix
	void CalcLaplace();
	// 4. Gravity term
	void CalcGravity();
	// 5. Strain coupling matrix
	void CalcStrainCoupling(int phase = 0);
	// 6. Thermal coupling
	void CalcRHS_by_ThermalDiffusion();
	// 7. Advection matrix
	void CalcAdvection();
	// 8. Storage matrix
	void CalcStorage();
	// 9. Content matrix
	void CalcContent();
	//
	void CalcSatution();  // WW
//
#ifdef E_NORM
	// 25.08.2008. WW
	void CalcEnergyNorm(double& err_norm0, double& err_normn);
	// 25.09.2008. WW
	void CalcEnergyNorm_Dual(double& err_norm0, double& err_normn);
//
#endif
	void CalcNodeMatParatemer();  // WW
	// Assembly
	void Assembly(bool updateA = true,
	              bool updateRHS = true,
	              const bool add2global = true);
	void Assembly(int option, int dimension);  // PCH for Fluid Momentum
	void Cal_Velocity();
	void Cal_Velocity_2();  // CB this is to provide velocity only at the
	                        // element center of gravity
	void Cal_GP_Velocity_FM(int* i_ind);  // SB 4900 interpolate node velocities
	                                      // to Gauss point velocities
	// BG, 04/2012: Provides the average element velocity over all gauss points
	double Get_Element_Velocity(int Index,
	                            CRFProcess* m_pcs,
	                            int phase_index,
	                            int dimension);
	// necessary for using precalculated density and viscosity BG, 11/2010
	double InterpolatePropertyToGausspoint(int GPIndex,
	                                       CRFProcess* m_pcs,
	                                       int Variableindex);
	//
	// OK
	void AssembleParabolicEquationRHSVector();


	//
	// CB added by CB: 090507
	void UpwindAlphaMass(double* alpha);
	// CB added by CB: 090507
	void UpwindUnitCoord(int p, int point, int ind);
	int UpwindElement(int option, int phase);  // PCH
	// CB added by CB: 090507
	void UpwindSummandMass(const int gp,
	                       int& gp_r,
	                       int& gp_s,
	                       int& gp_t,
	                       double* alpha,
	                       double* summand);
	// NW
	double CalcSUPGCoefficient(double* vel,
	                           int ip,
	                           const double* diff_tensor = NULL);
	// NW
	void CalcSUPGWeightingFunction(double* vel,
	                               int ip,
	                               double& tau,
	                               double* v_dN);
	// NW
	double CalcSUPGEffectiveElemenetLength(double* vel);
	// Gauss value
	void ExtropolateGauss(CRFProcess* m_pcs, const int idof);
	//

private:
	bool newton_raphson;  // 24.05.2007 WW
	long index;
	int dof_index;  // 24.02.2007 WW
	// Column index in the node value table
	int idx0, idx1, idxS, idxSn0, idxSn1, idx3;
	int idxp0, idxp1, idxp20, idxp21, idxT0, idxT1;
	int phase;
	int comp;        // Component
	int LocalShift;  // For RHS
	// Danymic
	int* idx_vel_disp, idx_pres;
	// Velocity
	int* idx_vel;  // WW
	// Material properties
	double* mat;
	double* eqs_rhs;  // For DDC WW
	bool heat_phase_change;

	CRFProcess* mfp_pcs;
	SolidProp::CSolidProperties* SolidProp;
	CFluidProperties* FluidProp;
	CFluidProperties* GasProp;
	CMediumProperties* MediaProp;
	CMediumProperties* MediaProp1;  // Matrix for the dual model. YD/WW
	SolidProp::CSolidProperties*
	    SolidProp1;  // Matrix for the dual model. YD/WW
	CRFProcess* pcs;
	::CRFProcess* cpl_pcs;  // Pointer to coupled process. WW
	process::CRFProcessDeformation* dm_pcs;
	bool flag_cpl_pcs;  // OK
	//-------------------------------------------------------
	// Auxillarary matrices
	Matrix* StiffMatrix;
	Matrix* AuxMatrix;
	Matrix* AuxMatrix1;
	// Gravity matrix;
	// 25.2.2007.WW  SymMatrix *GravityMatrix;
	// Gauss point value. Buffers. // Some changes. 27.2.2007 WW
	double TG, TG0, PG, PG0, PG2, PG20, drho_gw_dT;
	double Sw, rhow, poro, dSdp;
	double rho_gw, rho_ga, rho_g, p_gw, M_g, tort;
	//
	double* edlluse;  // WW edlluse[16]
	double* edttuse;  // WW edlluse[16]

	// Local matrices
	Matrix* Mass;  // CB symMatrix *Mass; // unsymmetric in case of upwinding
	Matrix* Mass2;
	Matrix* Laplace;
	Matrix* Advection;  // SB4200
	Matrix* Storage;    // SB4200
	Matrix* Content;    // SB4209
	Matrix* StrainCoupling;
	Vector* RHS;
	DiagonalMatrix* FCT_MassL;  // NW
	//-------------------------------------------------------
	void SetHighOrderNodes();  // 25.2.2007 WW
	// Primary as water head
	bool HEAD_Flag;
	//
	double CalCoefMass();
	// 25.2.2007 WW
	double CalCoefMass2(int dof_index);
	// 03.3.2009 PCH
	double CalCoefMassPSGLOBAL(int dof_index);
	void CalCoefLaplace(bool Gravity, int ip = 0);
	// 10 2008 PCH
	void CalCoefLaplaceMultiphase(int phase, int ip = 0);
	void CalCoefLaplace2(bool Gravity, int dof_index);
	// AKS/NB
	void CalCoefLaplacePSGLOBAL(bool Gravity, int dof_index);
	double CalCoefAdvection();  // SB4200 OK/CMCD
	double CalCoefStorage();    // SB4200
	double CalCoefContent();
	double CalCoefStrainCouping(const int phase = 0);

	double CalcCoefDualTransfer();
	// 27.2.2007 WW
	double CalCoef_RHS_T_MPhase(int dof_index);
	// 27.2.2007 WW
	double CalCoef_RHS_M_MPhase(int dof_index);
	double CalCoef_RHS_PSGLOBAL(int dof_index);
	//  NB
	double CalCoef_RHS_T_PSGlobal(int dof_index);
	// 03.2007 PCH
	void CalCoef_RHS_Pc(int dof_index);
	// AKS
	double CalCoef_RHS_AIR_FLOW(int dof_index);
	// AKS
	double CalCoef_RHS_HEAT_TRANSPORT(int dof_index);
	// AKS
	double CalCoef_RHS_HEAT_TRANSPORT2(int dof_index);
	void CalNodalEnthalpy();

	void ComputeAdditionalJacobi_H2();        // WW
	void ComputeAdditionalJacobi_Richards();  // WW

	//-----------------------------------------------------
	// Process type
	// L: Liquid flow
	// U: Unconfined flow
	// G: Groundwater flow
	// T: Two-phase flow
	// C: Componental flow
	// H: heat transport
	// M: Mass transport
	// O: Liquid flow
	// R: Richards flow
	// A: Gas flow
	// F: Fluid Momentum
	EnumProcessType PcsType;
	//-----------------------------------------------------
	// Local Assembly
	// Assembly of parabolic equation
	void AssembleParabolicEquation();  // OK4104
	void AssembleMixedHyperbolicParabolicEquation();
	void Assemble_strainCPL(const int phase = 0);  // Assembly of strain coupling
	void Assemble_strainCPL_Matrix(const double fac, const int phase = 0);
	void Assemble_totalStressCPL(const int phase = 0);

	void AssembleMassMatrix(int option);  // PCH
	// Assembly of RHS by Darcy's gravity term
	void Assemble_Gravity();
	void Assemble_Gravity_Multiphase();
	// Assembly of RHS by temperature for m-phase flow 27.2.2007 WW
	void Assemble_RHS_T_MPhaseFlow();
	// Assembly of RHS by deformation. 27.2.2007 WW
	void Assemble_RHS_M();
	void Assemble_RHS_Pc();               // 03.2009 PCH
	void Assemble_RHS_AIR_FLOW();         // AKS
	void Assemble_RHS_HEAT_TRANSPORT();   // AKS
	void Assemble_RHS_HEAT_TRANSPORT2();  // AKS
	void
	Assemble_RHS_T_PSGlobal();  // Assembly of RHS by temperature for PSGlobal
	void AssembleRHS(int dimension);  // PCH
	void Assemble_RHS_LIQUIDFLOW();   // NW
	void Assemble_DualTransfer();
	bool check_matrices;             // OK4104
	void AssembleRHSVector();        // OK
	void AssembleCapillaryEffect();  // PCH
	                                 // PCH for debugging
	void AssembleTHEquation(bool updateA = true, bool updateRHS = true);

#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
	void add2GlobalMatrixII(bool updateA = true, bool updateRHS = true);
	void add2GlobalMatrixII_Split(bool updateA = true, bool updateRHS = true);
#else
	void add2GlobalMatrixII(const int block_cols = 2);  // WW. 06.2011
#endif
	void PrintTheSetOfElementMatrices(std::string mark);
	// Friend classes, 01/07, WW
	friend class ::CMediumProperties;
	friend class SolidProp::CSolidProperties;
	friend class ::CFluidProperties;
	// Friend functions. WW
	friend double ::MFPCalcFluidsHeatCapacity(CFiniteElementStd* assem,
	                                          double* var);

	// Auxillarary vectors for node values
	// Vector of local node values, e.g. pressure, temperature.
	// Assume maximium element nodes is 20
	// double OldMatrix[64]; // For grid adapting
	double* NodalVal;
	double* NodalVal0;  //?? NodalValueSaturation, NodalValueTemperature; ...
	double* NodalVal1;
	double* NodalVal2;
	double* NodalVal3;
	double* NodalVal4;
	double* NodalValC;
	double* NodalValC1;
	double* NodalVal_Sat;
	double* NodalVal_SatNW;
	double* NodalVal_p0;
	double* NodalVal_p1;
	double* NodalVal_p2;
	double* NodalVal_p20;  // AKS
	double* NodalVal_T0;
	double* NodalVal_T1;  // AKS
	//
	double* weight_func;  // NW
	void CalcFEM_FCT();   // NW
	//
	friend class ::CRFProcess;

#if 0
	CFiniteElementStd* fem1;
	CFiniteElementStd* fem2;
#endif
	bool add2global;
};


}  // end namespace

/*------------------------------------------------------------------
   Finite element calculation for standard PDE.
   12.12.2004 WW
   ------------------------------------------------------------------*/
#endif
