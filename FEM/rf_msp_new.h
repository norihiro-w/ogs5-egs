/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef rf_msp_new_INC
#define rf_msp_new_INC

#include "fem_ele_std.h"
#include "rf_pcs.h"

#define MSP_FILE_EXTENSION ".msp"

namespace FiniteElement
{
class element;
class CFiniteElementVec;
class CFiniteElementStd;
class ElementValue_DM;
}

namespace Math_Group
{
class Matrix;
}
namespace process
{
class CRFProcessDeformation;
}

#if defined(WIN32)
class CMATGroupEditorDataEdit;  // WW
#endif
namespace SolidProp
{
using FiniteElement::CFiniteElementVec;
using FiniteElement::ElementValue_DM;
using FiniteElement::CFiniteElementStd;
using FiniteElement::ElementValue;
using Math_Group::Matrix;
using process::CRFProcessDeformation;
using ::CRFProcess;
/*---------------------------------------------------------------*/
class CSolidProperties
{
private:
	// Material parameters
	double PoissonRatio;
	int Youngs_mode;
	int excavation;  // 12.2009. WW
	bool excavated;  // 12.2009. To be ..... WW
	Matrix* data_Youngs;
	double ThermalExpansion;
	//
	double s_tol;  // 16.06.2008 WW
	double f_tol;  // 16.06.2008 WW
	double biot_const;
	int bishop_model;           // 05.2011 WX
	double bishop_model_value;  // 05.2011 WX
	double grav_const;          // WW
	Matrix* data_Density;
	//
	Matrix* data_Capacity;
	Matrix* data_Conductivity;
	//
	Matrix* data_Plasticity;
	Matrix* data_Creep;
	//
	int Density_mode;
	//
	int Capacity_mode;
	int Conductivity_mode;
	int Plasticity_type;
	double primary_variable[10];     // CMCD
	double primary_variable_t0[10];  // CMCD
	double primary_variable_t1[10];  // CMCD
	// Creep property
	// 1. Stationary Norton model
	int Creep_mode;
	//
	bool axisymmetry;

	int mode;  // CMCD
	// Swelling pressure
	int SwellingPressureType;
	double Max_SwellingPressure;
	//
	std::string CurveVariable_Conductivity;
	int CurveVariableType_Conductivity;
	// Secondary data
	// Elasticity
	double E;  // Youngs moduls calculated from data_Youngs
	double Lambda;
	double G;   // Shear stress modulus
	double K;   // Bulk modulus
	double Ks;  // WX solid Bulk modulus

	// Rotation matrices and their transpose: UJG 25.11.2009
	Matrix* Crotm;  // If this is needed by permaebility calculation, we keep
	                // it. Otherwise remove it. (To do, UJG/WW)
	Matrix* D_tran;

	// Plasticity
	double dl2;
	// 2. Single yield surface
	Matrix* d2G_dSdS;
	Matrix* d2G_dSdM;
	Matrix* LocalJacobi;  // To store local Jacobi matrix
	Matrix* inv_Jac;      // To store the inverse of the  Jacobi matrix
	Matrix* sumA_Matrix;
	double* rhs_l;  // To store local unknowns of 15
	double* x_l;    // To store local unknowns of 15
	int* Li;
	void AllocateMemoryforSYS();
	void ResizeMatricesSYS(const int Dim);

	// Direct stress integration for Drucker-Prager
	double* devS;
	double* dFds;
	double* dGds;
	double* D_dFds;
	double* D_dGds;
	double* dFtds;               // WX: 08.2010
	double* dGtds;               // WX: 08.2010
	Matrix* ConstitutiveMatrix;  // WX: 08.2010
	// Mini linear solver
	void Gauss_Elimination(const int DimE, Matrix& AA, int* L, double* xx);
	void Gauss_Back(const int DimE, Matrix& AA, double* rhs, int* L,
	                double* xx);
	// Thermal properties
	int thermal_conductivity_tensor_type;
	int thermal_conductivity_tensor_dim;
	double thermal_conductivity_tensor[9];
	std::string thermal_conductivity_tensor_type_name;

	//-------------------------------------------------------------
	// Numeric
	double CalulateValue(const Matrix* data, const double x) const;
	double Kronecker(const int ii, const int jj);

	// Friends that can access to this data explicitly
	friend bool MSPRead(std::string file_base_name);
	friend void MSPWrite(std::string);
	// WW
	friend class FiniteElement::CFiniteElementVec;
	friend class FiniteElement::CFiniteElementStd;
	friend class FiniteElement::ElementValue;
	friend class process::CRFProcessDeformation;
	friend class ::CRFProcess;
#if defined(WIN32)  // 15.03.2008 WW
	friend class ::CMATGroupEditorDataEdit;
#endif
	// WW
public:
	//
	CSolidProperties();
	~CSolidProperties();

	std::ios::pos_type Read(std::ifstream*);
	// CMCD
	FiniteElement::CFiniteElementStd* Fem_Ele_Std;
	std::string name;
	// IO
	std::string file_base_name;
	// Output
	void Write(std::fstream*);
	// CMCD
	void CalPrimaryVariable(std::vector<std::string>& pcs_name_vector);

	//-------------------------------------------------------------
	// Access to data
	//-------------------------------------------------------------
	// 1. Density
	double Density(double refence = 0.0);
	// 2. Thermal
	double Heat_Capacity(double refence = 0.0);
	// Boiling model
	double Heat_Capacity(double temperature, double porosity, double Sat);
	int GetCapacityModel() const { return Capacity_mode; }
	int GetConductModel() const { return Conductivity_mode; }
	bool CheckTemperature_in_PhaseChange(const double T0, const double T1);
	double Enthalpy(double temperature, const double latent_factor);
	double Heat_Conductivity(double refence = 0.0);
	void HeatConductivityTensor(const int dim, double* tensor, int group,
	                            double* variables = NULL);
//   int GetCapacityMode() {return Capacity_mode;};  ??
// 3. Elasticity
#ifdef RFW_FRACTURE
	double Youngs_Modulus(CElem* elem, double refence = 0.0);
	// RFW, for fracture calc
	double Get_Youngs_Min_Aperture(CElem* elem);
#endif
#ifndef RFW_FRACTURE
	double Youngs_Modulus(double refence = 0.0);
#endif
	void SetYoungsModulus(const double El) { (*data_Youngs)(0) = El; }
	double Poisson_Ratio() const { return PoissonRatio; }
	void CalcYoungs_SVV(const double strain_v);
	void CalcYoungs_Drained(const double stress_eff);
	double Thermal_Expansion() const { return ThermalExpansion; }
	double Solid_Bulk_Modulus() const { return Ks; }
	// 4. Plasticity
	int Plastictity() const { return Plasticity_type; }
	double GetPlasticParameter(const int index)
	{
		return (*data_Plasticity)(index);
	}
	// 5. Creep
	int CreepModel() const { return Creep_mode; }
	double GetCreepParameter(const int index) { return (*data_Creep)(index); }
	// Initilize density
	void NullDensity();

	//-------------------------------------------------------------
	// Manipulators of data
	//-------------------------------------------------------------
	// 1. Elasticity
	void Calculate_Lame_Constant();
	// For thermal elastic model
	void ElasticConsitutive(const int Dimension, Matrix* D_e) const;
	// For transverse isotropic linear elasticity: UJG 24.11.2009
	void ElasticConstitutiveTransverseIsotropic(const int Dimension);
	Matrix* getD_tran() const { return D_tran; }
	void CalculateTransformMatrixFromNormalVector(const int Dimension);
	// 2. Plasticity
	// 2.1 Drucker-Prager
	double GetAngleCoefficent_DP(const double Angle);
	double GetYieldCoefficent_DP(const double Angle);
	void CalulateCoefficent_DP();
	bool StressIntegrationDP(const int GPiGPj, const ElementValue_DM* ele_val,
	                         double* TryStress, double& dPhi, const int Update);
	void ConsistentTangentialDP(Matrix* Dep, const double dPhi, const int Dim);
	bool DirectStressIntegrationDP(const int GPiGPj,
	                               const ElementValue_DM* ele_val,
	                               double* TryStress,
	                               const int Update);
	int DirectStressIntegrationDPwithTension(const int GPiGPj,
	                                         Matrix* De,
	                                         const ElementValue_DM* ele_val,
	                                         double* TryStress,
	                                         const int Update,
	                                         double& mm);  // WX
	void TangentialDP(Matrix* Dep);
	void TangentialDP2(Matrix* Dep);                             // WX
	void TangentialDPwithTension(Matrix* Dep, double mm);        // WX
	void TangentialDPwithTensionCorner(Matrix* Dep, double mm);  // WX
	// 2.2 Single yield surface model
	void dF_dNStress(double* dFdS, const double* DevS,
	                 const double* S_Invariants, const double* MatN1,
	                 const int LengthStrs);
	void dF_dStress(double* dFdS, const double* RotV,
	                const double* S_Invariants, const double* MatN1,
	                const int LengthStrs);
	void dF_dMat(double* dFdM, const double* S_Invariants, const double* MatN1);
	void dG_dNStress(double* dGdS, const double* DevS,
	                 const double* S_Invariants, const double* MatN1,
	                 const int LengthStrs);
	void dG__dNStress_dNStress(const double* DevS, const double* S_Invariants,
	                           const double* MatN1, const int LengthStrs);
	void dG__dStress_dStress(const double* DevS, const double* RotV,
	                         const double* S_Invariants, const double* MatN1,
	                         const int LengthStrs);
	void dG_dSTress_dMat(const double* DevS, const double* S_Invariants,
	                     const double* MatN1, const int LengthStrs);
	void dfun2(const double* DevS, const double* RotV,
	           const double* S_Invariants, const double* MatN1,
	           const int LengthStrs);

	int CalStress_and_TangentialMatrix_SYS(const int GPiGPj,
	                                       const ElementValue_DM* ele_val,
	                                       const Matrix* De, Matrix* D_ep,
	                                       double* dStress, const int Update);
	// 2.2 Cam-clay model
	void CalStress_and_TangentialMatrix_CC(const int GPiGPj,
	                                       const ElementValue_DM* ele_val,
	                                       double* dStrain,
	                                       Matrix* Dep,
	                                       const int Update);
	// Substep integration. 16.05.2008 WW
	void CalStress_and_TangentialMatrix_CC_SubStep(
	    const int GPiGPj,
	    const ElementValue_DM* ele_val,
	    double* dStrain,
	    Matrix* Dep,
	    const int Update);
	// Parameter function for thermal elatic model. Last modifed on 15.03.2008
	// //WW
	double TEPSwellingParameter(const double mean_stress);
	void TEPSwellingParameter_kis(const double suction);
	// Strain inrement by creep
	void AddStain_by_Creep(const int ns,
	                       double* stress_n,
	                       double* dstrain,
	                       double temperature = 0.0);
	void AddStain_by_HL_ODS(const ElementValue_DM* ele_val,
	                        double* stress_n,
	                        double* dstrain,
	                        double temperature = 30);
	void CleanTrBuffer_HL_ODS();
	void AccumulateEtr_HL_ODS(const ElementValue_DM* ele_val, const int nGS);

	// Plasticity
	// 1. Drucker-Prager
	double Al;
	double Xi;
	double Y0;
	double BetaN;
	double Hard;
	double Hard_Loc;
	double tension;  // WX:08.2010 Tension strength

	// 4. Mohr-Coulomb	//WX: 11.2010. Mohr-Coulomb model
	double Ntheta;
	double Nphi;
	double csn;
	void CalculateCoefficent_MOHR(double ep);
	void CalPrinStrs(double* stresses, double* prin_stresses, int Size);
	void CalPrinDir(double* prin_str, double* stress, double* v, int Size);
	void CalTransMatrixA(double* v, Matrix* A, int Size);
	int DirectStressIntegrationMOHR(const int GPiGPj,
	                                const ElementValue_DM* ele_val,
	                                double* TryStress, const int Update,
	                                Matrix* Dep);
	int MohrCheckFailure(double* NormStr, int& failurestate, int Size);
	void TangentialMohrShear(Matrix* Dep);
	void TangentialMohrTension(Matrix* Dep);
	void Cal_Inv_Matrix(int Size, Matrix* MatrixA, Matrix* xx);
	double CalVarP(double* vec1, double* vec2, double* sigma_B,
	               double* sigma_l);
	double CalVar_t(double* vecl,
	                double* veclg,
	                Matrix* D,
	                double* sigma_B,
	                double* sigma_l,
	                int Size);
	void CalDep_l(double* vecl, double* veclg, Matrix* D, Matrix* Dep_l,
	              double fkt);
	void VecCrossProduct(double* vec1, double* vec2, double* result_vec);

	// 5. Hoek-Brown WX
	double HoekB_a;
	double HoekB_s;
	double HoekB_mb;
	double HoekB_sigci;
	double HoekB_tens;
	double HoekB_cohe;
	void CalculateCoefficent_HOEKBROWN();
	int StressIntegrationHoekBrown(const int GPiGPj,
	                               const ElementValue_DM* ele_val,
	                               double* TryStress, const int Update,
	                               Matrix* Dep);
	void TangentialHoekBrown(Matrix* Dep);
	void CalPrinStrDir(double* stress, double* prin_str, double* prin_dir,
	                   int Dim);
	// WW
	std::vector<std::string> capacity_pcs_name_vector;
	// WW
	std::vector<std::string> conductivity_pcs_name_vector;
};
}  // end namespace
extern std::vector<SolidProp::CSolidProperties*> msp_vector;
extern bool MSPRead(std::string file_base_name);
extern void MSPWrite(std::string);
extern void MSPDelete();
// OK
extern std::vector<std::string> msp_key_word_vector;
extern void MSPStandardKeywords();  // OK
// OK
extern SolidProp::CSolidProperties* MSPGet(std::string);

extern double StressNorm(const double* s, const int Dim);
extern double TensorMutiplication2(const double* s1, const double* s2,
                                   const int Dim);
extern double TensorMutiplication3(const double* s1,
                                   const double* s2,
                                   const double* s3,
                                   const int Dim);
extern double DeviatoricStress(double* Stress);
#endif
