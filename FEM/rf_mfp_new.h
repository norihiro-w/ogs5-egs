/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef rf_mfp_new_INC
#define rf_mfp_new_INC

#include <iostream>
#include <string>
#include <vector>

class CompProperties;
class CRFProcess;

extern double gravity_constant;

namespace FiniteElement
{
class CFiniteElementStd;
}

using FiniteElement::CFiniteElementStd;

class CFluidProperties
{
	friend bool MFPRead(std::string);

public:
	double getCriticalDensity() const { return rhoc; }
	double getCriticalTemperature() const { return Tc; }
	double getCriticalPressure() const { return pc; }
	double getSpecificGasConstant() const { return Rs; }
	double getMolarMass() const { return molar_mass; }
	double getUniversalGasConstant() const { return Ru; }
	/**
	 * get the acentric factor for Peng-Robinson EOS (omega)
	 * @return omega
	 */
	double getAzentricFactor() const { return omega; }
	double getDiffusion() const { return diffusion; }
	double getSpecificHeatCapacity() const { return specific_heat_capacity; }

	/**
	 * get the reference temperature
	 * @return
	 */
	double getReferenceTemperature() const { return T_0; }

private:
	double rhoc;   // critical_density; //NB
	double Tc;     // critical_temperature;
	double pc;     // critical_pressure;
	double Tt;     // triple_point_temperature;
	double pt;     // triple_point_pressure;
	double Rs;     // specific_gas_constant;
	double Ru;     // universal_gas_constant;
	double omega;  // azentric factor for Peng-Robinson EOS
	double molar_mass;
	double Vd, Zc, n0, m0, a, b, k1, k2, k3;  // constants are used in VTPR
	                                          /**
	                                           * density
	                                           */
	double rho_0;
	/**
	 * density deviated with respect to the pressure
	 */
	double drho_dp;
	/**
	 * density deviated with respect to the temperature
	 */
	double drho_dT;
	/**
	 * density deviated with respect to the concentration
	 */
	double drho_dC;
	double rho_p0;
	double rho_T0;
	double rho_C0;

	double diffusion; /*SB:2p */

	// Viscosity
	double viscosity;
	double viscosity_T_shift;  // JM in order to use some viscosity functions,
	                           // based on total temperature, within Richards
	double my_0;
	double my_rho0;
	double my_p0;
	double my_T0;
	double my_C0;
	double my_Tstar;
	double dmy_dp;
	double dmy_dT;
	double dmy_dC;

	// Thermal properties
	double specific_heat_capacity;
	double heat_conductivity;
	double temperature_buffer;  // YD, shifted to public JOD

	// State variables
	double p_0;
	/**
	 * state variable: reference temperature
	 */
	double T_0;
	double T_1;
	double C_0;
	double C_1;
	double Z;

	// Chemical properties
	double T_Latent1, T_Latent2, latent_heat;

public:
	int fluid_id;  // specification of substance (NB JUN 09)
	std::string name;

private:
	std::string fluid_name;  // NB4801

	// compressibility
	int compressibility_model_pressure;     // NB
	int compressibility_model_temperature;  // NB
	double compressibility_pressure = 0;           // NB
	double compressibility_temperature;        // NB
	int JTC;                                // NB

	int phase;

private:
	// Limits and coefficients for free Helmholtz Energy, NB JUN 09
	int limit[5];
	double k[2][8];
	double K[14][56];

public:
	// FEM
	CFiniteElementStd* Fem_Ele_Std;
	long node;

	// Density
	int density_model;

	// Viscosity
	int viscosity_model;
	int heat_conductivity_model;
	int heat_diffusion_model;
	int heat_capacity_model;

	// Electrical properties
	int electric_conductivity_model;
	int electric_conductivity_num_val;
	double* electric_conductivity_val;

	// Chemical properties
	int diffusion_model;
	int heat_phase_change_curve;

	// IO
	int mode;

	// PCS  YD
	std::vector<std::string> density_pcs_name_vector;
	std::vector<std::string> viscosity_pcs_name_vector;
	std::vector<std::string> specific_heat_capacity_pcs_name_vector;
	std::vector<std::string> heat_conductivity_pcs_name_vector;
	std::vector<std::string> enthalpy_pcs_name_vector;

	std::vector<CompProperties*> component_vector;

public:
	CFluidProperties(void);
	~CFluidProperties(void);
	std::ios::pos_type Read(std::ifstream*);
	void Write(std::ofstream*) const;
	void CalPrimaryVariable(std::vector<std::string>& pcs_name_vector);
	double Density(double* variables = NULL);
	double MixtureSubProperity(int properties,
	                           long idx_elem,
	                           double p,
	                           double T);
	double GetElementValueFromNodes(long ElementIndex,
	                                int GPIndex,
	                                int PhaseIndex,
									int Variable);
	double drhodP(double* variables);
	double drhodT(double* variables);
	double Viscosity(double* variables = NULL);
	double dViscositydP(double* variables);
	double dViscositydT(double* variables);

	double SpecificHeatCapacity(double* variables = NULL);
	void therm_prop(std::string caption);
	double PhaseChange();
	double HeatConductivity(double* variables = NULL);
	double CalcEnthalpy(double temperature);
	double vaporDensity(const double T);
	double vaporDensity_derivative(const double T);
	bool CheckGravityCalculation() const { return cal_gravity; }
	int GetHeatCapacityModel() const
	{
		return heat_capacity_model;
	}
	// Derivations of free Helmholtz energy, NB JUN 09
	double phi_r_d(double rho, double T) const;
	double phi_r_tt(double rho, double T) const;
	double phi_0_t(double T) const;
	double phi_r_t(double rho, double T) const;
	double phi_r_dt(double rho, double T) const;
	double phi_r_dd(double rho, double T) const;
	double phi_0_tt(double T) const;
	double SuperCompressibiltyFactor(int idx_elem, double p, double T);
	double dZ(int idx_elem, double p, double T, int nk);
	double MaxwellStefanDiffusionCoef(int idx_elem,
	                                  double p,
	                                  double T,
									  int CNm);
	bool drho_dT_unsaturated;

private:
	// State variables
	double primary_variable[10];
	double primary_variable_t0[10];
	double primary_variable_t1[10];
	bool cal_gravity;

	friend class FiniteElement::CFiniteElementStd;

	double GasViscosity_Reichenberg_1971(double, double);
	double MATCalcFluidDensityMethod8(double p, double T, double C);
	double MATCalcFluidDensityFabi(double p, double T);
	double LiquidViscosity_Yaws_1976(double);
	double LiquidViscosity_Marsily_1986(double);
	double LiquidViscosity_Fabi(double);
	double LiquidViscosity_NN(double, double);
	double LiquidViscosity_LJH_MP1(double c, double T);
	double LiquidViscosity_LJH_MP2(double c, double T, double rho);
	double LiquidViscosity_CMCD(double p, double T, double C);
	double LiquidViscosity_Ramey1974(double T);
	double MATCalcHeatConductivityMethod2(double p, double T, double C);
	double MATCalcFluidHeatCapacityMethod2(double p, double T, double C);
};

extern std::vector<CFluidProperties*> mfp_vector;
extern bool MFPRead(std::string);
extern void MFPWrite(std::string);
#define MFP_FILE_EXTENSION ".mfp"
extern double MFPCalcFluidsHeatCapacity(CFiniteElementStd* assem = NULL,
                                        double* var = NULL);
extern double MFPCalcFluidsHeatConductivity(long index,
                                            double* gp,
                                            double theta,
                                            CFiniteElementStd* assem = NULL);
extern void MFPDelete();
extern CFluidProperties* MFPGet(const std::string&);
extern CFluidProperties* MFPGet(int);
double MFPGetNodeValue(long, const std::string&, int);

#endif
