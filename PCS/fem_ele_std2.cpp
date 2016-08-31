/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_std.h"

#include <cfloat>

#include "Configure.h"
#include "memory.h"
#include "makros.h"

#include "mathlib.h"

#ifdef NEW_EQS
#include "equation_class.h"
#endif
#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif
#include "SparseMatrixDOK.h"

#include "ElementValue.h"
#include "eos.h"
#include "rfmat_cp.h"
#include "rf_mmp_new.h"
#include "rf_msp_new.h"
#include "rf_pcs_dm.h"

#ifdef NEW_EQS
using Math_Group::CSparseMatrix;
#endif

extern double gravity_constant;     // TEST, must be put in input file
#define COMP_MOL_MASS_AIR 28.96     // kg/kmol WW  28.96
#define COMP_MOL_MASS_WATER 18.016  // WW 18.016
#define GAS_CONSTANT 8314.41        // J/(kmol*K) WW
#define GAS_CONSTANT_V 461.5        // WW
#define T_KILVIN_ZERO 273.15        // WW

using namespace std;
using namespace MeshLib;

namespace FiniteElement
{

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   02/2007 WW Multi-phase flow
   05/2008 WW Generalization
**************************************************************************/
double CFiniteElementStd::CalCoefMass2(int dof_index)
{
	int Index = MeshElement->GetIndex();
	double val = 0.0;
	double expfactor = 0.0;
	double dens_arg[3];               // 08.05.2008 WW
	double pert = sqrt(DBL_EPSILON);  // 15.08.2011. WW

	bool diffusion = false;  // 08.05.2008 WW

	if (MediaProp->heat_diffusion_model == 1 && cpl_pcs) diffusion = true;
	dens_arg[1] = 293.15;
	//
	if (pcs->m_num->ele_mass_lumping) ComputeShapefct(1);
	// CB_merge_0513 in case of het K, store local K in permeability_tensor
	double* tensor = NULL;
	tensor = MediaProp->PermeabilityTensor(Index);
	MediaProp->local_permeability = tensor[0];

	switch (dof_index)
	{
		case 0:
			PG = interpolate(NodalVal1);  // Capillary pressure
			dens_arg[0] = PG;             // Should be P_w in some cases
			if (diffusion)
			{
				TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
				dens_arg[1] = TG;
			}
			Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
			rhow = FluidProp->Density(dens_arg);
			dSdp = MediaProp->PressureSaturationDependency(
			    Sw, true);  // JT: dSdp now returns actual sign (i.e. <0)
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			// Storativity   28.05.2008
			// val = MediaProp->StorageFunction(Index,unit,pcs->m_num->time_theta)
			// *Sw;
			// Fluid compressibility
			// val += poro  *Sw* FluidProp->drho_dp / rhow;
			if (SolidProp)
			{
				if (SolidProp->Ks > MKleinsteZahl)  // Storativity WX:28.05.2008
					val -= Sw * (SolidProp->biot_const - poro) / SolidProp->Ks *
					       Sw;
			}
			// Fluid compressibility
			if (fabs(FluidProp->drho_dp) > MKleinsteZahl)
				val -= poro * Sw * FluidProp->drho_dp;
			val += poro * dSdp;  // WX:04.2013 val = poro * dSdp;

			// Coupled (T)
			if (diffusion)
			{
				// Water vapour pressure
				expfactor = COMP_MOL_MASS_WATER / (rhow * GAS_CONSTANT * TG);
				rho_gw = FluidProp->vaporDensity(TG) * exp(-PG * expfactor);
				//
				val -= poro * dSdp * rho_gw / rhow;
				//
				val -= (1.0 - Sw) * poro * COMP_MOL_MASS_WATER * rho_gw /
				       (rhow * GAS_CONSTANT * TG * rhow);
			}
			break;
		case 1:  // 01
			val = 0.0;
			// WX:05.2012 Storgae
			PG = interpolate(NodalVal1);
			Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			if (SolidProp)
			{
				if (SolidProp->Ks > MKleinsteZahl)
					// WX:11.2012
					val += Sw * (SolidProp->biot_const - poro) / SolidProp->Ks;
			}
			// WX:05.2012 Compressibility
			if (fabs(FluidProp->drho_dp) > MKleinsteZahl)
				val += poro * Sw * FluidProp->drho_dp;

			break;
		case 2:  //
			// (1-S)n(d rhop_c/d p_c)
			PG2 = interpolate(NodalVal_p2);
			dens_arg[0] = PG2;  // 28.05.2008. WW
			val = 0.;           // 28.05.2008. WW
			if (diffusion)      // 28.05.2008. WW
			{
				val = (1.0 - Sw) * COMP_MOL_MASS_WATER * rho_gw /
				      (rhow * GAS_CONSTANT * TG * rhow);
				p_gw = rho_gw * GAS_CONSTANT * TG / COMP_MOL_MASS_WATER;
				dens_arg[0] -= p_gw;
				dens_arg[1] = TG;
			}
			rho_ga = GasProp->Density(dens_arg);  // 28.05.2008. WW
			val -= rho_ga * dSdp / rhow;
			val *= poro;
			// WX:11.2012.storage
			if (SolidProp)
			{
				if (SolidProp->Ks > MKleinsteZahl)
					val -= rho_ga / rhow *
					       ((1 - Sw) * (SolidProp->biot_const - poro) /
					        SolidProp->Ks * Sw);
			}

			break;
		case 3:  //
			// Approximation of d dens_g/dp_g 16.08.2011. WW
			dens_arg[0] = PG2 + pert;
			if (diffusion) dens_arg[1] = TG;
			/// d dens_g/dp_g:
			if (GasProp->density_model ==
			    2)  // dens_g/dp_g = drho_dp 02.2012. WX
				val = (1.0 - Sw) * poro * GasProp->rho_0 * GasProp->drho_dp /
				      rhow;
			else
				val = (1.0 - Sw) * poro *
				      (GasProp->Density(dens_arg) - rho_ga) / (pert * rhow);
			// Storage WX:11.2012
			if (SolidProp)
			{
				if (SolidProp->Ks > MKleinsteZahl)
					val += (SolidProp->biot_const - poro) / SolidProp->Ks *
					       (1 - Sw) * GasProp->Density(dens_arg) / rhow;
			}

			break;
	}
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for mass matrix
   Programing:
   03/2009 PCH Multi-phase flow
**************************************************************************/
double CFiniteElementStd::CalCoefMassPSGLOBAL(int dof_index)
{
	int Index = MeshElement->GetIndex();
	double val = 0.0, variables[3];
	double P, T;
	// OK411 double expfactor = 0.0;
	// WW bool diffusion = false;                     //08.05.2008 WW
	// WWif(MediaProp->heat_diffusion_model==273&&cpl_pcs)
	// WW  diffusion = true;
	//
	if (pcs->m_num->ele_mass_lumping) ComputeShapefct(1);
	switch (dof_index)
	{
		case 0:

			// compressibility also for the wetting phase NB
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			Sw = 1.0 - interpolate(NodalVal_SatNW);  // Sw = 1-Snw
			P = interpolate(NodalVal1);              // Pw
			T = interpolate(NodalValC1);
			variables[0] = P;
			variables[1] = T;
			val =
			    poro * (Sw)*FluidProp->drhodP(variables) / FluidProp->Density();
			//		cout << FluidProp->fluid_name << " Pressure: " << P << "
			//Temp:
			//" << ": drhodP: " << FluidProp->drhodP(P,T) << " density: " <<
			// FluidProp->Density() << endl;
			break;
		case 1:  // Snw in the wetting equation
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			val = -poro;
			break;
		case 2:  // Pw in the non-wetting equation
			Sw = 1.0 - interpolate(NodalVal_SatNW);  // Sw = 1 - Snw
			// Pnw = Pw + Pc(Sw)
			P = interpolate(NodalVal1) +
			    MediaProp->CapillaryPressureFunction(Sw);
			//      P = interpolate(NodalVal1);  // Pw
			T = interpolate(NodalValC1);
			variables[0] = P;
			variables[1] = T;
			val = poro * (1. - Sw) * GasProp->drhodP(variables) /
			      GasProp->Density();

			break;
		case 3:  // Snw in the non-wetting equation
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			val = poro;
			break;
	}

	return val;
}


/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix
   Programing:
   10/2008 PCH Implementation
**************************************************************************/
void CFiniteElementStd::CalCoefLaplaceMultiphase(int phase, int /*ip*/)
{
	int i = 0;
	double mat_fac = 1.0;
	double* tensor = NULL;
	// static double Hn[9],z[9];
	int Index = MeshElement->GetIndex();
	double k_rel;
	ComputeShapefct(1);  //  12.3.2007 WW

	// For nodal value interpolation
	//======================================================================
	switch (PcsType)
	{
		default:
			break;

		case T:  // Two-phase flow
			// PCH Rewriting...
			// PCH Laplace mat_fac is accounted for two phases here.
			// thought to be related to the reference pressure.
			tensor = MediaProp->PermeabilityTensor(Index);
			if (pcs->pcs_type_number == 0)
			{
				// PCH Laplace mat_fac is accounted for two phases here.
				// thought to be related to the reference pressure.
				double mat_fac = 0.0;

				idxS = cpl_pcs->GetNodeValueIndex("SATURATION2");

				for (i = 0; i < nnodes; i++)
					NodalVal_Sat[i] = cpl_pcs->GetNodeValue(nodes[i], idxS + 1);
				Sw = 1.0 - interpolate(NodalVal_Sat);

				k_rel = MediaProp->PermeabilitySaturationFunction(Sw, phase);

				// Note here mat_fac is += meaning adding two phases
				mat_fac =
				    time_unit_factor * k_rel / mfp_vector[phase]->Viscosity();

				for (size_t i = 0; i < dim * dim; i++)
					mat[i] = tensor[i] * mat_fac;
			}
			else if (pcs->pcs_type_number == 1)
			{
				int phase = pcs->pcs_type_number;

				idxS = pcs->GetNodeValueIndex("SATURATION2");
				for (i = 0; i < nnodes; i++)
					NodalVal_Sat[i] = pcs->GetNodeValue(nodes[i], idxS + 1);
				Sw = 1.0 - interpolate(NodalVal_Sat);
				k_rel = MediaProp->PermeabilitySaturationFunction(Sw, phase);

				// Here only the second phase accounted.
				mat_fac =
				    time_unit_factor * k_rel / mfp_vector[phase]->Viscosity();
				for (size_t i = 0; i < dim * dim; i++)
					mat[i] = tensor[i] * mat_fac;
			}
			break;
	}
}

///////
/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix of multi-phase
      flow
   Programing:
   02/2007 WW Implementation
   last modification:
**************************************************************************/
void CFiniteElementStd::CalCoefLaplace2(bool Gravity, int dof_index)
{
	double* tensor = NULL;
	double mat_fac = 1.0, m_fac = 0.;
	double fac_perm = 1.;  // WX: factor for Permeability as funktion of
	                       // pressure, strain, etc ... 05.2010
	double expfactor, D_gw, D_ga;
	expfactor = D_gw = D_ga = 0.0;
	double dens_arg[3];      // 08.05.2008 WW
	bool diffusion = false;  // 08.05.2008 WW
	if (MediaProp->heat_diffusion_model == 1 && cpl_pcs) diffusion = true;
	//
	dens_arg[1] = 293.15;
	//
	int Index = MeshElement->GetIndex();
	// CB_merge_0513 in case of het K, store local K in permeability_tensor
	tensor = MediaProp->PermeabilityTensor(Index);
	MediaProp->local_permeability = tensor[0];
	//
	ComputeShapefct(1);  //  12.3.2007 WW

	// WX: 11.05.2010
	PG = interpolate(NodalVal1);
	PG2 = interpolate(NodalVal_p2);
	// WX: cal factor for permeability 11.05.2010
	CFiniteElementStd* h_fem;
	h_fem = this;

	if (MediaProp->permeability_pressure_model > 0)  // 01.09.2011. WW
		fac_perm = MediaProp->PermeabilityFunctionPressure(Index, PG2);
	if (MediaProp->permeability_strain_model > 0)  // 01.09.2011. WW
		fac_perm *= MediaProp->PermeabilityFunctionStrain(Index, nnodes, h_fem);
	//======================================================================
	for (size_t i = 0; i < dim * dim; i++)
		mat[i] = 0.0;
	switch (dof_index)
	{
		case 0:
		{
			PG = interpolate(NodalVal1);
			Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
			//
			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = MediaProp->PermeabilitySaturationFunction(Sw, 0) /
			          FluidProp->Viscosity();
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = -tensor[i] * mat_fac * time_unit_factor *
				         fac_perm;  // WX:05.2010
			// For velocity caculation
			if (!Gravity)
			{
				dens_arg[0] = PG;  // Shdould be Pw in some cases
				if (diffusion)
				{
					TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
					dens_arg[1] = TG;
				}
				//
				rhow = FluidProp->Density(dens_arg);
				poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
				PG2 = interpolate(NodalVal_p2);
				dens_arg[0] = PG2;
				//
				if (diffusion)
				{
					tort = MediaProp->TortuosityFunction(Index, unit,
					                                     pcs->m_num->time_theta);
					tort *= MediaProp->base_heat_diffusion_coefficient *
					        (1 - Sw) * poro * pow(TG / T_KILVIN_ZERO, 1.8);
					expfactor =
					    COMP_MOL_MASS_WATER / (rhow * GAS_CONSTANT * TG);
					rho_gw = FluidProp->vaporDensity(TG) * exp(-PG * expfactor);
					p_gw = rho_gw * GAS_CONSTANT * TG / COMP_MOL_MASS_WATER;
					dens_arg[0] -= p_gw;
				}
				//
				rho_ga = GasProp->Density(dens_arg);
				//
				if (diffusion)
				{
					rho_g = rho_ga + rho_gw;
					// 1/Mg
					M_g = (rho_gw / COMP_MOL_MASS_WATER +
					       rho_ga / GasProp->molar_mass) /
					      rho_g;
					D_gw = tort * rho_g * COMP_MOL_MASS_WATER *
					       GasProp->molar_mass * M_g * M_g / rhow;
					D_gw *= rho_gw / (rhow * PG2);
					for (size_t i = 0; i < dim; i++)
						mat[i * dim + i] -= D_gw * time_unit_factor;
				}
			}
			break;
		}
		case 1:
			if (Gravity)
			{
				PG = interpolate(NodalVal1);
				PG2 = interpolate(NodalVal_p2);
				Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
				dens_arg[0] = PG;  // Shdould be Pw in some cases
				if (diffusion)
				{
					TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
					dens_arg[1] = TG;
				}
				// Liquid density
				rhow = FluidProp->Density(dens_arg);
				dens_arg[0] = PG2;
				rho_gw = 0.0;
				if (diffusion)
				{
					expfactor =
					    COMP_MOL_MASS_WATER / (rhow * GAS_CONSTANT * TG);
					rho_gw = FluidProp->vaporDensity(TG) * exp(-PG * expfactor);
					p_gw = rho_gw * GAS_CONSTANT * TG / COMP_MOL_MASS_WATER;
					dens_arg[0] -= p_gw;
				}
				rho_ga = GasProp->Density(dens_arg);
				rho_g = rho_ga + rho_gw;
			}
			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = MediaProp->PermeabilitySaturationFunction(Sw, 0) /
			          FluidProp->Viscosity();
			m_fac = 0.;
			if (diffusion)
				m_fac = rho_gw *
				        MediaProp->PermeabilitySaturationFunction(Sw, 1) /
				        (GasProp->Viscosity() * rhow);
			if (Gravity)
				mat_fac = mat_fac + m_fac * rho_g / rhow;
			else
				mat_fac += m_fac;
			//
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor *
				         fac_perm;  // WX:05.2010
			//
			if ((!Gravity) && diffusion)
			{
				D_gw = tort * COMP_MOL_MASS_WATER * GasProp->molar_mass * M_g *
				       M_g * rho_g / rhow;
				D_gw *= time_unit_factor * p_gw / (PG2 * PG2);
				for (size_t i = 0; i < dim; i++)
					mat[i * dim + i] -= D_gw;
			}
			break;
		case 2:
			if (diffusion)
			{
				D_ga = tort * COMP_MOL_MASS_WATER * GasProp->molar_mass * M_g *
				       M_g * rho_g / rhow;
				D_ga *= time_unit_factor * rho_gw / (PG2 * rhow);
			}
			else
				D_ga = 0.;
			for (size_t i = 0; i < dim; i++)
				mat[i * dim + i] = D_ga;
			break;
		case 3:
			// WX: for Cal_Velocity, rho_ga muss be calculated again before
			// used. 11.05.2010
			dens_arg[0] = PG2;
			rho_ga = GasProp->Density(dens_arg);
			//
			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = rho_ga *
			          MediaProp->PermeabilitySaturationFunction(Sw, 1) /
			          (GasProp->Viscosity() * rhow);
			//
			if (Gravity)
				//        mat_fac *= rhow/rho_ga;
				mat_fac *= rho_ga / rhow;  // 29.04.2009 WW
			//
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor *
				         fac_perm;  // WX:05.2010
			if ((!Gravity) && diffusion)
			{
				D_ga = tort * rho_g * COMP_MOL_MASS_WATER *
				       GasProp->molar_mass * M_g * M_g / rhow;
				D_ga *= p_gw / (PG2 * PG2);
				for (size_t i = 0; i < dim; i++)
					mat[i * dim + i] += D_ga * time_unit_factor;
			}
			break;
			//------------------------------------------------------------------
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate material coefficient for Laplacian matrix of PS multi-phase
      flow
   Programing:
   03/2009 PCH Implementation
   last modification:
**************************************************************************/
void CFiniteElementStd::CalCoefLaplacePSGLOBAL(bool Gravity, int dof_index)
{
	double* tensor = NULL;
	double mat_fac = 1.0;  // OK411 m_fac=0.;
	double k_rel = 0.0;
	double mfp_arg[2];
	double variables[3];

	int Index = MeshElement->GetIndex();
	//
	ComputeShapefct(1);  //  12.3.2007 WW
	//======================================================================
	for (size_t i = 0; i < dim * dim; i++)
		mat[i] = 0.0;
	switch (dof_index)
	{
		case 0:
			tensor = MediaProp->PermeabilityTensor(Index);
			if (pcs->m_num->ele_upwinding == 1)
			{
				// Doing Upwind elements for saturation by divergent of
				// pressure.
				// Pw upwind
				int WhichNode =
				    UpwindElement((int)(pcs->m_num->ele_upwind_method), 0);
				Sw = 1.0 - NodalVal_SatNW[WhichNode];
			}
			else
				Sw = 1.0 - interpolate(NodalVal_SatNW);
			k_rel = MediaProp->PermeabilitySaturationFunction(Sw, 0);

			// CB_merge_0513
			variables[0] = interpolate(NodalVal1);  // pressure
			variables[1] = interpolate(NodalValC);  // temperature

			mat_fac = k_rel / FluidProp->Viscosity(variables);
			// mat_fac = k_rel / FluidProp->Viscosity();
			// Since gravity for water phase is handled directly in
			// Assemble_Gravity,
			// no need of any code for water phase here.
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
		case 1:
			tensor = MediaProp->PermeabilityTensor(Index);
			mat_fac = 0.0;  // Snw has no laplace term
			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
		case 2:
			tensor = MediaProp->PermeabilityTensor(Index);
			if (pcs->m_num->ele_upwinding == 1)
			{
				// Doing Upwind elements for saturation by divergent of
				// pressure.
				// Pnw upwind
				int WhichNode =
				    UpwindElement((int)(pcs->m_num->ele_upwind_method), 0);
				Sw = 1.0 - NodalVal_SatNW[WhichNode];
			}
			else
				Sw = 1.0 - interpolate(NodalVal_SatNW);

			k_rel = MediaProp->PermeabilitySaturationFunction(Sw, 1);
			// Pnw = Pw + Pc(Sw) //TODO: could cause errors in some cases
			mfp_arg[0] = interpolate(NodalVal1) +
			             MediaProp->CapillaryPressureFunction(Sw);
			mfp_arg[1] = interpolate(NodalValC1);  // TEMPERATURE1 in most cases

			mat_fac = k_rel / GasProp->Viscosity(mfp_arg);

			// The density of the non-wetting phase fluid should be considered
			// here.
			// However, the default water phase density should be canceled out
			// simultaneously.
			if (Gravity) mat_fac *= GasProp->Density() / FluidProp->Density();

			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
		case 3:
			// Snw Laplace from Pc in Eqn 2
			tensor = MediaProp->PermeabilityTensor(Index);

			if (pcs->num_type_name.find("dPcdSwGradSnw") != string::npos)
			{
				//			double Snw = -1.0;
				//			if(pcs->m_num->ele_upwinding == 1)
				//			{
				// Doing Upwind elements for saturation by divergent of
				// pressure.
				// Pnw upwind
				//				int WhichNode =
				// UpwindElement((int)(pcs->m_num->ele_upwind_method), 1); //
				// TF:
				// set, but never used
				//				Snw = NodalVal_SatNW[WhichNode]; // TF: set, but
				//never
				// used
				//			}
				//			else
				//				Snw = interpolate(NodalVal_SatNW); // TF: set,
				//but
				// never used

				CMediumProperties* m_mmp = NULL;
				CElem* thisEle = pcs->m_msh->ele_vector[index];
				int matgrp = thisEle->GetPatchIndex();
				m_mmp = mmp_vector[matgrp];
				k_rel = MediaProp->PermeabilitySaturationFunction(Sw, 1);

				double dPcdSw = 0.0;
				dPcdSw = m_mmp->PressureSaturationDependency(Sw, false);

				// Pnw = Pw + Pc(Sw) // TODO: could cause errors in some cases
				mfp_arg[0] = interpolate(NodalVal1) +
				             MediaProp->CapillaryPressureFunction(Sw);
				mfp_arg[1] = interpolate(NodalValC1);
				mat_fac = k_rel / GasProp->Viscosity(mfp_arg) * (-dPcdSw);
			}
			else
				mat_fac = 0.0;

			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
		case 4:
			// For Vnw
			tensor = MediaProp->PermeabilityTensor(Index);

			//		double Snw = -1.0;
			//		if(pcs->m_num->ele_upwinding == 1)
			//		{
			// Doing Upwind elements for saturation by divergent of pressure.
			// Pnw upwind
			//			int WhichNode =
			// UpwindElement((int)(pcs->m_num->ele_upwind_method), 1); // TF:
			// set, but never used
			//			Snw = NodalVal_SatNW[WhichNode]; // TF: set, but never
			// used
			//		}
			//		else
			//			Snw = interpolate(NodalVal_SatNW); // TF: set, but never
			// used

			//		CElem* thisEle = pcs->m_msh->ele_vector[index]; // TF: set,
			// but never used
			//		int matgrp = thisEle->GetPatchIndex(); // TF: set, but never
			// used
			//		CMediumProperties* m_mmp = mmp_vector[matgrp];
			k_rel = MediaProp->PermeabilitySaturationFunction(Sw, 1);
			mat_fac = k_rel / GasProp->Viscosity();

			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac * time_unit_factor;
			break;
			//------------------------------------------------------------------
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcMass2
   Programming:
   02/2007   WW
 **************************************************************************/
void CFiniteElementStd::CalcMass2()
{
	int i, j, in, jn;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, mat_fac;
	// Material
	int dof_n = 2;
	mat_fac = 1.0;
	//----------------------------------------------------------------------
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeShapefct(1);  // Linear interpolation function
		for (in = 0; in < dof_n; in++)
			for (jn = 0; jn < dof_n; jn++)
			{
				// Material
				mat_fac = CalCoefMass2(in * dof_n + jn);
				mat_fac *= fkt;
				// Calculate mass matrix
				const int jsh = jn * nnodes;
//#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
//				for (i = 0; i < act_nodes; i++)
//				{
//					const int ia = local_idx[i];  // local_idx[ia]
//					const int ish = ia + in * nnodes;
//					for (j = 0; j < nnodes; j++)
//					{
//						(*Mass2)(ish, j + jsh) +=
//						    mat_fac * shapefct[ia] * shapefct[j];
//					}
//				}
//#else
				for (i = 0; i < nnodes; i++)
				{
					const int ish = i + in * nnodes;
					for (j = 0; j < nnodes; j++)
						(*Mass2)(ish, j + jsh) +=
						    mat_fac * shapefct[i] * shapefct[j];
				}
//#endif
			}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcMassPSGLOBAL
   Programming:
   03/2009   PCH
 **************************************************************************/
void CFiniteElementStd::CalcMassPSGLOBAL()
{
	int i, j, in, jn;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, mat_fac;
	// Material
	int dof_n = 2;
	mat_fac = 1.0;
	//----------------------------------------------------------------------
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeShapefct(1);  // Linear interpolation function
		for (in = 0; in < dof_n; in++)
			for (jn = 0; jn < dof_n; jn++)
			{
				// Material
				mat_fac = CalCoefMassPSGLOBAL(in * dof_n + jn);
				mat_fac *= fkt;
// Calculate mass matrix
//#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
//				for (i = 0; i < act_nodes; i++)
//				{
//					const int ia = local_idx[i];
//					for (j = 0; j < nnodes; j++)
//					{
//						(*Mass2)(ia + in * nnodes, j + jn * nnodes) +=
//						    mat_fac * shapefct[ia] * shapefct[j];
//					}
//				}

//#else
				for (i = 0; i < nnodes; i++)
					for (j = 0; j < nnodes; j++)
						(*Mass2)(i + in * nnodes, j + jn * nnodes) +=
						    mat_fac * shapefct[i] * shapefct[j];
//#endif
			}
	}
}

///
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcLumpedMass2
   Programming:
   02/2007   WW Multi-phase flow
 **************************************************************************/
void CFiniteElementStd::CalcLumpedMass2()
{
	int i, in, jn, gp_r, gp_s, gp_t;
	double factor, vol = 0.0;
	int dof_n = 2;
	//----------------------------------------------------------------------
	// Volume
	if (axisymmetry)
	{  // This calculation should be done in CompleteMesh.
		// However, in order not to destroy the concise of the code,
		// it is put here. Anyway it is computational cheap. WW
		vol = 0.0;
		for (gp = 0; gp < nGaussPoints; gp++)
			//---------------------------------------------------------
			//  Get local coordinates and weights
			//  Compute Jacobian matrix and its determinate
			//---------------------------------------------------------
			vol += GetGaussData(gp, gp_r, gp_s, gp_t);
	}
	else
		vol = MeshElement->GetVolume() * MeshElement->GetFluxArea();  // WW. 24.05.2012
	//----------------------------------------------------------------------
	// Initialize
	(*Mass2) = 0.0;
	// Center of the reference element
	SetCenterGP();
	for (in = 0; in < dof_n; in++)
	{
		const int ish = in * nnodes;
		for (jn = 0; jn < dof_n; jn++)
		{
			// Factor
			factor = CalCoefMass2(in * dof_n + jn);
			pcs->timebuffer = factor;  // Tim Control "Neumann"
			// Volume
			factor *= vol / (double)nnodes;
			const int jsh = jn * nnodes;  // WW
//#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
//			for (i = 0; i < act_nodes; i++)
//			{
//				const int ia = local_idx[i];
//				(*Mass2)(ia + ish, ia + jsh) = factor;
//			}
//#else
			for (i = 0; i < nnodes; i++)
				(*Mass2)(i + ish, i + jsh) = factor;
//#endif
		}
	}
	// TEST OUT
	// Mass2->Write();
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: CalcLumpedMassPSGLOBAL
   Programming:
   03/2009   PCH PS_GLOBAL for Multi-phase flow
 **************************************************************************/
void CFiniteElementStd::CalcLumpedMassPSGLOBAL()
{
	int i, in, jn, gp_r, gp_s, gp_t;
	double factor, vol = 0.0;
	int dof_n = 2;
	//----------------------------------------------------------------------
	// Volume
	if (axisymmetry)
	{  // This calculation should be done in CompleteMesh.
		// However, in order not to destroy the concise of the code,
		// it is put here. Anyway it is computational cheap. WW
		vol = 0.0;
		for (gp = 0; gp < nGaussPoints; gp++)
			//---------------------------------------------------------
			//  Get local coordinates and weights
			//  Compute Jacobian matrix and its determinate
			//---------------------------------------------------------
			vol += GetGaussData(gp, gp_r, gp_s, gp_t);
	}
	else
		vol = MeshElement->GetVolume();
	//----------------------------------------------------------------------
	// Initialize
	(*Mass2) = 0.0;
	// Center of the reference element
	SetCenterGP();
	for (in = 0; in < dof_n; in++)
	{
		const int ish = in * nnodes;  // WW
		for (jn = 0; jn < dof_n; jn++)
		{
			// Factor
			factor = CalCoefMassPSGLOBAL(in * dof_n + jn);
			pcs->timebuffer = factor;  // Tim Control "Neumann"
			// Volume
			factor *= vol / (double)nnodes;
			const int jsh = jn * nnodes;  // WW
//#if defined(USE_PETSC)  // || defined(other parallel libs)//03~04.3012. WW
//			for (i = 0; i < act_nodes; i++)
//			{
//				const int ia = local_idx[i];
//				(*Mass2)(ia + ish, ia + jsh) = factor;
//			}
//#else
			for (i = 0; i < nnodes; i++)
				(*Mass2)(i + ish, i + jsh) = factor;
//#endif
		}
	}
	// TEST OUT
	//  Mass2->Write();
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementStd:: Assemby_Gravity
   Aufgabe:
           Assemble the contribution of gravity to RHS in Darcy flow
           to the global system for the pressure equation of multiphase flow

   Programming:
   10/2008   PCH
 **************************************************************************/
// Local assembly
void CFiniteElementStd::Assemble_Gravity_Multiphase()
{
	if ((coordinate_system) % 10 != 2 && (!axisymmetry))
		// 27.2.2007 WW (*GravityMatrix) = 0.0;
		return;
	int i, ii;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp_t = 0;
	double fkt, rho;  //, rich_f;
	double k_rel_iteration;
	// GEO
	double geo_fac = MediaProp->geo_area;
	if (!FluidProp->CheckGravityCalculation()) return;
	long cshift = 0;  // WW
	//
	//
	int dof_n = 1;  // 27.2.2007 WW
	if (PcsType == V) dof_n = 2;

	// WW 05.01.07
	cshift = 0;
	if (pcs->dof > 1) cshift = NodeShift[pcs->continuum];

	// rich_f = 1.0;
	// if(PcsType==R) rich_f = -1.0; //WW

	k_rel_iteration = 1.0;

	for (i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;

	// (*GravityMatrix) = 0.0;
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determination
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Moved from CalCoefLaplace(). 12.3.2007 WW
		// Material
		// PCH
		// Pressure equation is the sum of all pressures for all the phases
		// so is gravity. Thus, this sumation should be considered depending on
		// solving the pressure or saturation equation.
		// Thus, the gravity term for the presure equations should cover
		// two phases.
		if (pcs->pcs_type_number == 0)  // if the pressure equation,
		{
			int numOfPhases = 2;
			for (int p = 0; p < numOfPhases; ++p)
			{
				rho = mfp_vector[p]->Density();
				if (gravity_constant < MKleinsteZahl)  // HEAD version
					rho = 1.0;
				else if (HEAD_Flag)
					rho = 1.0;
				else
					rho *= gravity_constant;

				// Initialization
				fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
				fkt *= rho;  //*rich_f;

				for (ii = 0; ii < dof_n; ii++)
				{
					// CalCoefLaplace does the sumation mat_fac for twophase
					// flow.
					// This is not right for two phase gravity terms, because
					// the equation
					// is not summable this way. It should be seperated.
					if (dof_n == 1) CalCoefLaplaceMultiphase(p);
					if (dof_n == 2) CalCoefLaplace2(false, ii * dof_n + 1);

					// Calculate mass matrix
					for (i = 0; i < nnodes; i++)
						for (size_t k = 0; k < dim; k++)
							NodalVal[i + ii * nnodes] -=
							    fkt * dshapefct[k * nnodes + i] *
							    mat[dim * k + dim - 1];
				}
			}
		}
		else
		{
			rho = mfp_vector[1]->Density();
			if (gravity_constant < MKleinsteZahl)  // HEAD version
				rho = 1.0;
			else if (HEAD_Flag)
				rho = 1.0;
			else
				rho *= gravity_constant;
			fkt *= rho;  //*rich_f;

			for (ii = 0; ii < dof_n; ii++)
			{
				// CalCoefLaplace does the sumation mat_fac for twophase flow.
				// This is not right for two phase gravity terms, because the
				// equation
				// is not summable this way. It should be seperated.
				if (dof_n == 1) CalCoefLaplace(false);
				if (dof_n == 2) CalCoefLaplace2(false, ii * dof_n + 1);

				// Calculate mass matrix
				for (i = 0; i < nnodes; i++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i + ii * nnodes] -= fkt *
						                             dshapefct[k * nnodes + i] *
						                             mat[dim * k + dim - 1];
			}
		}
	}

	//
	cshift += NodeShift[problem_dimension_dm];  // 05.01.07 WW
	int ii_sh = 0;
	for (ii = 0; ii < dof_n; ii++)  // 07.02.07 WW
	{
		cshift += NodeShift[ii];
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
			eqs_rhs[cshift + eqs_number[i]] +=
			    k_rel_iteration * geo_fac * NodalVal[i + ii_sh];
			(*RHS)(i + LocalShift + ii_sh) += NodalVal[i + ii_sh];
		}
	}
	// TEST OUTPUT
	// RHS->Write();
}
////////////////////////////////////////////////////////////////
/*
   void  CFiniteElementStd::Assemble_Gravity()
   {
   if((coordinate_system)%10!=2)
     return;
   int i, j, k, l;
   // ---- Gauss integral
   int gp, gp_r=0, gp_s=0, gp_t;
   gp_t = 0;
   double fkt, rho;
   double k_rel_iteration;
   // GEO
   double geo_fac = MediaProp->geo_area;

   k_rel_iteration = 1.0;

   (*GravityMatrix) = 0.0;
   // Loop over Gauss points
   for (gp = 0; gp < nGaussPoints; gp++)
   {
   //---------------------------------------------------------
   //  Get local coordinates and weights
   //  Compute Jacobian matrix and its determination
   //---------------------------------------------------------
   fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

   //---------------------------------------------------------
   // Compute geometry
   //---------------------------------------------------------
   ComputeGradShapefct(1); // Linear interpolation function

   // Material
   CalCoefLaplace(true);
   rho = FluidProp->Density(Index,unit,pcs->m_num->time_theta);
   if(gravity_constant<MKleinsteZahl) // HEAD version
   rho = 1.0;
   else if(HEAD_Flag) rho = 1.0;
   else
   rho *= gravity_constant;

   fkt *= rho;
   // Calculate mass matrix
   for (i = 0; i < nnodes; i++)
   for (j = 0; j < nnodes; j++)
   {
   if(j>i) continue;
   for (k = 0; k < dim; k++)
   {
   for(l=0; l<dim; l++)
   (*GravityMatrix)(i,j) += fkt*dshapefct[k*nnodes+i]
   *mat[dim*k+l]* dshapefct[l*nnodes+j];
   }
   }
   }

   //TEST OUTPUT
   //GravityMatrix->Write();

   double* G_coord = NULL;
   if((coordinate_system)/10==1)
   G_coord = X;
   else if((coordinate_system)/10==2)
   G_coord = Y;
   else if((coordinate_system)/10==3)
   G_coord = Z;

   for (i = 0; i < nnodes; i++)
   {
   NodalVal[i] = 0.0;
   for (j = 0; j < nnodes; j++)
   NodalVal[i] -= (*GravityMatrix)(i,j)* G_coord[j];
   }

   for (i=0;i<nnodes;i++)
   {
   pcs->eqs->b[NodeShift[problem_dimension_dm] + eqs_number[i]]
   += k_rel_iteration* geo_fac*NodalVal[i];
   (*RHS)(i+LocalShift) += NodalVal[i];
   }
   //TEST OUTPUT
   //RHS->Write();
   }
 */


/***************************************************************************
      GeoSys - Funktion:
              CFiniteElementStd:: Velocity calulation

      Programming:  WW
      08/2005
      03/2007   WW  Multi-phase flow
      11/2007   CB  this function was only introduced to allow the calculation
 of
                    the element center of gravity velocity for upwinding

 **************************************************************************/
// Local assembly
void CFiniteElementStd::Cal_Velocity_2()
{
	int k;
	static double vel[3], vel_g[3];
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	double coef = 0.0;
	int dof_n = 1;
	if (PcsType == V) dof_n = 2;
	//
	gp_t = 0;

	// Get room in the memory for local matrices
	SetMemory();
	// Set material
	SetMaterial();

	ElementValue* gp_ele = ele_gp_value[Index];

	// Loop over Gauss points
	k = (coordinate_system) % 10;
	if (PcsType == T)  // WW/CB
	{
		if (pcs->pcs_type_number == 0)
		{
			// gas pressure
			idx1 = pcs->GetNodeValueIndex("PRESSURE1") + 1;
			for (int i = 0; i < nnodes; i++)
				NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
		}
		else if (pcs->pcs_type_number == 1)
		{
			idxp21 = pcs->GetNodeValueIndex("PRESSURE_CAP");
			// gas pressure
			idx1 = cpl_pcs->GetNodeValueIndex("PRESSURE1") + 1;
			gp_ele = ele_gp_value[Index + (long)pcs->m_msh->ele_vector.size()];
			for (int i = 0; i < nnodes; i++)
				// P_l = P_g - P_cap
				NodalVal[i] = cpl_pcs->GetNodeValue(nodes[i], idx1) -
				              pcs->GetNodeValue(nodes[i], idxp21);
		}
	}
	else
		for (int i = 0; i < nnodes; i++)
			NodalVal[i] = pcs->GetNodeValue(nodes[i], idx1);
	//
	if (PcsType == V)
		for (int i = 0; i < nnodes; i++)
		{
			NodalVal[i] -= pcs->GetNodeValue(nodes[i], idxp21);
			NodalVal1[i] = pcs->GetNodeValue(nodes[i], idxp21);
		}
	//
	gp_ele->Velocity = 0.0;
	//

	gp = 0;

	// for (gp = 0; gp < nGaussPoints; gp++)
	//{
	//---------------------------------------------------------
	//  Get local coordinates and weights
	//  Compute Jacobian matrix and its determination
	//---------------------------------------------------------

	GetGaussData(gp, gp_r, gp_s, gp_t);
	// calculate the velocity at the element center of gravity
	if (PcsType == T) SetCenterGP();  // CB 11/2007

	//---------------------------------------------------------
	// Compute geometry
	//---------------------------------------------------------
	ComputeGradShapefct(1);  // Linear interpolation function
	ComputeShapefct(1);      // Moved from CalCoefLaplace(). 12.3.2007 WW
	if ((PcsType == T) && (pcs->pcs_type_number == 1))  // WW/CB
		flag_cpl_pcs = true;
	// Material
	if (dof_n == 1)
		CalCoefLaplace(true);
	else if (dof_n == 2)
		CalCoefLaplace2(true, 0);
	if ((PcsType == T) && (pcs->pcs_type_number == 1))  // WW/CB
		flag_cpl_pcs = false;

	// Velocity
	for (size_t i = 0; i < dim; i++)
	{
		vel[i] = 0.0;
		for (int j = 0; j < nnodes; j++)
			vel[i] += NodalVal[j] * dshapefct[i * nnodes + j];
		//			 vel[i] += fabs(NodalVal[j])*dshapefct[i*nnodes+j];
	}
	if (PcsType == V)
		for (size_t i = 0; i < dim; i++)
		{
			vel_g[i] = 0.0;
			for (int j = 0; j < nnodes; j++)
				vel_g[i] += NodalVal1[j] * dshapefct[i * nnodes + j];
		}
	// Gravity term
	if (k == 2 && (!HEAD_Flag))
	{
		coef = gravity_constant * FluidProp->Density();
		if (dim == 3 && ele_dim == 2)
		{
			auto &tensor(*MeshElement->getTransformTensor());
			for (size_t i = 0; i < dim; i++)
				for (size_t j = 0; j < ele_dim; j++)
				{
					vel[i] += coef * tensor(i, k) *
					          tensor(2, k);
					if (PcsType == V)
						vel_g[i] += rho_g * gravity_constant *
						            tensor(i, k) *
						            tensor(2, k);
				}
		}  // To be correctted
		else
		{
			if (PcsType == V)
			{
				vel[dim - 1] -= coef;
				vel_g[dim - 1] += gravity_constant * rho_g;
			}
			else
				vel[dim - 1] += coef;
		}
	}
	for (size_t i = 0; i < dim; i++)
		for (size_t j = 0; j < dim; j++)
			//            gp_ele->Velocity(i, gp) -= mat[dim*i+j]*vel[j];  //
			//            unit as that given in input file
			gp_ele->Velocity(i, gp) -=
			    mat[dim * i + j] * vel[j] / time_unit_factor;
	//
	if (PcsType == V)
	{
		CalCoefLaplace2(true, 3);
		coef = rhow / rho_ga;
		for (size_t i = 0; i < dim; i++)
			for (size_t j = 0; j < dim; j++)
				gp_ele->Velocity_g(i, gp) -=
				    coef * mat[dim * i + j] * vel_g[j] / time_unit_factor;
	}
	//
	//   cout << gp << " " << vel[0] << " " << vel[1] << " " << vel[2] << endl;
	//   //Test
	//} // for (gp = 0;...
	//
	if (pcs->Write_Matrix)
	{
		(*pcs->matrix_file) << "### Element: " << Index << "\n";
		(*pcs->matrix_file) << "---Velocity of water "
		                    << "\n";
		gp_ele->Velocity.Write(*pcs->matrix_file);
		if (gp_ele->Velocity_g.Size() > 0)
		{
			(*pcs->matrix_file) << "---Velocity of gas "
			                    << "\n";
			gp_ele->Velocity_g.Write(*pcs->matrix_file);
		}
	}
	// gp_ele->Velocity.Write();
}


/**************************************************************************
   FEMLib-Method:
   Task: Calculate  coefficient of displacement induced RHS of multi-phase
      flow
   Programing:
   02/2007 WW Implementation
   05/2008 WW Generalization
   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoef_RHS_M_MPhase(int dof_index)
{
	double val = 0.0;
	double expfactor = 0.0;
	double dens_aug[3];
	dens_aug[1] = 293.15;
	bool diffusion = false;  // 08.05.2008 WW
	if (MediaProp->heat_diffusion_model == 1 && cpl_pcs) diffusion = true;
	//======================================================================
	switch (dof_index)
	{
		case 0:
			PG = interpolate(NodalVal1);
			PG2 = interpolate(NodalVal_p2);  // JT
			Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
			if (diffusion)
			{
				TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
				dens_aug[1] = TG;
			}
			//
			dens_aug[0] =
			    PG2 -
			    PG;  // JT: this was wrong. density argument is water pressure.
			rhow = FluidProp->Density(dens_aug);
			val = Sw;
			//
			dens_aug[0] = PG2;
			//
			if (diffusion)
			{
				expfactor = COMP_MOL_MASS_WATER / (rhow * GAS_CONSTANT * TG);
				rho_gw = FluidProp->vaporDensity(TG) * exp(-PG * expfactor);
				p_gw = rho_gw * GAS_CONSTANT * TG / COMP_MOL_MASS_WATER;
				dens_aug[0] -= p_gw;
			}
			rho_ga = GasProp->Density(dens_aug);
			if (diffusion) val += (1.0 - Sw) * rho_gw / rhow;
			break;
		case 1:
			val = (1.0 - Sw) * rho_ga / rhow;
			break;
			//------------------------------------------------------------------
	}
	return val;
}


///////
/**************************************************************************
   FEMLib-Method:
   Task: Calculate  coefficient of temperature induced RHS of multi-phase
      flow
   Programing:
   02/2007 WW Implementation
   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoef_RHS_T_MPhase(int dof_index)
{
	double val = 0.0, D_gw = 0.0, D_ga = 0.0;
	double expfactor = 0.0, dens_arg[3];
	int Index = MeshElement->GetIndex();
	ComputeShapefct(1);
	//======================================================================
	switch (dof_index)
	{
		case 0:
			PG = interpolate(NodalVal1);
			Sw = MediaProp->SaturationCapillaryPressureFunction(PG);
			TG = interpolate(NodalValC1) + T_KILVIN_ZERO;
			TG0 = interpolate(NodalValC) + T_KILVIN_ZERO;
			PG2 = interpolate(NodalVal_p2);
			rhow = FluidProp->Density();
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			expfactor = COMP_MOL_MASS_WATER / (rhow * GAS_CONSTANT * TG);
			rho_gw = FluidProp->vaporDensity(TG) * exp(-PG * expfactor);
			//
			drho_gw_dT = (FluidProp->vaporDensity_derivative(TG) +
			              PG * expfactor * FluidProp->vaporDensity(TG) / TG) *
			             exp(-PG * expfactor);
			val = (1. - Sw) * poro * drho_gw_dT / rhow;
			//
			if (SolidProp)
				val -= (1.0 - poro) * ((1 - Sw) * rho_gw / rhow + Sw) *
				       SolidProp->Thermal_Expansion();
			//
			// val += n*(1.0-rho_gw/rhow)*(dSw/dT)
			val *= (TG - TG0);
			break;
		case 1:
			//
			val = -(1. - Sw) * poro * drho_gw_dT / rhow;
			//
			if (SolidProp)
				val -= (1.0 - poro) * (1 - Sw) * rho_ga *
				       SolidProp->Thermal_Expansion() / rhow;
			//
			// val -= n*rho_ga/rhow)*(dSw/dT)
			//---------------------------------------------------------------
			val *= (TG - TG0);
			break;
		case 2:
			//------------------------------------------------------------------------
			// From grad (p_gw/p_g)
			tort = MediaProp->TortuosityFunction(Index, unit,
			                                     pcs->m_num->time_theta);
			tort *= MediaProp->base_heat_diffusion_coefficient * (1 - Sw) *
			        poro * pow(TG / T_KILVIN_ZERO, 1.8);
			p_gw = rho_gw * GAS_CONSTANT * TG / COMP_MOL_MASS_WATER;
			dens_arg[0] = PG2 - p_gw;
			dens_arg[1] = TG;
			rho_ga = GasProp->Density(dens_arg);  // AKS SEP 2010
			// //(PG2-p_gw)*GasProp->molar_mass/(GAS_CONSTANT*TG);
			rho_g = rho_ga + rho_gw;
			// 1/Mg
			M_g =
			    (rho_gw / COMP_MOL_MASS_WATER + rho_ga / GasProp->molar_mass) /
			    rho_g;
			D_gw = tort * rho_g * COMP_MOL_MASS_WATER * GasProp->molar_mass *
			       M_g * M_g / rhow;
			val = D_gw * drho_gw_dT * GAS_CONSTANT * TG /
			      (COMP_MOL_MASS_WATER * PG2) * time_unit_factor;
			break;
		case 3:
			//---------------------------------------------------------------
			//
			D_ga = tort * rho_g * COMP_MOL_MASS_WATER * GasProp->molar_mass *
			       M_g * M_g / rhow;
			// From grad (p_gw/p_g)
			val = -D_ga * drho_gw_dT * GAS_CONSTANT * TG /
			      (COMP_MOL_MASS_WATER * PG2) * time_unit_factor;

			break;
			//------------------------------------------------------------------
	}
	return val;
}
/**************************************************************************
   FEMLib-Method:
   Task: Calculate coefficient of temperature induced RHS of PSGlobal scheme

   Programing:

   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoef_RHS_T_PSGlobal(int dof_index)
{
	double val = 0.0, variables[3];  // OK411 D_gw=0.0, D_ga=0.0;
	// OK411 double expfactor=0.0;
	double P, T;
	int Index = MeshElement->GetIndex();
	ComputeShapefct(1);
	//======================================================================
	switch (dof_index)
	{
		case 0:
			val = 0.;
			break;
		case 1:
			poro = MediaProp->Porosity(Index, pcs->m_num->time_theta);
			Sw = 1.0 - interpolate(NodalVal_SatNW);
			// Pnw = Pw + Pc(Sw)
			P = interpolate(NodalVal1) +
			    MediaProp->CapillaryPressureFunction(Sw);
			//      P  = interpolate(NodalVal1);  // Pw
			T = interpolate(NodalValC1);
			variables[0] = P;
			variables[1] = T;
			val = -(1. - Sw) * poro * GasProp->drhodT(variables) /
			      GasProp->Density();
			break;
		case 2:
			val = 0.;
			break;
		case 3:
			val = 0.;
			break;
			//------------------------------------------------------------------
	}
	return val;
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate  Capillary pressure for RHS in the global scheme
   Programing:
   03/2009 PCH Implementation
   last modification:
**************************************************************************/
void CFiniteElementStd::CalCoef_RHS_Pc(int dof_index)
{
	double* tensor = NULL;
	double mat_fac = 1.0;  // OK411 m_fac=0.;
	double k_rel = 0.0;

	int Index = MeshElement->GetIndex();
	//
	ComputeShapefct(1);  //  12.3.2007 WW
	//======================================================================
	for (size_t i = 0; i < dim * dim; i++)
		mat[i] = 0.0;

	switch (dof_index)
	{
		case 0:
			break;
		case 1:
			break;
		case 2:
			tensor = MediaProp->PermeabilityTensor(Index);
			Sw = 1.0 - interpolate(NodalVal_SatNW);
			k_rel = MediaProp->PermeabilitySaturationFunction(Sw, 1);
			mat_fac = k_rel / GasProp->Viscosity() * time_unit_factor;

			for (size_t i = 0; i < dim * dim; i++)
				mat[i] = tensor[i] * mat_fac;
			break;
		case 3:
			break;
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Calculate coefficient of Pc induced RHS of multi-phase
      flow
   Programing:
   03/2009 PCH Implementation
   last modification:
**************************************************************************/
double CFiniteElementStd::CalCoef_RHS_PSGLOBAL(int dof_index)
{
	ComputeShapefct(1);
	for (size_t i = 0; i < dim * dim; i++)
		mat[i] = 0.0;

	switch (dof_index)
	{
		case 0:
			return 0.0;
			break;
		case 1:
			return 0.0;
			break;
		case 2:
		{
			Sw = 1.0 - interpolate(NodalVal_SatNW);
			double k_rel(MediaProp->PermeabilitySaturationFunction(Sw, 1));
			return (k_rel / GasProp->Viscosity() * time_unit_factor);
			break;
		}
		case 3:
			return 0.0;
			break;
	}
	return 0.0;
}

/***************************************************************************
   GeoSys - Funktion:
          Assemble_RHS_T_MPhaseFlow
   Programming:
   02/2007   WW
 **************************************************************************/
void CFiniteElementStd::Assemble_RHS_T_MPhaseFlow()
{
	int i, j, ii;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, fac;
	// Material
	int dof_n = 2;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	// 02.2011 WW
	int dm_shift = 0;
	if (pcs->type / 10 == 4) dm_shift = problem_dimension_dm;
#endif
	//----------------------------------------------------------------------
	for (i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function
		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			fac = fkt * CalCoef_RHS_T_MPhase(ii) / dt;
			// Calculate THS
			for (i = 0; i < nnodes; i++)
				NodalVal[i + ii * nnodes] += fac * shapefct[i];
		}
		// grad T
		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			fac = fkt * CalCoef_RHS_T_MPhase(ii + dof_n);
			// Calculate THS
			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					for (size_t k = 0; k < dim; k++)
						NodalVal[i + ii * nnodes] +=
						    fac * dshapefct[k * nnodes + i] *
						    dshapefct[k * nnodes + j] *
						    (NodalValC1[j] + T_KILVIN_ZERO);
		}
	}
	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	for (ii = 0; ii < 2; ii++)
	{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
		i_sh = NodeShift[ii + dm_shift];
#endif
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i + ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) -= NodalVal[i + ii_sh];
		}
	}
	//
}
/***************************************************************************
   GeoSys - Function: Assemble_RHS_T_PSGlobal
   Programming:
   09/2009
 **************************************************************************/
void CFiniteElementStd::Assemble_RHS_T_PSGlobal()
{
	int i, ii;  // OK411 j, k,
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt, fac;
	// Material
	int dof_n = 2;
	double temp[10];

	for (i = 0; i < 10; i++)
		temp[i] = 0;  // remove

	//----------------------------------------------------------------------
	for (i = 0; i < dof_n * nnodes; i++)
		NodalVal[i] = 0.0;
	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function

		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			fac = fkt * CalCoef_RHS_T_PSGlobal(ii) / dt;
			// Calculate THS
			for (i = 0; i < nnodes; i++)
			{
				NodalVal[i + ii * nnodes] += fac * shapefct[i];
				// remove
				temp[i + ii * nnodes] += fac * shapefct[i];
			}
		}
	}
	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	for (ii = 0; ii < pcs->dof; ii++)
	{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
		i_sh = NodeShift[ii];
#endif
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[i_sh + eqs_number[i]] -= NodalVal[i + ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) -= NodalVal[i + ii_sh];
		}
	}
	//
}

/***************************************************************************
   GeoSys - Funktion:
          Assemble_RHS_Pc
   Programming:
   03/2009   PCH
 **************************************************************************/
void CFiniteElementStd::Assemble_RHS_Pc()
{
	int i, j, ii;
	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t = 0;
	double fkt;
	// Material
	int dof_n = 2;
	int ndx_p_cap = pcs->GetNodeValueIndex("PRESSURE_CAP");
	//----------------------------------------------------------------------

	//      double temp[20];

	for (i = 0; i < dof_n * nnodes; i++)
	{
		//         temp[i] = NodalVal[i] = 0.0;
		NodalVal[i] = 0.0;
		NodalVal1[i] = 0.0;
	}
	for (i = 0; i < nnodes; i++)
		//        temp[i+dof_n] = NodalVal1[i+dof_n] =
		//        -pcs->GetNodeValue(nodes[i],ndx_p_cap);
		NodalVal1[i + dof_n] = -pcs->GetNodeValue(nodes[i], ndx_p_cap);

	//======================================================================
	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);
		// Compute geometry
		ComputeGradShapefct(1);  // Linear interpolation function
		ComputeShapefct(1);      // Linear interpolation function
		// grad Pc
		for (ii = 0; ii < dof_n; ii++)
		{
			// Material
			CalCoef_RHS_Pc(ii + dof_n);
			// Calculate Pc
			for (i = 0; i < nnodes; i++)
				for (j = 0; j < nnodes; j++)
					for (size_t k = 0; k < dim; k++)
						for (size_t l = 0; l < dim; l++)
							NodalVal[dof_n + i] += fkt * mat[dim * k + l] *
							                       dshapefct[k * nnodes + i] *
							                       dshapefct[l * nnodes + j] *
							                       NodalVal1[j + dof_n];
		}
	}

	//      for(i=0; i<2*nnodes; ++i)
	//         temp[i]=NodalVal[i];

	int ii_sh;
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
	long i_sh;
#endif
	for (ii = 0; ii < pcs->dof; ii++)
	{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
		i_sh = NodeShift[ii];
#endif
		ii_sh = ii * nnodes;
		for (i = 0; i < nnodes; i++)
		{
#if !defined(USE_PETSC)  // && !defined(other parallel libs)//03~04.3012. WW
			eqs_rhs[i_sh + eqs_number[i]] += NodalVal[i + ii_sh];
#endif
			(*RHS)(i + LocalShift + ii_sh) += NodalVal[i + ii_sh];
		}
	}
	//
}

}  // end namespace

