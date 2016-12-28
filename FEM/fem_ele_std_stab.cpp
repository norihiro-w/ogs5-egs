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
#include <iostream>

#include "mathlib.h"

#include "ElementValue.h"
#include "rf_mmp_new.h"
#include "rf_pcs.h"

#define T_KILVIN_ZERO 273.15

using namespace std;

namespace FiniteElement
{

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2009 PCH Upwind Material Scheme
   Background:
   Now material property at the upstream node is taken to exclude future
   predition at the current due to abrupt change of material properties.
   This is conservative perticularly the nodes very close to the interface
   that divides two highly different materials. Thus,this is vector
   characterisitics. In our case, pressure gradient can be a determinant
   for permeability and saturation and so other material properties.

   Description:
   0: none yet
   1: Upwind element determined by div of pressure
**************************************************************************/
int CFiniteElementStd::UpwindElement(int option, int phase)
{
	int WhichNodeInTheElement = -1;  // Initialized to be none of elements
	double Pmin = 1.0 / DBL_MIN;     // Just set to be outrageously big.
	int GravityOn = 1;               // Initialized to be on

	// If no gravity, then set GravityOn to be zero.
	if ((coordinate_system) % 10 != 2 && (!axisymmetry)) GravityOn = 0;

	if (option == 1)  // If upwind by divergent of pressure
	{
		double PdivAtnode = -1.0 / DBL_MIN;  // Meaningless pressure.
		double Pc = 0.0;                     // Set to be no capillary at all.
		int idx_p1 = pcs->GetNodeValueIndex("PRESSURE1");
		int idx_pc = pcs->GetNodeValueIndex("PRESSURE_CAP");

		for (int i = 0; i < nnodes; ++i)
		{
			double Pw = pcs->GetNodeValue(nodes[i], idx_p1 + 1);
			if (phase == 0)
			{
				PdivAtnode = Pw;
				if (GravityOn)
					PdivAtnode -= FluidProp->Density() * gravity_constant;
			}
			else if (phase == 1)
			{
				Pc = pcs->GetNodeValue(nodes[i], idx_pc);
				PdivAtnode = Pw + Pc;
				if (GravityOn)
					PdivAtnode -= GasProp->Density() * gravity_constant;
			}
			else
			{
				std::cout << "Phase number is wrong in UpwindElement."
				          << "\n";
				abort();
			}

			if (Pmin > PdivAtnode)
			{
				Pmin = PdivAtnode;
				WhichNodeInTheElement = i;
			}
		}
	}

	if (WhichNodeInTheElement == -1)
	{
		std::cout << "UpwindElement is failed. Impossible node index!!!"
		          << "\n";
		std::cout << "Pmin = " << Pmin << "\n";
		abort();
	}
	return WhichNodeInTheElement;
}
// CB 090507
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 CB
   last modification:
**************************************************************************/
// void CFiniteElementStd::UpwindUnitCoord(int p, int point, int ind, double
// *rupw, double *supw, double *tupw)
void CFiniteElementStd::UpwindUnitCoord(int /*p*/, int point, int ind)
{
	double scale;
	double alpha[3] = {};
	int gp_r, gp_s, gp_t;
	bool mmp_element_integration_method_maximum = false;

	int upwind_meth;
	double upwind_para;
	double v[3] = {}, v_rst[3];

	//
	ElementValue* gp_ele = ele_gp_value[ind];
	if (pcs->pcs_type_number == 1)  // WW/CB
		gp_ele = ele_gp_value[ind + (long)pcs->m_msh->ele_vector.size()];

	// TF unused:  MshElemType::type eletyp = MeshElement->GetElementType();

	//
	upwind_para = pcs->m_num->ele_upwinding;
	upwind_meth = pcs->m_num->ele_upwind_method;

	// Numerik
	// *rupw = *supw = *tupw = 0.0;
	gp_r = gp_s = gp_t = 0;

	// get the velocities
	// gp_ele->GetEleVelocity(v);

	// CB: not sure if this is correct, as velocity
	// at each Gauss point is regarded here (within GP loop)
	// while in cel_mmp.cpp velocity is evaluated before the GP loop:
	// CalcVelo3Drst(phase, index, GetTimeCollocationupwind_MMP(), 0., 0., 0.,
	// v);
	// v[0] = gp_ele->Velocity(0, point);
	// v[1] = gp_ele->Velocity(1, point);
	// v[2] = gp_ele->Velocity(2, point);
	v[0] = gp_ele->Velocity(0, 0);
	v[1] = gp_ele->Velocity(1, 0);
	v[2] = gp_ele->Velocity(2, 0);
	// this would give v at GP 0
	// but: if(PcsType==T) SetCenterGP(); // CB 11/2007
	// was set in Cal_Velo_2(); which is calculated,
	// when mass matrix is assembled
	// hence V is at element center of gravity
	// otherwise use element averaged v?:
	// for(i=0; i<nGaussPoints; i++)
	//{
	//  v[0] += gp_ele->Velocity(0, i)/(double)nGaussPoints;
	//  v[1] += gp_ele->Velocity(1, i)/(double)nGaussPoints;
	//  v[2] += gp_ele->Velocity(2, i)/(double)nGaussPoints;
	//}

	for (size_t i = 0; i < 3; i++)
		v[i] *= time_unit_factor;

	// instead of r_upw, etc. we use unit[i], next function sets Gauss Integrals
	// unit[0] = MXPGaussPkt(nGauss, gp_r); -> r_upw ; etc. for unit[i]
	SetGaussPoint(point, gp_r, gp_s,
	              gp_t);  // this sets unit[] to standard coordinates

	// v transformation: Jacobi*v-->v_rst
	// computing the Jacobian at this point results in pressure difference
	// in comparison to the original model of CT
	// However, it seems not necessary, as Jacobian has already
	// been computed in function UpwindAlphaMass
	// computeJacobian(1); // order 1

	// multiply velocity vector with Jacobian matrix
	// Jacobi*v-->v_rst
	// This may need attention to tell different types of elements in the same
	// dimension	// PCH
	for (size_t i = 0; i < ele_dim; i++)
	{
		v_rst[i] = 0.0;
		for (size_t j = 0; j < ele_dim; j++)
			v_rst[i] += Jacobian[i * dim + j] * v[j];
	}
	//

	// These need to be rewritten according to different types of elements.  //
	// PCH
	if (MBtrgVec(v_rst, ele_dim) > MKleinsteZahl)
	{
		// Upwind-Faktoren
		for (size_t i = 0; i < ele_dim; i++)
			alpha[i] = -upwind_para * v_rst[i] /
			           (MBtrgVec(v_rst, ele_dim) + MKleinsteZahl);

		// moving the Gauss points
		if (upwind_meth == 1)  // limit Gauss point moving on Element domain
		{
			scale = 1.;
			for (size_t i = 0; i < ele_dim; i++)
			{
				// Integral over GaussPoints, not used
				if (mmp_element_integration_method_maximum)
				{
					if (fabs(unit[i] + alpha[i]) > 1.)
						scale =
						    MMin(scale, (1. - fabs(unit[i])) / fabs(alpha[i]));
				}
				else  // regard all quantities in the center of element
				    if (fabs(alpha[i]) > 1.)
					scale = MMin(scale, (1. / fabs(alpha[i])));
			}
			for (size_t i = 0; i < ele_dim; i++)
			{
				// Integral over GaussPoints, not used
				if (mmp_element_integration_method_maximum)
					unit[i] +=
					    scale *
					    alpha[i];  // scale is added to unit[i] (=Gaussintegral)
				else  // regard all quantities in the center of element
					unit[i] = scale * alpha[i];  // unit[i] (=Gaussintegral)
			}
		}
		else if (upwind_meth == 2)  // limit moving on -1<x<1

			// PCH this has never been used, but the code is only for line,
			// quad, and hex.
			for (size_t i = 0; i < ele_dim; i++)
			{
				// Integral ?er GaussPunkte
				if (mmp_element_integration_method_maximum)
					unit[i] = MRange(-1., unit[i] + alpha[i], 1.);
				else  // regard all quantities in the center of element
					unit[i] = MRange(-1., alpha[i], 1.);
			}
		// here only Methods 1 + 2; M3 = MaxMobilUW is done in CalcCoefLaplace
	}

#ifdef OLD_UPWINDING
	// test
	for (i = 0; i < ele_dim; i++)
		cout << unit[i] << " ";
	cout << endl;

	double ur, us, ut;

	switch (eletyp)
	{
		case 1:  // Line
		{
			// Elementgeometriedaten
			static double detjac, *invjac, jacobi[4];
			double l[3];
			// invjac = GetElementJacobiMatrix(index, &detjac);
			// Calc1DElementJacobiMatrix(ind, invjac, &detjac);  //ind = element
			// id number  // wird das irgendwo gebraucht?
			detjac = computeJacobian(1);  // order
			invjac = invJacobian;

			MNulleVec(l, 3);
			l[0] = X[1] - X[0];
			l[1] = Y[1] - Y[0];
			l[2] = Z[1] - Z[0];

			// hier nur Methoden 1 + 2; 3 wird in CalcCoefLaplace erledigt
			if ((upwind_meth == 1) || (upwind_meth == 2))
			{
				if (MBtrgVec(v, 3) > MKleinsteZahl)
				{
					if (MSkalarprodukt(v, l, 3) > 0.)
						// CB VZ ge?dert!!!
						*rupw = MRange(-1., -upwind_para, 1.);
					//*rupw = MRange(-1., upwind_para , 1.);
					else
						// CB VZ ge?dert!!!
						*rupw = MRange(-1., upwind_para, 1.);
					//*rupw = MRange(-1., -upwind_para , 1.); //
				}
				// else { // test
				//    cout << "-vau";
				//    if (MSkalarprodukt(v, l, 3) > 0.)
				//     *rupw = MRange(-1., upwind_para , 1.);
				//    else
				//     *rupw = MRange(-1., -upwind_para , 1.);
				//}
				if (aktueller_zeitschritt == 1)  // test
					*rupw = MRange(-1., upwind_para, 1.);
			}
			// Upwind-Faktor Fully upwinding
		}
		break;
		case 2:  // Quadrilateral
		{
			// Elementgeometriedaten
			static double detjac, *invjac, jacobi[4];
			// Elementdaten
			static double v_rs[2];
			// Initialisieren
			MNulleVec(v_rs, 2);

			// if (mmp_element_integration_method_maximum) ?? CB: was ist das
			if (1 > 0)
			{
				gp_r = (int)(point / nGauss);
				gp_s = point % nGauss;
				ur = MXPGaussPkt(nGauss, gp_r);
				us = MXPGaussPkt(nGauss, gp_s);
			}
			else
			{
				ur = 0.0;  // Alle Groessen nur in Elementmitte betrachten
				us = 0.0;  // Alle Groessen nur in Elementmitte betrachten
			}

			// Geschwindigkeitstransformation: a,b -> r,s
			// Calc2DElementJacobiMatrix(ind, 0., 0., invjac, &detjac);
			detjac = computeJacobian(1);  // order
			invjac = invJacobian;
			MKopierVec(invjac, jacobi, 4);
			M2Invertiere(jacobi);  // Jacobi-Matrix
			MMultMatVec(jacobi, 2, 2, v, 2, v_rs, 2);

			if (MBtrgVec(v_rs, 2) > MKleinsteZahl)
				// Upwind-Faktoren
				for (k = 0; k < 2; k++)
					alpha[k] = -upwind_para * v_rs[k] /
					           (MBtrgVec(v_rs, 2) + MKleinsteZahl);

			// hier nur Methoden 1 + 2; 3 wird in CalcCoefLaplace erledigt
			if (upwind_meth == 1)
			{
				// Verschiebungen der Gausspunkte auf Element begrenzen
				scale = 1.;
				if (fabs(ur + alpha[0]) > 1.)
					scale = MMin(scale, (1. - fabs(ur)) / fabs(alpha[0]));
				if (fabs(us + alpha[1]) > 1.)
					scale = MMin(scale, (1. - fabs(us)) / fabs(alpha[1]));
				*rupw = ur + scale* alpha[0];
				*supw = us + scale* alpha[1];
			}
			else if (upwind_meth == 2)
			{
				// Verschiebungen auf -1<x<1 begrenzen
				*rupw = MRange(-1., ur + alpha[0], 1.);
				*supw = MRange(-1., us + alpha[1], 1.);
			}
		}
		break;
		case 3:  // Hexahedra
		{
			/* Elementgeometriedaten */
			static double detjac, *invjac, jacobi[9];
			/* Elementdaten */
			static double v_rst[3];
			// Initialisieren
			MNulleVec(v_rst, 3);

			// if (mmp_element_integration_method_maximum) ?? CB: was ist das
			if (1 > 0)  // CB: to do ??
			{
				gp_r = (int)(point / (nGauss * nGauss));
				gp_s = (point % (nGauss * nGauss));
				gp_t = gp_s % nGauss;
				gp_s /= nGauss;
				ur = MXPGaussPkt(nGauss, gp_r);
				us = MXPGaussPkt(nGauss, gp_s);
				ut = MXPGaussPkt(nGauss, gp_t);
			}
			else
			{
				ur = 0.0;  // Alle Groessen nur in Elementmitte betrachten
				us = 0.0;  // Alle Groessen nur in Elementmitte betrachten
				ut = 0.0;  // Alle Groessen nur in Elementmitte betrachten
			}

			// Calc3DElementJacobiMatrix(ind, 0., 0., 0., invjac, &detjac);
			detjac = computeJacobian(1);  // order
			invjac = invJacobian;
			MKopierVec(invjac, jacobi, 9);
			M3Invertiere(jacobi); /* zurueck zur Jacobi-Matrix */
			MMultMatVec(jacobi, 3, 3, v, 3, v_rst, 3);

			if (MBtrgVec(v_rst, 3) > MKleinsteZahl) /* Upwind-Faktoren */
				for (l = 0; l < 3; l++)
					alpha[l] = -upwind_para * v_rst[l] / MBtrgVec(v_rst, 3) +
					           MKleinsteZahl;
			// hier nur Methoden 1 + 2; 3 wird in CalcCoefLaplace erledigt
			if (upwind_meth == 1)
			{
				// Verschiebungen der Gausspunkte auf Element begrenzen
				scale = 1.;
				if (fabs(ur + alpha[0]) > 1.)
					scale = MMin(scale, (1. - fabs(ur)) / fabs(alpha[0]));
				if (fabs(us + alpha[1]) > 1.)
					scale = MMin(scale, (1. - fabs(us)) / fabs(alpha[1]));
				if (fabs(ut + alpha[2]) > 1.)
					scale = MMin(scale, (1. - fabs(ut)) / fabs(alpha[2]));
				*rupw =
				    ur + scale* alpha[0];  // ist die reihenfolge hier richtig?
				*supw = us + scale* alpha[1];  // scale h?gt hier ja nur von dem
				                               // letzten if ab..
				*tupw = ut + scale* alpha[2];
			}
			else if (upwind_meth == 2)
			{
				// Verschiebungen auf -1<x<1 begrenzen
				*rupw = MRange(-1., ur + alpha[0], 1.);
				*supw = MRange(-1., us + alpha[1], 1.);
				*tupw = MRange(-1., ut + alpha[2], 1.);
			}
		}
		break;
		case 4:  // Triangle
			break;
		case 5:  // Tedrahedra
			break;
		case 6:  // Prism
			break;
	}

	// test
	for (i = 0; i < ele_dim; i++)
		cout << unit[i] << " ";
	cout << endl;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2009 NW
   last modification:
**************************************************************************/
double CFiniteElementStd::CalcSUPGCoefficient(double* vel, int ip,
                                              const double* diff_tensor)
{
	//--------------------------------------------------------------------
	// Collect following information to determine SUPG coefficient
	// + flow velocity
	// + diffusivity coefficient (scalar)
	// + characteristic element length
	// + (Peclet number)

	double v_mag = sqrt(vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2]);
	if (v_mag == 0.0) return 0.0;

	// Characteristic element length
	double ele_len = CalcSUPGEffectiveElemenetLength(vel);
	// Diffusivity = (effective heat conductivity) / (fluid heat capacity)
	const double* dispersion_tensor = diff_tensor;
	if (diff_tensor == NULL)
	{
		if (PcsType == H || PcsType == TH)  // heat
			dispersion_tensor =
			    MediaProp->HeatConductivityTensor(MeshElement->GetIndex());
		// mass
		else if (PcsType == M)
			// SB, BG
			dispersion_tensor = MediaProp->MassDispersionTensorNew(ip, 0);
	}

	double diff = .0;
	switch (pcs->m_num->ele_supg_method_diffusivity)
	{
		case 1:  // min
		{
			double min_diff = dispersion_tensor[0];
			for (size_t i = 1; i < dim * dim; i++)
				if (dispersion_tensor[i] < min_diff)
					min_diff = dispersion_tensor[i];
			diff = min_diff;
		}
		break;
		case 2:  // magnitude of diagonal
		{
			double tmp_diff = 0.0;
			for (size_t i = 0; i < dim; i++)
				tmp_diff = dispersion_tensor[i + i * dim] *
				           dispersion_tensor[i + i * dim];
			diff = sqrt(tmp_diff);
		}
		break;
		default:  // 0 or any invalid number: max. in dispersion coefficient
		{
			double max_diff = dispersion_tensor[0];
			for (size_t i = 1; i < dim * dim; i++)
				if (dispersion_tensor[i] > max_diff)
					max_diff = dispersion_tensor[i];
			diff = max_diff;
		}
		break;
	}

	if (diff_tensor == NULL)
	{
		if (PcsType == H || PcsType == TH)  // heat
		{
			if (FluidProp->density_model == 14)
			{
				double dens_arg[3];  // AKS
				int Index = MeshElement->GetIndex();
				dens_arg[0] = interpolate(NodalValC1);
				dens_arg[1] = interpolate(NodalVal1) + T_KILVIN_ZERO;
				dens_arg[2] = Index;
				diff /= (FluidProp->SpecificHeatCapacity(dens_arg) *
				         FluidProp->Density(dens_arg));
			}
			else
				diff /=
				    (FluidProp->SpecificHeatCapacity() * FluidProp->Density());
		}
	}

	//--------------------------------------------------------------------
	// Here calculates SUPG coefficient (tau)
	double tau = 0.0;
	switch (pcs->m_num->ele_supg_method)
	{
		case 1:
		{
			// this coefficient matches with the analytical solution in 1D
			// steady state case
			double alpha = 0.5 * v_mag * ele_len / diff;  // 0.5*Pe
			double func = MLangevin(alpha);
			tau = 0.5 * ele_len / v_mag * func;
		}
		break;
		case 2:
		{
			// taking into account time step
			//          tau = 1.0 / sqrt(pow(2.0/dt
			//          ,2.0)+pow(2.0*v_mag/ele_len,2.0));
			tau = 1.0 / sqrt((2.0 / dt) * (2.0 / dt) +
			                 (2.0 * v_mag / ele_len) * (2.0 * v_mag / ele_len) +
			                 (4.0 * diff / (ele_len * ele_len)) *
			                     (4.0 * diff / (ele_len * ele_len)));
		}
		break;
	}

	return tau;
}
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2009 NW
   last modification:
**************************************************************************/
void CFiniteElementStd::CalcSUPGWeightingFunction(double* vel, int ip,
                                                  double& tau, double* v_dN)
{
	if (pcs->m_num->ele_supg_method == 0)
	{
		cout << "***Warning in CFiniteElementStd::CalcSUPGWeightingFunction(): "
		        "SUPG option is not selected" << endl;
		return;
	}

	// tau
	tau = CalcSUPGCoefficient(vel, ip);

	// {v}[dN]
	for (int i = 0; i < nnodes; i++)
		v_dN[i] = 0.0;
	for (int i = 0; i < nnodes; i++)
		for (size_t k = 0; k < dim; k++)
			v_dN[i] += dshapefct[k * nnodes + i] * vel[k];
}
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   12/2009 NW
   last modification:
**************************************************************************/
double CFiniteElementStd::CalcSUPGEffectiveElemenetLength(double* /*vel*/)
{
	double L = 0.0;
	switch (this->ele_dim)
	{
		case 1:
		{
			L = this->MeshElement->GetVolume();
		}
		break;
		case 2:
		case 3:
		{
			switch (pcs->m_num->ele_supg_method_length)
			{
				case 1:  // min
				{
					double min = MeshElement->GetEdge(0)->getLength();
					for (size_t i = 1; i < MeshElement->GetEdgesNumber(); i++)
					{
						L = MeshElement->GetEdge(i)->getLength();
						if (L < min) min = L;
					}
					L = min;
				}
				break;
				case 2:  // average
				{
					double tmp_L = 0.0;
					for (size_t i = 1; i < MeshElement->GetEdgesNumber(); i++)
						tmp_L += MeshElement->GetEdge(i)->getLength();
					L = tmp_L / MeshElement->GetEdgesNumber();
				}
				break;
				case 3:  // stream line length
				{
					cout << "***Error: ele_supg_method_length <3> has not been "
					        "supported yet." << endl;
				}
				break;
				default:  // 0 or any invalid number: max edge length
				{
					double max = MeshElement->GetEdge(0)->getLength();
					for (size_t i = 1; i < MeshElement->GetEdgesNumber(); i++)
					{
						L = MeshElement->GetEdge(i)->getLength();
						if (L > max) max = L;
					}
					L = max;
				}
				break;
			}
		}
		break;
	}
	return L;
}

// CB 090507
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 CB
   last modification:
**************************************************************************/
void CFiniteElementStd::UpwindAlphaMass(double* alpha)
{
	// Laufvariablen
	static long i;
	// WW int no_phases;

	// static long *element_nodes;
	// WW double gp[3],
	double v_rst[3], v_tot[3] = {};
	// WW static double zeta;
	// static double *velovec, vg, v[2], vt[2], v_rs[2];
	// static double alpha_adv[3];
	double upwind_para;
	// WW double upwind_meth;
	// Numerik
	// WW zeta = 0.0;
	// WW gp[0]=0.0;   gp[1]=0.0;   gp[2]=0.0;

	int ind = MeshElement->GetIndex();
	ElementValue* gp_ele = ele_gp_value[ind];
	//  TF ununsed: MshElemType::type eletyp = MeshElement->GetElementType();

	// Elementdaten und globale Modellparameter
	// WW no_phases =(int)mfp_vector.size();
	//
	upwind_para = pcs->m_num->ele_upwinding;
	// WW upwind_meth = pcs->m_num->ele_upwind_method;

	// alpha initialisieren
	for (int i=0; i<3; i++)
		alpha[i] = 0;

	// get the velocities for phase 0
	v_tot[0] = gp_ele->Velocity(0, 0);
	v_tot[1] = gp_ele->Velocity(1, 0);
	v_tot[2] = gp_ele->Velocity(2, 0);
	// this would only give v at GP 0
	// but: if(PcsType==T) SetCenterGP(); // CB 11/2007
	// was set in Cal_Velo_2(); which is calculated,
	// when mass matrix is assembled
	// hence v is at element center of gravity
	// otherwise use following approximation:
	// for(i=0; i<nGaussPoints; i++) //element averaged v?
	//{
	//  v_tot[0] += gp_ele->Velocity(0, i)/nGaussPoints;
	//  v_tot[1] += gp_ele->Velocity(1, i)/nGaussPoints;
	//  v_tot[2] += gp_ele->Velocity(2, i)/nGaussPoints;
	//}

	// switch to next phase
	gp_ele = ele_gp_value[ind + (long)pcs->m_msh->ele_vector.size()];
	// get the velocities for phases 1 and add
	v_tot[0] += gp_ele->Velocity(0, 0);
	v_tot[1] += gp_ele->Velocity(1, 0);
	v_tot[2] += gp_ele->Velocity(2, 0);
	// for(i=0; i<nGaussPoints; i++) //element averaged v?
	//{
	//  v_tot[0] += gp_ele->Velocity(0, i)/nGaussPoints;
	//  v_tot[1] += gp_ele->Velocity(1, i)/nGaussPoints;
	//  v_tot[2] += gp_ele->Velocity(2, i)/nGaussPoints;
	//}
	for (i = 0; i < 3; i++)
		v_tot[i] *= time_unit_factor;

	// SetGaussPoint(point, gp_r, gp_s, gp_t);
	// velocity transformation a,b,c -> r,s,t
	computeJacobian(1);  // order 1
	// multiply velocity vector with Jacobian matrix
	// Jacobi*v-->v_rst
	for (size_t i = 0; i < ele_dim; i++)
	{
		v_rst[i] = 0.0;
		for (size_t j = 0; j < ele_dim; j++)
			v_rst[i] += Jacobian[i * dim + j] * v_tot[j];
	}

	// Upwind-Factors
	if (MBtrgVec(v_rst, ele_dim) >
	    MKleinsteZahl)  // if(lengthOftheVector > tolerance)

		for (size_t i = 0; i < ele_dim; i++)
			alpha[i] = -upwind_para * v_rst[i] /
			           (MBtrgVec(v_rst, ele_dim) + MKleinsteZahl);

#ifdef OLD_UPWINDING

	// test
	for (i = 0; i < ele_dim; i++)
		cout << alpha[i] << " ";
	cout << endl;

	switch (eletyp)
	{
		case 1:  // Line
		{
			// Elementgeometriedaten
			static double detjac, *invjac, jacobi[4];
			static double l[3];
			// invjac = GetElementJacobiMatrix(index, &detjac);
			// Calc1DElementJacobiMatrix(ind, invjac, &detjac);  //index =
			// element id number
			detjac = computeJacobian(1);  // order
			invjac = invJacobian;
			// element_nodes = ElGetElementNodes(ind);

			MNulleVec(l, 3);
			l[0] = X[1] - X[0];
			l[1] = Y[1] - Y[0];
			l[2] = Z[1] - Z[0];

			if (MBtrgVec(v_tot, 3) > MKleinsteZahl)
			{
				if (MSkalarprodukt(v_tot, l, 3) > 0.)
					zeta = 1.;  // upwind_para
				else
					zeta = -1.;  //-upwind_para
			}

			// aus RF 3.5.06 CT 1D elements: {
			//// detjac = A*L/2
			// vorfk = porosity * detjac * Mdrittel;
			//// Massenmatrix mit SUPG ohne Zeitanteile
			// mass[0] = (2.0 + 1.5 * mms_upwind_parameter * zeta) * vorfk;
			// mass[1] = (1.0 + 1.5 * mms_upwind_parameter * zeta) * vorfk;
			// mass[2] = (1.0 - 1.5 * mms_upwind_parameter * zeta) * vorfk;
			// mass[3] = (2.0 - 1.5 * mms_upwind_parameter * zeta) * vorfk; // }

			// Upwind-Faktor Fully upwinding
			// alpha[0]     = m_pcs->m_num->ele_upwinding * zeta;
			// alpha_adv[0] = m_pcs->m_num->ele_upwinding * zeta;
			alpha[0] = 1.0 * zeta;  //??
			                        // alpha_adv[0] = 1.0 * zeta;
			                        // Advection upwinding
			// if (MTM2_upwind_method == 2) alpha_adv[0] = ele_upwinding * zeta;
			// /
		}
		break;
		case 2:  // Quadrilateral
		{
			// Elementgeometriedaten
			static double detjac, *invjac, jacobi[4];
			// Elementdaten
			static double v_rs[3];

			// Geschwindigkeitstransformation: a,b -> r,s
			// Calc2DElementJacobiMatrix(ind, 0., 0., invjac, &detjac);
			detjac = computeJacobian(1);  // order
			invjac = invJacobian;
			MKopierVec(invjac, jacobi, 4);
			M2Invertiere(jacobi); /* Jacobi-Matrix */
			MMultMatVec(jacobi, 2, 2, v_tot, 2, v_rs, 2);

			if (MBtrgVec(v_rs, 2) > MKleinsteZahl)
				// Upwind-Faktoren
				for (k = 0; k < 2; k++)
					alpha[k] = -upwind_para * v_rs[k] /
					           (MBtrgVec(v_rs, 2) + MKleinsteZahl);
		}
		break;
		case 3:  // Hexahedra
		{
			/* Elementgeometriedaten */
			static double* invjac, jacobi[9], detjac;
			/* Elementdaten */
			// static double v_rst[3];

			if (MBtrgVec(v_tot, 3) > MKleinsteZahl)
			{
				/* Geschwindigkeitstransformation: x,y,z -> r,s,t */
				// Calc3DElementJacobiMatrix(ind, 0., 0., 0., invjac, &detjac);
				detjac = computeJacobian(1);  // order
				invjac = invJacobian;
				MKopierVec(invjac, jacobi, 9);
				M3Invertiere(jacobi); /* Jacobi-Matrix */
				MMultMatVec(jacobi, 3, 3, v_tot, 3, v_rst, 3);

				/* Upwind-Faktoren */
				for (l = 0; l < 3; l++)
					alpha[l] = -upwind_para * v_rst[l] /
					           (MBtrgVec(v_rst, 3) + MKleinsteZahl);
			}
		}
		break;
		case 4:  // Triangle
			break;
		case 5:  // Tedrahedra
			break;
		case 6:  // Prism
			break;
	}
	// test
	for (i = 0; i < ele_dim; i++)
		cout << alpha[i] << " ";
	cout << endl;
#endif
}

// CB 160507
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2007 CB
   last modification:
**************************************************************************/
void CFiniteElementStd::UpwindSummandMass(const int gp,
                                          int& gp_r,
                                          int& gp_s,
                                          int& gp_t,
                                          double* alpha,
                                          double* summand)

{
	int i;
	//
	GetGaussData(gp, gp_r, gp_s, gp_t);  // this sets unit[] to standard values
	GradShapeFunction(dshapefct, unit);
	for (i = 0; i < nnodes; i++)
	{
		summand[i] = 0.0;
		for (size_t k = 0; k < dim; k++)
			summand[i] += dshapefct[nnodes * k + i] * alpha[k];
		// summand[i] /= (double)nnodes;
	}

#ifdef OLD_UPWINDING

	double u1, u2, u3;
	u1 = u2 = u3 = 0;

	MshElemType::type eletyp = MeshElement->GetElementType();
	switch (eletyp)
	{
		case MshElemType::LINE:
		{
			// Line
			gp_r = gp;
			u1 = MXPGaussPkt(nGauss, gp_r);
			summand[0] = +alpha[0] * (1 + u1);  // CB: ?? hab ich mir so gedacht
			summand[1] = -alpha[0] * (1 - u1);  // CB: ?? hab ich mir so gedacht
			for (i = 0; i < 2; i++)
				summand[i] *= 0.5;
		}
		break;
		case MshElemType::QUAD:  // Quadrilateral
		{
			gp_r = (int)(gp / nGauss);
			gp_s = gp % nGauss;
			u1 = MXPGaussPkt(nGauss, gp_r);
			u2 = MXPGaussPkt(nGauss, gp_s);
			// derived from MPhi2D_SUPG
			summand[0] = +alpha[0] * (1 + u2) + alpha[1] * (1 + u1);
			summand[1] = -alpha[0] * (1 + u2) + alpha[1] * (1 - u1);
			summand[2] = -alpha[0] * (1 - u2) - alpha[1] * (1 - u1);
			summand[3] = +alpha[0] * (1 - u2) - alpha[1] * (1 + u1);
			for (i = 0; i < 4; i++)
				summand[i] *= 0.25;
		}
		break;
		case MshElemType::HEXAHEDRON:  // Hexahedra
		{
			gp_r = (int)(gp / (nGauss * nGauss));
			gp_s = (gp % (nGauss * nGauss));
			gp_t = gp_s % nGauss;
			gp_s /= nGauss;
			u1 = MXPGaussPkt(nGauss, gp_r);
			u2 = MXPGaussPkt(nGauss, gp_s);
			u3 = MXPGaussPkt(nGauss, gp_t);
			// derived from MPhi3D_SUPG
			summand[0] = +alpha[0] * (1 + u2) * (1 + u3) +
			             alpha[1] * (1 + u1) * (1 + u3) +
			             alpha[2] * (1 + u1) * (1 + u2);
			summand[1] = -alpha[0] * (1 + u2) * (1 + u3) +
			             alpha[1] * (1 - u1) * (1 + u3) +
			             alpha[2] * (1 - u1) * (1 + u2);
			summand[2] = -alpha[0] * (1 - u2) * (1 + u3) -
			             alpha[1] * (1 - u1) * (1 + u3) +
			             alpha[2] * (1 - u1) * (1 - u2);
			summand[3] = +alpha[0] * (1 - u2) * (1 + u3) -
			             alpha[1] * (1 + u1) * (1 + u3) +
			             alpha[2] * (1 + u1) * (1 - u2);
			summand[4] = +alpha[0] * (1 + u2) * (1 - u3) +
			             alpha[1] * (1 + u1) * (1 - u3) -
			             alpha[2] * (1 + u1) * (1 + u2);
			summand[5] = -alpha[0] * (1 + u2) * (1 - u3) +
			             alpha[1] * (1 - u1) * (1 - u3) -
			             alpha[2] * (1 - u1) * (1 + u2);
			summand[6] = -alpha[0] * (1 - u2) * (1 - u3) -
			             alpha[1] * (1 - u1) * (1 - u3) -
			             alpha[2] * (1 - u1) * (1 - u2);
			summand[7] = +alpha[0] * (1 - u2) * (1 - u3) -
			             alpha[1] * (1 + u1) * (1 - u3) -
			             alpha[2] * (1 + u1) * (1 - u2);
			for (i = 0; i < 8; i++)
				summand[i] *= 0.125;
		}
		break;
		case MshElemType::TRIANGLE:  // Triangle
		{
			// SamplePointTriHQ(gp, unit);
		}
		break;
		case MshElemType::TETRAHEDRON:  // Tedrahedra
		{
			// SamplePointTet5(gp, unit);
		}
		break;
		case MshElemType::PRISM:  // Prism
		{
			gp_r = gp % nGauss;
			gp_s = (int)(gp / nGauss);
			gp_t = (int)(nGaussPoints / nGauss);
			// u1 = MXPGaussPktTri(nGauss,gp_r,0); //femlib.cpp statt
			// mathlib.cpp, nicht verfügbar?
			// u2 = MXPGaussPktTri(nGauss,gp_r,1);
			// u3 = MXPGaussPkt(gp_t,gp_s);
		}
		break;
	}
#endif
}

}  // end namespace

