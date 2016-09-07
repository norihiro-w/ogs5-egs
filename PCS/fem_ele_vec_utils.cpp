/**
 * \copyright
 * Copyright (c) 2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_vec.h"

#include <valarray>

#include "makros.h"

#include "mathlib.h"
#include "matrix_class.h"

#include "ElementValueDM.h"
#include "rf_pcs_dm.h"

namespace  FiniteElement
{

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: setB_Matrix(const int LocalIndex)

   Aufgabe:
          Form B matric
   Formalparameter:
           E:
             const int LocalIndex   : Local node index

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::setB_Matrix(const int LocalIndex)
{
	switch (dim)
	{
		case 2:
			// B_11, dN/dx
			(*B_matrix)(0, 0) = dshapefctHQ[LocalIndex];
			// B_12, 0.0
			(*B_matrix)(0, 1) = 0.0;

			if (axisymmetry)  // Axisymmtry
			{
				// B_21, N/r
				(*B_matrix)(1, 0) = shapefctHQ[LocalIndex] / Radius;
				// B_22, 0.0
				(*B_matrix)(1, 1) = 0.0;
				// B_31, 0.0
				(*B_matrix)(2, 0) = 0.0;
				// B_32, dN/dz
				(*B_matrix)(2, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
			}
			else
			{
				// B_21, 0.0
				(*B_matrix)(1, 0) = 0.0;
				// B_22, dN/dy
				(*B_matrix)(1, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
				// B_31, 0.0
				(*B_matrix)(2, 0) = 0.0;
				// B_32, 0.0
				(*B_matrix)(2, 1) = 0.0;
			}
			// B_41, dN/dy
			(*B_matrix)(3, 0) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_42, dN/dx
			(*B_matrix)(3, 1) = dshapefctHQ[LocalIndex];

			break;
		case 3:
			// B_11, dN/dx
			(*B_matrix)(0, 0) = dshapefctHQ[LocalIndex];
			// B_22, dN/dy
			(*B_matrix)(1, 1) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_33, dN/dz
			(*B_matrix)(2, 2) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			//
			// B_41, dN/dy
			(*B_matrix)(3, 0) = dshapefctHQ[nnodesHQ + LocalIndex];
			// B_42, dN/dx
			(*B_matrix)(3, 1) = dshapefctHQ[LocalIndex];
			//
			// B_51, dN/dz
			(*B_matrix)(4, 0) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			// B_53, dN/dx
			(*B_matrix)(4, 2) = dshapefctHQ[LocalIndex];
			//
			// B_62, dN/dz
			(*B_matrix)(5, 1) = dshapefctHQ[2 * nnodesHQ + LocalIndex];
			// B_63, dN/dy
			(*B_matrix)(5, 2) = dshapefctHQ[nnodesHQ + LocalIndex];

			break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: setTransB_Matrix(const int LocalIndex)

   Aufgabe:
          Form the tanspose of B matric
   Formalparameter:
           E:
             const int LocalIndex   : Local node index

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::setTransB_Matrix(const int LocalIndex)
{
	setB_Matrix(LocalIndex);
	B_matrix->GetTranspose(*B_matrix_T);
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ComputeStrain()

   Aufgabe:
          Compute strains

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::ComputeStrain()
{
	int i, j = 0, k = 0;
	switch (dim)
	{
		case 2:
			dstrain = 0.0;
			if (axisymmetry)
			{
				for (i = 0; i < nnodesHQ; i++)
				{
					j = i + nnodesHQ;
					dstrain[0] += Disp[i] * dshapefctHQ[i];
					dstrain[1] += Disp[i] * shapefctHQ[i];
					dstrain[2] += Disp[j] * dshapefctHQ[j];
					dstrain[3] +=
					    Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				}
				dstrain[1] /= Radius;
			}
			else
				for (i = 0; i < nnodesHQ; i++)
				{
					j = i + nnodesHQ;
					dstrain[0] += Disp[i] * dshapefctHQ[i];
					dstrain[1] += Disp[j] * dshapefctHQ[j];
					dstrain[3] +=
					    Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				}
			break;
		case 3:
			dstrain = 0.0;
			for (i = 0; i < nnodesHQ; i++)
			{
				j = i + nnodesHQ;
				k = i + 2 * nnodesHQ;
				dstrain[0] += Disp[i] * dshapefctHQ[i];
				dstrain[1] += Disp[j] * dshapefctHQ[j];
				dstrain[2] += Disp[k] * dshapefctHQ[k];
				dstrain[3] +=
				    Disp[i] * dshapefctHQ[j] + Disp[j] * dshapefctHQ[i];
				dstrain[4] +=
				    Disp[i] * dshapefctHQ[k] + Disp[k] * dshapefctHQ[i];
				dstrain[5] +=
				    Disp[j] * dshapefctHQ[k] + Disp[k] * dshapefctHQ[j];
			}
			break;
	}
}

/***************************************************************************
   GeoSys - Funktion:
      CFiniteElementVec::ComputePrincipleStresses(double *dEStress)
      (2D only)
   Aufgabe:

   Formalparameter:
         E:
   const double *Stresses: Stresses

   Return: Angle of maxium principle stress component to x direction
   sig2<sig1

   Programming:
   06/2004     WW        Erste Version
 **************************************************************************/
double CFiniteElementVec::ComputePrincipleStresses(const double* Stresses, double* pr_stress)
{
	double prin_ang, var;
	// Angle of the principle plane
	if (fabs(Stresses[0] - Stresses[1]) < MKleinsteZahl) return 0.0;
	prin_ang = atan(2.0 * Stresses[3] / (Stresses[0] - Stresses[1]));
	// Principle stress 1
	pr_stress[0] = 0.5 * (Stresses[0] + Stresses[1]) +
	               0.5 * (Stresses[0] - Stresses[1]) * cos(prin_ang) +
	               Stresses[3] * sin(prin_ang);
	// Principle stress 2
	pr_stress[1] = 0.5 * (Stresses[0] + Stresses[1]) -
	               0.5 * (Stresses[0] - Stresses[1]) * cos(prin_ang) -
	               Stresses[3] * sin(prin_ang);
	pr_stress[2] = Stresses[2];

	prin_ang *= 0.5;

	// Jump direction
	if (pr_stress[1] >= pr_stress[0])
	{
		prin_ang += 0.5 * pai;
		var = pr_stress[0];
		pr_stress[0] = pr_stress[1];
		pr_stress[1] = var;
	}
	// if(pr_stress[0]<0.0) prin_ang += pai;
	return prin_ang;
}

/***************************************************************************
   GeoSys - Funktion:
            CFiniteElementVec:: CalcStrain_v()
   Aufgabe:
           Calculate effictive strain at Gauss points
   Formalparameter:
           E:

   Programming:
   01/2009   WW/UWG
 **************************************************************************/
double CFiniteElementVec::CalcStrain_v()
{
	int i;
	for (i = 0; i < ns; i++)
		dstrain[i] = 0.0;
	//
	for (i = 0; i < nnodesHQ; i++)
		dstrain[i] +=
		    pcs->GetNodeValue(nodes[i], Idx_Strain[0]) * shapefctHQ[i];
	double val = 0;
	for (i = 0; i < 3; i++)
		val += dstrain[i] * dstrain[i];
	for (i = 3; i < ns; i++)
		val += 0.5 * dstrain[i] * dstrain[i];

	return sqrt(2.0 * val / 3.);
}
/***************************************************************************
   GeoSys - Funktion:
            CFiniteElementVec:: CalcStress_eff()
   Aufgabe:
           Calculate effective stress at Gauss points
   Formalparameter:
           E:

   Programming:
   06/2014   AJ
 **************************************************************************/
double CFiniteElementVec::CalcStress_eff()
{
	size_t j;
	double val = 0.;
	ElementValue_DM* ele_value = ele_value_dm[MeshElement->GetIndex()];

	for (j = 0; j < dim; j++)
		val -= (*ele_value->Stress)(j, gp);

	val /= dim;

	return val;
}
} // namespace
