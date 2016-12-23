/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "fem_ele_vec.h"

#include <cfloat>

#include "Curve.h"
#if defined(NEW_EQS)
#include "equation_class.h"
#endif
#include "mathlib.h"
#if defined(USE_PETSC)
#include "PETSC/PETScLinearSolver.h"
#endif

#include "ElementMatrix.h"
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
#include "rf_pcs_dm.h"
#include "tools.h"

//
#define COMP_MOL_MASS_AIR 28.96  // kg/kmol WW  28.96
#define GAS_CONSTANT 8314.41     // J/(kmol*K) WW

namespace FiniteElement
{

CFiniteElementVec::CFiniteElementVec(process::CRFProcessDeformation* dm_pcs,
                                     const int C_Sys_Flad,
                                     const int order)
    : CElement(C_Sys_Flad, order),
	  pcs(dm_pcs)
{
	Tem = 273.15 + 23.0;
	dim = pcs->m_msh->GetMaxElementDim();  // overwrite dim in CElement
	ns = 4;
	if (dim == 3) ns = 6;

	AuxNodal.resize(8);
	AuxNodal1.resize(60);

	Idx_Stress = new int[ns];
	Idx_Strain = new int[ns];
	strain_ne.resize(ns);
	stress0.resize(ns);
	for (int i = 0; i < 4; i++)
		NodeShift[i] = pcs->Shift[i];

	// Indecex in nodal value table
	Idx_dm0[0] = pcs->GetNodeValueIndex("DISPLACEMENT_X1");
	Idx_dm0[1] = pcs->GetNodeValueIndex("DISPLACEMENT_Y1");
	Idx_dm1[0] = Idx_dm0[0] + 1;
	Idx_dm1[1] = Idx_dm0[1] + 1;
	if (dim == 3)
	{
		Idx_dm0[2] = pcs->GetNodeValueIndex("DISPLACEMENT_Z1");
		Idx_dm1[2] = Idx_dm0[2] + 1;
	}

	// idx_pls =  pcs->GetNodeValueIndex("STRAIN_PLS");
	// Strain
	Idx_Strain[0] = pcs->GetNodeValueIndex("STRAIN_XX");
	Idx_Strain[1] = pcs->GetNodeValueIndex("STRAIN_YY");
	Idx_Strain[2] = pcs->GetNodeValueIndex("STRAIN_ZZ");
	Idx_Strain[3] = pcs->GetNodeValueIndex("STRAIN_XY");
	// Stress
	Idx_Stress[0] = pcs->GetNodeValueIndex("STRESS_XX");
	Idx_Stress[1] = pcs->GetNodeValueIndex("STRESS_YY");
	Idx_Stress[2] = pcs->GetNodeValueIndex("STRESS_ZZ");
	Idx_Stress[3] = pcs->GetNodeValueIndex("STRESS_XY");

	// For cache
	vec_B_matrix.resize(20);
	vec_B_matrix_T.resize(20);
	for (int i = 0; i < (int)vec_B_matrix.size(); i++)
	{
		vec_B_matrix[i] = new Matrix(ns, dim);
		vec_B_matrix_T[i] = new Matrix(dim, ns);
	}

	//
	switch (dim)
	{
		case 2:
			Disp = new double[18];
			dT = new double[9];
			T1 = new double[9];

			Sxx = new double[9];
			Syy = new double[9];
			Szz = new double[9];
			Sxy = new double[9];
			pstr = new double[9];

			break;
		case 3:
			Disp = new double[60];
			dT = new double[20];
			T1 = new double[20];

			Sxx = new double[20];
			Syy = new double[20];
			Szz = new double[20];
			Sxy = new double[20];
			Sxz = new double[20];
			Syz = new double[20];
			pstr = new double[20];

			// Indecex in nodal value table
			Idx_Strain[4] = pcs->GetNodeValueIndex("STRAIN_XZ");
			Idx_Strain[5] = pcs->GetNodeValueIndex("STRAIN_YZ");
			//
			Idx_Stress[4] = pcs->GetNodeValueIndex("STRESS_XZ");
			Idx_Stress[5] = pcs->GetNodeValueIndex("STRESS_YZ");
			break;
	}

	dstress.resize(ns);
	dstrain.resize(ns);
	B_matrix = new Matrix(ns, dim);
	B_matrix_T = new Matrix(dim, ns);
	De = new Matrix(ns, ns);
	ConsistDep = new Matrix(ns, ns);
	AuxMatrix = new Matrix(dim, dim);
	AuxMatrix2 = new Matrix(dim, ns);

	*B_matrix = 0.0;
	*B_matrix_T = 0.0;

	if (pcs->Memory_Type == 0)  // Do not store local matrices
	{
		Stiffness = new Matrix(60, 60);
		RHS = new Vector(60);
		if (H_Process)
			PressureC = new Matrix(60, 20);
	}

	// Coupling
	idx_P0 = idx_pls = 0;
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		if (isFlowProcess(pcs_vector[i]->getProcessType()))
		{
			h_pcs = pcs_vector[i];
			if (h_pcs->type == 1 || h_pcs->type == 41)
			{
				if (h_pcs->getProcessType() == GROUNDWATER_FLOW)
					Flow_Type = 10;
				else
					Flow_Type = 0;
			}
			break;
		}
	}
	if (Flow_Type == 0)
	{
		idx_P1 = h_pcs->GetNodeValueIndex("PRESSURE1") + 1;
	}
	if (Flow_Type == 10)
		idx_P1 = h_pcs->GetNodeValueIndex("HEAD") + 1;

	if (T_Flag)
	{
		for (size_t i = 0; i < pcs_vector.size(); i++)
		{
			if (pcs_vector[i]->getProcessType() == HEAT_TRANSPORT)
			{
				t_pcs = pcs_vector[i];
				break;
			}
		}
		idx_T0 = t_pcs->GetNodeValueIndex("TEMPERATURE1");
		idx_T1 = idx_T0 + 1;
	}

	time_unit_factor = pcs->time_unit_factor;
#if defined(USE_PETSC)
	size_t size = 60;      // dim * nnodesHQ;
	row_ids = new int[size];  //> global indices of local matrix rows
	col_ids = new int[size];  //> global indices of local matrix columns
#endif
}


CFiniteElementVec::~CFiniteElementVec()
{
	delete B_matrix;
	delete B_matrix_T;
	delete De;
	delete ConsistDep;
	delete AuxMatrix;
	delete AuxMatrix2;
	delete[] Disp;
	delete[] dT;
	delete[] T1;
	delete[] Sxx;
	delete[] Syy;
	delete[] Szz;
	delete[] Sxy;
	delete[] pstr;
	delete[] Idx_Strain;
	delete[] Idx_Stress;
	delete[] Sxz;
	delete[] Syz;

	if (pcs->Memory_Type == 0)  // Do not store local matrices
	{
		delete Stiffness;
		delete RHS;
	}

	delete PressureC;

	for (int i = 0; i < (int)vec_B_matrix.size(); i++)
	{
		delete vec_B_matrix[i];
		delete vec_B_matrix_T[i];
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::SetMaterial(const int EleIndex)

   Aufgabe:
         Set material data for local assembly
   Formalparameter:
           E:

   Programming:
   11/2004     WW        Erste Version
 **************************************************************************/
void CFiniteElementVec::SetMaterial()
{
	int MatGroup = MeshElement->GetPatchIndex();

	m_msp = msp_vector[MatGroup];
	m_msp->axisymmetry = pcs->m_msh->isAxisymmetry();
	if (m_msp->Plasticity_type == 2)
		m_msp->ResizeMatricesSYS(ele_dim);

	if (F_Flag)
	{
		m_mfp = MFPGet("LIQUID");
		if (!m_mfp) m_mfp = mfp_vector[0];
	}
	m_mmp = mmp_vector[MatGroup];
}

void CFiniteElementVec::SetMemory()
{
	int size = 0;
	ElementMatrix* EleMat = NULL;

	// Prepare local matrices
	if (pcs->Memory_Type == 0)
	{
		// If local matrices are not stored, resize the matrix
		size = dim * nnodesHQ;
		Stiffness->LimitSize(size, size);
		if (PressureC) PressureC->LimitSize(size, nnodes);
		RHS->LimitSize(size);
	}
	else
	{
		EleMat = pcs->Ele_Matrices[Index];
		Stiffness = EleMat->GetStiffness();
		RHS = EleMat->GetRHS();
		if (PressureC) PressureC = EleMat->GetCouplingMatrixA();
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: CalDensity()

   Aufgabe:
          Caculate density of porous medium
   Programming:
   05/2005     WW        Erste Version
 **************************************************************************/
double CFiniteElementVec::CalDensity()
{
	double const solid_density = m_msp->Density();
	// If negative value is given in the .msp file, gravity by solid is
	// skipped
	if (solid_density <= 0.0)
		return 0.0;

	if (!F_Flag)
		return solid_density;

	double porosity = m_mmp->Porosity(this);
	double rho = (1. - porosity) * solid_density;

	// JT, should be 1.0, unless multiphase (calculate below)
	// (if unsaturated, fluid density would be negligible...
	// so still works)
	double Sw = 1.0;

	int const no_phases = (int)mfp_vector.size();
	int phase = 0; //TDODO
	if (no_phases > 0 && no_phases > phase)
	{
		rho += porosity * Sw * m_mfp->Density();
	}

	if (Flow_Type == 2 || Flow_Type == 3)
	{
		CFluidProperties* GasProp;
		GasProp = MFPGet("GAS");
		rho += porosity * (1.0 - Sw) * GasProp->Density();
	}

	return rho;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: ComputeMatrix_RHS(const double fkt)
   Aufgabe:
           as the name

   Programming:
   07/2004   WW
   08/2004   OK   MFP implementation
   01/2010   NW   use chache of B,B^T matrices
 **************************************************************************/
void CFiniteElementVec::ComputeMatrix_RHS(const double fkt, const Matrix* p_D)
{
	int i, j, k, l;
	double rho, fac, fac1, fac2, dN_dx, f_buff;
	fac = fac1 = fac2 = f_buff = 0.0;
	dN_dx = 0.0;
	rho = CalDensity();
	const int nnodesHQ = this->nnodesHQ;
	const int nnodes = this->nnodes;
	const int ele_dim = this->ele_dim;
	const int ns = this->ns;

	// NW cache B, B^T
	for (i = 0; i < nnodesHQ; i++)
	{
		setTransB_Matrix(i);
		(*this->vec_B_matrix[i]) = *B_matrix;
		(*this->vec_B_matrix_T[i]) = *B_matrix_T;
	}
	Matrix* old_B_matrix = B_matrix;
	Matrix* old_B_matrix_T = B_matrix_T;

	Matrix* tmp_B_matrix = NULL;
	Matrix* tmp_B_matrix_T = NULL;
	Matrix* tmp_AuxMatrix = AuxMatrix;
	Matrix* tmp_AuxMatrix2 = AuxMatrix2;
	Matrix* tmp_Stiffness = Stiffness;

	for (i = 0; i < nnodesHQ; i++)
	{
		tmp_B_matrix_T = this->vec_B_matrix_T[i];

		// Local assembly of A*u=int(B^t*sigma) for Newton-Raphson method
		for (j = 0; j < ele_dim; j++) {
			for (k = 0; k < ns; k++) {
				(*RHS)(j* nnodesHQ + i) += (*tmp_B_matrix_T)(j, k) * (dstress[k] - stress0[k]) * fkt;
			}
		}

		(*tmp_AuxMatrix2) = 0.0;
		tmp_B_matrix_T->multi(*p_D, *tmp_AuxMatrix2);
		for (j = 0; j < nnodesHQ; j++)
		{
			tmp_B_matrix = this->vec_B_matrix[j];
			// Compute stiffness matrix
			(*tmp_AuxMatrix) = 0.0;
			tmp_AuxMatrix2->multi(*tmp_B_matrix, *tmp_AuxMatrix);

			// Local assembly of stiffness matrix
			for (k = 0; k < ele_dim; k++) {
				for (l = 0; l < ele_dim; l++) {
					(*tmp_Stiffness)(i + k * nnodesHQ, j + l * nnodesHQ) += (*tmp_AuxMatrix)(k, l) * fkt;
				}
			}
		}  // loop j
	}      // loop i

	// should restore pointer NW
	B_matrix = old_B_matrix;
	B_matrix_T = old_B_matrix_T;

	//---------------------------------------------------------
	// Assemble coupling matrix
	//---------------------------------------------------------
	// LoadFactor: factor of incremental loading, prescibed in rf_pcs.cpp

	if (PressureC)
	{
		fac = fkt;

		if (axisymmetry)
		{
			for (k = 0; k < nnodesHQ; k++)
			{
				for (l = 0; l < nnodes; l++)
				{
					for (j = 0; j < ele_dim; j++)
					{
						dN_dx = dshapefctHQ[nnodesHQ * j + k];
						if (j == 0) dN_dx += shapefctHQ[k] / Radius;

						f_buff = fac * dN_dx * shapefct[l];
						(*PressureC)(nnodesHQ* j + k, l) += f_buff;
					}
				}
			}
		}
		else
		{
			for (k = 0; k < nnodesHQ; k++)
			{
				for (l = 0; l < nnodes; l++)
				{
					for (j = 0; j < ele_dim; j++)
					{
						f_buff = fac * dshapefctHQ[nnodesHQ * j + k] * shapefct[l];
						(*PressureC)(nnodesHQ* j + k, l) += f_buff;
					}
				}
			}
		}
	}
	//---------------------------------------------------------
	// Assemble gravity force vector
	//---------------------------------------------------------
	if (rho > 0.0)
	{
		// 2D, in y-direction
		// 3D, in z-direction
		i = (ele_dim - 1) * nnodesHQ;
		const double coeff = rho * m_msp->grav_const * fkt;
		for (k = 0; k < nnodesHQ; k++)
			(*RHS)(i + k) += coeff * shapefctHQ[k];
		//        (*RHS)(i+k) += LoadFactor * rho * smat->grav_const *
		//        shapefctHQ[k] * fkt;
	}
}
/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: LocalAssembly
   Aufgabe:
           Compute the local finite element matrices
   Formalparameter:
           E:
         const int update  : indicator to update stress and strain only
 **************************************************************************/
void CFiniteElementVec::LocalAssembly(const int update)
{
	int j;
	double* a_n = NULL;

	Index = MeshElement->GetIndex();
	SetMemory();
	SetMaterial();
	// 12.2009. WW
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];
	ns = 4;
	if (MeshElement->GetDimension() == 3) ns = 6;

#if defined(NEW_EQS)
	b_rhs = pcs->eqs_new->getRHS();
#endif
	(*RHS) = 0.0;
	(*Stiffness) = 0.0;
	// 07.2011. WW
	if (PressureC) (*PressureC) = 0.0;
	for (int i = 0; i < nnodesHQ; i++)
		eqs_number[i] = MeshElement->GetNode(i)->GetEquationIndex(true);

	// For strain and stress extropolation all element types
	// Number of elements associated to nodes
	for (int i = 0; i < nnodes; i++)
		dbuff[i] =
			(double)MeshElement->GetNode(i)->getConnectedElementIDs().size();

	// Get displacement_n
	for (size_t i = 0; i < dim; i++)
		for (j = 0; j < nnodesHQ; j++)
			Disp[j + i * nnodesHQ] =
				pcs->GetNodeValue(nodes[j], Idx_dm0[i]);

	LocalAssembly_continuum(update);

	if (update == 0)
	{
		GlobalAssembly();
		// Output matrices
		if (pcs->Write_Matrix)
		{
			(*pcs->matrix_file) << "### Element: " << Index << "\n";
			(*pcs->matrix_file) << "---Stiffness matrix: "
			                    << "\n";
			Stiffness->Write(*pcs->matrix_file);
			(*pcs->matrix_file) << "---RHS: "
			                    << "\n";
			RHS->Write(*pcs->matrix_file);
			if (PressureC)
			{
				(*pcs->matrix_file) << "Pressue coupling matrix: "
				                    << "\n";
				PressureC->Write(*pcs->matrix_file);
			}
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly()
   Aufgabe:
           Assemble local matrics and vectors to the global system
 **************************************************************************/
bool CFiniteElementVec::GlobalAssembly()
{
	GlobalAssembly_RHS();

	GlobalAssembly_Stiffness();

#ifdef USE_PETSC
	// ScreenMessage2("-> add2GlobalMatrixII()\n");
	add2GlobalMatrixII();
#endif

	return true;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_Stiffness()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
 **************************************************************************/
void CFiniteElementVec::GlobalAssembly_Stiffness()
{
	int i, j;
	double f1, f2;
	f1 = 1.0;
	f2 = -1.0;

	double biot = 1.0;
	biot = m_msp->biot_const;
#if defined(NEW_EQS)
	CSparseMatrix* A = NULL;
	A = pcs->eqs_new->getA();
#endif

#ifndef USE_PETSC
	// Assemble stiffness matrix
	for (i = 0; i < nnodesHQ; i++)
	{
		const long eqs_number_i = eqs_number[i];
		for (j = 0; j < nnodesHQ; j++)
		{
			const long eqs_number_j = eqs_number[j];
			// Local assembly of stiffness matrix
			for (size_t k = 0; k < ele_dim; k++)
			{
				const long globalRowId = eqs_number_i + NodeShift[k];
				const long localRowId = i + k * nnodesHQ;
				for (size_t l = 0; l < ele_dim; l++)
				{
					double globalColId = eqs_number_j + NodeShift[l];
					double val =
					    f1 * (*Stiffness)(localRowId, j + l * nnodesHQ);
#if defined(USE_PETSC)  // || defined(other parallel libs)//10.3012. WW
// TODO
#else
#ifdef NEW_EQS
					double& a = (*A)(globalRowId, globalColId);
					#pragma omp atomic
					a += val;
//(*A)(globalRowId,globalColId) += val;
#endif
#endif
				}
			}
		}  // loop j
	}      // loop i
#endif

	// TEST OUT
	// Stiffness->Write();
	if (pcs->type / 40 != 1)  // Not monolithic scheme
		return;

	if (PressureC)
	{
		i = 0;               // phase
		if (Flow_Type == 2)  // Multi-phase-flow
			i = 1;
		GlobalAssembly_PressureCoupling(PressureC, f2 * biot, i);
	}

	// TEST OUT
	// PressureC->Write();
}
//--------------------------------------------------------------------------
/*!
   \brief Assembe the pressure coupling matrix

    to the global stiffness matrix in the monolithic scheme

   \param pCMatrix: the matrix
   \param fct: factor
   \param phase: phasse index

    07.2011. WW
 */
#if defined(USE_PETSC)  // || defined(other parallel libs)//10.3012. WW
void CFiniteElementVec::GlobalAssembly_PressureCoupling(Matrix*,
                                                        double,
                                                        const int)
{
}
#else
void CFiniteElementVec::GlobalAssembly_PressureCoupling(Matrix* pCMatrix,
                                                        double fct,
                                                        const int phase)
{
#if defined(NEW_EQS)
	CSparseMatrix* A = pcs->eqs_new->getA();
#endif

	int dim_shift = dim + phase;
	// Add pressure coupling matrix to the stifness matrix
	for (int i = 0; i < nnodesHQ; i++)
	{
		for (int j = 0; j < nnodes; j++)
		{
			for (size_t k = 0; k < ele_dim; k++)
			{
#ifdef NEW_EQS
				(*A)(NodeShift[k] + eqs_number[i],
				     NodeShift[dim_shift] + eqs_number[j]) +=
				    fct * (*pCMatrix)(nnodesHQ* k + i, j);
#endif
			}
		}
	}
}
#endif

#ifdef USE_PETSC
//------------------------------------------------------
// Add the local stiff matrix to the global one
void CFiniteElementVec::add2GlobalMatrixII()
{
	int const ndof = pcs->pcs_number_of_primary_nvals;
	int const m_dim = nnodesHQ * ndof;
	int const n_dim = m_dim;
	double* local_matrix = Stiffness->getEntryArray();
	double* local_vec = RHS->getEntryArray();

	for (int i = 0; i < nnodesHQ; i++)
	{
		MeshLib::CNode const* node = MeshElement->GetNode(i);
		const int offset = node->GetEquationIndex(true) * ndof;
		const bool isGhost = !pcs->m_msh->isNodeLocal(node->GetIndex());
		for (int k = 0; k < ndof; k++)
		{
			const int ki = k * nnodesHQ + i;
			row_ids[ki] = isGhost ? -1 : (offset + k);
			col_ids[ki] = offset + k;
		}
	}

	pcs->eqs_new->addMatrixEntries(m_dim, row_ids, n_dim, col_ids, &local_matrix[0]);

	static double temp_vec[100];
	for (int i = 0; i < m_dim; i++)
		temp_vec[i] = -local_vec[i];  // r -= RHS
	pcs->eqs_new->setArrayValues(1, m_dim, row_ids, &temp_vec[0]);
}
#endif

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::  GlobalAssembly_RHS()
   Aufgabe:
           Assemble local matrics and vectors to the global system

   Programming:
   02/2005   WW
   05/2005   WW dyn
 **************************************************************************/
void CFiniteElementVec::GlobalAssembly_RHS()
{
	bool Residual = false;

	if (H_Process)
	{
		if (pcs->type / 10 == 4)  // Monolithic scheme
		{
			if (FiniteElement::isNewtonKind(pcs->m_num->nls_method))
				Residual = true;
		}
		else  // Partitioned scheme
			Residual = true;
	}

	// Assemble coupling matrix
	if (Residual)
	{
		// should calculate (p - p0) because OGS calculates (stress - stress0)
		double const* p0 = pcs->GetInitialFluidPressure();
		switch (Flow_Type)
		{
			case 0:  // Liquid flow
				// For monolithic scheme and liquid flow, the limit of positive
				// pressure must be removed
				for (int i = 0; i < nnodes; i++)
				{
					auto val_n = h_pcs->GetNodeValue(nodes[i], idx_P1);

					// Initial pressure should be subtracted, i.e. (p-p0)
					// because DEFORMATION
					// calculates stress balance of changes from the initial
					// stress.
					if (pcs->calcDiffFromStress0 && p0 != NULL)
					{
						val_n -= p0[nodes[i]];
					}
					AuxNodal[i] = val_n;
				}
				break;
			case 10:  // Ground_flow. Will be merged to case 0
				for (int i = 0; i < nnodes; i++)
					AuxNodal[i] = h_pcs->GetNodeValue(nodes[i], idx_P1);
				break;
		}

		const int dim_times_nnodesHQ(dim * nnodesHQ);
		// Coupling effect to RHS
		if (Flow_Type != 2)
		{
			for (int i = 0; i < dim_times_nnodesHQ; i++)
				AuxNodal1[i] = 0.0;
			PressureC->multi(AuxNodal, AuxNodal1);
		}
		for (int i = 0; i < dim_times_nnodesHQ; i++)
			(*RHS)(i) -= fabs(m_msp->biot_const) * AuxNodal1[i];
	}  // End if partioned

	// RHS->Write();

#ifndef USE_PETSC
	//	std::cerr << "e: " << Index << "\n";
	for (size_t i = 0; i < dim; i++)
	{
		for (int j = 0; j < nnodesHQ; j++)
		{
			//		    std::cerr << eqs_number[j] << "," << NodeShift[i] << ","
			//<< i * nnodesHQ + j << "\n";
			long global_id = eqs_number[j] + NodeShift[i];
			double val = (*RHS)(i* nnodesHQ + j);
#pragma omp atomic
			b_rhs[global_id] -= val;
		}
	}
#endif
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec:: LocalAssembly_continuum()
   Aufgabe:
           Compute the local finite element matrices with the framework
        of continuum assumption
   Formalparameter:
           E:
         const int update  : indicator to update stress and strain only

   Programming:
   02/2000   OK   Elastische Elemente
   09/2002   WW   Local assemby of stiffness matrix of elasto-plastic
   tangential model
   Local assemby of residual
   07/2003   WW   Quadratic triangle element
   06/2004   WW   Generalize for different element types as a member of class
   12/2005   WW   Creep
 **************************************************************************/
void CFiniteElementVec::LocalAssembly_continuum(const int update)
{
	long i;

	Matrix* p_D = NULL;
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];

	// ---- Gauss integral
	int gp_r = 0, gp_s = 0, gp_t;
	gp = 0;
	gp_t = 0;
	double fkt = 0.0;

	// WW double *DevStress ;
	const int PModel = m_msp->Plasticity_type;
	double dPhi = 0.0;  // Sclar factor for the plastic strain
	//  double J2=0.0;
	double dS = 0.0;

	double ThermalExpansion = 0.0;
	double t1 = 0.0;
	bool Strain_TCS = false;

	//
	ThermalExpansion = 0.0;
	// Thermal effect
	if (m_msp->Thermal_Expansion() > 0.0)
		ThermalExpansion = m_msp->Thermal_Expansion();

	// Get porosity model
	// ---- Material properties
	// For swelling pressure;
	double deporo = 0.0;
	// OK_MMP
	//--------------------------------------------------------------------
	// MMP medium properties
	int PoroModel = m_mmp->porosity_model;
	if (PoroModel == 4)
		// OK deporo =  PCSGetElementPorosityChangeRate(index)/(double)ele_dim;
		// MX  deporo =
		// PCSGetELEValue(index,NULL,1.0,"n_sw_Rate")/(double)ele_dim;
		deporo = h_pcs->GetElementValue(
		             Index, h_pcs->GetElementValueIndex("n_sw_rate")) /
		         (double)ele_dim;
	if (T_Flag)
		for (i = 0; i < nnodes; i++)
		{
			T1[i] = t_pcs->GetNodeValue(nodes[i], idx_T1);
			dT[i] = t_pcs->GetNodeValue(nodes[i], idx_T1) -
			          t_pcs->GetNodeValue(nodes[i], idx_T0);
		}
	//

	if (PModel == 1 || PModel == 10 ||
	    PModel == 11)  // WX:modified for DP with tension cutoff
		m_msp->CalulateCoefficent_DP();
	//
	if (PModel != 3 && m_msp->Youngs_mode != 2)  // modified due to transverse
	                                            // isotropic elasticity: UJG
	                                            // 24.11.2009
	{
		if (m_msp->Youngs_mode < 10 || m_msp->Youngs_mode > 13)
		{
			m_msp->Calculate_Lame_Constant();
			m_msp->ElasticConsitutive(ele_dim, De);
		}
		else
			*De = *(m_msp->getD_tran());  // UJG/WW
	}

	if (PModel == 5) m_msp->CalculateCoefficent_HOEKBROWN();

	//
	if (PoroModel == 4 || T_Flag || m_msp->Creep_mode > 0) Strain_TCS = true;
	//
	if (m_msp->CreepModel() == 1000)  // HL_ODS
		m_msp->CleanTrBuffer_HL_ODS();
	for (i = 0; i < ns; i++)
		stress0[i] = .0;

    auto &nodal_p1 = AuxNodal;
    auto &nodal_p0 = AuxNodal1;
    if (H_Process)
    {
        for (int i = 0; i < nnodes; i++)
        {
            nodal_p1[i] = h_pcs->GetNodeValue(nodes[i], idx_P1);
            nodal_p0[i] = h_pcs->GetNodeValue(nodes[i], idx_P1 - 1);
        }
    }

	// Loop over Gauss points
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		//---------------------------------------------------------
		//  Get local coordinates and weights
		//  Compute Jacobian matrix and its determinate
		//---------------------------------------------------------
		fkt = GetGaussData(gp, gp_r, gp_s, gp_t);

		//---------------------------------------------------------
		// Compute geometry
		//---------------------------------------------------------
		ComputeGradShapefct(2);
		ComputeShapefct(2);
		if (m_msp->Youngs_mode == 2)  // WW/UJG. 22.01.2009
		{
			m_msp->CalcYoungs_SVV(CalcStrain_v());
			m_msp->ElasticConsitutive(ele_dim, De);
		}
		if (m_msp->Youngs_mode == 3)  // AJ. 16.06.2014
		{
			m_msp->CalcYoungs_Drained(CalcStress_eff());
			m_msp->ElasticConsitutive(ele_dim, De);
		}

		ComputeStrain();
		if (update) RecordGuassStrain(gp, gp_r, gp_s, gp_t);
		if (F_Flag || T_Flag)
			ComputeShapefct(1);  // Linear order interpolation function
		//---------------------------------------------------------
		// Material properties (Integration of the stress)
		//---------------------------------------------------------
		// Initial the stress vector
		if (PModel != 3)
		{
			for (i = 0; i < ns; i++)
				dstress[i] = 0.0;
			De->multi(dstrain, dstress);
		}

		//---------------------------------------------------------
		// Integrate the stress by return mapping:
		//---------------------------------------------------------
		switch (PModel)
		{
			case -1:  // Pure elasticity
				// Non-linear elasticity: TE model. 10.03.2008. WW
				for (i = 0; i < ns; i++)
					dstress[i] += (*eleV_DM->Stress)(i, gp);
				break;
			case 1:  // Drucker-Prager model
				if (m_msp->StressIntegrationDP(gp, eleV_DM, &dstress[0], dPhi,
				                              update))

					// WW DevStress = smat->devS;
					m_msp->ConsistentTangentialDP(ConsistDep, dPhi, ele_dim);
				break;
			case 10:  // Drucker-Prager model, direct integration. 02/06 WW
				if (m_msp->DirectStressIntegrationDP(gp, eleV_DM, &dstress[0],
				                                    update))
				{
					*ConsistDep = *De;
					m_msp->TangentialDP(ConsistDep);
					dPhi = 1.0;
				}
				break;
			case 11:  // WX: 08.2010
			{
				double mm = 0.;  // WX:09.2010. for DP with Tension.
				switch (m_msp->DirectStressIntegrationDPwithTension(
					gp, De, eleV_DM, &dstress[0], update, mm))
				{
					case 1:
					{
						*ConsistDep = *De;
						m_msp->TangentialDP2(ConsistDep);
						dPhi = 1.0;
					}
					break;
					case 2:
					{
						*ConsistDep = *De;
						m_msp->TangentialDPwithTension(ConsistDep, mm);
						dPhi = 1.0;
					}
					break;
					case 3:
					{
						*ConsistDep = *De;
						m_msp->TangentialDPwithTensionCorner(ConsistDep, mm);
						dPhi = 1.0;
					}
					break;
					default:
						break;
				}
				break;
			}
			case 2:  // Rotational hardening model
				// Compute stesses and plastic multi-plier
				dPhi = 0.0;
				if (m_msp->CalStress_and_TangentialMatrix_SYS(
						gp, eleV_DM, De, ConsistDep, &dstress[0], update) > 0)
					dPhi = 1.0;
				break;
			case 3:  // Generalized Cam-Clay model
				for (i = 0; i < ns; i++)
					dstress[i] = dstrain[i];
				//
				if (m_msp->SwellingPressureType == 3)
				{
					double suc = interpolate(AuxNodal1);
					double dsuc = interpolate(AuxNodal);
					(*m_msp->data_Youngs)(7) = suc;
					(*m_msp->data_Youngs)(8) = dsuc;
					m_msp->CalStress_and_TangentialMatrix_CC_SubStep(
						gp, eleV_DM, &dstress[0], ConsistDep, update);
					//
					// double pn = -((*eleV_DM->Stress)(0,
					// gp)+(*eleV_DM->Stress)(1, gp)+
					//             (*eleV_DM->Stress)(2, gp))/3.0;
					// de_vsw =
					// -smat->TEPSwellingParameter(pn,suc)*dsuc/(suc+1.0e5);
				}
				/*
				   else if (smat->SwellingPressureType==4)  //Power. 07.05.2008
				   WW
				   {
				   double suc = interpolate(AuxNodal1);
				   double dsuc = interpolate(AuxNodal);
				   smat->TEPSwellingParameter_kis(suc);
				   S_Water = m_mmp->SaturationCapillaryPressureFunction(suc);
				   dS = S_Water -
				   m_mmp->SaturationCapillaryPressureFunction(suc-dsuc);
				   de_vsw = pow(S_Water, (*smat->data_Youngs)(2) )*dS;
				   }
				 */
				else
					m_msp->CalStress_and_TangentialMatrix_CC(
						gp, eleV_DM, &dstress[0], ConsistDep, update);

				dPhi = 1.0;
				break;
			case 4:  // Mohr-Coloumb	//WX:10.2010
				if (m_msp->DirectStressIntegrationMOHR(gp, eleV_DM, &dstress[0],
				                                      update, De))
				{
					*ConsistDep = *De;
					// also for tension
					m_msp->TangentialMohrShear(ConsistDep);
					// ConsistDep->Write();
					dPhi = 1.0;
				}
				break;
				/*case 5:
				   if(smat->StressIntegrationHoekBrown(gp, eleV_DM, dstress,
				   update, De))
				   {
				   *ConsistDep = *De;
				         smat->TangentialHoekBrown(ConsistDep);		//also for
				   tension
				         //ConsistDep->Write();
				         dPhi = 1.0;
				   }
				   break;*/
		}
		// --------------------------------------------------------------------
		// Stress increment by heat, swelling, or heat
		//
		if (Strain_TCS)
		{
			if (PModel == 3) m_msp->ElasticConsitutive(ele_dim, De);
			for (i = 0; i < ns; i++)
				strain_ne[i] = 0.0;
			if (PoroModel == 4)  // For swelling pressure
                for (i = 0; i < 3; i++)
					strain_ne[i] -= deporo;
			//
			if (T_Flag)  // Contribution by thermal expansion
			{
				Tem = 0.0;
				t1 = 0.0;
				for (i = 0; i < nnodes; i++)
				{
					Tem += shapefct[i] * dT[i];
					t1 += shapefct[i] * T1[i];
				}
                for (i = 0; i < 3; i++)
					strain_ne[i] -= ThermalExpansion * Tem;
			}
			// Strain increment by creep
			if (m_msp->Creep_mode == 1 || m_msp->Creep_mode == 2)
			{
				for (i = 0; i < ns; i++)
					stress_ne[i] = (*eleV_DM->Stress)(i, gp);
				m_msp->AddStain_by_Creep(ns, &stress_ne[0], &strain_ne[0], t1);
			}
			if (m_msp->Creep_mode == 1000)  // HL_ODS. Strain increment by creep
			{
				for (i = 0; i < ns; i++)
					stress_ne[i] = (*eleV_DM->Stress)(i, gp);
				m_msp->AddStain_by_HL_ODS(eleV_DM, &stress_ne[0], &strain_ne[0], t1);
			}
			// Stress deduced by thermal or swelling strain incremental:
			De->multi(strain_ne, dstress);
		}
		// Assemble matrices and RHS
		if (update < 1)
		{
			//---------------------------------------------------------
			// Assemble matrices and RHS
			//---------------------------------------------------------
			if (dPhi <= 0.0)
				p_D = De;
			else
				p_D = ConsistDep;
			if (pcs->calcDiffFromStress0)
			{
				for (i = 0; i < ns; i++)
					stress0[i] = (*eleV_DM->Stress0)(i, gp);
			}
			ComputeMatrix_RHS(fkt, p_D);
		}
		else  // Update stress
        {
			for (i = 0; i < ns; i++)
			{
				(*eleV_DM->Stress)(i, gp) = dstress[i];
                (*eleV_DM->dTotalStress)(i, gp) = (*eleV_DM->Stress)(i, gp) - (*eleV_DM->Stress_last_ts)(i, gp);
                (*eleV_DM->Strain)(i, gp) += dstrain[i];
			}
            if (H_Process)
            {
                const double gp_p1 = interpolate(nodal_p1);
                const double gp_p0 = interpolate(nodal_p0);
                const double d_bp = m_msp->biot_const * gp_p1 - m_msp->biot_const * gp_p0;
                for (long i = 0; i < 3; i++)
                    (*eleV_DM->dTotalStress)(i, gp) -= d_bp;
            }
        }
	}
	// The mapping of Gauss point strain to element nodes
	if (update)
        ExtropolateGaussStrain();
	else if (m_msp->Creep_mode == 1000)  // HL_ODS. Strain increment by creep
		m_msp->AccumulateEtr_HL_ODS(eleV_DM, nGaussPoints);

	/*
	   //TEST

	      //TEST
	        if(update&&Index==0)
	        {
	           oss<<" Stress "<<endl;
	           eleV_DM->Stress->Write(oss);
	           oss.close();
	        }
	 */
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::RecordGuassValues()
   Aufgabe:
           Accumulate stress at each nodes
   Formalparameter:
           E:

   Programming:
   06/2004   WW
 **************************************************************************/
bool CFiniteElementVec::RecordGuassStrain(const int gp, const int gp_r,
                                          const int gp_s, int gp_t)
{
	int LoIndex = 0;

	//---------------------------------------------------------
	// Accumulate strains
	//---------------------------------------------------------
	switch (MeshElement->GetElementType())
	{
		case MshElemType::QUAD:  // Quadralateral
			LoIndex = GetLocalIndex(gp_r, gp_s, gp_t);
			Sxx[LoIndex] = dstrain[0];
			Syy[LoIndex] = dstrain[1];
			Sxy[LoIndex] = dstrain[3];
			Szz[LoIndex] = dstrain[2];
			break;
		case MshElemType::TRIANGLE:  // Triangle
			Sxx[gp] = dstrain[0];
			Syy[gp] = dstrain[1];
			Szz[gp] = dstrain[2];
			Sxy[gp] = dstrain[3];
			break;
		case MshElemType::HEXAHEDRON:  // Hexahedra
			LoIndex = GetLocalIndex(gp_r, gp_s, gp_t);
			if (LoIndex < 0) return false;
			Sxx[LoIndex] = dstrain[0];
			Syy[LoIndex] = dstrain[1];
			Szz[LoIndex] = dstrain[2];
			Sxy[LoIndex] = dstrain[3];
			Sxz[LoIndex] = dstrain[4];
			Syz[LoIndex] = dstrain[5];
			break;
		case MshElemType::TETRAHEDRON:  // Tedrahedra
			Sxx[gp] = dstrain[0];
			Syy[gp] = dstrain[1];
			Szz[gp] = dstrain[2];
			Sxy[gp] = dstrain[3];
			Sxz[gp] = dstrain[4];
			Syz[gp] = dstrain[5];
			break;
		case MshElemType::PYRAMID:
			Sxx[gp] = dstrain[0];
			Syy[gp] = dstrain[1];
			Szz[gp] = dstrain[2];
			Sxy[gp] = dstrain[3];
			Sxz[gp] = dstrain[4];
			Syz[gp] = dstrain[5];
			break;
		default:
			break;
			// 3D
	}
	return false;
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ExtropolateGuassStrain()
   Aufgabe:
           Extropolate the Gauss point strains to nodes
   Formalparameter:
           E:

   Programming:
   06/2004   WW
   02/2007   Make it work for all 2nd variables
 **************************************************************************/
void CFiniteElementVec::ExtropolateGaussStrain()
{
	int i, j;
	//  int l1,l2,l3,l4; //, counter;
	double ESxx, ESyy, ESzz, ESxy, ESxz, ESyz;
	double avgESxx, avgESyy, avgESzz, avgESxy, avgESxz, avgESyz;
	int i_s, i_e, ish, k = 0;
	gp = 0;
	// double Area1, Area2, Tol=10e-9;
	// For strain and stress extropolation all element types
	for (int i = 0; i < nnodes; i++)
		dbuff[i] = (double)MeshElement->GetNode(i)->getConnectedElementIDs().size();


	MshElemType::type ElementType = MeshElement->GetElementType();
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		if (ElementType == MshElemType::QUAD ||
		    ElementType == MshElemType::HEXAHEDRON)
		{
			int gp_r, gp_s, gp_t;
			SetGaussPoint(gp, gp_r, gp_s, gp_t);
			i = GetLocalIndex(gp_r, gp_s, gp_t);
			if (i == -1) continue;
		}
		else
			i = gp;

		Sxx[i] = (*eleV_DM->Strain)(0, gp);
		Syy[i] = (*eleV_DM->Strain)(1, gp);
		Szz[i] = (*eleV_DM->Strain)(2, gp);
		Sxy[i] = (*eleV_DM->Strain)(3, gp);
		if (ele_dim == 3)
		{
			Sxz[i] = (*eleV_DM->Strain)(4, gp);
			Syz[i] = (*eleV_DM->Strain)(5, gp);
		}
	}


	// l1=l2=l3=l4=0;
	if (ElementType == MshElemType::QUAD ||
	    ElementType == MshElemType::HEXAHEDRON)
		Xi_p = CalcXi_p();

	//
	i_s = 0;
	i_e = nnodes;
	ish = 0;
	if (ElementType == MshElemType::TETRAHEDRON)  // tet
	{
		i_s = 1;
		i_e = nnodes + 1;
		ish = 1;
	}
	//---------------------------------------------------------
	// Mapping Gauss point strains to nodes and update nodes
	// strains:
	//---------------------------------------------------------
	avgESxx = avgESyy = avgESzz = avgESxy = avgESxz = avgESyz = 0.0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
	{
		// average
		avgESxx = CalcAverageGaussPointValues(Sxx);
		avgESyy = CalcAverageGaussPointValues(Syy);
		avgESzz = CalcAverageGaussPointValues(Szz);
		avgESxy = CalcAverageGaussPointValues(Sxy);
		avgESxz = CalcAverageGaussPointValues(Sxz);
		avgESyz = CalcAverageGaussPointValues(Syz);
	}

	for (i = 0; i < nnodes; i++)
	{
		ESxx = ESyy = ESzz = ESxy = ESxz = ESyz = 0.0;

		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			SetExtropoGaussPoints(i);
			ComputeShapefct(1);  // Linear interpolation function
			//
			for (j = i_s; j < i_e; j++)
			{
				k = j - ish;
				ESxx += Sxx[j] * shapefct[k];
				ESyy += Syy[j] * shapefct[k];
				ESxy += Sxy[j] * shapefct[k];
				ESzz += Szz[j] * shapefct[k];
				if (ele_dim == 3)
				{
					ESxz += Sxz[j] * shapefct[k];
					ESyz += Syz[j] * shapefct[k];
				}
			}
		}
		else if (this->GetExtrapoMethod() ==
		         ExtrapolationMethod::EXTRAPO_AVERAGE)
		{
			// average
			ESxx = avgESxx;
			ESyy = avgESyy;
			ESxy = avgESxy;
			ESzz = avgESzz;
			if (ele_dim == 3)
			{
				ESxz = avgESxz;
				ESyz = avgESyz;
			}
		}

		// Average value of the contribution of all neighbor elements
		ESxx /= dbuff[i];
		ESyy /= dbuff[i];
		ESxy /= dbuff[i];
		ESzz /= dbuff[i];

		ESxx += pcs->GetNodeValue(nodes[i], Idx_Strain[0]);
		ESyy += pcs->GetNodeValue(nodes[i], Idx_Strain[1]);
		ESzz += pcs->GetNodeValue(nodes[i], Idx_Strain[2]);
		ESxy += pcs->GetNodeValue(nodes[i], Idx_Strain[3]);

		pcs->SetNodeValue(nodes[i], Idx_Strain[0], ESxx);
		pcs->SetNodeValue(nodes[i], Idx_Strain[1], ESyy);
		pcs->SetNodeValue(nodes[i], Idx_Strain[2], ESzz);
		pcs->SetNodeValue(nodes[i], Idx_Strain[3], ESxy);

		if (ele_dim == 3)
		{
			ESxz /= dbuff[i];
			ESyz /= dbuff[i];
			ESxz += pcs->GetNodeValue(nodes[i], Idx_Strain[4]);
			ESyz += pcs->GetNodeValue(nodes[i], Idx_Strain[5]);
			//
			pcs->SetNodeValue(nodes[i], Idx_Strain[4], ESxz);
			pcs->SetNodeValue(nodes[i], Idx_Strain[5], ESyz);
		}
	}
}

/***************************************************************************
   GeoSys - Funktion:
           CFiniteElementVec::ExtropolateGuassStress()
   Aufgabe:
           Extropolate the Gauss point strains to nodes
   Formalparameter:
           E:

   Programming:
   06/2004   WW
   03/2007   WW  Generize for all 2nd variables
 **************************************************************************/
void CFiniteElementVec::ExtropolateGaussStress()
{
	int i, j, gp_r, gp_s, gp_t;
	// int l1,l2,l3,l4; //, counter;
	double ESxx, ESyy, ESzz, ESxy, ESxz, ESyz, Pls;
	double avgESxx, avgESyy, avgESzz, avgESxy, avgESxz, avgESyz, avgPls;
	int i_s, i_e, ish, k = 0;
	MshElemType::type ElementType = MeshElement->GetElementType();
	long node_i = 0;
	// For strain and stress extropolation all element types
	// Number of elements associated to nodes
	for (i = 0; i < nnodes; i++)
		dbuff[i] =
			(double)MeshElement->GetNode(i)->getConnectedElementIDs().size();
	//
	gp = gp_r = gp_s = gp_t = 0;
	eleV_DM = ele_value_dm[MeshElement->GetIndex()];
	if (eleV_DM->pStrain)  // 08.02.2008 WW
		idx_pls = pcs->GetNodeValueIndex("STRAIN_PLS");
	//
	for (gp = 0; gp < nGaussPoints; gp++)
	{
		if (ElementType == MshElemType::QUAD ||
		    ElementType == MshElemType::HEXAHEDRON)
		{
			SetGaussPoint(gp, gp_r, gp_s, gp_t);
			i = GetLocalIndex(gp_r, gp_s, gp_t);
			if (i == -1) continue;
		}
		else
			i = gp;

		Sxx[i] = (*eleV_DM->Stress)(0, gp);
		Syy[i] = (*eleV_DM->Stress)(1, gp);
		Szz[i] = (*eleV_DM->Stress)(2, gp);
		Sxy[i] = (*eleV_DM->Stress)(3, gp);
		if (eleV_DM->pStrain)
			pstr[i] = (*eleV_DM->pStrain)(gp);
		else
			pstr[i] = 0.0;  // 08.02.2008 WW
		if (ele_dim == 3)
		{
			Sxz[i] = (*eleV_DM->Stress)(4, gp);
			Syz[i] = (*eleV_DM->Stress)(5, gp);
		}
	}
	//
	if (ElementType == MshElemType::QUAD ||
	    ElementType == MshElemType::HEXAHEDRON)
		Xi_p = CalcXi_p();

	//
	i_s = 0;
	i_e = nnodes;
	ish = 0;
	if (ElementType == MshElemType::TETRAHEDRON)  // tet
	{
		i_s = 1;
		i_e = nnodes + 1;
		ish = 1;
	}
	//---------------------------------------------------------
	// Mapping Gauss point strains to nodes and update nodes
	// strains:
	//---------------------------------------------------------
	avgESxx = avgESyy = avgESzz = avgESxy = avgESxz = avgESyz = avgPls = 0.0;
	if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_AVERAGE)
	{
		// average
		avgESxx = CalcAverageGaussPointValues(Sxx);
		avgESyy = CalcAverageGaussPointValues(Syy);
		avgESzz = CalcAverageGaussPointValues(Szz);
		avgESxy = CalcAverageGaussPointValues(Sxy);
		avgESxz = CalcAverageGaussPointValues(Sxz);
		avgESyz = CalcAverageGaussPointValues(Syz);
		avgPls = CalcAverageGaussPointValues(pstr);
	}

	for (i = 0; i < nnodes; i++)
	{
		ESxx = ESyy = ESzz = ESxy = ESxz = ESyz = Pls = 0.0;

		// Calculate values at nodes
		if (this->GetExtrapoMethod() == ExtrapolationMethod::EXTRAPO_LINEAR)
		{
			//
			SetExtropoGaussPoints(i);
			//
			ComputeShapefct(1);  // Linear interpolation function
			//
			for (j = i_s; j < i_e; j++)
			{
				k = j - ish;
				ESxx += Sxx[j] * shapefct[k];
				ESyy += Syy[j] * shapefct[k];
				ESxy += Sxy[j] * shapefct[k];
				ESzz += Szz[j] * shapefct[k];
				Pls += pstr[j] * shapefct[k];
				if (ele_dim == 3)
				{
					ESxz += Sxz[j] * shapefct[k];
					ESyz += Syz[j] * shapefct[k];
				}
			}
		}
		else if (this->GetExtrapoMethod() ==
		         ExtrapolationMethod::EXTRAPO_AVERAGE)
		{
			// average
			ESxx = avgESxx;
			ESyy = avgESyy;
			ESxy = avgESxy;
			ESzz = avgESzz;
			Pls = avgPls;
			if (ele_dim == 3)
			{
				ESxz = avgESxz;
				ESyz = avgESyz;
			}
		}

		// Average value of the contribution of ell neighbor elements
		ESxx /= dbuff[i];
		ESyy /= dbuff[i];
		ESxy /= dbuff[i];
		ESzz /= dbuff[i];
		Pls /= dbuff[i];
		//
		node_i = nodes[i];
		ESxx += pcs->GetNodeValue(node_i, Idx_Stress[0]);
		ESyy += pcs->GetNodeValue(node_i, Idx_Stress[1]);
		ESzz += pcs->GetNodeValue(node_i, Idx_Stress[2]);
		ESxy += pcs->GetNodeValue(node_i, Idx_Stress[3]);
		if (eleV_DM->pStrain)  // 08.02.2008 WW
			Pls += pcs->GetNodeValue(node_i, idx_pls);

		pcs->SetNodeValue(node_i, Idx_Stress[0], ESxx);
		pcs->SetNodeValue(node_i, Idx_Stress[1], ESyy);
		pcs->SetNodeValue(node_i, Idx_Stress[2], ESzz);
		pcs->SetNodeValue(node_i, Idx_Stress[3], ESxy);
		if (eleV_DM->pStrain)  // 08.02.2008 WW
			pcs->SetNodeValue(node_i, idx_pls, fabs(Pls));

		if (ele_dim == 3)
		{
			ESxz /= dbuff[i];
			ESyz /= dbuff[i];

			ESxz += pcs->GetNodeValue(node_i, Idx_Stress[4]);
			ESyz += pcs->GetNodeValue(node_i, Idx_Stress[5]);

			pcs->SetNodeValue(node_i, Idx_Stress[4], ESxz);
			pcs->SetNodeValue(node_i, Idx_Stress[5], ESyz);
		}
	}
}

} // end namespace FiniteElement
