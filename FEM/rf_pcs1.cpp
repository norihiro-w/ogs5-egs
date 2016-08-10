/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#if defined(USE_PETSC)
#include "rf_pcs.h"

#include <cmath>
#include <cfloat>

#include "display.h"

#include "PETSC/PETScLinearSolver.h"
#include "FEMEnums.h"
#include "msh_lib.h"

using namespace petsc_group;
using namespace std;
using namespace MeshLib;
using namespace FiniteElement;
/*!
   Initialize the RHS array of the system of equations with the previous
   solution.

   03.2012. WW
*/
void CRFProcess::InitializeRHS_with_u0(const bool quad)
{
	int j, ish;
	vector<int> ix;
	vector<double> val;
	int nidx1;

	int g_nnodes = m_msh->getNumNodesLocal();  // GetNodesNumber(false);
	if (quad) g_nnodes = m_msh->getNumNodesLocal_Q();
	const int size = g_nnodes * pcs_number_of_primary_nvals;
	ix.resize(size);
	val.resize(size);

	for (int i = 0; i < pcs_number_of_primary_nvals; i++)
	{
		// new time
		nidx1 = GetNodeValueIndex(pcs_primary_function_name[i]) + 1;
		for (j = 0; j < g_nnodes; j++)
		{
			ish = pcs_number_of_primary_nvals * j + i;
			ix[ish] =
			    pcs_number_of_primary_nvals * m_msh->Eqs2Global_NodeIndex[j] +
			    i;
			val[ish] = GetNodeValue(j, nidx1);
		}
	}
	// Assign u0 to x
	eqs_new->setArrayValues(0, size, &ix[0], &val[0], INSERT_VALUES);
	eqs_new->AssembleUnkowns_PETSc();
}

#if 0
/*!
    Parallel definition of this function for DDC by node

    WW 03.2012
*/
double CRFProcess::CalcIterationNODError(int method)
{
  long i, g_nnodes;
  double error_l, change_l, max_cl, min_cl;
  double error, change, max_c, min_c;
  //-----------------------------------------------------
  int nidx1;
  int ii;
  g_nnodes = m_msh->getNumNodesLocal(); 
  //-----------------------------------------------------
  //double* eqs_x = NULL;     // local solution

  error_l = 0.;
  change_l = 0.;
  
  max_cl = 0.;
  min_cl = 1.e99;
  
  double val = 0.;
  switch (method)
    {
    default:
    case 0:
      return 0.;
      // Maximum error
    case 1:
      //
      for(ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
	  //new time
	  nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
	  for (i = 0l; i < g_nnodes; i++) 
	    {
	      val = eqs_x[pcs_number_of_primary_nvals*m_msh->Eqs2Global_NodeIndex[i] + ii];
	      error_l = max(error_l, fabs(GetNodeValue(i, nidx1) - val));
	    }
	}
      MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      break;
    case 2:
      for(ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
	  //new time
	  nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
	  for (i = 0l; i < g_nnodes; i++) 
	    {
	      val = eqs_x[pcs_number_of_primary_nvals*m_msh->Eqs2Global_NodeIndex[i] + ii];
	      error_l =  max(error_l,  fabs(val - GetNodeValue(i, nidx1))
			   / (fabs(val) + fabs(GetNodeValue(i, nidx1)) + MKleinsteZahl) );
	    }
	}
      MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      break;
    case 3:
      for(ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
	  //new time
	  nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
	  for (i = 0l; i < g_nnodes; i++) 
	    {
	      val =  eqs_x[pcs_number_of_primary_nvals*m_msh->Eqs2Global_NodeIndex[i] + ii];
	      error_l = max(error_l, fabs(val - GetNodeValue(i, nidx1)));
	      max_cl = max(max(max_cl, fabs(fabs( val ))), fabs(GetNodeValue(i, nidx1)));
	    }
	}
      MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&max_cl, &max_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      error /= (max_c + MKleinsteZahl);
      break;
    case 4:
      for(ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
	  //new time
	  nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
	  for (i = 0l; i < g_nnodes; i++) 
	    {
	      val = eqs_x[pcs_number_of_primary_nvals*m_msh->Eqs2Global_NodeIndex[i] + ii]; 
	      error_l =
		max(error_l, fabs(val - GetNodeValue(i, nidx1)));
	      min_cl = min(min_cl, fabs(val));
	      max_cl = max(max_cl, fabs(val));
	    }
	}
      MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&max_cl, &max_c, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&min_cl, &min_c, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
      error /= (max_c - min_c + MKleinsteZahl);
      break;
    case 5:
      //
      for(ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
	  //new time
	  nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
	  for (i = 0l; i < g_nnodes; i++) 
	    {
	      val =  eqs_x[pcs_number_of_primary_nvals*m_msh->Eqs2Global_NodeIndex[i] + ii];
	      error_l = max(error_l, fabs(val -  GetNodeValue(i, nidx1))
		    / (fabs(val - GetNodeValue(i, nidx1 - 1)) + MKleinsteZahl));
	    }
	}
      MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      break;
    case 6:
      for(ii = 0; ii < pcs_number_of_primary_nvals; ii++)
	{
	  //new time
	  nidx1 = GetNodeValueIndex(pcs_primary_function_name[ii]) + 1;
	  for (i = 0l; i < g_nnodes; i++) 
	    {
	      val =  eqs_x[pcs_number_of_primary_nvals*m_msh->Eqs2Global_NodeIndex[i] + ii];
	      error_l = max(error_l, fabs(val -  GetNodeValue(i, nidx1)));
	      change_l = max(change_l, fabs(val - GetNodeValue(i, nidx1 - 1)));
	    }
	}
      MPI_Allreduce(&error_l, &error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      MPI_Allreduce(&change_l, &change, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
      error /= (change + MKleinsteZahl);
      break;
    }
  
  return error; 
}
#endif

/// Set PETSc solver. 10.2012. WW
void CRFProcess::setSolver(petsc_group::PETScLinearSolver* petsc_solver)
{
	eqs_new = petsc_solver;
	eqs_new->Config(m_num->ls_error_tolerance, m_num->ls_max_iterations,
	                m_num->getLinearSolverName(),
	                m_num->getPreconditionerName(), m_num->ls_extra_arg,
	                convertProcessTypeToString(this->getProcessType()) + "_");
}
//------------------------------------------------------------
/*!
     PETSc version of CreateEQS_LinearSolver()
     03.2012 WW
     10.2012 WW
*/
void CreateEQS_LinearSolver()
{
	int rank_p;
	int size_p;
	MPI_Comm_rank(PETSC_COMM_WORLD, &rank_p);
	MPI_Comm_size(PETSC_COMM_WORLD, &size_p);

	MeshLib::CFEMesh* mesh = fem_msh_vector[0];

	const int n_linear_nodes = mesh->getNumNodesGlobal();
	const int n_quadratic_nodes = mesh->getNumNodesGlobal_Q();
	// const int nnl = mesh->getNumNodesLocal();
	// const int nnl_q = mesh->getNumNodesLocal_Q();
	const int dim = mesh->GetMaxElementDim();

	int sparse_info[4] = {}; // d_nz, o_nz, nz, m_size_loc
	for (size_t i = 0; i < pcs_vector.size(); i++)
	{
		CRFProcess* a_pcs = pcs_vector[i];
		FiniteElement::ProcessType pcs_type = a_pcs->getProcessType();

		if (pcs_type == FLUID_MOMENTUM) continue;

		ScreenMessage(
		    "-> Process: %s\n",
		    FiniteElement::convertProcessTypeToString(pcs_type).c_str());

		const bool quadratic = isDeformationProcess(pcs_type);
		a_pcs->CheckMarkedElement();
		mesh->ConnectedNodes(quadratic);
		const int max_connected_nodes = mesh->calMaximumConnectedNodes();
		ScreenMessage("-> max. number of connected nodes = %d\n", max_connected_nodes);

//		const int max_connected_eles = mesh->getMaxNumConnectedElements(); // (dim < 3) ? 4 : 8;
//		const int max_ele_nodes = mesh->getMaxNumNodesOfElement(quadratic); //(dim < 3) ? 9 : 20;
//		ScreenMessage("-> max. number of connected elements = %d\n", max_connected_eles);
//		ScreenMessage("-> max. number of element nodes = %d\n", max_ele_nodes);

		PETScLinearSolver* eqs = NULL;
		if (isDeformationProcess(pcs_type))
		{
			int eqs_dim = n_quadratic_nodes * dim;
			//vector<int> global_n_id;
			if (pcs_type == DEFORMATION_H2)
				eqs_dim += 2 * n_linear_nodes;
			else if (pcs_type == DEFORMATION_FLOW)
				eqs_dim += n_linear_nodes;

			// number of nonzeros per row in DIAGONAL portion
			sparse_info[0] = max_connected_nodes * dim;
			//sparse_info[0] = max_connected_eles * max_ele_nodes * dim;
			// number of nonzeros per row in the OFF-DIAGONAL portion
			sparse_info[1] = sparse_info[0];
			sparse_info[2] = max_connected_nodes * dim;
			sparse_info[3] = mesh->getNumNodesLocal_Q() * dim;
			eqs = new PETScLinearSolver(eqs_dim);
			eqs->Init(sparse_info);
			eqs->set_rank_size(rank_p, size_p);
		}
		else if ((pcs_type == MULTI_PHASE_FLOW) || (pcs_type == TWO_PHASE_FLOW))
		{
			sparse_info[0] = max_connected_nodes * 2;
			sparse_info[1] = 0;  // max_cnct_nodes * 2;
			sparse_info[2] = max_connected_nodes * 2;
			sparse_info[3] = mesh->getNumNodesLocal() * 2;

			eqs = new PETScLinearSolver(2 * n_linear_nodes);
			eqs->Init(sparse_info);
			eqs->set_rank_size(rank_p, size_p);
		}
		else
		{
			sparse_info[0] = max_connected_nodes * a_pcs->GetPrimaryVNumber();
			sparse_info[1] = max_connected_nodes * a_pcs->GetPrimaryVNumber();  // 0;
			sparse_info[2] = max_connected_nodes * a_pcs->GetPrimaryVNumber();
			sparse_info[3] =
			    mesh->getNumNodesLocal() * a_pcs->GetPrimaryVNumber();

			eqs = new PETScLinearSolver(n_linear_nodes * a_pcs->GetPrimaryVNumber());
			eqs->Init(sparse_info);
			eqs->set_rank_size(rank_p, size_p);
		}
		a_pcs->setSolver(eqs);
	}
}

#endif // PETSC
