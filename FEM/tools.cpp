/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "tools.h"

#include <cfloat>

#include "display.h"
#include "FileToolsRF.h"
#include "makros.h"
#include "memory.h"

#include "msh_elem.h"
#include "msh_lib.h"

#include "femlib.h"
#include "mathlib.h"
#include "rf_fct.h"
#include "rf_mmp_new.h"
#include "rf_num_new.h"
#include "rf_pcs.h"
#include "rf_tim_new.h"
#include "sparse_table.h"

using namespace std;


/**********************************************************************
   Function for interpolation between two points P1(x1|z1) and P2(x2|z2)
   The function returns the value zn at the position xn.
   function is used in GetMatrixValue();

   Programming:
   08/2008 NB
 ***********************************************************************/
double interpol(double x1, double x2, double zx1, double zx2, double xn)
{
	if (x1 == x2)
		return zx1;
	else
		return zx1 + (xn - x1) * (zx2 - zx1) / (x2 - x1);
}

/**********************************************************************
Function GetMatrixValue (double var1, double var2, int *gueltig)

This function reads a data matrix with two independent arguments and
function values in an FCT-File and returns the value corresponding to
var1 and var2. If var1 and var2 are not in the matrix, the function
returns a value interpolated between both arguments.

Programming:
11-2011 NB/TF
***********************************************************************/
double GetMatrixValue(double var1,
                      double var2,
                      std::string const& caption,
                      int* gueltig)
{
	CFunction* matrix;
	// WW int anz_variables, anz_data;
	int dim_x, dim_y;
	int i1 = 0;
	int i2 = 0;
	int j1 = 0;
	int j2 = 0;

	matrix = FCTGet(caption);
	dim_x = matrix->matrix_dimension[0];  // NB 4.8.01
	dim_y = matrix->matrix_dimension[1];  // NB

	if (var1 < *matrix->variable_data_vector[0])  // is var1 smaller then the
	                                              // smallest argument?
	{
		*gueltig = 0;
		i1 = i2 = 0;
	}
	else if (var1 > *matrix->variable_data_vector[dim_x - 1])  // is var1 larger
	                                                           // then largest
	                                                           // argument?
	{
		*gueltig = 0;
		i1 = i2 = dim_x - 1;
	}
	else
	{
		i1 = searchElement(var1, 0, dim_x - 1, matrix->variable_data_vector);
		i2 = i1 + 1;
	}

	if (var2 < *matrix->variable_data_vector[dim_x])  // is var1 smaller then
	                                                  // the smallest argument?
	{
		*gueltig = 0;
		j1 = j2 = dim_x;
	}
	else if (var2 >
	         *matrix->variable_data_vector[dim_y + dim_x - 1])  // is var1
	                                                            // larger then
	                                                            // largest
	                                                            // argument?
	{
		*gueltig = 0;
		j1 = j2 = dim_y + dim_x - 1;
	}
	else
	{
		j1 = searchElement(
		    var2, dim_x, dim_y + dim_x - 1, matrix->variable_data_vector);
		j2 = j1 + 1;
	}

	if (fabs(var1 - *matrix->variable_data_vector[i1]) <
	    std::numeric_limits<double>::epsilon())  // var 1 is in the matrix
	{
		if (fabs(var2 - *matrix->variable_data_vector[j1]) <
		    std::numeric_limits<double>::epsilon())  // var 2 is in the matrix
		{
			return *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                     (j1 - dim_x) * dim_x];
		}
		else  // only v1 is in the matrix
		{
			double zx1y1, zx1y2;
			double y1, y2;
			zx1y1 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j1 - dim_x) * dim_x];
			zx1y2 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j2 - dim_x) * dim_x];
			y1 = *matrix->variable_data_vector[j1];
			y2 = *matrix->variable_data_vector[j2];

			return interpol(y1, y2, zx1y1, zx1y2, var2);
		}
	}
	else  // v1 is not in the matrix
	{
		if (fabs(var2 - *matrix->variable_data_vector[dim_x + j1]) <
		    std::numeric_limits<double>::epsilon())  // only var 2 is in the
		                                             // matrix
		{
			double zx1y1, zx2y1;
			double x1, x2;
			zx1y1 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j1 - dim_x) * dim_x];
			zx2y1 = *matrix->variable_data_vector[dim_x + dim_y + i2 +
			                                      (j1 - dim_x) * dim_x];
			x1 = *matrix->variable_data_vector[i1];
			x2 = *matrix->variable_data_vector[i2];

			return interpol(x1, x2, zx1y1, zx2y1, var1);
		}

		else  // neither var1 nor var2 are in the matrix
		{
			double interp1, interp2;
			double zx1y1, zx2y1, zx1y2, zx2y2;
			double x1, x2, y1, y2;
			zx1y1 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j1 - dim_x) * dim_x];
			zx2y1 = *matrix->variable_data_vector[dim_x + dim_y + i2 +
			                                      (j1 - dim_x) * dim_x];
			zx1y2 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j2 - dim_x) * dim_x];
			zx2y2 = *matrix->variable_data_vector[dim_x + dim_y + i2 +
			                                      (j2 - dim_x) * dim_x];

			x1 = *matrix->variable_data_vector[i1];
			x2 = *matrix->variable_data_vector[i2];
			y1 = *matrix->variable_data_vector[j1];
			y2 = *matrix->variable_data_vector[j2];

			interp1 = interpol(x1, x2, zx1y1, zx2y1, var1);
			interp2 = interpol(x1, x2, zx1y2, zx2y2, var1);

			return interpol(y1, y2, interp1, interp2, var2);
		}
	}
}

/**********************************************************************
   Function GetMatrixValue (double var1, double var2, int *gueltig)

   This function reads a data matrix with two independent arguments and
   function values in an FCT-File and returns the value corresponding to
   var1 and var2. If var1 and var2 are not in the matrix, the function
   returns a value interpolated between both arguments.

   Programming:
   08/2008 NB
 ***********************************************************************/
// double GetMatrixValue(double var1, double var2, std::string caption, int*
// gueltig)
//{
//	CFunction* matrix;
//	//WW int anz_variables, anz_data;
//	int dim_x, dim_y;
//	int i1 = 0;
//	int i2 = 0;
//	int j1 = 0;
//	int j2 = 0;
//	int counter;
//	double x1 = 0.0,x2 = 0.0,y1 = 0.0,y2 = 0.0; //OK411
//	double zx1y1,zx2y1,zx1y2,zx2y2;
//
//	matrix = FCTGet(caption);
//	//WW anz_variables = (int)matrix->variable_names_vector.size();
//	//dim_x = matrix->matrix_dimension_x;
//	//dim_y = matrix->matrix_dimension_y;
//	dim_x = matrix->matrix_dimension[0];  //NB 4.8.01
//	dim_y = matrix->matrix_dimension[1];  //NB
//	//WW anz_data = (int)matrix->variable_data_vector.size()-dim_x-dim_y;
//	//----------------------------------------------------------------------
//	if (var1 < *matrix->variable_data_vector[0]) //is var1 smaller then the
// smallest argument?
//	{
//		x1 = x2 = *matrix->variable_data_vector[0];
//		*gueltig = 0;
//		i1 = i2 = 0;
//	}
//	else
//	//is var1 larger then largest argument?
//	if (var1 > *matrix->variable_data_vector[dim_x - 1])
//	{
//		x1 = x2 = *matrix->variable_data_vector[dim_x - 1];
//		*gueltig = 0;
//		i1 = i2 = dim_x - 1;
//	}
//	else
//		for (counter = 0; counter < dim_x; counter++)
//		{
//			//does var1 fit an argument in the matrix exactly?
//			if (var1 == *matrix->variable_data_vector[counter])
//			{
//				x1 = x2 = *matrix->variable_data_vector[counter];
//				i1 = i2 = counter;
//				break;
//			}
//			else
//			//var1 is between two arguments in the matrix
//			if (var1 < *matrix->variable_data_vector[counter])
//			{
//				x1 = *matrix->variable_data_vector[counter - 1];
//				x2 = *matrix->variable_data_vector[counter];
//				i2 = counter;
//				i1 = i2 - 1;
//				break;
//			}
//		}
//	//same procedure for var2:
//	if (var2 < *matrix->variable_data_vector[dim_x])
//	{
//		y1 = y2 = *matrix->variable_data_vector[dim_x];
//		*gueltig = 0;
//		j1 = j2 = dim_x;
//	}
//	else
//	if (var2 > *matrix->variable_data_vector[dim_x + dim_y - 1])
//	{
//		y1 = y2 = *matrix->variable_data_vector[dim_x + dim_y - 1];
//		*gueltig = 0;
//		j1 = j2 = dim_x + dim_y - 1;
//	}
//	else
//		for (counter = dim_x; counter < dim_x + dim_y; counter++)
//		{
//			if (var2 == *matrix->variable_data_vector[counter])
//			{
//				y1 = y2 = *matrix->variable_data_vector[counter];
//				j1 = j2 = counter;
//				break;
//			}
//			else
//			if (var2 < *matrix->variable_data_vector[counter])
//			{
//				y1 = *matrix->variable_data_vector[counter - 1];
//				y2 = *matrix->variable_data_vector[counter];
//				j2 = counter;
//				j1 = j2 - 1;
//				break;
//			}
//		}
//	//getting the corresponding Z values for the arguments from the data vector
//	zx1y1 = *matrix->variable_data_vector[(j1 - dim_x) * dim_x + (i1 + dim_x +
// dim_y)];
//	zx2y1 = *matrix->variable_data_vector[(j1 - dim_x) * dim_x + (i2 + dim_x +
// dim_y)];
//	zx1y2 = *matrix->variable_data_vector[(j2 - dim_x) * dim_x + (i1 + dim_x +
// dim_y)];
//	zx2y2 = *matrix->variable_data_vector[(j2 - dim_x) * dim_x + (i2 + dim_x +
// dim_y)];
//	return interpol (y1,y2,
//	                 interpol (x1,x2,zx1y1,zx2y1,  var1),interpol (x1,
//	                                                               y1,
//	                                                               zx1y2,
//	                                                               zx2y2,
//	                                                               var1),var2);
//}



double getNodalValueFromElementValue(
    CRFProcess& pcs,
    const std::map<long, double>& map_eleId_val,
    EleToNodeInterpolationMethod::type interpolation_type,
    long node_id)
{
	MeshLib::CFEMesh* m_msh = pcs.m_msh;
	MeshLib::CNode* m_node = m_msh->nod_vector[node_id];
	double nodal_val = 0.0;

	switch (interpolation_type)
	{
		case EleToNodeInterpolationMethod::VOLUME_WEIGHTED:
		{
			const size_t n_conn_ele = m_node->getConnectedElementIDs().size();
			double sum_weights = 0.0;

			// sum the distance weighted data from each connected block
			for (size_t j = 0; j < n_conn_ele; j++)
			{
				const size_t e_id = m_node->getConnectedElementIDs()[j];
				MeshLib::CElem* e = m_msh->ele_vector[e_id];

				double ele_v = .0;
				std::map<long, double>::const_iterator itr =
				    map_eleId_val.find(e_id);

				if (itr != map_eleId_val.end())
				{
					ele_v = itr->second;
				}
				else
				{
					// skip elements which are not in the given list
					continue;
				}

				// calculate representative volume of the considered node in
				// each connected element for weighting
				double volume =
				    e->GetVolume() / ((double)e->GetNodesNumber(false));
				// weight
				double weight = 1.0 / volume;
				sum_weights += weight;

				nodal_val += ele_v * weight;
			}

			nodal_val /= sum_weights;
		}
		break;
		case EleToNodeInterpolationMethod::GAUSS_EXTRAPOLATION:
		{
			const size_t n_conn_ele = m_node->getConnectedElementIDs().size();
			if (n_conn_ele == 0) return 0;

			MeshLib::CElem* e0 =
			    m_msh->ele_vector[m_node->getConnectedElementIDs()[0]];
			if (e0->GetElementType() != MshElemType::QUAD &&
			    e0->GetElementType() != MshElemType::HEXAHEDRON)
				return 0;

			if (e0->GetElementType() == MshElemType::QUAD)
			{
				//size_t mid_nod_id = node_id;
				if (n_conn_ele < 4)
				{
					std::map<long, int> set_mid_node;
					for (size_t i = 0; i < e0->GetNodesNumber(false); i++)
					{
						long tmp_node_id = e0->GetNodeIndex(i);
						if (tmp_node_id == node_id) continue;
						if (m_msh->nod_vector[tmp_node_id]
						        ->getConnectedElementIDs()
						        .size() < 4)
							continue;
						if (set_mid_node.count(tmp_node_id) == 0)
							set_mid_node[tmp_node_id] = 0;
						else
							set_mid_node[tmp_node_id]++;
					}

					for (std::map<long, int>::iterator itr =
					         set_mid_node.begin();
					     itr != set_mid_node.end();
					     ++itr)
					{
						if (itr->second > 1)
						{
							//mid_nod_id = itr->first;
							break;
						}
					}
				}
			}

			std::vector<double> ele_values(n_conn_ele);

			for (size_t j = 0; j < n_conn_ele; j++)
			{
				const size_t e_id = m_node->getConnectedElementIDs()[j];
				MeshLib::CElem* e = m_msh->ele_vector[e_id];

				if (e->GetElementType() != MshElemType::QUAD &&
				    e->GetElementType() != MshElemType::HEXAHEDRON)
					break;

				double ele_v = .0;
				std::map<long, double>::const_iterator itr =
				    map_eleId_val.find(e_id);

				if (itr != map_eleId_val.end())
				{
					ele_v = itr->second;
				}

				ele_values[j] = ele_v;
			}
		}
		break;
	}

	return nodal_val;
}

void convertElementDataToNodalData(
    CRFProcess& pcs,
    const std::vector<double>& vec_ele_data,
    EleToNodeInterpolationMethod::type interpolation_type,
    std::vector<double>& vec_nod_values)
{
	MeshLib::CFEMesh* m_msh = pcs.m_msh;
	const size_t n_nodes = m_msh->nod_vector.size();
	vec_nod_values.resize(n_nodes);

	if (interpolation_type == EleToNodeInterpolationMethod::VOLUME_WEIGHTED)
	{
		for (size_t i = 0; i < n_nodes; i++)
		{
			MeshLib::CNode* m_node = m_msh->nod_vector[i];
			const size_t n_conn_ele = m_node->getConnectedElementIDs().size();

			double nodal_val = 0.0;
			double sum_weights = 0.0;

			// sum the distance weighted data from each connected block
			for (size_t j = 0; j < n_conn_ele; j++)
			{
				const size_t e_id = m_node->getConnectedElementIDs()[j];
				MeshLib::CElem* e = m_msh->ele_vector[e_id];

				double weight = 0.0;

				// calculate representative volume of the considered node in
				// each connected element for weighting
				double volume =
				    e->GetVolume() / ((double)e->GetNodesNumber(false));

				// weight
				weight = 1.0 / volume;
				sum_weights += weight;

				nodal_val += vec_ele_data[e_id] * weight;
			}

			nodal_val /= sum_weights;
			vec_nod_values[i] = nodal_val;
		}
	}
}

#ifdef NEW_EQS
void CreateSparseTable(MeshLib::CFEMesh* msh, Math_Group::SparseTable* &sparse_graph, Math_Group::SparseTable* &sparse_graph_H)
{
	Math_Group::StorageType stype;
	stype = Math_Group::JDS;
	for (int i = 0; i < (int)num_vector.size(); i++)
		if (num_vector[i]->ls_storage_method == 100)
		{
			stype = Math_Group::CRS;
			break;
		}

	// Symmetry case is skipped.
	// 1. Sparse_graph_H for high order interpolation. Up to now, deformation
	if (msh->GetNodesNumber(false) != msh->GetNodesNumber(true))
		sparse_graph_H = new Math_Group::SparseTable(msh, true, false, stype);
	// 2. M coupled with other processes with linear element
	if (sparse_graph_H)
	{
		if ((int)pcs_vector.size() > 1)
			sparse_graph = new Math_Group::SparseTable(msh, false, false, stype);
	}
	// 3. For process with linear elements
	else
		sparse_graph = new Math_Group::SparseTable(msh, false, false, stype);

	//  sparse_graph->Write();
	//  sparse_graph_H->Write();
	//
	// ofstream Dum("sparse.txt", ios::out);
	// sparse_graph_H->Write(Dum);
}
#endif  //#ifndef NON_PROCESS

