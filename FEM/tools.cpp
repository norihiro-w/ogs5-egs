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

#include "Curve.h"
#include "mathlib.h"

#include "msh_elem.h"
#include "msh_lib.h"

#include "files0.h"
#include "rf_fct.h"
#include "rf_mmp_new.h"
#include "rf_num_new.h"
#include "rf_tim_new.h"

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


/****************************************************************************
 * Finds and returns the positive minimum of a vector.
 * Programming: NB Dec 08
 *****************************************************************************/
double FindMin(vector<double> Vec)
{
	double x = DBL_MAX;
	int unsigned i;

	for (i = 0; i < Vec.size(); i++)
		if ((Vec[i] >= 0) && (Vec[i] < x)) x = Vec[i];
	return x;
}

/****************************************************************************
 * Finds and returns the maximum of a vector.
 * Programming: NB Jan 08
 *****************************************************************************/
double FindMax(vector<double> Vec)
{
	double x = DBL_MIN;
	int unsigned i;

	for (i = 0; i < Vec.size(); i++)
		if (Vec[i] > x) x = Vec[i];
	return x;
}

/****************************************************************************
 * Finds all real roots of a third grade polynomial in the form:
 * P(x) = x^3 + px^2 + qx + r
 * roots are returned in a vector
 *
 * Programming: NB, Dec08
 *****************************************************************************/
void NsPol3(double p, double q, double r, vector<double>* roots)
{
	double eps = 7E-15;
	double a, b, h, phi, D, z[3];
	double pi = 3.1415926535897;
	double nz;
	int i;

	b = (p / 3) * (p / 3);
	a = q / 3 - b;
	b = b * p / 3 + 0.5 * (r - p / 3 * q);
	h = sqrt(fabs(a));

	if (b < 0) h = -h;

	D = MathLib::fastpow(a, 3) + b * b;

	if (D <= (-eps))
	{
		nz = 3;
		phi = acos(b / MathLib::fastpow(h, 3)) / 3;
		z[0] = 2 * h * cos(pi / 3 - phi) - p / 3;
		z[1] = 2 * h * cos(pi / 3 + phi) - p / 3;
		z[2] = -2 * h * cos(phi) - p / 3;
	}
	else if (D < eps)
	{
		nz = 3;
		z[0] = -2 * h - p / 3;
		z[1] = h - p / 3;
		z[2] = z[1];
	}
	else
	{
		nz = 1;
		if (a >= eps)
		{
			b = b / MathLib::fastpow(h, 3);
			phi = log(b + sqrt(b * b + 1)) / 3;
			z[0] = -2 * h * sinh(phi) - p / 3;
		}
		else if (a > (-eps))
		{
			z[0] = pow((2 * abs(b)), 1. / 3.);
			if (b > 0) z[0] = -z[0];
			z[0] = z[0] - p / 3;
		}
		else
		{
			b = b / MathLib::fastpow(h, 3);
			phi = log(b + sqrt(b * b - 1)) / 3;
			z[0] = -2 * h * cosh(phi) - p / 3;
		}
	}

	for (i = 0; i < nz; i++)
		roots->push_back(z[i]);
}

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
				size_t mid_nod_id = node_id;
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
							mid_nod_id = itr->first;
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
Math_Group::SparseTable* createSparseTable(MeshLib::CFEMesh* a_mesh,
                         bool quadratic,
                         bool symm)
{
	Math_Group::SparseTable* st = new Math_Group::SparseTable();
	st->symmetry = symm;

	//
	// In sparse table, = number of nodes
	st->rows = a_mesh->GetNodesNumber(quadratic);
	st->size_entry_column = 0;
	st->diag_entry = new long[st->rows];

	long** larraybuffer = nullptr;
	if (st->symmetry)
	{
		larraybuffer = new long* [st->rows];
		for (long i = 0; i < st->rows; i++)
		{
			MeshLib::CNode* node = a_mesh->nod_vector[i];
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			const long lbuff1 = (long)node->getConnectedNodes().size();
			larraybuffer[i] = new long[lbuff1 + 1];
			//
			larraybuffer[i][0] = lbuff1;
			for (long j = 0; j < lbuff1; j++)
				larraybuffer[i][j + 1] = node->getConnectedNodes()[j];
			node->getConnectedNodes().clear();
			for (long j = 0; j < lbuff1; j++)
			{
				long jj = larraybuffer[i][j + 1];
				if (i <= jj)
					node->getConnectedNodes().push_back(jj);
			}
		}
	}

	// num_column_entries saves vector ptr of CRS
	st->num_column_entries = new long[st->rows + 1];

	std::vector<long> A_index;

	for (long i = 0; i < st->rows; i++)
	{
		st->num_column_entries[i] = (long)A_index.size();

		MeshLib::CNode const* node = a_mesh->nod_vector[i];
		std::vector<size_t> vec_node_ids(node->getConnectedNodes());
		vec_node_ids.push_back(i);
		std::sort(vec_node_ids.begin(), vec_node_ids.end());
		for (size_t j = 0; j < vec_node_ids.size(); j++)
		{
			long col_index = vec_node_ids[j];

			if ((!quadratic) && (col_index >= st->rows))
				continue;

			if (i == col_index)
				st->diag_entry[i] = (long)A_index.size();
			A_index.push_back(col_index);
		}
	}

	st->size_entry_column = (long)A_index.size();
	st->num_column_entries[st->rows] = st->size_entry_column;

	st->entry_column = new long[st->size_entry_column];
	for (long i = 0; i < st->size_entry_column; i++)
		st->entry_column[i] = A_index[i];

	// For the case of symmetry matrix
	if (st->symmetry)
	{
		for (long i = 0; i < st->rows; i++)
		{
			long lbuff0 = larraybuffer[i][0];
			a_mesh->nod_vector[i]->getConnectedNodes().resize(lbuff0);
			//
			for (long j = 0; j < lbuff0; j++)
				a_mesh->nod_vector[i]->getConnectedNodes()[j] =
				    larraybuffer[i][j + 1];
		}
		for (long i = 0; i < st->rows; i++)
		{
			delete[] larraybuffer[i];
			larraybuffer[i] = nullptr;
		}
		delete[] larraybuffer;
		larraybuffer = nullptr;
	}

	return st;
}

void CreateSparseTable(MeshLib::CFEMesh* msh, Math_Group::SparseTable* &sparse_graph, Math_Group::SparseTable* &sparse_graph_H)
{
	// Symmetry case is skipped.
	// 1. Sparse_graph_H for high order interpolation. Up to now, deformation
	if (msh->GetNodesNumber(false) != msh->GetNodesNumber(true))
		sparse_graph_H = createSparseTable(msh, true, false);
	// 2. M coupled with other processes with linear element
	if (sparse_graph_H)
	{
		if ((int)pcs_vector.size() > 1)
			sparse_graph = createSparseTable(msh, false, false);
	}
	// 3. For process with linear elements
	else
		sparse_graph = createSparseTable(msh, false, false);

	//  sparse_graph->Write();
	//  sparse_graph_H->Write();
	//
	// ofstream Dum("sparse.txt", ios::out);
	// sparse_graph_H->Write(Dum);
}
#endif  //#ifndef NON_PROCESS

