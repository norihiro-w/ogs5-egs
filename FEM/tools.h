/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef tools_INC
#define tools_INC

#include "rf_pcs.h"
#include "sparse_matrix.h"

double GetMatrixValue(double var1,
                      double var2,
                      std::string const& caption,
                      int* gueltig);

double FindMin(std::vector<double> Vec);
double FindMax(std::vector<double> Vec);

/**
 * Author: NB 4.9.05
 * Finds all real roots of a third grade polynomial in the form:
 * P(x) = x^3 + px^2 + qx + r
 * employing the <a
 * href="http://de.wikipedia.org/wiki/Kubische_Gleichung#Analytische_Bestimmung_reeller_L.C3.B6sungen">Cardano
 * method</a>
 * @param p coefficient of cubic equation
 * @param q coefficient of cubic equation
 * @param r coefficient of cubic equation
 * @param t vector contains the roots
 */
void NsPol3(double p, double q, double r, std::vector<double>* t);

struct EleToNodeInterpolationMethod
{
	enum type
	{
		VOLUME_WEIGHTED,
		GAUSS_EXTRAPOLATION
	};
};

double getNodalValueFromElementValue(
    CRFProcess& pcs,
    const std::map<long, double>& map_eleId_val,
    EleToNodeInterpolationMethod::type interpolation_type,
    long node_id);
void convertElementDataToNodalData(
    CRFProcess& pcs,
    const std::vector<double>& vec_ele_data,
    EleToNodeInterpolationMethod::type interpolation_type,
    std::vector<double>& vec_nod_data);

#ifdef NEW_EQS
void CreateSparseTable(MeshLib::CFEMesh* msh, Math_Group::SparseTable* &sparse_graph, Math_Group::SparseTable* &sparse_graph_H);
#endif

#endif
