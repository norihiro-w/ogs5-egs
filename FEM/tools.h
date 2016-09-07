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

#include <map>
#include <string>
#include <vector>

namespace Math_Group {struct SparseTable; }
namespace MeshLib {class CFEMesh; }
class CRFProcess;

double GetMatrixValue(double var1,
                      double var2,
                      std::string const& caption,
                      int* gueltig);

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
