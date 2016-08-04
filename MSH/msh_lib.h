/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   MSHLib - Object:
   Task:
   Programing:
   08/2005 WW/OK Encapsulation from rf_ele_msh
   last modified
**************************************************************************/
#ifndef msh_lib_INC
#define msh_lib_INC

#include "msh_mesh.h"

extern std::vector<MeshLib::CFEMesh*> fem_msh_vector;

extern MeshLib::CFEMesh* FEMGet(const std::string& msh_name);

/**
 * reads a mesh. The following formats are possible:
 * <ol>
 * <li>native OGS meshes</li>
 * <li>TetGen</li>
 * <li>GMSH</li>
 * </ol>
 * @param mesh_fname file name of the mesh
 * @param mesh_vec a vector, the new mesh will be put in this vector
 * @param geo_obj object, that manages the geometric entities
 * @param unique_name the name of geometric data
 */
void FEMRead(const std::string& mesh_fname,
             std::vector<MeshLib::CFEMesh*>& mesh_vec,
             GEOLIB::GEOObjects* geo_obj = NULL,
             std::string* unique_name = NULL);

void Read_RFI(std::istream& msh_file, MeshLib::CFEMesh* m_msh);

extern void CompleteMesh();
extern void FEMDeleteAll();

extern int MSHSetMaxMMPGroups();
extern bool MSHTestMATGroups(size_t mmp_vector_size);

extern void GEOGetNodesInMaterialDomain(MeshLib::CFEMesh const* const,
										int,
										std::vector<long>&,
										bool);

#endif
