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

namespace MeshLib
{

MeshLib::CFEMesh* FEMGet(const std::string& msh_name);

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

void MSHWrite(std::string);
void CompleteMesh();
void FEMDeleteAll();

void Read_RFI(std::istream& msh_file, MeshLib::CFEMesh* m_msh);

int MSHSetMaxMMPGroups();

void MSHDelete(std::string);
//void DATWriteParticleFile(int);

MeshLib::CFEMesh* MSHGet(const std::string& mat_type_name);

long* GetPointsIn(Surface*, long*);
void GEOGetNodesInMaterialDomain(MeshLib::CFEMesh const* const,
                                        int,
                                        std::vector<long>&,
                                        bool);
void GEOGetNodesInMaterialDomain(const int MatIndex, std::vector<long>& Nodes);

long* MSHGetNodesClose(long*, CGLPolyline*);
void MSHGetNodesClose(std::vector<long>&, CGLPoint*);
long MSHGetNextNode(long startnode, MeshLib::CFEMesh* m_msh);
long* MSHGetNodesInColumn(long nextnode, int anz_zeilen,
						  MeshLib::CFEMesh* m_msh);

} // MeshLib

#endif
