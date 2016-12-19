/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef msh_tools_INC
#define msh_tools_INC

class CRFProcess;

namespace  MeshLib
{
class CFEMesh;

void MSHDefineMobile(CRFProcess*);
void MSHMoveNODUcFlow(CRFProcess*);
void DefineMobileNodes(CFEMesh*msh, CRFProcess* m_pcs);
bool MSHTestMATGroups();

} // MeshLib

#endif
