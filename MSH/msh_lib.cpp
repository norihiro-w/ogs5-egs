/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "msh_lib.h"

#include <set>
#include <string>
#include <vector>

#include "display.h"
#include "memory.h"

#include "mathlib.h"

#include "geo_lib.h"

std::vector<MeshLib::CFEMesh*> fem_msh_vector;

#define FEM_FILE_EXTENSION ".msh"

using namespace Math_Group;

namespace  MeshLib
{

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   04/2005 OK Implementation
**************************************************************************/
void MSHDelete(std::string m_msh_name)
{
	CFEMesh* m_fem_msh = NULL;
	size_t fem_msh_vector_size = fem_msh_vector.size();
	for (size_t i = 0; i < fem_msh_vector_size; i++)
	{
		m_fem_msh = fem_msh_vector[i];
		if (m_fem_msh->pcs_name.compare(m_msh_name) == 0)
		{
			delete m_fem_msh;
			fem_msh_vector.erase((fem_msh_vector.begin() + i));
		}
	}
}

/**************************************************************************
   FEMLib-Method:
   03/2005 OK Implementation
   05/2005 TK modified
   05/2006 TK modified
**************************************************************************/
void FEMDeleteAll()
{
	for (int i = 0; i < (int)fem_msh_vector.size(); i++)
	{
		delete fem_msh_vector[i];
		fem_msh_vector[i] = NULL;
	}
	fem_msh_vector.clear();
}

#ifndef USE_PETSC  // && not defined(other parallel method with ddc)
/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   08/2005 WW Topology construction and rfi compatible
   10/2005 OK BINARY
   08/2010 KR deleted binary mesh read
   03/2011 KR cleaned up code
   08/2011 WW Recovery multi-mesh
   09/2011 TF changed signature of function in order to read more than one mesh
**************************************************************************/
void FEMRead(const std::string& file_base_name, std::vector<CFEMesh*>& mesh_vec,
             GEOLIB::GEOObjects* geo_obj, std::string* unique_name)
{
	ScreenMessage("MSHRead:  ASCII file\n");

	CFEMesh* mesh(NULL);
	std::string msh_file_name(file_base_name + FEM_FILE_EXTENSION);

	std::ifstream msh_file_ascii(msh_file_name.data(), std::ios::in);
	if (!msh_file_ascii.is_open()) {
		ScreenMessage("*** error: Could not open file %s\n", msh_file_name.c_str());
		return;
	}

	std::string line_string("");
	getline(msh_file_ascii, line_string);

	bool more_mesh = false;                                 // 12.08.2011. WW
	if (line_string.find("#FEM_MSH") != std::string::npos)  // OGS mesh file
	{
		mesh = new CFEMesh(geo_obj, unique_name);
		more_mesh = mesh->Read(&msh_file_ascii);
		mesh_vec.push_back(mesh);  // TF

		// Multi-mesh 12.08.2011 WW
		if (more_mesh)
			while (!msh_file_ascii.eof())
			{
				// getline(msh_file_ascii, line_string);
				// if(line_string.find("#FEM_MSH")!=std::string::npos)
				mesh = new CFEMesh(geo_obj, unique_name);
				more_mesh = mesh->Read(&msh_file_ascii);
				mesh_vec.push_back(mesh);  // TF
				if (!more_mesh) break;
			}
		//  if(line_string.find("#STOP")!=std::string::npos)
		//     break;
	}
	else  // RFI mesh file
	{
		msh_file_ascii.seekg(0L, std::ios::beg);
		mesh = new CFEMesh(geo_obj, unique_name);
		Read_RFI(msh_file_ascii, mesh);
		mesh_vec.push_back(mesh);  // 12.08.2011 WW
	}

	msh_file_ascii.close();

	for (size_t i = 0; i < mesh_vec.size(); i++)
	{
		CFEMesh* msh = mesh_vec[i];
		ScreenMessage("-> Mesh %d: found %d nodes, %d elements\n",
		              i, msh->nod_vector.size(), msh->ele_vector.size());
	}
}
#endif

/**************************************************************************
   MSHLib-Method: Read rfi file ()
   Task:
   Programing:
   08/2005 WW Re-implememtation
**************************************************************************/
void Read_RFI(std::istream& msh_file, CFEMesh* m_msh)
{
	long id;
	long i = 0;
	int NumNodes = 0;
	int NumElements = 0;
	int End = 1;
	double x, y, z;
	std::string strbuffer;

	MeshLib::CNode* node = NULL;
	MeshLib::CElem* elem = NULL;
	//----------------------------------------------------------------------
	while (End)
	{
		getline(msh_file, strbuffer);  // The first line
		msh_file >> i >> NumNodes >> NumElements >> std::ws;
		//....................................................................
		// Node data
		for (i = 0; i < NumNodes; i++)
		{
			msh_file >> id >> x >> y >> z >> std::ws;
			node = new MeshLib::CNode(id, x, y, z);
			m_msh->nod_vector.push_back(node);
		}
		for (i = 0; i < NumElements; i++)
		{
			elem = new MeshLib::CElem(i);
			elem->Read(msh_file, 1);
			m_msh->ele_vector.push_back(elem);
		}
		End = 0;
	}
}

/**************************************************************************
   MSHLib-Method:
   02/2006 WW Implementation
**************************************************************************/
void CompleteMesh()
{
	for (int i = 0; i < (int)fem_msh_vector.size(); i++)
	{
		fem_msh_vector[i]->ConstructGrid();
		fem_msh_vector[i]->FillTransformMatrix();
	}
}

/**************************************************************************
   FEMLib-Method:
   Task: Master write functionn
   Programing:
   03/2005 OK Implementation
   10/2005 OK BINARY
   last modification:
   08/2010	KR binary case deleted
**************************************************************************/
void MSHWrite(std::string file_base_name)
{
	// File handling
	std::string fem_msh_file_name = file_base_name + FEM_FILE_EXTENSION;
	std::ofstream fem_msh_file(fem_msh_file_name.c_str(),
	                           std::ios::trunc | std::ios::out);

	if (!fem_msh_file.good()) return;

	for (size_t i = 0; i < fem_msh_vector.size(); i++)
	{
		FileIO::OGSMeshIO meshIO;
		meshIO.setPrecision(12);
		meshIO.setFormat(std::ios::scientific);
		meshIO.setMesh(fem_msh_vector[i]);
		fem_msh_file << meshIO.writeToString();
	}
	fem_msh_file << "#STOP";
	fem_msh_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   03/2005 OK Implementation
   last modification:
**************************************************************************/
MeshLib::CFEMesh* FEMGet(const std::string& msh_name)
{
	size_t no_msh = fem_msh_vector.size();
	// If there is only one msh file available, use it for all process. WW
	if (no_msh == 1) return fem_msh_vector[0];  // WW
	for (size_t i = 0; i < no_msh; i++)
		if (fem_msh_vector[i]->pcs_name.compare(msh_name) == 0)
			return fem_msh_vector[i];
	return NULL;
}

/**************************************************************************
   MSHLib-Method:
   12/2005 OK Implementation
   07/2007 OK PCS
**************************************************************************/
MeshLib::CFEMesh* MSHGet(const std::string& geo_name)
{
	MeshLib::CFEMesh* m_msh = NULL;
	for (int i = 0; i < (int)fem_msh_vector.size(); i++)
	{
		m_msh = fem_msh_vector[i];
		if (m_msh->geo_name.compare(geo_name) == 0) return m_msh;
		if (m_msh->pcs_name.compare(geo_name) == 0) return m_msh;
	}
	return NULL;
}

/**************************************************************************/
/* ROCKFLOW - Function: MSHGetNextNode
 */
/* Task:
   Find the next node to the starting node in user defined direction
 */
/* Parameter: (I: Input; R: Return; X: Both)
   I: startnode, direction
 */
/* Return:
   nextnode
 */
/* Programming:
   09/2002     MB         First Version
   08/2005     MB
 */
/**************************************************************************/
long MSHGetNextNode(long startnode, MeshLib::CFEMesh* m_msh)
{
	size_t NumberOfNodes(m_msh->nod_vector.size());
	long NumberOfNodesPerLayer =
	    NumberOfNodes / (m_msh->getNumberOfMeshLayers() + 1);
	return startnode + NumberOfNodesPerLayer;
}


/**************************************************************************
   ROCKFLOW - Function: MSHGetNodesInColumn

   Task:
   Gets nodes of a column searching downward from startnode.

   Parameter: (I: Input; R: Return; X: Both)
           I: long node, int anz_zeilen

   Return:
  *long strang

   Programming:
   09/2002   MB   First Version
   08/2005   MB   m_msh
**************************************************************************/
long* MSHGetNodesInColumn(long nextnode, int anz_zeilen,
                          MeshLib::CFEMesh* m_msh)
{
	int i;
	long startnode;
	long* strang = NULL;

	for (i = 0; i < anz_zeilen + 1; i++)
	{
		strang = (long*)Realloc(strang, (i + 1) * sizeof(long));
		strang[i] = nextnode;
		startnode = nextnode;
		// nextnode = MSHGetNextNode (startnode, direction);
		nextnode = MSHGetNextNode(startnode, m_msh);
	}
	return strang;
}


/**************************************************************************
   MSHLib-Method:
   Task:
   Programing:
   09/2004 OK Implementation
   ToDo evtl. vector<CGLPoint>
   08/2005 CC Modification: CGLPoint* e_pnt - Move from GeoLib to MshLib
**************************************************************************/
void MSHGetNodesClose(std::vector<long>& msh_point_vector, CGLPoint* e_pnt)
{
	e_pnt = e_pnt;
	msh_point_vector.size();
	/*OK411
	   long i;
	   CGLPoint m_pnt;
	   // Node loop
	   for (i=0;i<NodeListSize();i++) {
	    if (GetNode(i)==NULL) continue;
	    m_pnt.x = GetNodeX(i);
	    m_pnt.y = GetNodeY(i);
	    m_pnt.z = GetNodeZ(i);
	    if(e_pnt->PointDis(&m_pnt)<=(e_pnt->epsilon+MKleinsteZahl))
	      msh_point_vector.push_back(i);
	   }
	 */
}

/*************************************************************************
   ROCKFLOW - Function: MSHGetNodesClose
   Task: Searching grid points which are close to a polyline
   Programming:
   10/2002 OK Encapsulated from ExecuteSourceSinkMethod11 (CT)
   01/2003 OK Test
   last modified: 20.01.2003 OK
   08/2005 CC Modification Move from GeoLib to MSHLib
 **************************************************************************/
long* MSHGetNodesClose(long* number_of_nodes, CGLPolyline* m_ply)
{
	long* nodes_all = NULL;
	m_ply = m_ply;
	number_of_nodes = number_of_nodes;
	/*OK411
	   long j,k,l;
	   double pt1[3],line1[3],line2[3],pt0[3];
	   double mult_eps = 1.0;
	   double dist1p,dist2p,*length,laenge;
	   long anz_relevant = 0;
	   typedef struct {
	     long knoten;
	     long abschnitt;
	     double laenge;
	   } INFO;
	   INFO *relevant=NULL;
	   int weiter;
	   double w1,w2;
	   long knoten_help;
	   double laenge_help;
	   double gesamte_laenge = 0.;
	   long polyline_point_vector_size;

	   m_ply->sbuffer.clear();
	   m_ply->ibuffer.clear();

	   if (m_ply) {

	   length = (double*) Malloc(sizeof(double)
	   *(long)m_ply->point_vector.size());

	   pt0[0] = m_ply->point_vector[0]->x;
	   pt0[1] = m_ply->point_vector[0]->y;
	   pt0[2] = m_ply->point_vector[0]->z;

	   polyline_point_vector_size =(long)m_ply->point_vector.size();
	   for (k=0;k<polyline_point_vector_size-1;k++) {
	   line1[0] = m_ply->point_vector[k]->x;
	   line1[1] = m_ply->point_vector[k]->y;
	   line1[2] = m_ply->point_vector[k]->z;
	   line2[0] = m_ply->point_vector[k+1]->x;
	   line2[1] = m_ply->point_vector[k+1]->y;
	   line2[2] = m_ply->point_vector[k+1]->z;
	   length[k] = MCalcDistancePointToPoint(line2, line1);
	   gesamte_laenge += length[k];
	   }

	   // Wiederholen bis zumindest ein Knoten gefunden wurde
	   while(anz_relevant==0) {

	   for (j=0;j<NodeListSize();j++) {
	   if (GetNode(j)==NULL) continue;

	   polyline_point_vector_size =(long)m_ply->point_vector.size();
	   for (k=0;k<polyline_point_vector_size-1;k++) {

	   pt1[0] = GetNodeX(j);
	   pt1[1] = GetNodeY(j);
	   pt1[2] = GetNodeZ(j);

	   line1[0] = m_ply->point_vector[k]->x;
	   line1[1] = m_ply->point_vector[k]->y;
	   line1[2] = m_ply->point_vector[k]->z;
	   line2[0] = m_ply->point_vector[k+1]->x;
	   line2[1] = m_ply->point_vector[k+1]->y;
	   line2[2] = m_ply->point_vector[k+1]->z;

	   if ( MCalcDistancePointToLine(pt1,line1,line2) <= mult_eps*m_ply->epsilon
	   ) {
	   MCalcProjectionOfPointOnLine(pt1,line1,line2,pt1);
	   dist1p = MCalcDistancePointToPoint(line1, pt1);
	   dist2p = MCalcDistancePointToPoint(line2, pt1);
	   if ((dist1p+dist2p-length[k]) <=  mult_eps*m_ply->epsilon ) {

	   // For boundary conditions. WW
	   m_ply->sbuffer.push_back(dist1p);
	   m_ply->ibuffer.push_back(k);
	   // ---------------------------

	   anz_relevant++;
	   nodes_all = (long *) Realloc(nodes_all,sizeof(long)*anz_relevant);
	   relevant = (INFO *) Realloc(relevant, sizeof(INFO) * anz_relevant);
	   nodes_all[anz_relevant-1] = j;
	   laenge = 0.;
	   for (l=0; l < k; l++)
	   laenge += length[l];
	   relevant[anz_relevant-1].knoten = j;
	   relevant[anz_relevant-1].laenge = laenge + dist1p;
	   k =(long)m_ply->point_vector.size();
	   }
	   }
	   }
	   }
	   if(anz_relevant==0) mult_eps *=2.;
	   }

	   if (mult_eps > 1.)
	   std::cout << "!!! Epsilon increased in sources!" << endl;

	   do {
	   weiter = 0;
	   for (k=0;k<anz_relevant-1;k++) {
	   w1=relevant[k].laenge;
	   w2=relevant[k+1].laenge;
	   if (w1>w2) { // Die Eintraege vertauschen
	   knoten_help = relevant[k].knoten;
	   laenge_help = relevant[k].laenge;
	   relevant[k].knoten = relevant[k+1].knoten;
	   relevant[k].laenge = relevant[k+1].laenge;
	   relevant[k+1].knoten = knoten_help;
	   relevant[k+1].laenge = laenge_help;
	   weiter=1;
	   }
	   }
	   } while (weiter);

	   relevant = (INFO*) Free(relevant);
	   *number_of_nodes = anz_relevant;
	   }
	 */
	return nodes_all;
}

/**************************************************************************
   GeoLib-Method: GetPointsIn
   Task:
   Programing:
   01/2004 OK Implementation
   08/2005 CC Modification Move from Geolib to Mshlib
**************************************************************************/
long* GetPointsIn(Surface* m_sfc, long* number_of_nodes)
{
	long* nodes = NULL;
	number_of_nodes = number_of_nodes;
	m_sfc = m_sfc;
	/*OK411
	   long i;
	   double *xp=NULL,*yp=NULL,*zp=NULL;
	   long anz_relevant = 0;
	   CGLPoint m_pnt;
	   // Inside polygon
	   if(!m_sfc->polygon_point_vector.empty()) {
	    xp = (double*)
	   Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
	    yp = (double*)
	   Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
	    zp = (double*)
	   Malloc(((long)m_sfc->polygon_point_vector.size())*sizeof(double));
	    long polygon_point_vector_length =
	   (long)m_sfc->polygon_point_vector.size();
	   for(i=0;i<polygon_point_vector_length;i++) {
	   xp[i] = m_sfc->polygon_point_vector[i]->x;
	   yp[i] = m_sfc->polygon_point_vector[i]->y;
	   zp[i] = m_sfc->polygon_point_vector[i]->z;
	   }

	   //-----------------------------------------------------------------
	   for(i=0;i<NodeListSize();i++) {
	   if (GetNode(i)==NULL) continue;
	   m_pnt.x = GetNodeX(i);
	   m_pnt.y = GetNodeY(i);
	   m_pnt.z = GetNodeZ(i);
	   if(m_pnt.IsInsidePolygonPlain(
	   xp,yp,zp,\
	   (long)m_sfc->polygon_point_vector.size())) {
	   anz_relevant++;
	   nodes = (long *) Realloc(nodes,sizeof(long)*anz_relevant);
	   nodes[anz_relevant-1] = i;
	   }
	   }
	   }
	   // Destructions
	   // nodes extern
	   xp = (double*) Free(xp);
	   yp = (double*) Free(yp);
	   zp = (double*) Free(zp);
	   //
	   *number_of_nodes = anz_relevant;
	 */
	return nodes;
}

/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: Get element nodes in a material domain
   Programing:
   10/2004 WW Implementation
**************************************************************************/
void GEOGetNodesInMaterialDomain(const int MatIndex, std::vector<long>& Nodes)
{
	(void)MatIndex;
	(void)Nodes;
	/*OK411
	   MatIndex;
	   Nodes.size();
	   long index, *element_nodes;
	   int i, j, Size, nn, order = 2;
	   const int L_Nodes = GetLowOrderNodeNumber();
	   bool exist;
	   if(L_Nodes==NodeListSize()) order = 1;
	   if(L_Nodes==0) order = 1;

	   Nodes.resize(0);
	   nn = 0;
	   for (index=0;index<ElListSize();index++)
	   {
	   if (ElGetElement(index)!=NULL)
	   {  // Eelement exist
	   if (ElGetElementActiveState(index))
	   {  // Element active
	   if(order==1) nn = NumbersOfElementNode(index);
	   if(order==2) nn = NumbersOfElementNodeHQ(index);
	   if(ElGetElementGroupNumber(index)==MatIndex)
	   {
	   Size = (int)Nodes.size();
	   element_nodes = ElGetElementNodes(index);
	   for(i=0; i<nn; i++)
	   {
	   exist = false;
	   for(j=0; j<Size; j++)
	   {
	   if(element_nodes[i]==Nodes[j])
	   {
	   exist = true;
	   break;
	   }
	   }
	   if(!exist) Nodes.push_back(element_nodes[i]);
	   }
	   }
	   }
	   }
	   }
	 */
}

/**************************************************************************
   GeoLib-Method: GEOGetVolume
   Task: Get element nodes in a material domain
   Programing:
   10/2004 WW Implementation
   06/2012 NW Made this function faster using std::set for large data set
**************************************************************************/
void GEOGetNodesInMaterialDomain(CFEMesh const* const msh,
                                 int MatIndex,
                                 std::vector<long>& Nodes,
                                 bool Order)
{
	Nodes.resize(0);
	std::set<long> set_nodes;
	const size_t n_ele(msh->ele_vector.size());
	for (size_t e = 0; e < n_ele; e++)
	{
		MeshLib::CElem const* const elem(msh->ele_vector[e]);
		if (elem->GetMark())  // Marked for use
		{
			const int nn(elem->GetNodesNumber(Order));
			if (elem->GetPatchIndex() == static_cast<size_t>(MatIndex))
			{
				for (int i = 0; i < nn; i++)
				{
					set_nodes.insert(elem->GetNodeIndex(i));
				}
			}
		}  // if
	}      // For
	Nodes.assign(set_nodes.begin(), set_nodes.end());
}

/**************************************************************************
   MSHLib-Method:
   01/2006 OK Implementation
   06/2009 OK Bug fix
**************************************************************************/
int MSHSetMaxMMPGroups()
{
	int i;
	long j;
	CFEMesh* m_msh = NULL;
	//----------------------------------------------------------------------
	size_t msh_max_mmp_groups;
	for (i = 0; i < (int)fem_msh_vector.size(); i++)
	{
		m_msh = fem_msh_vector[i];
		m_msh->max_mmp_groups = 0;
		msh_max_mmp_groups = 0;
		for (j = 0; j < (long)m_msh->ele_vector.size(); j++)
			if ((m_msh->ele_vector[j]->GetPatchIndex() + 1) >
			    msh_max_mmp_groups)
				msh_max_mmp_groups++;
		m_msh->max_mmp_groups = msh_max_mmp_groups;
	}
	//----------------------------------------------------------------------
	size_t g_msh_max_mmp_groups = 0;
	for (i = 0; i < (int)fem_msh_vector.size(); i++)
		if (m_msh->max_mmp_groups > g_msh_max_mmp_groups)
			g_msh_max_mmp_groups++;
	//----------------------------------------------------------------------
	return g_msh_max_mmp_groups;
}

} // MeshLib
