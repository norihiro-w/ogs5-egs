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

std::vector<MeshLib::CFEMesh*> fem_msh_vector;

#define FEM_FILE_EXTENSION ".msh"


using namespace Math_Group;
using namespace MeshLib;

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
   07/2007 OK Implementation
**************************************************************************/
bool CompleteMesh(std::string pcs_name)
{
	bool succeed = false;
	for (int i = 0; i < (int)fem_msh_vector.size(); i++)
		if (fem_msh_vector[i]->pcs_name.compare(pcs_name) == 0)
		{
			fem_msh_vector[i]->ConstructGrid();
			fem_msh_vector[i]->FillTransformMatrix();
			succeed = true;
		}
	return succeed;
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

/**************************************************************************
   MSHLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool MSHTestMATGroups(size_t mmp_vector_size)
{
	int g_max_mmp_groups = MSHSetMaxMMPGroups();
	if (g_max_mmp_groups > (int)mmp_vector_size)
	{
		std::cout << "Error: not enough MMP data";
		return false;  // abort();
	}
	return true;
}
