/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

//**************************************************************************
/*!
   \file msh_lib_ext.cpp

   Define the member function of class CFEMesh for reading
   subdomain mesh

   02/2012 WW/
   last modified
*/
//**************************************************************************
#include "msh_lib.h"

//#ifdef USET_PETSC
// kg44 included for memcpy needed by petsc
#include <stdio.h>
#include <string.h>
#include "petscksp.h"
//#elif USE_MPI
#include <mpi.h>
//#endif
#include <sstream>

#include "StringTools.h"
#include "display.h"

using namespace std;
using namespace MeshLib;

extern int myrank;
extern int mysize;

void BuildNodeStruc(MeshNodes* anode, MPI_Datatype* MPI_Node_ptr);

void FEMRead(const string& file_base_name, vector<MeshLib::CFEMesh*>& mesh_vec,
             GEOLIB::GEOObjects* geo_obj, string* unique_name)
{
	ScreenMessage("MSHRead\n");

	// 0 long size_sbd_nodes = 0;
	// 1 long size_sbd_nodes_l = 0;
	// 2 long size_sbd_nodes_h = 0;
	// 3 long size_sbd_elems = 0;
	// 4 long size_g_elems = 0;
	static const int nheaders = 10;
	int mesh_header[nheaders];

	MeshNodes* s_nodes = (MeshNodes*)malloc(sizeof(MeshNodes));
	int* elem_info = (int*)malloc(1);

	MPI_Datatype MPI_node;
	int tag[] = {0, 1, 2};
	MPI_Status status;

	ifstream is;
#ifdef MULTI_MESH_FILE
	stringstream ss(stringstream::in | stringstream::out);
#endif

	if (myrank == 0)
	{
		std::string str_var = file_base_name + "_partitioned.msh";
		is.open(str_var.c_str());
		if (!is.good())
		{
			is.close();
			is.clear();
			str_var =
			    file_base_name + "_partitioned_" + number2str(mysize) + ".msh";
			is.open(str_var.c_str());
			if (is.good())
			{
				ScreenMessage("-> found %s\n", str_var.c_str());
			}
			else
			{
				ScreenMessage("-> cannot find a partitioned mesh file\n");
				PetscFinalize();
				exit(1);
			}
		}
		getline(is, str_var);
		int num_parts;
		is >> num_parts >> ws;
		if (num_parts != mysize)
		{
			string str_m =
			    "Sorry, I have to quit the simulation now because that "
			    " the number of the requested computer cores "
			    "is not identical to the number of subdomains.";
			cout << str_m << endl;
			ScreenMessage("(MPI size = %d, the number of subdomains = %d)\n",
			              mysize, num_parts);
			PetscFinalize();
			exit(1);
		}
	}

	CFEMesh* mesh = new CFEMesh(geo_obj, unique_name);
	mesh_vec.push_back(mesh);

	MPI_Barrier(MPI_COMM_WORLD);
	ScreenMessage("-> Parallel reading the partitioned mesh\n");
	for (int i = 0; i < mysize; i++)
	{
		if (myrank == 0)
		{
#ifdef MULTI_MESH_FILE
			ss.clear();
			ss.str("");
			ss << i;
			str_var = file_base_name + "_" + ss.str() + ".msh";

			is.open(str_var.c_str());
			getline(is, str_var);
#endif

			// cout<<"-->Parallel reading the partitioned mesh: "<<i<<endl;

			for (int j = 0; j < nheaders; j++)
				is >> mesh_header[j];
			is >> ws;
		}
		MPI_Bcast(mesh_header, nheaders, MPI_INT, 0, MPI_COMM_WORLD);
		if (i == myrank)
		{
			std::stringstream ss;
			for (int j = 0; j < nheaders; j++)
				ss << mesh_header[j] << " ";
			ScreenMessage2d("-> header: %s\n", ss.str().c_str());
		}
		// cout<<"\ncccccccccc "<<mesh_header[0]<<"    "<<mesh_header[1]
		//	<<"   " <<mesh_header[2]<<"  "<<mesh_header[3]<<endl;

		//-------------------------------------------------------------------------
		// Node
		s_nodes =
		    (MeshNodes*)realloc(s_nodes, sizeof(MeshNodes) * mesh_header[0]);
		if (i > 0) BuildNodeStruc(s_nodes, &MPI_node);

		if (myrank == 0)
		{
			// Nodes
			for (int k = 0; k < mesh_header[0]; k++)
			{
				MeshNodes* anode = &s_nodes[k];
				is >> anode->index;
				is >> anode->x >> anode->y >> anode->z >> ws;
			}
			/*
			//TEST
			for(k=0; k<mesh_header[0]; k++)
			{
			MeshNodes *anode = &s_nodes[k];
			cout<<anode->index<<" ";

			cout<<anode->x<<" "<<anode->y<<" "<<anode->z<<endl;

			}
			*/
			if (i == 0)
			{
				mesh->setSubdomainNodes(&mesh_header[0], s_nodes);
			}
			else
			{
				MPI_Send(s_nodes, mesh_header[0], MPI_node, i, tag[0],
				         MPI_COMM_WORLD);
			}
		}

		if (i > 0)
		{
			if (myrank == i)
			{
				MPI_Recv(s_nodes, mesh_header[0], MPI_node, 0, tag[0],
				         MPI_COMM_WORLD, &status);
				mesh->setSubdomainNodes(mesh_header, s_nodes);

				/*
				//TEST
				for(k=0; k<mesh_header[0]; k++)
				{
				    MeshNodes *anode = &s_nodes[k];
				    cout<<anode->index<<" ";
				    cout<<anode->x<<" "<<anode->y<<" "<<anode->z<<endl;
				}
				*/
			}
		}

		//-------------------------------------------------------------------------
		// Element
		const int size_elem_info = mesh_header[2] + mesh_header[8];
		elem_info = (int*)realloc(elem_info, sizeof(int) * size_elem_info);
		for (int ii = 0; ii < size_elem_info; ii++)
			elem_info[ii] = -1;
		if (myrank == 0)
		{
			int counter = mesh_header[2];
			for (int j = 0; j < mesh_header[2]; j++)
			{
				elem_info[j] = counter;
				is >> elem_info[counter];  // mat. idx
				counter++;
				is >> elem_info[counter];  // type
				counter++;
				is >> elem_info[counter];  // nnodes
				const int nn_e = elem_info[counter];
				counter++;
				for (int k = 0; k < nn_e; k++)
				{
					is >> elem_info[counter];
					counter++;
				}
			}

			if (i == 0)
				mesh->setSubdomainElements(mesh_header, elem_info, true);
			else
				MPI_Send(elem_info, size_elem_info, MPI_INT, i, tag[1],
				         MPI_COMM_WORLD);
		}
		if (i > 0)
		{
			if (myrank == i)
			{
				MPI_Recv(elem_info, size_elem_info, MPI_INT, 0, tag[1],
				         MPI_COMM_WORLD, &status);
				mesh->setSubdomainElements(mesh_header, elem_info, true);
			}
		}

		// if(elem_info)
		//  {
		// free(elem_info);
		// elem_info = NULL;
		//  }

		//-------------------------------------------------------------------------
		// Ghost element
		const int size_elem_g_info = mesh_header[3] + mesh_header[9];
		elem_info = (int*)realloc(elem_info, sizeof(int) * size_elem_g_info);
		if (myrank == 0)
		{
			int counter = mesh_header[3];
			for (int j = 0; j < mesh_header[3]; j++)
			{
				elem_info[j] = counter;
				is >> elem_info[counter];  // mat. idx
				counter++;
				is >> elem_info[counter];  // type
				counter++;
				is >> elem_info[counter];  // nnodes
				const int nn_e = elem_info[counter];
				counter++;
				for (int k = 0; k < nn_e; k++)
				{
					is >> elem_info[counter];
					counter++;
				}
				is >> elem_info[counter];
				// const int nn_e_g =  elem_info[counter];
				counter++;
				// ghost nodes for linear element
				is >> elem_info[counter];
				const int nn_e_g_quad = elem_info[counter];
				counter++;
				for (int k = 0; k < nn_e_g_quad;
				     k++)  // NW use nn_e_g_quad instead of nn_e_g
				{
					is >> elem_info[counter];
					counter++;
				}
			}

			if (i == 0)
				mesh->setSubdomainElements(mesh_header, elem_info, false);
			else
				MPI_Send(elem_info, size_elem_g_info, MPI_INT, i, tag[2],
				         MPI_COMM_WORLD);
		}
		if (i > 0)
		{
			if (myrank == i)
			{
				MPI_Recv(elem_info, size_elem_g_info, MPI_INT, 0, tag[2],
				         MPI_COMM_WORLD, &status);
				mesh->setSubdomainElements(mesh_header, elem_info, false);
			}
		}
	}

	if (s_nodes)
	{
		free(s_nodes);
		s_nodes = NULL;
	}
	if (elem_info)
	{
		free(elem_info);
		elem_info = NULL;
	}

	if (myrank == 0)
	{
		is.clear();
		is.close();
	}

	if (mysize > 1) MPI_Type_free(&MPI_node);

	MPI_Barrier(MPI_COMM_WORLD);
	ScreenMessage2(
	    "-> nelements_g=%d, nnodes_gl=%d, nnodes_gq=%d, nnodes_ll=%d, "
	    "nnodes_lq=%d\n",
	    mesh->getElementVector().size(), mesh->GetNodesNumber(false),
	    mesh->GetNodesNumber(true), mesh->getNumNodesLocal(),
	    mesh->getNumNodesLocal_Q());
	MPI_Barrier(MPI_COMM_WORLD);

	mesh->ConstructGrid();
	mesh->FillTransformMatrix();
	// mesh->calMaximumConnectedNodes();
	MPI_Barrier(MPI_COMM_WORLD);
}

namespace MeshLib
{
/*!
   Fill data for subdomain mesh

   @param header  : mesh header
   @param s_nodes : mesh nodes

   WW. 02~03.2012
*/
void CFEMesh::setSubdomainNodes(int* header, const MeshNodes* s_nodes)
{
	int k;

	NodesNumber_Quadratic = header[0];
	NodesNumber_Linear = header[1];

	loc_NodesNumber_Linear = header[4];
	loc_NodesNumber_Quadratic = header[5];
	glb_NodesNumber_Linear = header[6];
	glb_NodesNumber_Quadratic = header[7];

	nod_vector.resize(header[0]);
	_global_local_nodeids.reserve(header[0]);
	for (k = 0; k < header[0]; k++)
	{
		const MeshNodes* anode = &s_nodes[k];
		CNode* new_node = new CNode(k, anode->x, anode->y, anode->z);
		new_node->SetEquationIndex(anode->index);
		nod_vector[k] = new_node;
		_global_local_nodeids.push_back(
		    std::make_pair((std::size_t)anode->index, (std::size_t)k));
	}
}

/*!
   Fill data for subdomain mesh

   @param header    : mesh header
   @param elem_info : element information
   @param inside :    indicator for inside domain
   WW. 02~03.2012
*/
void CFEMesh::setSubdomainElements(int* header, const int* elem_info,
                                   const bool inside)
{
	int i;
	int k;
	int counter;
	int mat_idx;
	int e_type;
	int nnodes;

	int ne = 0;
	if (inside)
		ne = header[2];
	else
		ne = header[3];

	// Element
	for (i = 0; i < ne; i++)
	{
		CElem* new_elem = new CElem(ele_vector.size());
		ele_vector.push_back(new_elem);

		counter = elem_info[i];

		mat_idx = elem_info[counter];
		counter++;
		e_type = elem_info[counter];
		counter++;
		nnodes = elem_info[counter];
		counter++;

		new_elem->nnodesHQ = nnodes;
		new_elem->nodes_index.resize(new_elem->nnodesHQ);
		for (k = 0; k < new_elem->nnodesHQ; k++)
		{
			new_elem->nodes_index[k] = elem_info[counter];
			counter++;
		}

		if (!inside)
		{
			const int nn_gl = elem_info[counter];
			counter++;
			const int nn_g = elem_info[counter];
			counter++;

			new_elem->g_index = new int[nn_g + 2];
			int* ele_gnidx = new_elem->g_index;
			ele_gnidx[0] = nn_gl;
			ele_gnidx[1] = nn_g;
			for (k = 2; k < nn_g + 2; k++)
			{
				ele_gnidx[k] = elem_info[counter];
				counter++;
			}
		}

		new_elem->patch_index = static_cast<size_t>(mat_idx);

		// element type
		switch (e_type)
		{
			case 1:
				new_elem->geo_type = MshElemType::LINE;
				new_elem->nnodes = 2;
				new_elem->ele_dim = 1;
				new_elem->nfaces = 2;
				new_elem->nedges = 0;
				break;
			case 2:
				new_elem->geo_type = MshElemType::QUAD;
				new_elem->nnodes = 4;
				new_elem->ele_dim = 2;
				new_elem->nfaces = 4;
				new_elem->nedges = 4;
				break;
			case 3:
				new_elem->geo_type = MshElemType::HEXAHEDRON;
				new_elem->nnodes = 8;
				new_elem->ele_dim = 3;
				new_elem->nfaces = 6;
				new_elem->nedges = 12;
				break;
			case 4:
				new_elem->geo_type = MshElemType::TRIANGLE;
				new_elem->nnodes = 3;
				new_elem->ele_dim = 2;
				new_elem->nfaces = 3;
				new_elem->nedges = 3;
				break;
			case 5:
				new_elem->geo_type = MshElemType::TETRAHEDRON;
				new_elem->nnodes = 4;
				new_elem->ele_dim = 3;
				new_elem->nfaces = 4;
				new_elem->nedges = 6;
				break;
			case 6:
				new_elem->geo_type = MshElemType::PRISM;
				new_elem->nnodes = 6;
				new_elem->ele_dim = 3;
				new_elem->nfaces = 5;
				new_elem->nedges = 9;
				break;
			case 7:  //:PYRAMID:
				new_elem->geo_type = MshElemType::PYRAMID;
				new_elem->nnodes = 5;
				new_elem->ele_dim = 3;
				new_elem->nfaces = 5;
				new_elem->nedges = 8;
				break;
			default:
				break;
		}

		new_elem->InitializeMembers();

		// new_elem->WriteIndex();
	}
}

/*!
    Configure the high order elements for parallel computing

    03.2012. WW
*/
void CFEMesh::ConfigHighOrderElements()
{
	int k;
	size_t i, kk;

	CEdge* edge = NULL;
	CElem* elem = NULL;

	int egde_line[] = {2};
	int egde_tri[] = {3, 4, 5};
	int egde_quad[] = {4, 5, 6};
	int egde_tet[] = {4, 5, 6, 7, 8, 9};
	int egde_hex[] = {8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19};
	int egde_pri[] = {6, 7, 8, 9, 10, 11, 12, 13, 14};

	int* middle_node = NULL;

	for (i = 0; i < edge_vector.size(); i++)
	{
		edge_vector[i]->SetMark(false);
	}
	for (i = 0; i < ele_vector.size(); i++)
	{
		elem = ele_vector[i];

		// Nodes
		elem->nodes.resize(elem->nnodesHQ);
		for (k = 0; k < elem->nnodesHQ; k++)
			elem->nodes[k] = nod_vector[elem->nodes_index[k]];

		switch (elem->geo_type)
		{
			case MshElemType::LINE:
				middle_node = egde_line;
				break;
			case MshElemType::QUAD:
				middle_node = egde_quad;
				break;
			case MshElemType::HEXAHEDRON:
				middle_node = egde_hex;
				break;
			case MshElemType::TRIANGLE:
				middle_node = egde_tri;
				break;
			case MshElemType::TETRAHEDRON:
				middle_node = egde_tet;
				break;
			case MshElemType::PRISM:
				middle_node = egde_pri;
				break;
			case MshElemType::PYRAMID:
				break;
			default:
				break;
		}

		// Edges
		for (kk = 0; kk < elem->nedges; kk++)
		{
			edge = elem->edges[kk];

			if (edge->GetMark()) continue;

			edge->SetNode(kk, elem->nodes[middle_node[kk]]);

			edge->SetMark(true);
		}
	}
}
//
int CFEMesh::calMaximumConnectedNodes()
{
	size_t max_connected_nodes = 0;
	for (auto node : nod_vector)
	{
		max_connected_nodes = std::max(max_connected_nodes, node->getNumConnectedNodes());
	}
	ScreenMessage2d("-> max. connected nodes = %d\n", max_connected_nodes);

	int local_max = max_connected_nodes;
	int global_max = 0;
	MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	return global_max;
}

}  // end namespace

void BuildNodeStruc(MeshNodes* anode, MPI_Datatype* MPI_Node_ptr)
{
	MPI_Datatype my_comp_type[4];
	int nblocklen[4];
	MPI_Aint disp[4], base;
	int j;

	my_comp_type[0] = MPI_INT;
	my_comp_type[1] = MPI_DOUBLE;
	my_comp_type[2] = MPI_DOUBLE;
	my_comp_type[3] = MPI_DOUBLE;
	nblocklen[0] = 1;
	nblocklen[1] = 1;
	nblocklen[2] = 1;
	nblocklen[3] = 1;

// Compute displacement of struct MeshNodes

/*
 disp[0] = 0;
 MPI_Get_address(&(anode.nnodes),&base);
 MPI_Get_address(&(anode.index), disp+1);
 disp[1] -= base;
 MPI_Get_address(&(anode.coord),disp+2);
 disp[2] -= base;

for (j=0; j <3; j++)
{

    cout<<"j=" <<j<<" "<<disp[j]<<endl;
}
*/

#ifndef __MPIUNI_H
	MPI_Get_address(anode, disp);
	MPI_Get_address(&(anode[0].x), disp + 1);
	MPI_Get_address(&(anode[0].y), disp + 2);
	MPI_Get_address(&(anode[0].z), disp + 3);
	base = disp[0];
	for (j = 0; j < 4; j++)
	{
		disp[j] -= base;

		//      cout<<"j=" <<j<<" "<<disp[j]<<endl;
	}

#endif

	// build datatype describing structure
	MPI_Type_create_struct(4, nblocklen, disp, my_comp_type, MPI_Node_ptr);
	MPI_Type_commit(MPI_Node_ptr);
}
