/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "msh_lib.h"

#include <cstdio>
#include <cstring>
#include <sstream>

#include <mpi.h>
#include <petscksp.h>

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
	static const int nheaders = 11;
	int mesh_header[nheaders] = {-1};
	MeshHeader meshHeader;

	MeshNodes* s_nodes = (MeshNodes*)malloc(sizeof(MeshNodes));
	int* elem_info = (int*)malloc(1);

	MPI_Datatype MPI_node;
	int tag[] = {0, 1, 2};
	MPI_Status status;

	ifstream is; // only root opens the file

	const bool isRoot = (myrank == 0);
	if (isRoot)
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
				MPI_Abort(MPI_COMM_WORLD, 1);
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
	ScreenMessage("-> reading the partitioned mesh\n");
	for (int i = 0; i < mysize; i++)
	{
		// read a subdomain for rank i
		//-------------------------------------------------------------------------
		// root reads the header and boradcasts
		if (isRoot)
		{
			for (int j = 0; j < nheaders; j++)
				is >> mesh_header[j];
			is >> ws;

			//if (mesh->_vec_globalNodeID2domID.empty())
			//{
			//	mesh->_vec_globalNodeID2domID.resize(mesh_header[7], -1); // global nnodes
			//}
		}
		MPI_Bcast(mesh_header, nheaders, MPI_INT, 0, MPI_COMM_WORLD);
		meshHeader.set(mesh_header);

		// debug output
		if (i == myrank)
		{
			std::stringstream ss;
			for (int j = 0; j < nheaders; j++)
				ss << mesh_header[j] << " ";
			ScreenMessage2d("-> header: %s\n", ss.str().c_str());
		}

		//-------------------------------------------------------------------------
		// Node
		const bool hasQuadraticNodes = (meshHeader.n_dom_nodes_L != meshHeader.n_dom_nodes_Q);
		s_nodes = (MeshNodes*)realloc(s_nodes, sizeof(MeshNodes) * meshHeader.n_dom_nodes_Q);

		if (i > 0) // no need to broadcast if rank i is root
			BuildNodeStruc(s_nodes, &MPI_node);

		// root reads node data and broadcasts to rank-i
		if (isRoot)
		{
			for (int k = 0; k < meshHeader.n_dom_nodes_Q; k++)
			{
				MeshNodes* anode = &s_nodes[k];
				is >> anode->global_id >> anode->dom_id >> anode->eqs_id;
				if (hasQuadraticNodes)
					is >> anode->eqs_id_Q;
				is >> anode->x >> anode->y >> anode->z >> ws;

				//mesh->_vec_globalNodeID2domID[anode->global_id] = anode->dom_id;
			}
			if (i == 0)
			{
				mesh->setSubdomainNodes(meshHeader, s_nodes);
			}
			else
			{
				MPI_Send(s_nodes, meshHeader.n_dom_nodes_Q, MPI_node, i, tag[0], MPI_COMM_WORLD);
			}
		}

		// rank-i recieve the node data
		if (i > 0)
		{
			if (myrank == i)
			{
				MPI_Recv(s_nodes, meshHeader.n_dom_nodes_Q, MPI_node, 0, tag[0],
				         MPI_COMM_WORLD, &status);
				mesh->setSubdomainNodes(meshHeader, s_nodes);
			}
		}

		//-------------------------------------------------------------------------
		// Element
		const int size_elem_info = meshHeader.n_inner_elements + meshHeader.n_element_integers;
		elem_info = (int*)realloc(elem_info, sizeof(int) * size_elem_info);
		for (int ii = 0; ii < size_elem_info; ii++)
			elem_info[ii] = -1;
		if (isRoot)
		{
			int counter = meshHeader.n_inner_elements;
			for (int j = 0; j < meshHeader.n_inner_elements; j++)
			{
				elem_info[j] = counter;
				is >> elem_info[counter++];  // global ele. ID
				is >> elem_info[counter++];  // mat. idx
				is >> elem_info[counter++];  // type
				is >> elem_info[counter++];    // nnodes
				const int nn_e = elem_info[counter-1];
				for (int k = 0; k < nn_e; k++)
					is >> elem_info[counter++];
			}

			if (i == 0)
				mesh->setSubdomainElements(meshHeader, elem_info, true);
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
				mesh->setSubdomainElements(meshHeader, elem_info, true);
			}
		}

		// if(elem_info)
		//  {
		// free(elem_info);
		// elem_info = NULL;
		//  }

		//-------------------------------------------------------------------------
		// Ghost element
		const int size_elem_g_info = meshHeader.n_ghost_elements + meshHeader.n_ghost_element_integers;
		elem_info = (int*)realloc(elem_info, sizeof(int) * size_elem_g_info);
		if (myrank == 0)
		{
			int counter = meshHeader.n_ghost_elements;
			for (int j = 0; j < meshHeader.n_ghost_elements; j++)
			{
				elem_info[j] = counter;
				is >> elem_info[counter++];  // ele. idx
				is >> elem_info[counter++];  // mat. idx
				is >> elem_info[counter++];  // type
				is >> elem_info[counter];  // nnodes
				const int nn_e = elem_info[counter];
				counter++;
				for (int k = 0; k < nn_e; k++)
				{
					is >> elem_info[counter++];
				}
				is >> elem_info[counter++];
				// const int nn_e_g =  elem_info[counter++];
				// ghost nodes for linear element
				is >> elem_info[counter];
				const int nn_e_g_quad = elem_info[counter++];
				// NW use nn_e_g_quad instead of nn_e_g
				for (int k = 0; k < nn_e_g_quad; k++)
				{
					is >> elem_info[counter++];
				}
			}

			if (i == 0)
				mesh->setSubdomainElements(meshHeader, elem_info, false);
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
				mesh->setSubdomainElements(meshHeader, elem_info, false);
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

	//-------------------------------------------------------------------------------------
	ScreenMessage("-> global: nnodes_l=%d, nnodes_g=%d\n", mesh->getNumNodesGlobal(), mesh->getNumNodesGlobal_Q());
	MPI_Barrier(MPI_COMM_WORLD);
	ScreenMessage2(
	    "-> nelements=%d, nnodes=%d, nnodes_q=%d, nnodes_ll=%d, "
	    "nnodes_lq=%d\n",
	    mesh->getElementVector().size(), mesh->GetNodesNumber(false),
	    mesh->GetNodesNumber(true), mesh->getNumNodesLocal(),
	    mesh->getNumNodesLocal_Q());
	MPI_Barrier(MPI_COMM_WORLD);

	//-------------------------------------------------------------------------------------
#if 0
	if (mesh->hasHigherOrderNodes())
	{
		// assign global equation IDs to local nodes. IDs should be continuous even for quadratic nodes
		ScreenMessage("-> assgin global equation IDs to local nodes\n");
		int global_eqs_id_Q = 0;
		for (int i=0; i<mysize; i++)
		{
			if (i==myrank)
			{
				//ScreenMessage2("-> assgin global equation index for quadratic from id %d\n", global_eqs_id_Q);
				for (CNode* node : mesh->getNodeVector())
				{
					if (mesh->isNodeLocal(node->GetIndex()))
						node->SetEquationIndex(global_eqs_id_Q++, true);
					else
						node->SetEquationIndex(-1, true);
				}
				//ScreenMessage2("-> set global_eqs_id_Q=%d\n", global_eqs_id_Q);
			}
			MPI_Bcast(&global_eqs_id_Q, 1, MPI_INT, i, MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
		}
		if (myrank == mysize-1)
		{
			//ScreenMessage2("-> global_eqs_id_Q=%d, mesh->getNumNodesGlobal_Q=%d\n", global_eqs_id_Q, mesh->getNumNodesGlobal_Q());
			if (global_eqs_id_Q != mesh->getNumNodesGlobal_Q())
				ScreenMessage2("*** error: global_eqs_id_Q (%d) != mesh->getNumNodesGlobal_Q (%d)\n", global_eqs_id_Q, mesh->getNumNodesGlobal_Q());
		}

		// get global equation IDs of ghost nodes from other ranks
		ScreenMessage("-> assgin global equation IDs to ghost nodes\n");
		for (int i=0; i<mysize; i++)
		{
			// send a list of ghost nodes in this rank
			std::vector<int> vec_ghost_node_localIDs;
			std::vector<int> vec_ghost_node_globalIDs;
			if (i==myrank)
			{
				//ScreenMessage2("-> send a list of ghost node IDs\n");
				for (CNode* node : mesh->getNodeVector())
				{
					if (mesh->isNodeLocal(node->GetIndex()))
						continue;
					vec_ghost_node_localIDs.push_back(node->GetIndex());
					vec_ghost_node_globalIDs.push_back(node->GetEquationIndex());
				}
				ScreenMessage2("-> looking for equations ids for %d ghost nodes\n", vec_ghost_node_globalIDs.size());
			}
			int vec_ghost_node_globalIDs_size = vec_ghost_node_globalIDs.size();
			MPI_Bcast(&vec_ghost_node_globalIDs_size, 1, MPI_INT, i, MPI_COMM_WORLD);
			if (i!=myrank)
				vec_ghost_node_globalIDs.resize(vec_ghost_node_globalIDs_size);
			MPI_Bcast(vec_ghost_node_globalIDs.data(), vec_ghost_node_globalIDs.size(), MPI_INT, i, MPI_COMM_WORLD);
			//if (i!=myrank)
			//	ScreenMessage2("-> received %d ghost node IDs from rank %d\n", vec_ghost_node_globalIDs.size(), i);

			// for each ghost node
			std::vector<int> vec_rank_data(mysize);
			for (size_t j=0; j<vec_ghost_node_globalIDs.size(); j++)
			{
				int node_globalID = vec_ghost_node_globalIDs[j];
				//ScreenMessage("-> look for a node with global ID %d\n", node_globalID);
				MPI_Barrier(MPI_COMM_WORLD);
				// who owns this?
				int global_eqs_id_Q = -1;
				if (i != myrank)
				{
					CNode* node = mesh->findNodeByGlobalID(node_globalID);
					if (node) {
						global_eqs_id_Q = node->GetEquationIndex(true);
						//ScreenMessage2("-> found a node (local id=%d, global id=%d, global id(Q)=%d)\n", node->GetIndex(), node_globalID, global_eqs_id_Q);
					}
				}
				MPI_Gather(&global_eqs_id_Q, 1, MPI_INT, vec_rank_data.data(), 1, MPI_INT, i, MPI_COMM_WORLD);
				//
				if (i == myrank)
				{
					//ScreenMessage2("-> gathered a result of node look up\n");
					int reccv_eqsid = -1;
					for (size_t ii=0; ii<vec_rank_data.size(); ii++)
					{
						if (vec_rank_data[ii]!=-1) {
							reccv_eqsid = vec_rank_data[ii];
							break;
						}
					}
					//ScreenMessage2("-> global node id=%d, global equation ID (Q)=%d\n", node_globalID, reccv_eqsid);
					if (reccv_eqsid<0)
						ScreenMessage2("*** error: not found global equation ID (Q) for node %d\n", vec_ghost_node_localIDs[j]);
					mesh->nod_vector[vec_ghost_node_localIDs[j]]->SetEquationIndex(reccv_eqsid, true);
				}
			}
		}

	}
#endif

	//-------------------------------------------------------------------------------------
	mesh->ConstructGrid();
	mesh->FillTransformMatrix();
	// mesh->calMaximumConnectedNodes();
	MPI_Barrier(MPI_COMM_WORLD);
}

namespace MeshLib
{

CNode* CFEMesh::findNodeByGlobalID(long global_node_id) const
{
	for (CNode* node : nod_vector)
	{
		if (node->GetEquationIndex() != global_node_id)
			continue;
		if (!isNodeLocal(node->GetIndex()))
			continue;
		return node;
	}
	return nullptr;
}

/*!
   Fill data for subdomain mesh

   @param header  : mesh header
   @param s_nodes : mesh nodes

   WW. 02~03.2012
*/
void CFEMesh::setSubdomainNodes(MeshHeader const& header, const MeshNodes* s_nodes)
{
	NodesNumber_Quadratic = header.n_dom_nodes_Q;
	NodesNumber_Linear = header.n_dom_nodes_L;

	loc_NodesNumber_Linear = header.n_inner_nodes_L;
	loc_NodesNumber_Quadratic = header.n_inner_nodes_Q;
	glb_NodesNumber_Linear = header.n_global_nodes_L;
	glb_NodesNumber_Quadratic = header.n_global_nodes_Q;
	glb_ElementsNumber = header.n_global_elements;

	nod_vector.resize(NodesNumber_Quadratic);
	_global_local_nodeids.reserve(NodesNumber_Quadratic);
	_vec_node_dom_ids.reserve(NodesNumber_Quadratic);
	for (size_t k = 0; k < NodesNumber_Quadratic; k++)
	{
		const MeshNodes* anode = &s_nodes[k];
		CNode* new_node = new CNode(k, anode->x, anode->y, anode->z);
		new_node->SetGlobalIndex(anode->global_id);
		new_node->SetEquationIndex(anode->eqs_id, false);
		new_node->SetEquationIndex(anode->eqs_id_Q, true);
		nod_vector[k] = new_node;
		_global_local_nodeids.push_back(
		    std::make_pair((std::size_t)anode->global_id, (std::size_t)k));
		_vec_node_dom_ids.push_back(anode->dom_id);
	}
}

/*!
   Fill data for subdomain mesh

   @param header    : mesh header
   @param elem_info : element information
   @param inside :    indicator for inside domain
   WW. 02~03.2012
*/
void CFEMesh::setSubdomainElements(MeshHeader const& header, const int* elem_info,
                                   const bool inside)
{
	int ne = inside ? header.n_inner_elements : header.n_ghost_elements;

	// Element
	std::stringstream ss;
	for (int i = 0; i < ne; i++)
	{
		CElem* new_elem = new CElem(ele_vector.size());
		ele_vector.push_back(new_elem);

		int counter = elem_info[i];

		new_elem->global_index = elem_info[counter++];
		int mat_idx = elem_info[counter++];
		int e_type = elem_info[counter++];
		int nnodes = elem_info[counter++];

		new_elem->nnodesHQ = nnodes;
		new_elem->nodes_index.resize(new_elem->nnodesHQ);
		for (int k = 0; k < new_elem->nnodesHQ; k++)
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
			for (int k = 2; k < nn_g + 2; k++)
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

		new_elem->WriteIndex(ss);
	}
	ScreenMessage2("\n%s\n", ss.str().data());
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

int CFEMesh::calMaximumConnectedLocalNodes(bool quadratic, std::vector<int> &d_nnz)
{
	d_nnz.resize(quadratic ? this->getNumNodesLocal_Q() : this->getNumNodesLocal());
	size_t max_connected_nodes = 0;
	const int node0_eqsId = nod_vector[0]->GetEquationIndex();
	for (size_t i=0; i<this->GetNodesNumber(quadratic); i++)
	{
		CNode* node = nod_vector[i];
		if (!this->isNodeLocal(node->GetIndex()))
			continue;
		size_t cnt_local = 0;
		for (auto node_id : node->getConnectedNodes())
			if (this->isNodeLocal(node_id))
				cnt_local++;
		max_connected_nodes = std::max(max_connected_nodes, cnt_local);
		int eqs_id = node->GetEquationIndex() - node0_eqsId;
		d_nnz[eqs_id] = cnt_local;
	}
	ScreenMessage2d("-> max. connected local nodes = %d\n", max_connected_nodes);
	return max_connected_nodes;
}

int CFEMesh::calMaximumConnectedGhostNodes(bool quadratic, std::vector<int> &o_nnz)
{
	o_nnz.resize(quadratic ? this->getNumNodesLocal_Q() : this->getNumNodesLocal());
	size_t max_connected_nodes = 0;
	const int node0_eqsId = nod_vector[0]->GetEquationIndex();
	for (size_t i=0; i<this->GetNodesNumber(quadratic); i++)
	{
		CNode* node = nod_vector[i];
		if (!this->isNodeLocal(node->GetIndex()))
			continue;
		size_t cnt_ghost = 0;
		for (auto node_id : node->getConnectedNodes())
			if (!this->isNodeLocal(node_id))
				cnt_ghost++;
		max_connected_nodes = std::max(max_connected_nodes, cnt_ghost);
		int eqs_id = node->GetEquationIndex() - node0_eqsId;
		o_nnz[eqs_id] = cnt_ghost;
	}
	ScreenMessage2d("-> max. connected ghost nodes = %d\n", max_connected_nodes);
	return max_connected_nodes;
}

int CFEMesh::getMaxNumNodesOfElement(bool quadratic) const
{
	int max_nnodes = 0;
	if (this->getNumberOfLines() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 3 : 2);
	if (this->getNumberOfTris() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 6 : 3);
	if (this->getNumberOfQuads() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 9 : 4);
	if (this->getNumberOfHexs() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 20 : 8);
	if (this->getNumberOfTets() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 10 : 4);
	if (this->getNumberOfPrisms() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 15 : 6);
	if (this->getNumberOfPyramids() > 0)
		max_nnodes = std::max(max_nnodes, quadratic ? 13 : 5);

	int local_max = max_nnodes;
	int global_max = 0;
	MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	return global_max;
}

int CFEMesh::getMaxNumConnectedElements() const
{
	size_t max_n_conn_eles = 0;
	for (CNode const* node : nod_vector)
	{
		max_n_conn_eles = std::max(max_n_conn_eles, node->getConnectedElementIDs().size());
	}

	int local_max = max_n_conn_eles;
	int global_max = 0;
	MPI_Allreduce(&local_max, &global_max, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

	return global_max;
}

}  // end namespace

void BuildNodeStruc(MeshNodes* anode, MPI_Datatype* MPI_Node_ptr)
{
	static const int ncomp = 7;
	MPI_Datatype my_comp_type[ncomp];
	int nblocklen[ncomp];
	MPI_Aint disp[ncomp], base;

	int counter = 0;
	my_comp_type[counter++] = MPI_INT;
	my_comp_type[counter++] = MPI_INT;
	my_comp_type[counter++] = MPI_INT;
	my_comp_type[counter++] = MPI_INT;
	my_comp_type[counter++] = MPI_DOUBLE;
	my_comp_type[counter++] = MPI_DOUBLE;
	my_comp_type[counter++] = MPI_DOUBLE;
	counter = 0;
	for (int i=0; i<ncomp; i++)
		nblocklen[counter++] = 1;

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
	MPI_Get_address(&anode->dom_id, disp + 1);
	MPI_Get_address(&anode->eqs_id, disp + 2);
	MPI_Get_address(&anode->eqs_id_Q, disp + 3);
	MPI_Get_address(&(anode[0].x), disp + 4);
	MPI_Get_address(&(anode[0].y), disp + 5);
	MPI_Get_address(&(anode[0].z), disp + 6);
	base = disp[0];
	for (int j = 0; j < ncomp; j++)
	{
		disp[j] -= base;

		//      cout<<"j=" <<j<<" "<<disp[j]<<endl;
	}

#endif

	// build datatype describing structure
	MPI_Type_create_struct(ncomp, nblocklen, disp, my_comp_type, MPI_Node_ptr);
	MPI_Type_commit(MPI_Node_ptr);
}
