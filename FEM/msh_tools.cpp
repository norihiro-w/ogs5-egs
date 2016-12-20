/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "msh_tools.h"

#include <string>
#include <vector>

#include "display.h"
#include "memory.h"

#include "mathlib.h"

#include "geo_lib.h"

#include "msh_mesh.h"
#include "msh_lib.h"

#include "rf_mmp_new.h"
#include "rf_pcs.h"

using namespace Math_Group;

namespace  MeshLib
{


/**************************************************************************/
/* ROCKFLOW - Function: MSHSelectFreeSurfaceNodes
 */
/* Task:
   Selection of free surface nodes, i.e. Setting free surface node flag = 1
   for the uppermost row and free surface node flag = 2 for the lowermost
   row. (Lowermost row is not moving)
 */
/* Parameter: (I: Input; R: Return; X: Both)   - void -
 */
/* Return:
   - void -
 */
/* Programming:
   03/2003     MB   First Version
   09/2004     MB   PCS
   08/2005     MB msh
 */
/**************************************************************************/
void MSHSelectFreeSurfaceNodes(MeshLib::CFEMesh* m_msh)
{
#ifndef OGS_ONLY_TH
	// Number of nodes per node layer
	size_t NumberOfNodesPerLayer =
		m_msh->nod_vector.size() / (m_msh->getNumberOfMeshLayers() + 1);
	size_t no_unconfined_layer = 0;
	// create array with nodes in vertical column
	size_t* strang(new size_t[m_msh->getNumberOfMeshLayers()]);
	for (size_t i = 0; i < NumberOfNodesPerLayer; i++)
	{
		if (m_msh->nod_vector[i]->free_surface == 4)
		{
			size_t nextnode = i;
			no_unconfined_layer = 0;
			for (size_t j = 0; j < m_msh->getNumberOfMeshLayers(); j++)
			{
				//				strang = (long*)
				// Realloc(strang,(j+1)*sizeof(long));
				strang[j] = nextnode;
				size_t startnode = nextnode;
				nextnode = MSHGetNextNode(startnode, m_msh);
				if (m_msh->nod_vector[nextnode]->free_surface == 4)
				{
					strang[j + 1] = nextnode;
					no_unconfined_layer++;
				}
				else
					continue;
			}
		}  // endif free_surface==4

		// mark start of vertical column with 1 - end of column with 2
		// this is than used in MSHMoveNODUcFlow
		m_msh->nod_vector[strang[0]]->free_surface = 1;
		m_msh->nod_vector[strang[no_unconfined_layer]]->free_surface = 2;
	} /*endfor*/
	delete[] strang;
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: Searches mobile nodes and sets node->free_surface = 4
   Programing:
   09/2004 OK / MB Implementation
   05/2005 OK Bugfix
   07/2005 MB MMP keyword
   08/2005 MB m_pcs
**************************************************************************/
void MSHDefineMobile(CRFProcess* m_pcs)
{
	long* mobile_nodes = NULL;
	long no_mobile_nodes = -1;
	//----------------------------------------------------------------------
	// Define mobile MSH nodes
	//----------------------------------------------------------------------
	// MMP Groups
	if (mmp_vector.size() == 0) return;
	////Schleife �ber alle Gruppen
	for (size_t i = 0; i < mmp_vector.size(); i++)
	{
		CMediumProperties const* const m_mat_mp(mmp_vector[i]);

		// WW    int test = m_pcs->m_msh->GetMaxElementDim();
		// m_pcs->m_msh->cross_section

		// if (m_mat_mp->unconfined_flow_group ==1 &&
		// m_pcs->m_msh->GetMaxElementDim() == 3){
		if ((m_mat_mp->unconfined_flow_group == 1 &&
		     m_pcs->m_msh->GetMaxElementDim() == 3) ||
		    m_pcs->m_msh->hasCrossSection())
		{
			// if (m_mat_mp->unconfined_flow_group ==1){
			// if (m_pcs->m_msh->cross_section){
			//....................................................................
			// DOMAIN
			//			if(m_mat_mp->geo_type_name.find("DOMAIN") !=
			// std::string::npos)
			if (m_mat_mp->getGeoType() == GEOLIB::GEODOMAIN)
			{
				// CGLDomain *m_domain = NULL;
				// m_domain = m_domain->Get(m_mat_mp->geo_name);
				// mobile_nodes = m_domain->GetPointsIn(&no_mobile_nodes);
				// ToDo einlesen von domains ????
				for (size_t nodes = 0; nodes < m_pcs->m_msh->nod_vector.size();
				     nodes++)
				{
					mobile_nodes = (long*)Realloc(mobile_nodes,
					                              sizeof(long) * (nodes + 1));
					mobile_nodes[nodes] = nodes;
				}
				i = m_pcs->m_msh->nod_vector.size();
				no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
			}

			// SURFACE
			//			if(m_mat_mp->geo_type_name.find("SURFACE") !=
			// std::string::npos)
			if (m_mat_mp->getGeoType() == GEOLIB::SURFACE)
			{
				Surface* m_surface = NULL;
				// CC
				m_surface = GEOGetSFCByName(m_mat_mp->geo_name);
				// CC
				mobile_nodes = GetPointsIn(m_surface, &no_mobile_nodes);
			}

			// VOLUME
			//			if(m_mat_mp->geo_type_name.find("VOLUME") !=
			// std::string::npos)
			if (m_mat_mp->getGeoType() == GEOLIB::VOLUME)
			{
				// WW CGLVolume *m_volume = NULL;
				// CC 10/05
				// WW  m_volume = GEOGetVOL(m_mat_mp->geo_name);
				// ToDo TK
				// OK411 mobile_nodes
				// =GetPointsInVolume(m_volume,&no_mobile_nodes);//CC 10/05
			}
		}  // end if unconfined flow group
	}      // end for mmp vector

#ifndef OGS_ONLY_TH
	// Set mobile MSH nodes flag
	for (long i = 0; i < no_mobile_nodes; i++)
		m_pcs->m_msh->nod_vector[i]->free_surface = 4;
#endif

	if (no_mobile_nodes > 0)
	{
		m_pcs->mobile_nodes_flag = 1;
		MSHSelectFreeSurfaceNodes(m_pcs->m_msh);
	}
}

/**************************************************************************/
/* ROCKFLOW - Function: MSHMoveNODUcFlow
 */
/* Task:
   Moves free surface nodes according to the pressure distribution
 */
/* Parameter: (I: Input; R: Return; X: Both)   - void -
 */
/* Return:
   - void -
 */
/* Programming:
   09/2002     MB       First Version
   05/2003     MB       verallgemeinert f�r Prismen und Vierecke
   09/2004     MB       Methode vereinfacht
   09/2004     MB       PCS
   08/2005      MB       m_msh */
/**************************************************************************/
void MSHMoveNODUcFlow(CRFProcess* m_pcs)
{
#ifndef OGS_ONLY_TH
	long nextnode = -1;
	long startnode;
	int anz_zeilen = 0;
	int i;
	double spanne_ges;
	double spanne_rel;
	long* strang = NULL;
	double head = 0.0;
	int xxflag;
	int nidy;
	// Number of nodes per node layer
	const size_t NumberOfNodesPerLayer(
	    m_pcs->m_msh->nod_vector.size() /
	    (m_pcs->m_msh->getNumberOfMeshLayers() + 1));
	double MinThickness = 1e-1;  // OKMB
	double z_bottom;             // OKMB

	for (size_t node = 0; node < NumberOfNodesPerLayer; node++)

		if (m_pcs->m_msh->nod_vector[node]->free_surface == 1)
		{
			/* Z�hlen der Zeilen (-> anz_zeilen) */
			anz_zeilen = 0;
			xxflag = 0;
			nextnode = node;
			do
			{
				startnode = nextnode;
				nextnode = MSHGetNextNode(startnode, m_pcs->m_msh);

				/* Test2: Geh�rt der n�chste Knoten zu unterer Reihe ==> Abbruch
				 */
				if (m_pcs->m_msh->nod_vector[nextnode]->free_surface == 2)
					xxflag = 1;
				anz_zeilen++; /* Anzahl der beweglichen Zeilen (ohne die feste
				                 untere Zeile) */
			} while (xxflag != 1);
			/** Ende Z�hlen der Zeilen zwischen den oberen free surface node
			 * etc... und den Unteren **/

			/* Die Knoten unterhalb eines Free Surface Knotens bilden einen
			 * Strang */
			/* Die Knoten eines Stranges werden zwischengespeichert */
			strang = MSHGetNodesInColumn(node, anz_zeilen, m_pcs->m_msh);

			/* Die Knoten eines Stranges werden entsprechend der neuen
			 * Druckverteilung  verformt */
			/* Standrohrspiegelh�he bestimmen */
			nidy = m_pcs->GetNodeValueIndex("HEAD") + 1;
			if (GetRFProcessDensityFlow()) /* mit Dichteunterschiede */
			{
				// OK_MOD     head = MODCalcHeadInColumn_MB(strang, anz_zeilen);
			}
			else /* ohne Dichteunterschiede */
				head = m_pcs->GetNodeValue(strang[0], nidy);

			/* Set minimum thickness */
			z_bottom =
			    m_pcs->m_msh->nod_vector[strang[anz_zeilen]]->getData()[2];
			if (head - z_bottom < MinThickness) head = z_bottom + MinThickness;

			/* Berechnung der Differenz */
			spanne_ges = head - z_bottom;
			spanne_rel = spanne_ges / anz_zeilen;
			m_pcs->m_msh->nod_vector[strang[0]]->SetZ(head);

			if (spanne_ges != 0)
				/* Setzen der neuen Z-Werte entlang eines Stranges */
				for (i = 1; i < anz_zeilen;
				     i++) /* Schleife �ber Anzahl der Zeilen */
					m_pcs->m_msh->nod_vector[strang[i]]->SetZ(head -
					                                          i * spanne_rel);

			strang = (long*)Free(strang);
		} /*endif index ==1 */
          /* end for Schleife �ber alle Knoten */
#endif
}

/**************************************************************************
   FEMLib-Method:
   Task: Searches mobile nodes and sets node->free_surface = 4
   Programing:
   09/2004 OK / MB Implementation
   05/2005 OK Bugfix
   07/2005 MB MMP keyword
   08/2005 MB m_pcs
   01/2006 OK LAYER
**************************************************************************/
void DefineMobileNodes(MeshLib::CFEMesh*msh, CRFProcess* m_pcs)
{
#ifndef OGS_ONLY_TH
	long* mobile_nodes = NULL;
	long no_mobile_nodes = -1;
	long i, j;
	//----------------------------------------------------------------------
	// Define mobile MSH nodes
	//----------------------------------------------------------------------
	//......................................................................
	// DOMAIN
	if (m_pcs->geo_type.find("DOMAIN") != std::string::npos)
	{
		for (i = 0; i < (long)msh->nod_vector.size(); i++)
		{
			mobile_nodes = (long*)Realloc(mobile_nodes, sizeof(long) * (i + 1));
			mobile_nodes[i] = i;
		}
		no_mobile_nodes = (long)m_pcs->m_msh->nod_vector.size();
	}
	//......................................................................
	// LAYER
	if (m_pcs->geo_type.find("LAYER") != std::string::npos)
	{
		std::string m_string;
		long no_nodes_per_layer =
			(long)msh->nod_vector.size() / (msh->getNumberOfMeshLayers() + 1);
		int pos = 0;
		int layer_start = 0, layer_end = 0;
		if (m_pcs->geo_type_name.find("-") != std::string::npos)
		{
			pos = m_pcs->geo_type_name.find("-") != std::string::npos;
			m_string = m_pcs->geo_type_name.substr(0, pos);
			layer_start = strtol(m_string.c_str(), NULL, 0);
			m_string = m_pcs->geo_type_name.substr(pos + 1, std::string::npos);
			layer_end = strtol(m_string.c_str(), NULL, 0);
		}
		else
		{
			layer_start = strtol(m_pcs->geo_type_name.c_str(), NULL, 0);
			layer_end = layer_start;
		}
		int no_layers = layer_end - layer_start + 1;
		no_mobile_nodes = (no_layers + 1) * no_nodes_per_layer;
		mobile_nodes = new long[no_mobile_nodes];
		for (i = 0; i < no_layers + 1; i++)
			for (j = 0; j < no_nodes_per_layer; j++)
				mobile_nodes[i * no_nodes_per_layer + j] =
				    j + (layer_start - 1 + i) * no_nodes_per_layer;
	}
	//......................................................................
	// SURFACE
	if (m_pcs->geo_type.find("SURFACE") != std::string::npos)
	{
		Surface* m_sfc = NULL;
		// CC
		m_sfc = GEOGetSFCByName(m_pcs->geo_type_name);
		if (m_sfc)
			// CC
			mobile_nodes = GetPointsIn(m_sfc, &no_mobile_nodes);
		else
			std::cout << "Warning in CFEMesh::DefineMobileNodes - no GEO data"
			          << "\n";
	}
	//......................................................................
	// VOLUME
	/*OK411
	   if(m_pcs->geo_type.find("VOLUME")!=std::string::npos)
	   {
	    CGLVolume *m_vol = NULL;
	    m_vol = GEOGetVOL(m_pcs->geo_type_name);//CC 10/05
	    if(m_vol)
	      mobile_nodes = GetPointsInVolume(m_vol,&no_mobile_nodes);//CC 10/05
	    else
	      std::cout << "Warning in CFEMesh::DefineMobileNodes - no GEO data" <<
	   endl;
	   }
	 */
	//----------------------------------------------------------------------
	// Set mobile MSH nodes flag
	//----------------------------------------------------------------------
	for (i = 0; i < (long)msh->nod_vector.size(); i++)
		msh->nod_vector[i]->free_surface = -1;
	for (i = 0; i < no_mobile_nodes; i++)
		msh->nod_vector[i]->free_surface = 4;
	// nod_vector[mobile_nodes[i]]->free_surface = 4;
	//----------------------------------------------------------------------
	if (no_mobile_nodes > 0)
	{
		m_pcs->mobile_nodes_flag = 1;
		MSHSelectFreeSurfaceNodes(msh);
	}
	//----------------------------------------------------------------------
	delete[] mobile_nodes;
	mobile_nodes = NULL;
#endif
}

/**************************************************************************
   MSHLib-Method:
   07/2007 OK Implementation
**************************************************************************/
bool MSHTestMATGroups()
{
	int g_max_mmp_groups = MSHSetMaxMMPGroups();
	if (g_max_mmp_groups > (int)mmp_vector.size())
	{
		std::cout << "Error: not enough MMP data";
		return false;  // abort();
	}
	return true;
}

} // MeshLib
