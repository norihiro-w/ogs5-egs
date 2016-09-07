/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "DistributionTools.h"

#include <algorithm>
#include <set>

#include "display.h"
#include "InterpolationAlgorithms/PiecewiseLinearInterpolation.h"
#include "mathlib.h"
#include "msh_mesh.h"
#include "LinearFunctionData.h"

void getNodesOnDistribution(DistributionData& dis_data,
                            MeshLib::CFEMesh& msh,
                            std::vector<long>& nodes_vector)
{
	if (dis_data.geo_type == GEOLIB::POINT)
	{
		// ScreenMessage2("-> looking for nodes on POINT %s\n",
		// dis_data.geo_name.c_str());
		long node_id = msh.GetNODOnPNT(
		    static_cast<const GEOLIB::Point*>(dis_data.geo_obj));
		if (node_id >= 0)
		{
			nodes_vector.push_back(node_id);
			ScreenMessage2("-> node ID %d is found for POINT %s\n",
			               node_id,
			               dis_data.geo_name.c_str());
		}
	}
	else if (dis_data.geo_type == GEOLIB::POLYLINE)
	{
		ScreenMessage2d("-> looking for nodes on POLYLINE %s\n",
		                dis_data.geo_name.c_str());
		CGLPolyline* m_polyline = GEOGetPLYByName(dis_data.geo_name);
		if (m_polyline)
		{
			GEOLIB::Polyline const* ply(
			    static_cast<const GEOLIB::Polyline*>(dis_data.geo_obj));
			double msh_min_edge_length = msh.getMinEdgeLength();
			msh.setMinEdgeLength(m_polyline->epsilon);
			std::vector<size_t> my_nodes_vector;
			msh.GetNODOnPLY(ply, my_nodes_vector);
			msh.setMinEdgeLength(msh_min_edge_length);
			nodes_vector.resize(my_nodes_vector.size());
			for (size_t k(0); k < my_nodes_vector.size(); k++)
				nodes_vector[k] = my_nodes_vector[k];
			// for some benchmarks we need the vector entries sorted by index
			std::sort(nodes_vector.begin(), nodes_vector.end());
		}
	}
	else if (dis_data.geo_type == GEOLIB::SURFACE)
	{
		ScreenMessage2d("-> looking for nodes on SURFACE %s\n",
		                dis_data.geo_name.c_str());
		GEOLIB::Surface const* sfc(
		    static_cast<const GEOLIB::Surface*>(dis_data.geo_obj));

		Surface* m_surface = GEOGetSFCByName(dis_data.geo_name);
		if (m_surface)
		{
			nodes_vector.clear();

			//					m_msh->GetNODOnSFC(m_surface, nodes_vector);
			//#ifndef NDEBUG
			//					GEOLIB::GEOObjects const& geo_obj(*
			// m_msh->getGEOObjects());
			//					std::string const& geo_project_name (*
			// m_msh->getProjectName());
			//					std::string sfc_name;
			//					geo_obj.getSurfaceVecObj(geo_project_name)->getNameOfElement(sfc,
			// sfc_name);
			//					std::string debug_fname("MeshNodesOld-BC-" + sfc_name
			//+
			//".gli");
			//					std::ofstream debug_out (debug_fname.c_str());
			//					debug_out << "#POINTS" << "\n";
			//					for (size_t k(0); k<nodes_vector.size(); k++) {
			//						debug_out << k << " " <<
			//							GEOLIB::Point((m_msh->getNodeVector())[nodes_vector[k]]->getData())
			//<<
			//							" $NAME " << nodes_vector[k] << "\n";
			//					}
			//					debug_out << "#STOP" << "\n";
			//					debug_out.close();
			//#endif
			std::vector<size_t> msh_nod_vec;
			msh.GetNODOnSFC(sfc, msh_nod_vec);
			ScreenMessage2("-> %d nodes found for SURFACE %s\n", msh_nod_vec.size(), dis_data.geo_name.c_str());
			//#ifndef NDEBUG
			//					debug_fname = "MeshNodesNew-BC-" + sfc_name +
			//".gli";
			//					debug_out.open (debug_fname.c_str());
			//					debug_out << "#POINTS" << "\n";
			//					for (size_t k(0); k<msh_nod_vec.size(); k++) {
			//						debug_out << k << " " <<
			//							GEOLIB::Point((m_msh->getNodeVector())[msh_nod_vec[k]]->getData())
			//<<
			//							" $NAME " << msh_nod_vec[k] << "\n";
			//					}
			//					debug_out << "#STOP" << "\n";
			//					debug_out.close();
			//#endif
			//					nodes_vector.clear();
			for (size_t k(0); k < msh_nod_vec.size(); k++)
			{
				//						std::cout << "\t" << k << "\t" <<
				// nodes_vector_old[k] << "\t" << msh_nod_vec[k] << "\n";
				nodes_vector.push_back(msh_nod_vec[k]);
			}
		}
	}
	else if (dis_data.geo_type == GEOLIB::GEODOMAIN)
	{
		nodes_vector.resize(msh.nod_vector.size());
		for (size_t i = 0; i < msh.nod_vector.size(); i++)
			nodes_vector[i] = i;
	}
	else if (!dis_data.mesh_type_name.empty())
	{
		if (dis_data.mesh_type_name == "NODE")
			nodes_vector.push_back(dis_data.mesh_node_id);
	}
	else if (dis_data.mat_id >= 0)
	{
		std::cout << "-> looking for nodes with material group "
		          << dis_data.mat_id << "\n";
		std::set<long> set_node_ids;
		for (size_t ii = 0; ii < msh.ele_vector.size(); ii++)
		{
			MeshLib::CElem* elem = msh.ele_vector[ii];
			if (elem->GetPatchIndex() != static_cast<size_t>(dis_data.mat_id))
				continue;
			size_t nn = elem->GetNodesNumber(msh.getOrder());
			for (size_t i = 0; i < nn; i++)
				set_node_ids.insert(elem->GetNodeIndex(i));
		}
		nodes_vector.assign(set_node_ids.begin(), set_node_ids.end());
	}
	if (!nodes_vector.empty() && dis_data.geo_type != GEOLIB::POINT)
		ScreenMessage2d("-> %d nodes are found on the given distribution\n",
		                nodes_vector.size());
}

void setDistributionConstant(MeshLib::CFEMesh& /*msh*/,
                             std::vector<long>& vec_node_ids,
                             std::vector<double>& vec_node_values,
                             double val)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++)
		vec_node_values[i] = val;
}

void setDistributionConstantGeo(MeshLib::CFEMesh& /*msh*/,
                                std::vector<long>& vec_node_ids,
                                std::vector<double>& vec_node_values,
                                double val)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++)
		vec_node_values[i] = val / (double)nodes_vector_length;
}

void setDistributionGradient(MeshLib::CFEMesh& msh,
                             std::vector<long>& vec_node_ids,
                             std::vector<double>& vec_node_values,
                             double ref_depth,
                             double ref_value,
                             double gradient)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++)
	{
		long node_id = vec_node_ids[i];
		vec_node_values[i] =
		    gradient * (ref_depth - msh.nod_vector[node_id]->getData()[2]) +
		    ref_value;
	}
}

void setDistributionFunction(MeshLib::CFEMesh& msh,
                             std::vector<long>& vec_node_ids,
                             std::vector<double>& vec_node_values,
                             LinearFunctionData& linear_f)
{
	const size_t nodes_vector_length = vec_node_ids.size();
	for (size_t i = 0; i < nodes_vector_length; i++)
	{
		long node_id = vec_node_ids[i];
		double const* const coords(msh.nod_vector[node_id]->getData());
		vec_node_values[i] = linear_f.getValue(coords[0], coords[1], coords[2]);
	}
}

// Interpolation of polygon values to nodes_on_sfc
void SurfaceInterpolation(Surface* m_surface,
                          MeshLib::CFEMesh* m_msh,
                          std::vector<long>& nodes_on_sfc,
                          std::vector<double>& node_value_vector)
{
	//----------------------------------------------------------------------
	const double Tol = m_msh->getMinEdgeLength();
	double gC[3], p1[3], p2[3], vn[3], unit[3], NTri[3];
	//
	std::vector<CGLPolyline*>::iterator p =
	    m_surface->polyline_of_surface_vector.begin();

	for (size_t j = 0; j < nodes_on_sfc.size(); j++)
	{
		double const* const pn(m_msh->nod_vector[nodes_on_sfc[j]]->getData());
		node_value_vector[j] = 0.0;
		bool Passed = false;
		// nodes close to first polyline
		p = m_surface->polyline_of_surface_vector.begin();
		while (p != m_surface->polyline_of_surface_vector.end())
		{
			CGLPolyline* m_polyline = *p;
			// Gravity center of this polygon
			for (int i = 0; i < 3; i++)
				gC[i] = 0.0;
			vn[2] = 0.0;
			size_t nPointsPly = m_polyline->point_vector.size();
			for (size_t i = 0; i < nPointsPly; i++)
			{
				gC[0] += m_polyline->point_vector[i]->x;
				gC[1] += m_polyline->point_vector[i]->y;
				gC[2] += m_polyline->point_vector[i]->z;
				vn[2] += m_polyline->point_vector[i]->getPropert();
			}
			for (size_t i = 0; i < 3; i++)
				gC[i] /= (double)nPointsPly;
			// BC value at center is an average of all point values of polygon
			vn[2] /= (double)nPointsPly;
			// Area of this polygon by the gravity center
			for (size_t i = 0; i < nPointsPly; i++)
			{
				p1[0] = m_polyline->point_vector[i]->x;
				p1[1] = m_polyline->point_vector[i]->y;
				p1[2] = m_polyline->point_vector[i]->z;
				size_t k = i + 1;
				if (i == nPointsPly - 1) k = 0;
				p2[0] = m_polyline->point_vector[k]->x;
				p2[1] = m_polyline->point_vector[k]->y;
				p2[2] = m_polyline->point_vector[k]->z;
				vn[0] = m_polyline->point_vector[i]->getPropert();
				vn[1] = m_polyline->point_vector[k]->getPropert();

				double Area1 = fabs(ComputeDetTri(p1, gC, p2));

				// Check if pn is in the triangle by points (p1, gC, p2)
				double Area2 = fabs(ComputeDetTri(p2, gC, pn));
				unit[0] = fabs(ComputeDetTri(gC, p1, pn));
				unit[1] = fabs(ComputeDetTri(p1, p2, pn));
				Area2 += unit[0] + unit[1];
				if (fabs(Area1 - Area2) < Tol)
				{
					// Interpolation within a triangle (p1,p2,gC)
					// Shape function
					for (size_t l = 0; l < 2; l++)
						unit[l] /= Area1;
					ShapeFunctionTri(NTri, unit);
					for (size_t l = 0; l < 3; l++)
						node_value_vector[j] += vn[l] * NTri[l];
					Passed = true;
					break;
				}
			}
			//
			++p;
			if (Passed) break;
		}  // while
	}      // j
}

void setDistributionLinearPolyline(MeshLib::CFEMesh& msh,
                                   std::vector<long>& /*vec_node_ids*/,
                                   std::vector<double>& vec_node_values,
                                   GEOLIB::Polyline const* ply,
                                   std::vector<double>& DistribedBC,
                                   std::vector<int>& PointsHaveDistribedBC)
{
	std::vector<double> nodes_as_interpol_points;
	msh.getPointsForInterpolationAlongPolyline(ply, nodes_as_interpol_points);
	double msh_min_edge_length = msh.getMinEdgeLength();
	msh.setMinEdgeLength(msh_min_edge_length);

	std::vector<double> interpolation_points;
	std::vector<double> interpolation_values;
	for (size_t i(0); i < DistribedBC.size(); i++)
	{
		for (size_t j = 0; j < ply->getNumberOfPoints(); j++)
		{
			if (PointsHaveDistribedBC[i] == (int)ply->getPointID(j))
			{
				if (std::abs(DistribedBC[i]) <
				    std::numeric_limits<double>::epsilon())
					DistribedBC[i] = 1.0e-20;
				interpolation_points.push_back(ply->getLength(j));
				interpolation_values.push_back(DistribedBC[i]);
				break;
			}
		}
	}
	MathLib::PiecewiseLinearInterpolation(interpolation_points,
	                                      interpolation_values,
	                                      nodes_as_interpol_points,
	                                      vec_node_values);
}

void setDistributionLinearSurface(MeshLib::CFEMesh& msh,
                                  std::vector<long>& vec_node_ids,
                                  std::vector<double>& vec_node_values,
                                  Surface* m_surface,
                                  std::vector<double>& DistribedBC,
                                  std::vector<int>& PointsHaveDistribedBC)
{
	std::vector<CGLPolyline*>::iterator p =
	    m_surface->polyline_of_surface_vector.begin();
	p = m_surface->polyline_of_surface_vector.begin();
	while (p != m_surface->polyline_of_surface_vector.end())
	{
		CGLPolyline* m_polyline = *p;
		for (size_t i(0); i < DistribedBC.size(); i++)
		{
			for (size_t j = 0; j < m_polyline->point_vector.size(); j++)
				if (PointsHaveDistribedBC[i] == m_polyline->point_vector[j]->id)
				{
					if (std::abs(DistribedBC[i]) <
					    std::numeric_limits<double>::epsilon())
						DistribedBC[i] = 1.0e-20;
					m_polyline->point_vector[j]->setPropert(DistribedBC[i]);
					break;
				}
		}
		++p;
	}
	SurfaceInterpolation(m_surface, &msh, vec_node_ids, vec_node_values);
}

void setDistribution(DistributionData& dis_data,
                     MeshLib::CFEMesh& msh,
                     std::vector<long>& vec_node_ids,
                     std::vector<double>& vec_node_values)
{
	if (dis_data.dis_type == FiniteElement::CONSTANT ||
	    dis_data.dis_type == FiniteElement::CONSTANT_NEUMANN)
	{
		setDistributionConstant(
		    msh, vec_node_ids, vec_node_values, dis_data.dis_parameters[0]);
	}
	else if (dis_data.dis_type == FiniteElement::CONSTANT_GEO)
	{
		setDistributionConstantGeo(
		    msh, vec_node_ids, vec_node_values, dis_data.dis_parameters[0]);
	}
	else if (dis_data.dis_type == FiniteElement::LINEAR ||
	         dis_data.dis_type == FiniteElement::LINEAR_NEUMANN)
	{
		if (dis_data.geo_type == GEOLIB::POLYLINE)
		{
			GEOLIB::Polyline const* ply(
			    static_cast<const GEOLIB::Polyline*>(dis_data.geo_obj));
			setDistributionLinearPolyline(msh,
			                              vec_node_ids,
			                              vec_node_values,
			                              ply,
			                              dis_data._DistribedBC,
			                              dis_data._PointsHaveDistribedBC);
		}
		else if (dis_data.geo_type == GEOLIB::SURFACE)
		{
			Surface* m_surface = GEOGetSFCByName(dis_data.geo_name);
			setDistributionLinearSurface(msh,
			                             vec_node_ids,
			                             vec_node_values,
			                             m_surface,
			                             dis_data._DistribedBC,
			                             dis_data._PointsHaveDistribedBC);
		}
	}
	else if (dis_data.dis_type == FiniteElement::GRADIENT)
	{
		setDistributionGradient(msh,
		                        vec_node_ids,
		                        vec_node_values,
		                        dis_data.dis_parameters[0],
		                        dis_data.dis_parameters[1],
		                        dis_data.dis_parameters[2]);
	}
	else if (dis_data.dis_type == FiniteElement::FUNCTION)
	{
		setDistributionFunction(
		    msh, vec_node_ids, vec_node_values, *dis_data.linear_f);
	}
}
