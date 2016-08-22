/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CSVOutput.h"

#include <algorithm>
#include <iomanip>
#include <iterator>

#include "display.h"

#include "msh_lib.h"
#include "msh_mesh.h"
#include "msh_node.h"

#include "FEMEnums.h"
#include "fem_ele_std.h"
#include "matrix_class.h"
#include "Output.h"
#include "ProcessInfo.h"
#include "rf_mmp_new.h"
#include "rf_pcs.h"

using namespace std;

namespace
{

template <typename T>
void writeLine(std::ostream &os, std::vector<T> const& vec, std::string const& delim)
{
	if (vec.empty()) return;
	auto itr = vec.begin();
	os << *itr;
	++itr;
	for (;itr!=vec.end(); ++itr)
		os << delim << *itr;
	os << "\n";
}

void writeNodeData(COutput* output, std::string const& filename)
{
	if (output->_nod_value_vector.empty())
		return;

	std::ofstream os(filename);
	if (!os) {
		ScreenMessage2("***Error: cannot open %s for writing\n", filename.data());
		return;
	}
	os.setf(ios::scientific, std::ios::floatfield);
	os.precision(12);


	//
	std::vector<std::string> vec_output_val_name;
	vec_output_val_name.push_back("NodeID");
	std::vector<CRFProcess*> vec_pcs;
	std::vector<int> vec_pcs_value_index;
	for (size_t i=0; i<output->_nod_value_vector.size(); i++)
	{
		auto& value_name = output->_nod_value_vector[i];
		auto pcs = PCSGet(value_name, true);
		if (!pcs)
			continue;
		auto value_id = pcs->GetNodeValueIndex(value_name);
		if (value_id<0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_pcs.push_back(pcs);
		vec_pcs_value_index.push_back(value_id);
	}

	auto const n_nodal_values = vec_pcs_value_index.size();

	std::string const delim = ", ";
	// header
	writeLine(os, vec_output_val_name, delim);

	// values
	MeshLib::CFEMesh* msh = output->getMesh();
	std::vector<double> tmp_values(n_nodal_values);
	for (long i=0; i<(long)msh->GetNodesNumber(false); i++)
	{
		MeshLib::CNode* node = msh->getNodeVector()[i];
#ifdef USE_PETSC
		os << node->GetEquationIndex() << delim;
#else
		os << node->GetIndex() << delim;
#endif
		for (size_t j=0; j<vec_pcs_value_index.size(); j++)
		{
			auto pcs = vec_pcs[j];
			tmp_values[j] = pcs->GetNodeValue(i, vec_pcs_value_index[j]);
		}
		writeLine(os, tmp_values, delim);
	}
}

struct ELEMENT_MMP_VALUES
{
	static double getValue(CMediumProperties* mmp, int mmp_id, long i_e,
						   double* gp, double theta)
	{
		double mat_value = .0;
		switch (mmp_id)
		{
			case 0:
				mat_value = mmp->Porosity(i_e, theta);
				break;
			case 1:
				mat_value = mmp->PermeabilityTensor(i_e)[0];
				break;
			case 2:
				mat_value = mmp->StorageFunction(i_e, gp, theta);
				break;
			default:
				cout << "CVTK::WriteElementValues: no MMP values "
						"specified"
					 << "\n";
				break;
		}
		return mat_value;
	}

	static int getMMPIndex(const std::string& mmp_name)
	{
		int mmp_id = -1;
		if (mmp_name.compare("POROSITY") == 0)
		{
			mmp_id = 0;
		}
		else if (mmp_name.compare("PERMEABILITY") == 0)
		{
			mmp_id = 1;
		}
		else if (mmp_name.compare("STORAGE") == 0)
		{
			mmp_id = 2;
		}
		else
		{
			cout << "CVTK::WriteElementValues: no valid MMP values "
					"specified. " << mmp_name << "\n";
		}
		return mmp_id;
	}
};

struct ELEMENT_MFP_VALUES
{
	static double getValue(CFluidProperties* mfp, int mfp_id)
	{
		double mat_value = .0;
		switch (mfp_id)
		{
			case 0:
				mat_value = mfp->Density();
				break;
			case 1:
				mat_value = mfp->Viscosity();
				break;
			default:
				cout << "CVTK::WriteElementValues: no MFP values "
						"specified"
					 << "\n";
				break;
		}
		return mat_value;
	}

	static int getMFPIndex(const std::string& mfp_name)
	{
		int mfp_id = -1;
		if (mfp_name.compare("DENSITY") == 0)
		{
			mfp_id = 0;
		}
		else if (mfp_name.compare("VISCOSITY") == 0)
		{
			mfp_id = 1;
		}
		else
		{
			cout << "CVTK::WriteElementValues: no valid MFP values "
					"specified. " << mfp_name << "\n";
		}
		return mfp_id;
	}
};

void writeElementData(COutput* output, std::string const& filename)
{
	if (output->getElementValueVector().empty())
		return;

	std::ofstream os(filename);
	if (!os) {
		ScreenMessage2("***Error: cannot open %s for writing\n", filename.data());
		return;
	}
	os.setf(ios::scientific, std::ios::floatfield);
	os.precision(12);

	//
	std::vector<std::string> vec_output_val_name;
	vec_output_val_name.push_back("ElementID");
	std::vector<CRFProcess*> vec_pcs;
	std::vector<int> vec_pcs_value_index;
	for (auto const& value_name : output->getElementValueVector())
	{
		auto pcs = output->GetPCS_ELE(value_name);
		if (!pcs)
			continue;
		auto value_id = pcs->GetElementValueIndex(value_name);
		if (value_id<0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_pcs.push_back(pcs);
		vec_pcs_value_index.push_back(value_id);
	}


	std::vector<int> vec_mmp_id;
	for (auto const& value_name : output->mmp_value_vector)
	{
		int mmp_id = ELEMENT_MMP_VALUES::getMMPIndex(value_name);
		if (mmp_id < 0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_mmp_id.push_back(mmp_id);
	}

	std::vector<int> vec_mfp_id;
	for (auto const& value_name : output->mfp_value_vector)
	{
		int mfp_id = ELEMENT_MFP_VALUES::getMFPIndex(value_name);
		if (mfp_id < 0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_mfp_id.push_back(mfp_id);
	}

	auto const n_output_values = vec_output_val_name.size() - 1; // exclude elemnt ID

	std::string const delim = ", ";
	// header
	writeLine(os, vec_output_val_name, delim);

	// values
	MeshLib::CFEMesh* msh = output->getMesh();
	std::vector<double> tmp_values(n_output_values);
	double gp[3] = {.0, .0, .0};
	int gp_r=0, gp_s=0, gp_t=0;
	double const theta = 1.0;
	auto mat_pcs = pcs_vector[0];
	for (long i=0; i<(long)msh->getElementVector().size(); i++)
	{
		auto ele = msh->ele_vector[i];
		ele->SetOrder(false);
		// pcs element values
		for (size_t j=0; j<vec_pcs_value_index.size(); j++)
		{
			auto pcs = vec_pcs[j];
			tmp_values[j] = pcs->GetElementValue(i, vec_pcs_value_index[j]);
		}
		auto shift = vec_pcs_value_index.size();
		// mmp values
		for (size_t j=0; j<vec_mmp_id.size(); j++)
		{
			CFiniteElementStd* fem = mat_pcs->GetAssember();
			fem->ConfigElement(ele, false);
			fem->Config();
			fem->SetGaussPoint(0, gp_r, gp_s, gp_t);
			fem->ComputeShapefct(1);
			CMediumProperties* mmp = mmp_vector[ele->GetPatchIndex()];
			double mat_value = ELEMENT_MMP_VALUES::getValue(
				mmp, vec_mmp_id[j], i, gp, theta);
			tmp_values[j + shift] = mat_value;
		}
		shift += vec_mmp_id.size();
		// mfp values
		for (size_t j=0; j<vec_mfp_id.size(); j++)
		{
			CFiniteElementStd* fem = mat_pcs->GetAssember();
			fem->ConfigElement(ele, false);
			fem->Config();
			fem->SetGaussPoint(0, gp_r, gp_s, gp_t);
			fem->ComputeShapefct(1);
			CFluidProperties* mfp = mfp_vector[0];
			mfp->Fem_Ele_Std = fem;
			double mat_value = ELEMENT_MFP_VALUES::getValue(mfp, vec_mfp_id[j]);
			tmp_values[j + shift] = mat_value;
		}
		// output
		os << i << delim;
		writeLine(os, tmp_values, delim);
	}
}

#ifdef USE_PETSC
void writeNodeDataMPI(COutput* output, std::string const& filename)
{
	if (output->_nod_value_vector.empty())
		return;

	std::ofstream* p_os = nullptr;
	if (myrank == 0)
	{
		p_os = new std::ofstream(filename);
		if (!*p_os) {
			ScreenMessage2("***Error: cannot open %s for writing\n", filename.data());
			return;
		}
		p_os->setf(ios::scientific, std::ios::floatfield);
		p_os->precision(12);
	}
	std::ofstream &os = *p_os;


	//
	std::vector<std::string> vec_output_val_name;
	vec_output_val_name.push_back("NodeID");
	std::vector<CRFProcess*> vec_pcs;
	std::vector<int> vec_pcs_value_index;
	for (size_t i=0; i<output->_nod_value_vector.size(); i++)
	{
		auto& value_name = output->_nod_value_vector[i];
		auto pcs = PCSGet(value_name, true);
		if (!pcs)
			continue;
		auto value_id = pcs->GetNodeValueIndex(value_name);
		if (value_id<0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_pcs.push_back(pcs);
		vec_pcs_value_index.push_back(value_id);
	}

	auto const n_nodal_values = vec_pcs_value_index.size();

	std::string const delim = ", ";
	// header
	if (myrank == 0)
		writeLine(os, vec_output_val_name, delim);

	// values
	MeshLib::CFEMesh* msh = output->getMesh();
	std::vector<double> global_node_id_values;

	if (myrank == 0)
	{
		global_node_id_values.resize(msh->getNumNodesGlobal() * n_nodal_values);
		for (size_t i=0; i<msh->GetNodesNumber(false); i++)
		{
			if (!msh->isNodeLocal(i))
				continue;
			MeshLib::CNode* node = msh->getNodeVector()[i];
			for (size_t j=0; j<n_nodal_values; j++)
			{
				auto pcs = vec_pcs[j];
				global_node_id_values[node->GetGlobalIndex()*n_nodal_values + j] = pcs->GetNodeValue(i, vec_pcs_value_index[j]);
			}
		}

		for (int i=1; i<mysize; i++)
		{
			int n_dom_nodes = 0;
			MPI_Recv(&n_dom_nodes, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::vector<int> dom_global_node_ids(n_dom_nodes);
			MPI_Recv(&dom_global_node_ids[0], dom_global_node_ids.size(), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::vector<double> tmp_values(n_nodal_values * n_dom_nodes);
			MPI_Recv(&tmp_values[0], tmp_values.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (size_t j=0; j<dom_global_node_ids.size(); j++)
			{
				for (size_t k=0; k<n_nodal_values; k++)
				{
					global_node_id_values[dom_global_node_ids[j]*n_nodal_values + k] = tmp_values[j*n_nodal_values + k];
				}
			}
		}
	} else {
		int n_local_nodes = msh->getNumNodesLocal();
		MPI_Send(&n_local_nodes, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		std::vector<int> vec_local_node_ids;
		vec_local_node_ids.reserve(n_local_nodes);
		for (size_t i=0; i<msh->GetNodesNumber(false); i++)
		{
			if (!msh->isNodeLocal(i))
				continue;
			vec_local_node_ids.push_back(msh->getNodeVector()[i]->GetGlobalIndex());
		}
		MPI_Send(&vec_local_node_ids[0], vec_local_node_ids.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);

		std::vector<double> tmp_values(n_nodal_values * n_local_nodes);
		int counter = 0;
		for (size_t i=0; i<msh->GetNodesNumber(false); i++)
		{
			if (!msh->isNodeLocal(i))
				continue;
			for (size_t j=0; j<vec_pcs_value_index.size(); j++)
			{
				auto pcs = vec_pcs[j];
				tmp_values[counter*n_nodal_values + j] = pcs->GetNodeValue(i, vec_pcs_value_index[j]);
			}
			counter++;
		}
		MPI_Send(&tmp_values[0], tmp_values.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	if (myrank == 0)
	{
		std::vector<double> tmp_values(n_nodal_values);
		for (long i=0; i<msh->getNumNodesGlobal(); i++)
		{
			for (size_t j=0; j<n_nodal_values; j++)
				tmp_values[j] = global_node_id_values[i*n_nodal_values + j];

			os << i << delim;
			if (myrank == 0)
				writeLine(os, tmp_values, delim);
		}
	}

	delete p_os;
}

void getElementOutputValues(MeshLib::CElem* ele,
							std::vector<CRFProcess*> const& vec_pcs,
							std::vector<int> const& vec_pcs_value_index,
							std::vector<int> const& vec_mmp_id,
							std::vector<int> const& vec_mfp_id,
							std::vector<double> &tmp_values)
{
	double gp[3] = {.0, .0, .0};
	int gp_r=0, gp_s=0, gp_t=0;
	double const theta = 1.0;
	auto mat_pcs = pcs_vector[0];

	ele->SetOrder(false);
	// pcs element values
	for (size_t j=0; j<vec_pcs_value_index.size(); j++)
	{
		auto pcs = vec_pcs[j];
		tmp_values[j] = pcs->GetElementValue(ele->GetIndex(), vec_pcs_value_index[j]);
	}
	auto shift = vec_pcs_value_index.size();
	// mmp values
	for (size_t j=0; j<vec_mmp_id.size(); j++)
	{
		CFiniteElementStd* fem = mat_pcs->GetAssember();
		fem->ConfigElement(ele, false);
		fem->Config();
		fem->SetGaussPoint(0, gp_r, gp_s, gp_t);
		fem->ComputeShapefct(1);
		CMediumProperties* mmp = mmp_vector[ele->GetPatchIndex()];
		double mat_value = ELEMENT_MMP_VALUES::getValue(
			mmp, vec_mmp_id[j], ele->GetIndex(), gp, theta);
		tmp_values[j + shift] = mat_value;
	}
	shift += vec_mmp_id.size();
	// mfp values
	for (size_t j=0; j<vec_mfp_id.size(); j++)
	{
		CFiniteElementStd* fem = mat_pcs->GetAssember();
		fem->ConfigElement(ele, false);
		fem->Config();
		fem->SetGaussPoint(0, gp_r, gp_s, gp_t);
		fem->ComputeShapefct(1);
		CFluidProperties* mfp = mfp_vector[0];
		mfp->Fem_Ele_Std = fem;
		double mat_value = ELEMENT_MFP_VALUES::getValue(mfp, vec_mfp_id[j]);
		tmp_values[j + shift] = mat_value;
	}
}

void writeElementDataMPI(COutput* output, std::string const& filename)
{
	if (output->getElementValueVector().empty())
		return;

	std::ofstream* p_os = nullptr;
	if (myrank == 0)
	{
		p_os = new std::ofstream(filename);
		if (!*p_os) {
			ScreenMessage2("***Error: cannot open %s for writing\n", filename.data());
			return;
		}
		p_os->setf(ios::scientific, std::ios::floatfield);
		p_os->precision(12);
	}
	std::ofstream &os = *p_os;


	//
	std::vector<std::string> vec_output_val_name;
	vec_output_val_name.push_back("ElementID");
	std::vector<CRFProcess*> vec_pcs;
	std::vector<int> vec_pcs_value_index;
	for (auto const& value_name : output->getElementValueVector())
	{
		auto pcs = output->GetPCS_ELE(value_name);
		if (!pcs)
			continue;
		auto value_id = pcs->GetElementValueIndex(value_name);
		if (value_id<0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_pcs.push_back(pcs);
		vec_pcs_value_index.push_back(value_id);
	}

	std::vector<int> vec_mmp_id;
	for (auto const& value_name : output->mmp_value_vector)
	{
		int mmp_id = ELEMENT_MMP_VALUES::getMMPIndex(value_name);
		if (mmp_id < 0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_mmp_id.push_back(mmp_id);
	}

	std::vector<int> vec_mfp_id;
	for (auto const& value_name : output->mfp_value_vector)
	{
		int mfp_id = ELEMENT_MFP_VALUES::getMFPIndex(value_name);
		if (mfp_id < 0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_mfp_id.push_back(mfp_id);
	}

	auto const n_output_values = vec_output_val_name.size() - 1; // exclude elemnt ID

	std::string const delim = ", ";
	// header
	if (myrank == 0)
		writeLine(os, vec_output_val_name, delim);

	// values
	MeshLib::CFEMesh* msh = output->getMesh();

	std::vector<double> global_element_values; // nr. global elements x nr. values

	if (myrank == 0)
	{
		global_element_values.resize(msh->getNumElementsGlobal() * n_output_values);
		std::vector<double> tmp_values(n_output_values);
		for (long i=0; i<(long)msh->getElementVector().size(); i++)
		{
			auto ele = msh->ele_vector[i];
			getElementOutputValues(ele, vec_pcs, vec_pcs_value_index, vec_mmp_id, vec_mfp_id, tmp_values);
			for (size_t j=0; j<tmp_values.size(); j++)
				global_element_values[ele->GetGlobalIndex()*n_output_values + j] = tmp_values[j];
		}

		for (int i=1; i<mysize; i++)
		{
			int n_dom_eles = 0;
			MPI_Recv(&n_dom_eles, 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::vector<int> dom_global_ele_ids(n_dom_eles);
			MPI_Recv(&dom_global_ele_ids[0], dom_global_ele_ids.size(), MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::vector<double> tmp_values(n_output_values * n_dom_eles);
			MPI_Recv(&tmp_values[0], tmp_values.size(), MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			for (size_t j=0; j<dom_global_ele_ids.size(); j++)
			{
				for (size_t k=0; k<n_output_values; k++)
				{
					global_element_values[dom_global_ele_ids[j]*n_output_values + k] = tmp_values[j*n_output_values + k];
				}
			}
		}
	} else {
		int n_local_ele = msh->getElementVector().size();
		MPI_Send(&n_local_ele, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
		std::vector<int> vec_local_ele_ids;
		vec_local_ele_ids.reserve(n_local_ele);
		for (auto ele : msh->getElementVector())
			vec_local_ele_ids.push_back(ele->GetGlobalIndex());
		MPI_Send(&vec_local_ele_ids[0], vec_local_ele_ids.size(), MPI_INT, 0, 0, MPI_COMM_WORLD);

		std::vector<double> local_values(n_output_values * n_local_ele);
		std::vector<double> tmp_values(n_output_values);
		for (long i=0; i<(long)msh->getElementVector().size(); i++)
		{
			auto ele = msh->ele_vector[i];
			getElementOutputValues(ele, vec_pcs, vec_pcs_value_index, vec_mmp_id, vec_mfp_id, tmp_values);
			for (size_t j=0; j<tmp_values.size(); j++)
				local_values[i*n_output_values + j] = tmp_values[j];
		}
		MPI_Send(&local_values[0], local_values.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}

	if (myrank == 0)
	{
		std::vector<double> tmp_values(n_output_values);
		for (long i=0; i<msh->getNumElementsGlobal(); i++)
		{
			for (size_t j=0; j<n_output_values; j++)
				tmp_values[j] = global_element_values[i*n_output_values + j];

			os << i << delim;
			if (myrank == 0)
				writeLine(os, tmp_values, delim);
		}
	}

	delete p_os;
}
#endif

} // namespace

#ifdef USE_PETSC

void CSVOutput::writeDomain(COutput* output,
                      int timestep_number,
                      double /*current_time*/,
                      std::string const& baseFilename)
{
	std::string filename_base = baseFilename;
	if (output->getProcessType() != FiniteElement::INVALID_PROCESS)
		filename_base += "_" + FiniteElement::convertProcessTypeToString(output->getProcessType());

	std::string filename_node_data = filename_base + "_node";
	filename_node_data += "_" + std::to_string(timestep_number);
	filename_node_data += ".csv";

	writeNodeDataMPI(output, filename_node_data);

	std::string filename_ele_data = filename_base + "_ele";
	filename_ele_data += "_" + std::to_string(timestep_number);
	filename_ele_data += ".csv";
	writeElementDataMPI(output, filename_ele_data);

}

#else
void CSVOutput::writeDomain(COutput* output,
                      int timestep_number,
                      double /*current_time*/,
                      std::string const& baseFilename)
{
	std::string filename_base = baseFilename;
	if (output->getProcessType() != FiniteElement::INVALID_PROCESS)
		filename_base += "_" + FiniteElement::convertProcessTypeToString(output->getProcessType());

	std::string filename_node_data = filename_base + "_node";
#ifdef USE_PETSC
	filename_node_data += "_part" + std::to_string(myrank);
#endif
	filename_node_data += "_" + std::to_string(timestep_number);
	filename_node_data += ".csv";
	writeNodeData(output, filename_node_data);

	std::string filename_ele_data = filename_base + "_ele";
#ifdef USE_PETSC
	filename_ele_data += "_part" + std::to_string(myrank);
#endif
	filename_ele_data += "_" + std::to_string(timestep_number);
	filename_ele_data += ".csv";
	writeElementData(output, filename_ele_data);

}

#endif

void CSVOutput::writePoyline(COutput* output,
							   int timestep_number,
							   double /*current_time*/,
							   std::string const& baseFilename)
{
	if (output->_nod_value_vector.empty())
		return;

	//----------------------------------------------------------------------------------
	GEOLIB::Polyline const* const ply(
	    dynamic_cast<GEOLIB::Polyline const* const>(output->getGeoObj()));
	if (output->getGeoType() != GEOLIB::POLYLINE || ply == NULL)
	{
		std::cerr
		    << "COutput::NODWritePLYDataTEC geometric object is not a polyline"
		    << "\n";
		return;
	}

	MeshLib::CFEMesh* msh = output->getMesh();
	msh->SwitchOnQuadraticNodes(false);
	std::vector<long> nodes_vector;

	CGLPolyline* m_ply = GEOGetPLYByName(output->getGeoName());
	double tmp_min_edge_length(msh->getMinEdgeLength());
	msh->setMinEdgeLength(m_ply->epsilon);
	msh->GetNODOnPLY(ply, nodes_vector);
	if (nodes_vector.empty()) {
		ScreenMessage2("-> Found no nodes on polyline %s\n", output->getGeoName().data());
		return;
	}

	std::vector<double> interpolation_points;
	msh->getPointsForInterpolationAlongPolyline(ply, interpolation_points);
	msh->setMinEdgeLength(tmp_min_edge_length);



	//----------------------------------------------------------------------------------
	std::string filename_base = baseFilename;
	if (output->getProcessType() != FiniteElement::INVALID_PROCESS)
		filename_base += "_" + FiniteElement::convertProcessTypeToString(output->getProcessType());

	std::string filename_node_data = filename_base + "_ply_" + output->getGeoName();
	filename_node_data += "_" + std::to_string(timestep_number);
	filename_node_data += ".csv";


	std::ofstream os(filename_node_data);
	if (!os) {
		ScreenMessage2("***Error: cannot open %s for writing\n", filename_node_data.data());
		return;
	}
	os.setf(ios::scientific, std::ios::floatfield);
	os.precision(12);

	//----------------------------------------------------------------------------------
	std::vector<std::string> vec_output_val_name;
	vec_output_val_name.push_back("DIST");
	std::vector<CRFProcess*> vec_pcs;
	std::vector<int> vec_pcs_value_index;
	for (size_t i=0; i<output->_nod_value_vector.size(); i++)
	{
		auto& value_name = output->_nod_value_vector[i];
		auto pcs = PCSGet(value_name, true);
		if (!pcs)
			continue;
		auto value_id = pcs->GetNodeValueIndex(value_name);
		if (value_id<0)
			continue;
		vec_output_val_name.push_back(value_name);
		vec_pcs.push_back(pcs);
		vec_pcs_value_index.push_back(value_id);
	}

	auto const n_nodal_values = vec_pcs_value_index.size();

	std::string const delim = ", ";
	// header
	writeLine(os, vec_output_val_name, delim);

	// values
	std::vector<double> tmp_values(n_nodal_values);
	for (size_t i=0; i<nodes_vector.size(); i++)
	{
		auto nodeID = nodes_vector[i];
		os << interpolation_points[i] << delim;
		for (size_t j=0; j<vec_pcs_value_index.size(); j++)
		{
			auto pcs = vec_pcs[j];
			tmp_values[j] = pcs->GetNodeValue(nodeID, vec_pcs_value_index[j]);
		}
		writeLine(os, tmp_values, delim);
	}
}
