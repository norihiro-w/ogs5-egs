/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "rf_out_new.h"

#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <list>
#include <string>
#include <vector>

#if defined(USE_PETSC) || defined(USE_MPI)
#include <mpi.h>
#endif

#include "Configure.h"
#include "display.h"
#include "makros.h"
#include "StringTools.h"

#include "GeoIO.h"
#include "GEOObjects.h"
#include "geo_ply.h"
#include "geo_sfc.h"

#include "msh_lib.h"

#include "CSVOutput.h"
#include "fem_ele_std.h"
#include "files0.h"
#include "LegacyVtkInterface.h"
#include "mathlib.h"
#include "Output.h"
#include "problem.h"
#include "rf_msp_new.h"
#include "rf_pcs.h"
#include "rf_random_walk.h"
#include "rf_tim_new.h"
#include "vtk.h"


#ifdef SUPERCOMPUTER
// kg44 this is usefull for io-buffering as endl flushes the buffer
#define endl '\n'
#define MY_IO_BUFSIZE 4096
#endif

#define OUT_FILE_EXTENSION ".out"

std::vector<COutput*> out_vector;
extern size_t max_dim;

using namespace std;
using MeshLib::CFEMesh;

/**************************************************************************
   FEMLib-Method:
   Task: OUT read function
   Programing:
   06/2004 OK Implementation
   08/2004 WW Remove the old files
   01/2005 OK Boolean type
   01/2005 OK Destruct before read
   06/2006 WW Remove the old files by new way
   06/2010 TF reformated, restructured, signature changed, use new GEOLIB data
structures
**************************************************************************/
bool OUTRead(const std::string& file_base_name,
             const GEOLIB::GEOObjects& geo_obj,
             const std::string& unique_name)
{
	char line[MAX_ZEILE];
	std::string line_string;
	ios::pos_type position;
	bool output_version = false;  // 02.2011. WW

#if defined(USE_PETSC) || \
    defined(USE_MPI)  //|| defined(other parallel libs)//03.3012. WW
	int rank, msize;
	string rank_str;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &msize);
	std::ifstream is;
	stringstream ss(stringstream::in | stringstream::out);
	ss.clear();
	ss.str("");
	ss << rank;
	rank_str = ss.str();
	ss.clear();
#endif
	// File handling
	std::string out_file_name = file_base_name + OUT_FILE_EXTENSION;
	std::ifstream out_file(out_file_name.data(), ios::in);
	if (!out_file.good()) return false;
	out_file.seekg(0L, ios::beg);

	// Keyword loop
	ScreenMessage("OUTRead\n");
	while (!out_file.eof())
	{
		out_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != string::npos) return true;

		COutput* out = NULL;
		// Give version in file name
		// 15.01.2008. WW
		if (line_string.find("#VERSION") != string::npos)
			output_version = true;  // 02.2011. WW
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#OUTPUT") != string::npos)
		{
			out = new COutput(out_vector.size());
#if defined(USE_PETSC) || \
    defined(USE_MPI)  //|| defined(other parallel libs)//03.3012. WW
			out->setMPI_Info(rank, msize, rank_str);
#endif
			out->getFileBaseName() = file_base_name;
			position = out->Read(out_file, geo_obj, unique_name);

			if (output_version)  //// 02.2011. WW
			{
				std::string VersionStr = OGS_VERSION;  // 02.2011 WX
				int curPos = 0;
				int pos = 0;
				while ((pos = VersionStr.find("/", curPos)) != -1)
				{
					VersionStr.replace(pos, 1, "_");
					curPos = pos + 1;
				}
				out->getFileBaseName().append("(V");
				out->getFileBaseName().append(VersionStr);
				out->getFileBaseName().append(")");
			}

			out_vector.push_back(out);

			//			char number_char[3]; //OK4709
			//			sprintf(number_char, "%i", (int) out_vector.size() - 1);
			////OK4709
			//			out->ID = number_char; //OK4709
			//			out->setID (out_vector.size() - 1);

			out_file.seekg(position, ios::beg);
		}  // keyword found
	}      // eof
	return true;
}

/**************************************************************************
   FEMLib-Method: OUTWrite
   Task: master write function
   Programing:
   06/2004 OK Implementation
   last modification:
**************************************************************************/
void OUTWrite(string base_file_name)
{
	//========================================================================
	// File handling
	string out_file_name = base_file_name + OUT_FILE_EXTENSION;
	fstream out_file(out_file_name.data(), ios::trunc | ios::out);
	out_file.setf(ios::scientific, ios::floatfield);
	out_file.precision(12);
	if (!out_file.good()) return;
	out_file.seekg(0L, ios::beg);
#ifdef SUPERCOMPUTER
	// kg44 buffer the output
	char mybuffer[MY_IO_BUFSIZE * MY_IO_BUFSIZE];
	out_file.rdbuf()->pubsetbuf(mybuffer, MY_IO_BUFSIZE * MY_IO_BUFSIZE);
//
#endif
	//========================================================================
	out_file << "GeoSys-OUT: Output "
	            "------------------------------------------------\n";
	//========================================================================
	// OUT vector
	size_t out_vector_size(out_vector.size());
	for (size_t i = 0; i < out_vector_size; i++)
		out_vector[i]->Write(&out_file);
	out_file << "#STOP";
	out_file.close();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   08/2004 OK Implementation
   08/2004 WW Output by steps given in .out file
   03/2005 OK RFO format
   05/2005 OK MSH
   05/2005 OK Profiles at surfaces
   12/2005 OK VAR,MSH,PCS concept
   03/2006 WW Flag to remove existing files
   08/2006 OK FLX calculations
   08/2007 WW Output initial values of variables
**************************************************************************/
void OUTData(double time_current, int time_step_number, bool force_output)
{
	for (size_t i = 0; i < out_vector.size(); i++)
	{
		COutput* m_out = out_vector[i];
		m_out->setTime(time_current);

		CFEMesh* m_msh = m_out->getMesh();
		if (!m_msh)
			cout << "Warning in OUTData - no MSH data" << endl;

		bool doOutput = false;
		if (time_step_number == 0 || force_output)
		{
			doOutput = true;
		}
		else if (m_out->nSteps > 0 && (time_step_number % m_out->nSteps == 0))
		{
			doOutput = true;
		}
		else if (!m_out->time_vector.empty())
		{
			auto it = std::find_if(m_out->time_vector.begin(), m_out->time_vector.end(), [&](double t){return fabs(time_current-t)<MKleinsteZahl;});
			if (it!=m_out->time_vector.end())
				doOutput = true;
		}

		if (!doOutput)
			continue;

		//======================================================================
		// TECPLOT
		if (m_out->dat_type_name.compare("TECPLOT") == 0 ||
		    m_out->dat_type_name.compare("MATLAB")  == 0 ||
		    m_out->dat_type_name.compare("GNUPLOT") == 0)
		{

			switch (m_out->getGeoType())
			{
				case GEOLIB::GEODOMAIN:
					ScreenMessage("-> Data output: Domain");
					if (!m_out->_pcon_value_vector.empty())
						m_out->PCONWriteDOMDataTEC();
					else
					{
						m_out->NODWriteDOMDataTEC();
						m_out->ELEWriteDOMDataTEC();
					}
					if (!m_out->_new_file_opened)
						m_out->_new_file_opened = true;
					break;
				//------------------------------------------------------------------
				case GEOLIB::POLYLINE:
				{
					ScreenMessage("-> Data output: Polyline profile - %s\n", m_out->getGeoName().c_str());
					double tim_value = m_out->NODWritePLYDataTEC(time_step_number);
					if (tim_value > 0.0)
						m_out->TIMValue_TEC(tim_value);
					if (!m_out->_new_file_opened)
						m_out->_new_file_opened = true;
					//..............................................................
					break;
				}
				//------------------------------------------------------------------
				case GEOLIB::POINT:
					ScreenMessage("-> Data output: Breakthrough curves - %s\n", m_out->getGeoName().c_str());
					m_out->NODWritePNTDataTEC(time_current, time_step_number);
					if (!m_out->_new_file_opened)
						m_out->_new_file_opened = true;
					break;
				//------------------------------------------------------------------
				case GEOLIB::SURFACE:
					ScreenMessage("-> Data output: Surface profile\n");
					if (m_out->getProcessDistributionType() == FiniteElement::AVERAGE)
					{
						m_out->NODWriteSFCAverageDataTEC(time_current,
														 time_step_number);
						if (!m_out->_new_file_opened)
							m_out->_new_file_opened = true;
					}
					//..............................................................
					else
					{
						m_out->NODWriteSFCDataTEC(time_step_number);
						if (!m_out->_new_file_opened)
							m_out->_new_file_opened = true;
					}
					//..............................................................
					// ELE data
					if (m_out->getElementValueVector().size() > 0)
						m_out->ELEWriteSFC_TEC();
					//..............................................................
					break;
				default:
					break;
			}
		}
		//--------------------------------------------------------------------
		// CSV
		else if (m_out->dat_type_name.compare("CSV") == 0)
		{
			switch (m_out->getGeoType())
			{
				case GEOLIB::GEODOMAIN:
				{
					ScreenMessage("-> Data output: Domain - CSV\n");
					CSVOutput::writeDomain(m_out, time_step_number, m_out->_time, m_out->file_base_name);
					break;
				}
				case GEOLIB::POINT:
				{
					ScreenMessage("-> Data output: Point - %s\n", m_out->getGeoName().c_str());
					m_out->NODWritePNTDataTEC(time_current, time_step_number);
					if (!m_out->_new_file_opened)
						m_out->_new_file_opened = true;
					break;
				}
				case GEOLIB::POLYLINE:
				{
					ScreenMessage("-> Data output: Polyline - %s\n", m_out->getGeoName().c_str());
					CSVOutput::writePoyline(m_out, time_step_number, m_out->_time, m_out->file_base_name);
					break;
				}
				default:
					ScreenMessage("Warning - the given GEO_TYPE is not supported in CSV output\n");
					break;
			}
		}
		//--------------------------------------------------------------------
		// vtk
		else if (m_out->dat_type_name.compare("VTK") == 0)
		{
			LegacyVtkInterface vtkOutput(m_msh,
										 m_out->_nod_value_vector,
										 m_out->_ele_value_vector,
										 m_out->mmp_value_vector,
										 m_out->msh_type_name,
										 m_out);
#if defined(USE_PETSC)
			ScreenMessage2("-> output a VTK file\n");
			vtkOutput.WriteDataVTKPETSC(time_step_number,
										m_out->_time,
										m_out->file_base_name);
#else
			vtkOutput.WriteDataVTK(time_step_number,
								   m_out->_time,
								   m_out->file_base_name);
#endif
			if (!m_out->_new_file_opened)
				m_out->_new_file_opened = true;
		}  // PVD (ParaView)
		else if (m_out->dat_type_name.find("PVD") != string::npos)
		{
			if (m_out->vtk == NULL)
				m_out->CreateVTKInstance();  // WW m_out->vtk = new CVTK();
			CVTK* vtk = m_out->vtk;

			bool vtk_appended = false;
			if (m_out->dat_type_name.find("PVD_A") != string::npos)
				vtk_appended = true;

			if (m_out->getGeoType() == GEOLIB::GEODOMAIN)
			{
				ScreenMessage("-> Data output: Domain - PVD\n");
				if (time_step_number == 0)
				{
					std::string pcs_type("");
					if (m_out->getProcessType() !=
					    FiniteElement::INVALID_PROCESS)
						pcs_type = FiniteElement::convertProcessTypeToString(
						    m_out->getProcessType());
					vtk->InitializePVD(
					    m_out->file_base_name, pcs_type, vtk_appended);
				}
				// Set VTU file name and path
				const std::string vtk_file_path_base =
				    vtk->pvd_vtk_file_path_base + vtk->pvd_vtk_file_name_base;
#ifdef USE_PETSC
				std::string vtu_file_name =
				    vtk_file_path_base + "_part" + number2str(myrank) + "_" +
				    number2str(time_step_number) + ".vtu";
				std::string pvd_vtk_file_name = vtk_file_path_base + "_" +
				                                number2str(time_step_number) +
				                                ".pvtu";
#else
				std::string vtu_file_name = vtk_file_path_base + "_" +
				                            number2str(time_step_number) +
				                            ".vtu";
				std::string pvd_vtk_file_name = vtu_file_name;
#endif
				vtk->WriteXMLUnstructuredGrid(
				    vtu_file_name, m_out, time_step_number);
#ifdef USE_PETSC
				if (myrank == 0)
					vtk->WriteXMLPUnstructuredGrid(
					    vtk_file_path_base, m_out, time_step_number);
#endif
				VTK_Info dat;
				dat.timestep = m_out->getTime();
				dat.vtk_file = pvd_vtk_file_name;
				vtk->vec_dataset.push_back(dat);
#ifdef USE_PETSC
				if (myrank == 0)
				{
#endif
					vtk->UpdatePVD(vtk->pvd_file_name, vtk->vec_dataset);
#ifdef USE_PETSC
				}
#endif
			}
		}
		else if (m_out->dat_type_name.find("PETREL") != string::npos)
		{
			m_out->WritePetrelElementData(time_step_number);
			if (!m_out->_new_file_opened)
				m_out->_new_file_opened = true;
		}
		else if (m_out->dat_type_name.compare("WATER_BALANCE") == 0)
			m_out->NODWriteWaterBalance(time_current);
		else if (m_out->dat_type_name.compare("COMBINE_POINTS") == 0)
			m_out->NODWritePointsCombined(time_current);

	}  // OUT loop
	//======================================================================
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   01/2005 OK Implementation
   last modified:
**************************************************************************/
void OUTDelete()
{
	const size_t no_out = out_vector.size();
	for (size_t i = 0; i < no_out; i++)
		delete out_vector[i];
	out_vector.clear();
}

/**************************************************************************
   FEMLib-Method:
   Task:
   Programing:
   05/2005 OK Implementation
   last modification: 03/2010 JT
   09/2010 TF
**************************************************************************/
COutput* OUTGet(const std::string& out_name)
{
	FiniteElement::ProcessType pcs_type(
	    FiniteElement::convertProcessType(out_name));
	for (size_t i = 0; i < out_vector.size(); i++)
		if (out_vector[i]->getProcessType() == pcs_type) return out_vector[i];
	return NULL;
}

/**************************************************************************
   FEMLib-Method:
   Task: Return output object for variable name (based on OUTGet)
   Programing:
   03/2010 JT Implementation
   last modification:
**************************************************************************/
COutput* OUTGetRWPT(const std::string& out_name)
{
	for (size_t i = 0; i < out_vector.size(); i++)
	{
		COutput* out(out_vector[i]);
		for (size_t j = 0;
		     j < out->getRandomWalkParticleTracingValueVector().size();
		     j++)
			if (out->getRandomWalkParticleTracingValueVector()[j].compare(
			        out_name) == 0)
				return out;
	}
	return NULL;
}

/***************************************************************************************
   Function checking the consistency of the output data as specified in the
 input file
   This means up to now, that data for missing processes is not written
   05/09	SB
 *****************************************************************************************/
void OUTCheck()
{
	ScreenMessage("Checking output data ...\n");
	// Go through all out objects (#OUTPUT-section in input file)
	for (size_t i = 0; i < out_vector.size(); i++)
		out_vector[i]->checkConsistency();
	ScreenMessage("done\n");
}
