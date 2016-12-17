/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: files0.c
 */
/* Aufgabe:
   Enthaelt die uebergeordneten Datei- Ein- und Ausgaberoutinen, sowie
   das Speichern der Durchbruchskurven.
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   02/1999     CT         Bugfix: Anpassung fuer mit 0 beginnende Elemente.
                          Kennung fuer Version im Kopf. Versionsabhaengiges
                          lesen.
   02/1999     CT         Korrekturen fuer herausgeloeste Materialgruppen
   03/1999     AH         Korrekturen nicht noetig, da die Materialzuweisung
                          nicht mehr auf dieser Ebene erfolgt. Die Abfrage ist
                          jetzt auskommentiert. Sie war hier nur um
   Kompatibilitaeten
                          zum Rockflow aufzubewahren.
   03/1999     CT         anz_matxx, start_mat entfernt
   02/2000     CT         Restart wieder hergestellt
   07/2000     AH         Vorbereitungen zum HGM.
   9/2000     CT         Neu: RefreshNodeOutputData, Warnungen beseitigt
   10/2001     AH         Inverse Modellierung
   Trennung und Anpassung (CreateFileData)
   In DestroyFileData Datenfeld auskommentiert.
   Neues Konzept fuer Datenbank-Verwaltung ah inv
   01/2002     MK         DisplayMsgX-Umleitung in *.msg-Datei: OpenMsgFile
   08/2002     MK         GetPathRFDFile
   ConfigFileData aus CreateFileData herausgeloest
   03/2003     RK         Quellcode bereinigt, Globalvariablen entfernt
 */
/**************************************************************************/

#include "files0.h"

#include "Configure.h"
//#include <iostream>
#include "MemWatch.h"
#include "display.h"
#include "memory.h"
#include "FileTools.h"
#include "FileToolsRF.h"
#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"
#include "rfmat_cp.h"
#include "tools.h"
//#include "rf_pcs.h"
#include "rf_out_new.h"
//#include "rf_tim_new.h"
#include "rf_mfp_new.h"
#include "rf_msp_new.h"
//#include "rf_num_new.h"
#include "rf_fct.h"             //OK
#include "rf_fluid_momentum.h"  // PCH
#include "rf_kinreact.h"
#include "rf_random_walk.h"  // PCH
#include "rf_react.h"
#ifdef CHEMAPP
#include "eqlink.h"  //MX
#endif
#include "fct_mpi.h"
/* Tools */
//#include "mathlib.h"
//#include "femlib.h"
/* GeoLib */
//#include "geo_lib.h"
// MSHLib
//#include "msh_lib.h"
//#include "gs_project.h"
/* Dateinamen */
char* crdat = NULL;     /*MX*/
char* file_name = NULL; /* dateiname */
static char* msgdat = NULL;

#define RFD_FILE_EXTENSION ".rfd"  // OK
#ifndef MFC                        // WW
void CURRead(std::string);         // OK
#endif
std::ios::pos_type CURReadCurve(std::ifstream*);  // OK
void CURWrite();                                  // OK

#define KEYWORD '#'
#define SUBKEYWORD '$'

// GEOLIB
//#include "GEOObjects.h"

// FileIO
#include "OGSIOVer4.h"
#include "readNonBlankLineFromInputStream.h"
//#include "FEMIO.h"

using namespace std;

/**************************************************************************/
/* ROCKFLOW - Funktion: ReadString
 */
/* Aufgabe:
   Liest Zeichenkette von Standardeingabe
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   - void -
 */
/* Ergebnis:
   gelesene Zeichenkette
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
 */
/**************************************************************************/
char* ReadString(void)
{
	char* s = (char*)malloc(256);
	// char *s = new char[256];//CC
	scanf(" %s%*[^\n]%*c", s);
	//  int a = (int)strlen(s);
	//  delete[] s;
	// s = new char[a+1];//CC
	s = (char*)realloc(s, ((int)strlen(s) + 1));
	return s;
}

static bool isValidTextFileFormat(const std::string& basename,
                                  const std::string& fext)
{
	const std::string fname(basename + fext);
	if (!IsFileExisting(fname)) return true;
#ifdef _WIN32
	const bool is_win32 = true;
#else
	const bool is_win32 = false;
#endif
	if (is_win32 == HasCRInLineEnding(fname))
	{
		return true;
	}
	else
	{
		if (is_win32)
			ScreenMessage("*** ERROR: Detect UNIX file format %s\n",
			              fname.data());
		else
			ScreenMessage("*** ERROR: Detect Windows file format %s\n",
			              fname.data());
		return false;
	}
}

static bool checkFormatOfInputFiles(const std::string& basename)
{
	bool valid = true;
	valid &= isValidTextFileFormat(basename, ".gli");
	valid &= isValidTextFileFormat(basename, ".msh");
	valid &= isValidTextFileFormat(basename, ".pcs");
	valid &= isValidTextFileFormat(basename, ".ic");
	valid &= isValidTextFileFormat(basename, ".bc");
	valid &= isValidTextFileFormat(basename, ".st");
	valid &= isValidTextFileFormat(basename, ".mfp");
	valid &= isValidTextFileFormat(basename, ".msp");
	valid &= isValidTextFileFormat(basename, ".mmp");
	valid &= isValidTextFileFormat(basename, ".mcp");
	valid &= isValidTextFileFormat(basename, ".out");
	valid &= isValidTextFileFormat(basename, ".tim");

	return valid;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: ReadData
 */
/* Aufgabe:
   Liest Daten aus den Eingabedateien ein
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *dateiname: Dateiname ohne Extension
 */
/* Ergebnis:
   0 bei Fehler oder Ende aufgrund Dateitest, sonst 1
 */
/* Programmaenderungen:
   07/1996    MSR     Erste Version
   07/2000    AH      Eingabe fuer Netzgenerator
   10/2001    AH      Trennung und Anpassung (CreateFileData)
                      Inverse Modellierung
                      Neues Konzept fuer Datenbank-Verwaltung
   10/2002   OK   DATCreateFileNames(dateiname)
                  DATDestroyFileNames()

   last modified: OK 16.10.2002
 */
/**************************************************************************/
int ReadData(char* dateiname,
             GEOLIB::GEOObjects& geo_obj,
             std::string& unique_name)
{
	ScreenMessage("\n---------------------------------------------\n");
	ScreenMessage("Data input:\n");

	/* Dateinamen generieren */
	// OK  DATCreateFileNames(dateiname);
	static int datlen;
	datlen = (int)strlen(dateiname) + 5;
	crdat = (char*)Malloc(datlen); /*MX*/
	                               /*MX*/
	crdat = strcat(strcpy(crdat, dateiname), CHEM_REACTION_EXTENSION);
	msgdat = (char*)Malloc(datlen);
	msgdat = strcat(strcpy(msgdat, dateiname), RF_MESSAGE_EXTENSION);
	FILE* f = NULL;
	if ((f = fopen(msgdat, "r")) == NULL) /* MSG-Datei existiert nicht */
		msgdat = (char*)Free(msgdat);
	else
	{
		fclose(f);
		if ((f = fopen(msgdat, "a")) ==
		    NULL) /* MSG-Schreibzugriff nicht moeglich */
			msgdat = (char*)Free(msgdat);
		else
			fclose(f);
	}
	if (!checkFormatOfInputFiles(dateiname))
	{
		ScreenMessage("terminate this program");
		exit(0);
	}
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	ScreenMessage("\tcurrent mem: %d MB\n",
	              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
	//----------------------------------------------------------------------
	// Read GEO data
	ScreenMessage("GLIRead\n");
	GEOLIB_Read_GeoLib(dateiname);

	std::string geo_file_name(dateiname);
	geo_file_name += ".gli";
	std::vector<std::string> file_read_errors;
	FileIO::readGLIFileV4(
	    geo_file_name, &geo_obj, unique_name, file_read_errors);
#ifndef WIN32
	ScreenMessage("\tcurrent mem: %d MB\n",
	              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif

	//----------------------------------------------------------------------
	// Read object data
	PCSRead(dateiname);
	MFPRead(dateiname);
	// HS PCS immediately followed by the MCP read
	CPRead(dateiname);  // SB:GS4
	BCRead(dateiname, geo_obj, unique_name);
	STRead(dateiname, geo_obj, unique_name);
	ICRead(dateiname, geo_obj, unique_name);
	OUTRead(dateiname, geo_obj, unique_name);
	TIMRead(dateiname);

	MSPRead(dateiname);
	MMPRead(dateiname);
	RCRead(dateiname);
	KRRead(dateiname, geo_obj, unique_name);
	KRWrite(dateiname);
#ifdef CHEMAPP
	CHMRead(dateiname);  // MX for CHEMAPP
#endif
	NUMRead(dateiname);

	FEMDeleteAll();  // KR moved from FEMRead()
	std::vector<CFEMesh*> mesh_vec;
	FEMRead(dateiname, mesh_vec, &geo_obj, &unique_name);
#ifndef WIN32
	ScreenMessage("\tcurrent mem: %d MB\n",
	              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
	if (mesh_vec.empty())  // KR
	{
		abort();
	}
	fem_msh_vector.insert(fem_msh_vector.end(),
						  mesh_vec.begin(),
						  mesh_vec.end());  // re-inserted by KR
#ifndef USE_PETSC
	CompleteMesh();  // WW
#else
#if 0
	ScreenMessage("Optimize geometric objects\n");
	geo_obj.optimiseObjects(unique_name, *fem_msh_vector[0]->getGrid());
	ScreenMessage("\tcurrent mem: %d MB\n", mem_watch.getVirtMemUsage() / (1024*1024) );
#endif
#endif

// SBOK4209 MSHWrite(dateiname);
// PCTRead is bounded by msh
#ifndef OGS_ONLY_TH
	PCTRead(dateiname);  // PCH
#endif
	FMRead(dateiname);   // PCH
	FCTRead(dateiname);  // OK
	CURRead(dateiname);  // OK
// CURWrite(); //OK
#ifdef USE_PETSC
	FCT_MPI::FCTCommRead(dateiname);
#endif
	//----------------------------------------------------------------------
	// Read Excel/CVS data
	// PNTPropertiesRead(dateiname);

	msgdat = (char*)Free(msgdat);

	return 1;
}

/**************************************************************************
   GeoSys-Method: FEMOpen->RFDOpen
   Task:
   Programing:
   11/2003 OK Implementation
   08/2004 OK PCS2
   01/2005 OK Boolean type
**************************************************************************/
bool RFDOpen(std::string file_name_base)
{
	(void)file_name_base;
	return false;
}



/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
void CURRead(std::string base_file_name)
{
	char line[MAX_ZEILE];
	std::string sub_line;
	std::string line_string;
	std::ios::pos_type position;
	//----------------------------------------------------------------------
	StuetzStellen* stuetz = NULL;
	anz_kurven = 1;
	stuetz = (StuetzStellen*)Malloc(sizeof(StuetzStellen));
	stuetz[0].punkt = 1.0;
	stuetz[0].wert = 1.0;
	kurven = (Kurven*)Malloc(sizeof(Kurven));
	kurven[anz_kurven - 1].anz_stuetzstellen = 1;
	kurven[anz_kurven - 1].stuetzstellen = stuetz;
	//----------------------------------------------------------------------
	// file handling
	std::string cur_file_name;
	cur_file_name = base_file_name + RFD_FILE_EXTENSION;
	std::ifstream cur_file(cur_file_name.data(), std::ios::in);
	if (!cur_file.good()) return;
	cur_file.seekg(0L, std::ios::beg);
	//========================================================================
	// keyword loop
	std::cout << "CURRead"
	          << "\n";
	while (!cur_file.eof())
	{
		cur_file.getline(line, MAX_ZEILE);
		line_string = line;
		if (line_string.find("#STOP") != std::string::npos) return;
		//----------------------------------------------------------------------
		// keyword found
		if (line_string.find("#CURVE") != std::string::npos)
		{
			position = CURReadCurve(&cur_file);
			cur_file.seekg(position, std::ios::beg);
		}  // keyword found
	}      // eof
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
std::ios::pos_type CURReadCurve(std::ifstream* cur_file)
{
	bool new_keyword = false;
	std::string hash("#");
	std::string line_string;
	ios::pos_type position;
	std::stringstream line_stream;
	int anz = 0;
	double d1, d2;
	StuetzStellen* stuetz = NULL;
	//----------------------------------------------------------------------
	while (!new_keyword)
	{
		position = cur_file->tellg();
		// OK    cur_file->getline(buffer,MAX_ZEILE);
		// OK    line_string = buffer;
		line_string = GetLineFromFile1(cur_file);
		if (line_string.size() < 1) continue;
		//....................................................................
		// Test next keyword
		if (line_string.find(hash) != string::npos)
		{
			new_keyword = true;
			continue;
		}
		//--------------------------------------------------------------------
		if (line_string.find(";") != string::npos) continue;
		//--------------------------------------------------------------------
		// DATA
		// OK    cur_file->seekg(position,ios::beg);
		// OK    *cur_file >> d1 >> d2;
		line_stream.str(line_string);
		line_stream >> d1 >> d2;
		anz++;
		stuetz = (StuetzStellen*)Realloc(stuetz, (anz * sizeof(StuetzStellen)));
		stuetz[anz - 1].punkt = d1;
		stuetz[anz - 1].wert = d2;
		line_stream.clear();
		//--------------------------------------------------------------------
	}
	//----------------------------------------------------------------------
	if (anz >= 1l)
	{
		anz_kurven++;
		kurven = (Kurven*)Realloc(kurven, (anz_kurven * sizeof(Kurven)));
		kurven[anz_kurven - 1].anz_stuetzstellen = anz;
		kurven[anz_kurven - 1].stuetzstellen = stuetz;
	}
	return position;
}

/**************************************************************************
   FEMLib-Method:
   04/2007 OK Implementation
**************************************************************************/
void CURWrite()
{
	//========================================================================
	// File handling
	std::string fct_file_name = "test.cur";
	std::fstream fct_file(fct_file_name.c_str(), ios::trunc | ios::out);
	fct_file.setf(ios::scientific, ios::floatfield);
	fct_file.precision(12);
	if (!fct_file.good()) return;
	fct_file << "GeoSys-CUR: Functions "
	            "------------------------------------------------"
	         << "\n";
	//========================================================================
	int j;
	StuetzStellen stuetz;
	for (int i = 0; i < anz_kurven; i++)
	{
		fct_file << "#CURVES" << endl;
		for (j = 0; j < kurven[i].anz_stuetzstellen; j++)
		{
			stuetz = kurven[i].stuetzstellen[j];
			fct_file << stuetz.punkt << " " << stuetz.wert << endl;
		}
	}
	fct_file << "#STOP";
}

