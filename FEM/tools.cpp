/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "tools.h"

#include <cfloat>

#include "makros.h"
#include "display.h"
#include "FileToolsRF.h"
#include "memory.h"

#include "Curve.h"
#include "mathlib.h"

#include "msh_elem.h"
#include "msh_lib.h"

#include "femlib.h"
#include "files0.h"
#include "rf_fct.h"
#include "rf_mmp_new.h"
#include "rf_num_new.h"
#include "rf_tim_new.h"

using namespace std;

double* fracture_aperture_array = NULL;
hetfields* hf = NULL;
long fracture_aperture_anz = 0l;

#define TIMER_CEN_LIST "CEN_LIST"

/**************************************************************************/
/* ROCKFLOW - Funktion: FctCurves
 */
/* Aufgabe:
   Liest die zu dem Schluesselwort CURVES gehoerigen Daten ein und erstellt
   den zugehoerigen Protokollabschnitt.
   CURVES: Kurven fuer funktionale Zusammenhaenge
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *data: Zeiger auf die aus der Datei eingelesenen Zeichen
   E int found: Schluesselwort gefunden: 1, sonst 0
   E FILE *f: Dateizeiger der Protokolldatei
 */
/* Ergebnis:
   0 bei Fehler, sonst 1
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   07/1996     RK         Rockflow-Anpassung
   12/1999     OK         Protokoll für RF-Shell

   letzte Aenderung: OK 05.05.2000 Bugfix
 */
/**************************************************************************/
int FctCurves(char* data, int found, FILE* f)
{
	int ok = 1;
	int pos = 0;
	int p = 0;
	StuetzStellen* stuetz;
	long anz;
	double d1, d2;
	int curve_counter = 0;
	int i;

#ifdef TESTFILES
	DisplayMsgLn("Eingabedatenbank, Schluesselwort CURVES");
#endif

	LineFeed(f);
	FilePrintString(f,
	                "; 9 Functions "
	                "---------------------------------------------------------"
	                "-");
	LineFeed(f);

	/* Erste Konstante 1-Zeitkurve einfuegen (Index 0) - muss vorhanden sein
	 * TODO */
	if ((found == 0) || (found == 1))
	{
		anz_kurven = 1;
		stuetz = (StuetzStellen*)Malloc(sizeof(StuetzStellen));
		stuetz[0].punkt = 1.0;
		stuetz[0].wert = 1.0;
		kurven = (Kurven*)Malloc(sizeof(Kurven));
		kurven[anz_kurven - 1].anz_stuetzstellen = 1;
		kurven[anz_kurven - 1].stuetzstellen = stuetz;
	}
	/* einlesen */
	if (!found)
	{
#ifdef EXT_RFD
		FilePrintString(f, "; Schluesselwort: #CURVES (z.B. Zeitkurven)");
		LineFeed(f);
#endif
	}
	else
	{
		/* CURVES gefunden */
		while (StrTestHash(&data[p], &pos) || (found == 2))
		{
			/* Pruef-Funktion*/
			if ((found == 2) && (curve_counter >= anz_kurven)) break;
			FilePrintString(f, "#CURVES");
			LineFeed(f);
#ifdef EXT_RFD
			FilePrintString(f,
			                "; Schluesselwort: #CURVES (z.B. fuer Zeitkurven)");
			LineFeed(f);
			FilePrintString(f,
			                "; Das Schluesselwort muss nicht vorhanden sein.");
			LineFeed(f);
			FilePrintString(f,
			                "; Es folgen beliebig viele Abschnitte, die "
			                "jeweils eine eigene Kurve");
			LineFeed(f);
			FilePrintString(f,
			                "; darstellen und jeweils mit dem Schluesselwort "
			                "#CURVES eingeleitet werden");
			LineFeed(f);
			FilePrintString(
			    f, "; muessen. Jede Dieser Kurven hat folgendes Aussehen:");
			LineFeed(f);
			FilePrintString(f,
			                ";   Eine Kurve setzt sich aus Stuetzstellen und "
			                "Werten an den Stuetzstellen");
			LineFeed(f);
			FilePrintString(f, ";   zusammen:");
			LineFeed(f);
			FilePrintString(
			    f, ";   Es folgen beliebig viele Bloecke mit je 2 Werten.");
			LineFeed(f);
			FilePrintString(f, ";   - Stuetzstelle [double]");
			LineFeed(f);
			FilePrintString(f,
			                ";   - Funktionswert an der Stuetzstelle [double]");
			LineFeed(f);
#endif
			if (found == 1)
			{
				stuetz = NULL;
				anz = 0l;
				while (StrTestDouble(&data[p += pos]))
				{
					ok = (StrReadDouble(&d1, &data[p], f, &pos) && ok);
					ok = (StrReadDouble(&d2, & data[p += pos], f, &pos) && ok);
					LineFeed(f);
					anz++;
					stuetz = (StuetzStellen*)Realloc(
					    stuetz, (anz * sizeof(StuetzStellen)));
					stuetz[anz - 1].punkt = d1;
					stuetz[anz - 1].wert = d2;
				}
				if (anz >= 1l)
				{
					/* gueltige Zeitkurve, d.h. mind. 1 gueltige Stuetzstelle */
					anz_kurven++;
					kurven =
					    (Kurven*)Realloc(kurven, (anz_kurven * sizeof(Kurven)));
					kurven[anz_kurven - 1].anz_stuetzstellen = anz;
					kurven[anz_kurven - 1].stuetzstellen = stuetz;
				}
				else
				{
					FilePrintString(
					    f,
					    "* vorhergehende Zeitkurve unzulaessig, ignoriert !");
					LineFeed(f);
					stuetz = (StuetzStellen*)Free(stuetz);
					/* stuetz = Free(stuetz); MFC */
				}
			}
			else if (found == 2)
				if (curve_counter > 0) /* Dummy-Kurve nicht ausgeben */
				{
					stuetz = kurven[curve_counter].stuetzstellen;
					for (i = 0; i < kurven[curve_counter].anz_stuetzstellen;
					     i++)
					{
						fprintf(f, " %e ", stuetz[i].punkt);
						fprintf(f, " %e ", stuetz[i].wert);
						LineFeed(f);
					}
					stuetz = NULL;
				}
			curve_counter++;
		}
	}
	return ok;
}

///**************************************************************************/
///* ROCKFLOW - Funktion: FctReadHeterogeneousFields
// */
///* Aufgabe:
//   Liest zu jedem Knoten einen Wert der Permeabilität ein.
//   Identifikation über Koordinaten
// */
///* Ergebnis:
//    0 bei Fehler, sonst 1
// */
///* Programmaenderungen:
//    09/2003     SB  First Version
//   01/2004     SB  Verallgemeinert auf beliebig viele Parameter
//    06/2005     MB  msh / layer
//    08/2005     MB $NUM_TYPE NEW
// */
///**************************************************************************/
// int FctReadHeterogeneousFields(char* name_file, CMediumProperties* m_mat_mp)
//{
//	int ok = 0, method;
//	double* convertfact, * values, * defaultvalues, ** invals;
//	long i, j, no_values, nof, ihet;
//	int* help;
//	char zeile[MAX_ZEILE], outname[80];
//	string line;
//	int interpolation_option = 1;
//	ifstream ein;
//	ofstream out, out1;
//	std::stringstream in;
//	long NumberOfElements;
//	CFEMesh* m_msh = NULL;
//	int layer = 1;
//	int material_properties_index = -1;
//	int EleStart = -1;
//	int EleEnd = -1;
//	long NumberOfElementsPerLayer = -1;
//	MeshLib::CElem* m_ele = NULL;
//	//------------------------------------------------------------------------
//	DisplayMsgLn("Input file Heterogeneous Fields ");
//	//------------------------------------------------------------------------
//	// File handling
//	ein.open(name_file);
//	if(!ein)
//	{
//		//DisplayMsgLn(" ERROR opening file with heterogeneous fields!");
//		//DisplayMsgLn(" File does not exist.");
//		cout << " FctReadHeterogeneousFields" << endl;
//		cout << " Cannot find " << name_file << endl;
//		exit(1);
//	}
//	//------------------------------------------------------------------------
//	// Read MSH data
//	string line_string;
//	GetLineFromFile(zeile,&ein);
//	in.str((string ) zeile);
//	in >> line_string;
//	in.clear();
//	if(line_string.find("$MSH_TYPE") != string::npos)
//	{
//		GetLineFromFile(zeile,&ein);
//		in.str((string)zeile);
//		in >> line_string;
//		in.clear();
//		m_msh = MeshLib::FEMGet(line_string);
//		if(!m_msh)
//			cout << "FctReadHeterogeneousFields: no MSH data" << endl;
//	}
//	//------------------------------------------------------------------------
//	// Read Interpolation option
//	GetLineFromFile(zeile,&ein);
//	in.str((string ) zeile);
//	in >> line_string;
//	in.clear();
//	if(line_string.find("$INTERPOLATION") != string::npos)
//	{
//		GetLineFromFile(zeile,&ein);
//		in.str((string)zeile);
//		in >> line_string;
//		in.clear();
//		if(line_string.find("NEAREST_VALUE") != string::npos)
//			interpolation_option = 1;
//		if(line_string.find("GEOMETRIC_MEAN") != string::npos)
//			interpolation_option = 2;
//	}
//	//------------------------------------------------------------------------
//	GetLineFromFile(zeile,&ein);
//	in.str((string ) zeile);
//	in >> nof;
//	in.clear();
//	hf = Createhetfields(nof,name_file);
//	//------------------------------------------------------------------------
//	convertfact = (double*)Malloc(nof * sizeof(double));
//	values = (double*) Malloc(nof * sizeof(double));
//	defaultvalues = (double*) Malloc(nof * sizeof(double));
//	for(i = 0; i < nof; i++)
//	{
//		convertfact[i] = 0.0;
//		defaultvalues[i] = 0.0;
//	}
//	//------------------------------------------------------------------------
//	GetLineFromFile(zeile,&ein);
//	line = (string ) zeile;
//	in.str(line);
//	for(i = 0; i < nof; i++)
//	{
//		in >> outname;
//		set_hetfields_name(hf,i,outname);
//
//		//set material_properties_index
//		if(line.find("permeability") == 0)
//			material_properties_index = 0;
//		if(line.find("porosity") == 0)
//			material_properties_index = 1;
//	}
//	in.clear();
//	//------------------------------------------------------------------------
//	GetLineFromFile(zeile,&ein);
//	in.str((string ) zeile);
//	for(i = 0; i < nof; i++)
//		in >> convertfact[i];
//	in.clear();
//	//------------------------------------------------------------------------
//	// read default
//	GetLineFromFile(zeile,&ein);
//	in.str((string ) zeile);
//	for(i = 0; i < nof; i++)
//		in >> defaultvalues[i];
//	in.clear();
//	// for(i=0;i<nof;i++)
//	//  defaultvalues[i] *= convertfact[i]; MB->SB finde ich eher verwirrend ?
//	//------------------------------------------------------------------------
//	GetLineFromFile(zeile,&ein);
//	in.str((string ) zeile);
//	in >> no_values >> method;
//	in.clear();
//	//------------------------------------------------------------------------
//	NumberOfElements = (long)m_msh->ele_vector.size();
//	//------------------------------------------------------------------------
//	NumberOfElementsPerLayer = NumberOfElements /
// m_msh->getNumberOfMeshLayers();
//
//	//layers
//	if(m_mat_mp->geo_type_name.compare("LAYER") == 0)
//	{
//		char* temp = strdup(m_mat_mp->geo_name.c_str());
//		layer = atoi(temp);
//		EleStart = (layer - 1) * NumberOfElementsPerLayer;
//		EleEnd = layer * NumberOfElementsPerLayer;
//	}
//	//complete mesh
//	if(m_mat_mp->geo_type_name.compare("DOMAIN") == 0)
//	{
//		layer = 1;
//		EleStart = 0;
//		EleEnd = NumberOfElementsPerLayer;
//	}
//	//Warning
//	if(no_values < NumberOfElementsPerLayer)
//		DisplayMsgLn(
//		        "Warning! Fewer element values in File for heterogeneous
// permeability field than elements in element list");
//	//------------------------------------------------------------------------
//	/* field (int) for helping sort */
//	help = (int*) Malloc(NumberOfElements * sizeof(int));
//	for(i = 0; i < NumberOfElements; i++)
//		help[i] = 0;
//
//	/* initialize element values in element list; this is for the case, if not
// for all
//	   elements values are given in the input file */
//
//	//WW double test1;
//	//WW double test2;
//	//WW double test;
//	for(i = EleStart; i < EleEnd; i++)
//	{
//		m_ele = m_msh->ele_vector[i];
//		if (m_ele->mat_vector.Size() == 0)
//			m_ele->mat_vector.resize(2);
//		m_ele->mat_vector(material_properties_index) = defaultvalues[0];
//		// test1 = m_ele->mat_vector(0);
//	}
//	m_ele = m_msh->ele_vector[0];
//	//WW test1 = m_ele->mat_vector(0);
//	//WW test2 = m_ele->mat_vector(1);
//	//------------------------------------------------------------------------
//	//METHOD = 0:  read in unsorted values for coordinates and distribute to
// corresponding elements */
//	if(method == 0)
//	{
//		// allocate storage to read in file with het values and initialize
//		invals = (double**) Malloc((no_values) * sizeof(double*));
//		for(i = 0; i < no_values; i++)
//			invals[i] = (double*) Malloc((3 + nof + 1) * sizeof(double));
//		// initialize
//		for(i = 0; i < no_values; i++)
//			for(j = 0; j < (3 + nof + 1); j++) //+1 wegen GetAverageHetVal
//				invals[i][j] = 0.0;
//
//		//------------------------------------------------------------------------
//		// read values
//		for(i = 0; i < no_values; i++)
//		{
//			GetLineFromFile(zeile,&ein);
//			in.str((string) zeile);
//			//		in >> x >> y >> z ;
//			in >> invals[i][0] >> invals[i][1] >> invals[i][2];
//			for(j = 0; j < nof; j++)
//				in >> invals[i][j + 3];
//			in.clear();
//
//			// convert values by convertfact
//			for(j = 0; j < nof; j++)
//				invals[i][j + 3] = invals[i][j + 3] * convertfact[j];
//		}                         // end for read values
//
//		//------------------------------------------------------------------------
//		// element loop
//		for(i = EleStart; i < EleEnd; i++)
//		{
//			//.....................................................................
//			//Get the values that are nearest to element mid point
//			if(interpolation_option == 1)
//			{
//				ihet = GetNearestHetVal(i, m_msh, no_values, invals);
//				if(ihet < 0)
//					DisplayMsgLn(" Error getting nearest het_value location");
//				else
//					for(j = 0; j < nof; j++)
//						values[j] = invals[ihet][j + 3];
//				//DisplayMsg(" Het Val for element: "); DisplayLong(i);
//DisplayMsg("
// with coordinates ");
//				//DisplayDouble(x,0,0); DisplayMsg(", "); DisplayDouble(y,0,0);
// DisplayMsg(", "); DisplayDouble(z,0,0); DisplayMsg("       found at: ");
//				//DisplayDouble(invals[ihet][0],0,0); DisplayMsg(", ");
// DisplayDouble(invals[ihet][1],0,0); DisplayMsg(", ");
// DisplayDouble(invals[ihet][2],0,0); DisplayMsgLn(". ");
//			}
//			//.....................................................................
//			//Get all values in Element and calculate the geometric mean
//			if(interpolation_option == 2)
//			{
//				values[0] = GetAverageHetVal(i, m_msh, no_values, invals);
//
//				DisplayMsgLn(" AverageHetVal ");
//				DisplayMsg(" Element ");
//				DisplayDouble(i,0,0);
//				DisplayMsg("  Value: ");
//				DisplayDouble(values[0],0,0);
//			}
//			// save values
//			m_ele = m_msh->ele_vector[i];
//			m_ele->mat_vector(material_properties_index) = values[0];
//		}                         //end for element loop
//		ein.close();
//		// free storage for input values
//		invals = (double**) Free(invals);
//		//------------------------------------------------------------------------
//		//------------------------------------------------------------------------
//		//  OUT   write out fields sorted by element number
//		//------------------------------------------------------------------------
//		// Header
//		sprintf(outname,"%s%i",name_file,1);
//		out.open(outname);
//		out << "$MSH_TYPE" << endl;
//		out << "  GROUNDWATER_FLOW" << endl; //ToDo as variable
//		//out << "$LAYER" << endl;
//		//out << "  " << layer << endl;
//		out << "$INTERPOLATION" << endl;
//		out << "  GEOMETRIC_MEAN" << endl; //ToDo as variable
//		/* Field name */
//		out << nof << endl;
//		for(i = 0; i < nof; i++)
//			out << get_hetfields_name(hf,i) << ' ';
//		out << endl;
//		/* conversion factor is one in this case */
//		for(i = 0; i < nof; i++)
//			out << 1.0 << ' ';
//		out << endl;
//		// default values
//		for(i = 0; i < nof; i++)
//			out << defaultvalues[i] << ' ';
//		out << endl;
//		//out << NumberOfElements << ' ' << 1 << endl;
//		out << NumberOfElementsPerLayer << ' ' << 1 << endl;
//
//		out.setf(ios::scientific);
//		out.precision(5);
//
//		//------------------------------------------------------------------------
//		// Element data
//		for(i = EleStart; i < EleEnd; i++)
//		{
//			m_ele = m_msh->ele_vector[i];
//			for(j = 0; j < nof; j++)
//				out << m_ele->mat_vector(material_properties_index) << ' ';
//			out << endl;
//		}
//		out.close();
//	}                                     /* end if (method == 0) */
//	//------------------------------------------------------------------------
//	//METHOD = 1:  read in one sorted column, index is element number no
// conversion, no sorting
//	if(method == 1)
//	{
//		for(i = EleStart; i < EleEnd; i++)
//		{
//			GetLineFromFile(zeile,&ein);
//			in.str((string)zeile);
//			for(j = 0; j < nof; j++)
//				in >> values[j];
//			in.clear();
//			m_ele = m_msh->ele_vector[i];
//			m_ele->mat_vector(material_properties_index) = values[0];
//			//WW  test = m_ele->mat_vector(material_properties_index);
//		}
//		ein.close();
//	}                                     /* end if (method == 1) */
//
//	/* data structures deallocate */
//	convertfact = (double*) Free(convertfact);
//	values = (double*) Free(values);
//	defaultvalues = (double*) Free(defaultvalues);
//	help = (int*) Free(help);
//
//	return ok;
//}

/**************************************************************************
   MSHLib-Method: GetAverageHetVal
   Task:
   Programing:
   06/2005 MB Implementation
**************************************************************************/
double GetAverageHetVal(long EleIndex,
                        CFEMesh* m_msh,
                        long no_values,
                        double** invals)
{
	long i, j, ihet;
	double average;
	double xp[3], yp[3];
	double value;
	double NumberOfValues;
	//  double InvNumberOfValues;
	CGLPoint* m_point = NULL;
	MeshLib::CElem* m_ele = NULL;
	j = 0;  // only for 1 value
	//-----------------------------------------------------------------------
	// Get element data
	m_ele = m_msh->ele_vector[EleIndex];
	for (j = 0; j < 3; j++)
	{
		double const* const pnt(m_ele->GetNode(j)->getData());
		xp[j] = pnt[0];
		yp[j] = pnt[1];
		// zp[j] = 0.0;
	}
	//-----------------------------------------------------------------------
	// Find data points in the element
	NumberOfValues = 0;
	// WW  InvNumberOfValues = 0;
	m_point = new CGLPoint;
	average = -1;
	value = 0;
	for (i = 0; i < no_values; i++)
		if (invals[i][4] != -999999.0)  // Data point not within an element yet
		{
			m_point->x = invals[i][0];
			m_point->y = invals[i][1];
			m_point->z = 0.0;
			//....................................................................
			// Calculate the product of values in element
			// CC 10/05
			if (m_point->IsInTriangleXYProjection(xp, yp))
			{
				value = value + invals[i][3];
				NumberOfValues++;
				invals[i][4] = -999999.0;  // used as marker
			}
		}
	// end for
	//........................................................................
	if (NumberOfValues ==
	    0)  // if no data points in element --> get neares value
	{
		ihet = GetNearestHetVal(EleIndex, m_msh, no_values, invals);
		if (ihet < 0)
			DisplayMsgLn(" Error getting nearest het_value location");
		else
			average = invals[ihet][j + 3];
	}
	//........................................................................
	else  // if data points in element --> Calculate arithmetic mean

		average = value / NumberOfValues;
	delete m_point;
	return average;
}

/**************************************************************************
   MSHLib-Method: GetAverageHetVal
   Task:
   Programing:
   0?/2004 SB Implementation
   09/2005 MB EleClass
**************************************************************************/
long GetNearestHetVal(long EleIndex,
                      CFEMesh* m_msh,
                      long no_values,
                      double** invals)
{
	long i, nextele;
	double ex, ey, ez, dist, dist1;  // WW, dist2;
	double x, y, z;
	MeshLib::CElem* m_ele = NULL;
	//----------------------------------------------------------------------
	// MB ToDo
	// EleIndex = -1;
	// m_msh = NULL;
	//----------------------------------------------------------------------
	x = 0.0;
	y = 0.0;
	z = 0.0;
	dist = 10000000.0;  // Startwert
	// WW  dist2 = 0.01;	    // Abstand zwischen eingelesenen Knoten und
	// Geometrieknoten-RF;
	// Achtung, doppelbelegung möglich bei kleinen Gitterabständen
	nextele = -1;
	// Get element data
	m_ele = m_msh->ele_vector[EleIndex];
	double const* center = m_ele->GetGravityCenter();
	x = center[0];
	y = center[1];
	z = center[2];
	// Calculate distances
	for (i = 0; i < no_values; i++)
	{
		ex = invals[i][0];
		ey = invals[i][1];
		ez = invals[i][2];
		dist1 = (ex - x) * (ex - x) + (ey - y) * (ey - y) + (ez - z) * (ez - z);
		if (dist1 < dist)
		{
			dist = dist1;
			nextele = i;
		}
	}
	return nextele;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: Createhetfields
 */
/* Aufgabe:
   Erstellt die Datenstruktur zu den heterogenen Feldern
   und initialisiert diese                                   */
/* Ergebnis:
    -
 */
/* Programmaenderungen:
    01/2004     SB  First Version
 */
/**************************************************************************/
hetfields* Createhetfields(int n, char* name_file)
{
	int i;
	hetfields* hf;

	hf = (hetfields*)Malloc(sizeof(hetfields));
	hf->nof = n;
	hf->filename = name_file;
	hf->convertfact = (double*)Malloc(n * sizeof(double));
	hf->names = (char**)Malloc(n * sizeof(char*));
	for (i = 0; i < n; i++)
		hf->names[i] = (char*)Malloc(256 * sizeof(char));

	return hf;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: set_hetfields_name
 */
/* Aufgabe:
   set Funktion zur Struktur hetfields
   gesetzt wird ein Name aus dem Eingabefile
 */
/* Ergebnis:
    -
 */
/* Programmaenderungen:
    01/2004     SB  First Version
 */
/**************************************************************************/
void set_hetfields_name(hetfields* hf, int i, char* name)
{
	strcpy(hf->names[i], name);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: get_hetfields_name
 */
/* Programmaenderungen:
    01/2004     SB  First Version
 */
/**************************************************************************/
char* get_hetfields_name(hetfields* hf, int i)
{
	return hf->names[i];
}

/**************************************************************************/
/* ROCKFLOW - Funktion: get_hetfields_number
 */
/* Programmaenderungen:
    01/2004     SB  First Version
 */
/**************************************************************************/
int get_hetfields_number(hetfields* hf)
{
	return hf->nof;
}

/**************************************************************************/
/* ROCKFLOW - Funktion: GetHetValue
 */
/* Aufgabe:
   Liest zu einen übergebenen Namen den enstrechenden Wert aus
   der elementstruktur het_fields aus.

 */
/* Ergebnis:
    value = Ausgelesener Wert
 */
/* Programmaenderungen:
    01/2004     SB  First Version
 */
/**************************************************************************/
double GetHetValue(int ele_no, char* inname)
{
	int i, n;
	char* name;
	double value;
	int material_properties_index = 0;
	CFEMesh* m_msh = NULL;
	n = get_hetfields_number(hf);
	m_msh = MeshLib::FEMGet("GROUNDWATER_FLOW");
	for (i = 0; i < n; i++)
	{
		name = get_hetfields_name(hf, i);
		if (strstr(inname, name) && (strcmp(inname, name) == 0))
		{
			/* found */
			//  value = ELEGetHetFieldValue(ele_no,i);
			value = m_msh->ele_vector[ele_no]->mat_vector(
			    material_properties_index);
			return value;
		}
	}
	DisplayMsgLn(" Error - No het. Field Value found");
	return -1.0;
}

/**********************************************************************
   Function for interpolation between two points P1(x1|z1) and P2(x2|z2)
   The function returns the value zn at the position xn.
   function is used in GetMatrixValue();

   Programming:
   08/2008 NB
 ***********************************************************************/
double interpol(double x1, double x2, double zx1, double zx2, double xn)
{
	if (x1 == x2)
		return zx1;
	else
		return zx1 + (xn - x1) * (zx2 - zx1) / (x2 - x1);
}

/**********************************************************************
Function GetMatrixValue (double var1, double var2, int *gueltig)

This function reads a data matrix with two independent arguments and
function values in an FCT-File and returns the value corresponding to
var1 and var2. If var1 and var2 are not in the matrix, the function
returns a value interpolated between both arguments.

Programming:
11-2011 NB/TF
***********************************************************************/
double GetMatrixValue(double var1,
                      double var2,
                      std::string caption,
                      int* gueltig)
{
	CFunction* matrix;
	// WW int anz_variables, anz_data;
	int dim_x, dim_y;
	int i1 = 0;
	int i2 = 0;
	int j1 = 0;
	int j2 = 0;

	matrix = FCTGet(caption);
	dim_x = matrix->matrix_dimension[0];  // NB 4.8.01
	dim_y = matrix->matrix_dimension[1];  // NB

	if (var1 < *matrix->variable_data_vector[0])  // is var1 smaller then the
	                                              // smallest argument?
	{
		*gueltig = 0;
		i1 = i2 = 0;
	}
	else if (var1 > *matrix->variable_data_vector[dim_x - 1])  // is var1 larger
	                                                           // then largest
	                                                           // argument?
	{
		*gueltig = 0;
		i1 = i2 = dim_x - 1;
	}
	else
	{
		i1 = searchElement(var1, 0, dim_x - 1, matrix->variable_data_vector);
		i2 = i1 + 1;
	}

	if (var2 < *matrix->variable_data_vector[dim_x])  // is var1 smaller then
	                                                  // the smallest argument?
	{
		*gueltig = 0;
		j1 = j2 = dim_x;
	}
	else if (var2 >
	         *matrix->variable_data_vector[dim_y + dim_x - 1])  // is var1
	                                                            // larger then
	                                                            // largest
	                                                            // argument?
	{
		*gueltig = 0;
		j1 = j2 = dim_y + dim_x - 1;
	}
	else
	{
		j1 = searchElement(
		    var2, dim_x, dim_y + dim_x - 1, matrix->variable_data_vector);
		j2 = j1 + 1;
	}

	if (fabs(var1 - *matrix->variable_data_vector[i1]) <
	    std::numeric_limits<double>::epsilon())  // var 1 is in the matrix
	{
		if (fabs(var2 - *matrix->variable_data_vector[j1]) <
		    std::numeric_limits<double>::epsilon())  // var 2 is in the matrix
		{
			return *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                     (j1 - dim_x) * dim_x];
		}
		else  // only v1 is in the matrix
		{
			double zx1y1, zx1y2;
			double y1, y2;
			zx1y1 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j1 - dim_x) * dim_x];
			zx1y2 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j2 - dim_x) * dim_x];
			y1 = *matrix->variable_data_vector[j1];
			y2 = *matrix->variable_data_vector[j2];

			return interpol(y1, y2, zx1y1, zx1y2, var2);
		}
	}
	else  // v1 is not in the matrix
	{
		if (fabs(var2 - *matrix->variable_data_vector[dim_x + j1]) <
		    std::numeric_limits<double>::epsilon())  // only var 2 is in the
		                                             // matrix
		{
			double zx1y1, zx2y1;
			double x1, x2;
			zx1y1 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j1 - dim_x) * dim_x];
			zx2y1 = *matrix->variable_data_vector[dim_x + dim_y + i2 +
			                                      (j1 - dim_x) * dim_x];
			x1 = *matrix->variable_data_vector[i1];
			x2 = *matrix->variable_data_vector[i2];

			return interpol(x1, x2, zx1y1, zx2y1, var1);
		}

		else  // neither var1 nor var2 are in the matrix
		{
			double interp1, interp2;
			double zx1y1, zx2y1, zx1y2, zx2y2;
			double x1, x2, y1, y2;
			zx1y1 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j1 - dim_x) * dim_x];
			zx2y1 = *matrix->variable_data_vector[dim_x + dim_y + i2 +
			                                      (j1 - dim_x) * dim_x];
			zx1y2 = *matrix->variable_data_vector[dim_x + dim_y + i1 +
			                                      (j2 - dim_x) * dim_x];
			zx2y2 = *matrix->variable_data_vector[dim_x + dim_y + i2 +
			                                      (j2 - dim_x) * dim_x];

			x1 = *matrix->variable_data_vector[i1];
			x2 = *matrix->variable_data_vector[i2];
			y1 = *matrix->variable_data_vector[j1];
			y2 = *matrix->variable_data_vector[j2];

			interp1 = interpol(x1, x2, zx1y1, zx2y1, var1);
			interp2 = interpol(x1, x2, zx1y2, zx2y2, var1);

			return interpol(y1, y2, interp1, interp2, var2);
		}
	}
}

/**********************************************************************
   Function GetMatrixValue (double var1, double var2, int *gueltig)

   This function reads a data matrix with two independent arguments and
   function values in an FCT-File and returns the value corresponding to
   var1 and var2. If var1 and var2 are not in the matrix, the function
   returns a value interpolated between both arguments.

   Programming:
   08/2008 NB
 ***********************************************************************/
// double GetMatrixValue(double var1, double var2, std::string caption, int*
// gueltig)
//{
//	CFunction* matrix;
//	//WW int anz_variables, anz_data;
//	int dim_x, dim_y;
//	int i1 = 0;
//	int i2 = 0;
//	int j1 = 0;
//	int j2 = 0;
//	int counter;
//	double x1 = 0.0,x2 = 0.0,y1 = 0.0,y2 = 0.0; //OK411
//	double zx1y1,zx2y1,zx1y2,zx2y2;
//
//	matrix = FCTGet(caption);
//	//WW anz_variables = (int)matrix->variable_names_vector.size();
//	//dim_x = matrix->matrix_dimension_x;
//	//dim_y = matrix->matrix_dimension_y;
//	dim_x = matrix->matrix_dimension[0];  //NB 4.8.01
//	dim_y = matrix->matrix_dimension[1];  //NB
//	//WW anz_data = (int)matrix->variable_data_vector.size()-dim_x-dim_y;
//	//----------------------------------------------------------------------
//	if (var1 < *matrix->variable_data_vector[0]) //is var1 smaller then the
// smallest argument?
//	{
//		x1 = x2 = *matrix->variable_data_vector[0];
//		*gueltig = 0;
//		i1 = i2 = 0;
//	}
//	else
//	//is var1 larger then largest argument?
//	if (var1 > *matrix->variable_data_vector[dim_x - 1])
//	{
//		x1 = x2 = *matrix->variable_data_vector[dim_x - 1];
//		*gueltig = 0;
//		i1 = i2 = dim_x - 1;
//	}
//	else
//		for (counter = 0; counter < dim_x; counter++)
//		{
//			//does var1 fit an argument in the matrix exactly?
//			if (var1 == *matrix->variable_data_vector[counter])
//			{
//				x1 = x2 = *matrix->variable_data_vector[counter];
//				i1 = i2 = counter;
//				break;
//			}
//			else
//			//var1 is between two arguments in the matrix
//			if (var1 < *matrix->variable_data_vector[counter])
//			{
//				x1 = *matrix->variable_data_vector[counter - 1];
//				x2 = *matrix->variable_data_vector[counter];
//				i2 = counter;
//				i1 = i2 - 1;
//				break;
//			}
//		}
//	//same procedure for var2:
//	if (var2 < *matrix->variable_data_vector[dim_x])
//	{
//		y1 = y2 = *matrix->variable_data_vector[dim_x];
//		*gueltig = 0;
//		j1 = j2 = dim_x;
//	}
//	else
//	if (var2 > *matrix->variable_data_vector[dim_x + dim_y - 1])
//	{
//		y1 = y2 = *matrix->variable_data_vector[dim_x + dim_y - 1];
//		*gueltig = 0;
//		j1 = j2 = dim_x + dim_y - 1;
//	}
//	else
//		for (counter = dim_x; counter < dim_x + dim_y; counter++)
//		{
//			if (var2 == *matrix->variable_data_vector[counter])
//			{
//				y1 = y2 = *matrix->variable_data_vector[counter];
//				j1 = j2 = counter;
//				break;
//			}
//			else
//			if (var2 < *matrix->variable_data_vector[counter])
//			{
//				y1 = *matrix->variable_data_vector[counter - 1];
//				y2 = *matrix->variable_data_vector[counter];
//				j2 = counter;
//				j1 = j2 - 1;
//				break;
//			}
//		}
//	//getting the corresponding Z values for the arguments from the data vector
//	zx1y1 = *matrix->variable_data_vector[(j1 - dim_x) * dim_x + (i1 + dim_x +
// dim_y)];
//	zx2y1 = *matrix->variable_data_vector[(j1 - dim_x) * dim_x + (i2 + dim_x +
// dim_y)];
//	zx1y2 = *matrix->variable_data_vector[(j2 - dim_x) * dim_x + (i1 + dim_x +
// dim_y)];
//	zx2y2 = *matrix->variable_data_vector[(j2 - dim_x) * dim_x + (i2 + dim_x +
// dim_y)];
//	return interpol (y1,y2,
//	                 interpol (x1,x2,zx1y1,zx2y1,  var1),interpol (x1,
//	                                                               y1,
//	                                                               zx1y2,
//	                                                               zx2y2,
//	                                                               var1),var2);
//}

/****************************************************************************
 * Finds and returns the positive minimum of a vector.
 * Programming: NB Dec 08
 *****************************************************************************/
double FindMin(vector<double> Vec)
{
	double x = DBL_MAX;
	int unsigned i;

	for (i = 0; i < Vec.size(); i++)
		if ((Vec[i] >= 0) && (Vec[i] < x)) x = Vec[i];
	return x;
}

/****************************************************************************
 * Finds and returns the maximum of a vector.
 * Programming: NB Jan 08
 *****************************************************************************/
double FindMax(vector<double> Vec)
{
	double x = DBL_MIN;
	int unsigned i;

	for (i = 0; i < Vec.size(); i++)
		if (Vec[i] > x) x = Vec[i];
	return x;
}

/****************************************************************************
 * Finds all real roots of a third grade polynomial in the form:
 * P(x) = x^3 + px^2 + qx + r
 * roots are returned in a vector
 *
 * Programming: NB, Dec08
 *****************************************************************************/
void NsPol3(double p, double q, double r, vector<double>* roots)
{
	double eps = 7E-15;
	double a, b, h, phi, D, z[3];
	double pi = 3.1415926535897;
	double nz;
	int i;

	b = (p / 3) * (p / 3);
	a = q / 3 - b;
	b = b * p / 3 + 0.5 * (r - p / 3 * q);
	h = sqrt(fabs(a));

	if (b < 0) h = -h;

	D = MathLib::fastpow(a, 3) + b * b;

	if (D <= (-eps))
	{
		nz = 3;
		phi = acos(b / MathLib::fastpow(h, 3)) / 3;
		z[0] = 2 * h * cos(pi / 3 - phi) - p / 3;
		z[1] = 2 * h * cos(pi / 3 + phi) - p / 3;
		z[2] = -2 * h * cos(phi) - p / 3;
	}
	else if (D < eps)
	{
		nz = 3;
		z[0] = -2 * h - p / 3;
		z[1] = h - p / 3;
		z[2] = z[1];
	}
	else
	{
		nz = 1;
		if (a >= eps)
		{
			b = b / MathLib::fastpow(h, 3);
			phi = log(b + sqrt(b * b + 1)) / 3;
			z[0] = -2 * h * sinh(phi) - p / 3;
		}
		else if (a > (-eps))
		{
			z[0] = pow((2 * abs(b)), 1. / 3.);
			if (b > 0) z[0] = -z[0];
			z[0] = z[0] - p / 3;
		}
		else
		{
			b = b / MathLib::fastpow(h, 3);
			phi = log(b + sqrt(b * b - 1)) / 3;
			z[0] = -2 * h * cosh(phi) - p / 3;
		}
	}

	for (i = 0; i < nz; i++)
		roots->push_back(z[i]);
}

double getNodalValueFromElementValue(
    CRFProcess& pcs,
    const std::map<long, double>& map_eleId_val,
    EleToNodeInterpolationMethod::type interpolation_type,
    long node_id)
{
	MeshLib::CFEMesh* m_msh = pcs.m_msh;
	MeshLib::CNode* m_node = m_msh->nod_vector[node_id];
	double nodal_val = 0.0;

	switch (interpolation_type)
	{
		case EleToNodeInterpolationMethod::VOLUME_WEIGHTED:
		{
			const size_t n_conn_ele = m_node->getConnectedElementIDs().size();
			double sum_weights = 0.0;

			// sum the distance weighted data from each connected block
			for (size_t j = 0; j < n_conn_ele; j++)
			{
				const size_t e_id = m_node->getConnectedElementIDs()[j];
				MeshLib::CElem* e = m_msh->ele_vector[e_id];

				double ele_v = .0;
				std::map<long, double>::const_iterator itr =
				    map_eleId_val.find(e_id);

				if (itr != map_eleId_val.end())
				{
					ele_v = itr->second;
				}
				else
				{
					// skip elements which are not in the given list
					continue;
				}

				// calculate representative volume of the considered node in
				// each connected element for weighting
				double volume =
				    e->GetVolume() / ((double)e->GetNodesNumber(false));
				// weight
				double weight = 1.0 / volume;
				sum_weights += weight;

				nodal_val += ele_v * weight;
			}

			nodal_val /= sum_weights;
		}
		break;
		case EleToNodeInterpolationMethod::GAUSS_EXTRAPOLATION:
		{
			const size_t n_conn_ele = m_node->getConnectedElementIDs().size();
			if (n_conn_ele == 0) return 0;

			MeshLib::CElem* e0 =
			    m_msh->ele_vector[m_node->getConnectedElementIDs()[0]];
			if (e0->GetElementType() != MshElemType::QUAD &&
			    e0->GetElementType() != MshElemType::HEXAHEDRON)
				return 0;

			if (e0->GetElementType() == MshElemType::QUAD)
			{
				size_t mid_nod_id = node_id;
				if (n_conn_ele < 4)
				{
					std::map<long, int> set_mid_node;
					for (size_t i = 0; i < e0->GetNodesNumber(false); i++)
					{
						long tmp_node_id = e0->GetNodeIndex(i);
						if (tmp_node_id == node_id) continue;
						if (m_msh->nod_vector[tmp_node_id]
						        ->getConnectedElementIDs()
						        .size() < 4)
							continue;
						if (set_mid_node.count(tmp_node_id) == 0)
							set_mid_node[tmp_node_id] = 0;
						else
							set_mid_node[tmp_node_id]++;
					}

					for (std::map<long, int>::iterator itr =
					         set_mid_node.begin();
					     itr != set_mid_node.end();
					     ++itr)
					{
						if (itr->second > 1)
						{
							mid_nod_id = itr->first;
							break;
						}
					}
				}
			}

			std::vector<double> ele_values(n_conn_ele);

			for (size_t j = 0; j < n_conn_ele; j++)
			{
				const size_t e_id = m_node->getConnectedElementIDs()[j];
				MeshLib::CElem* e = m_msh->ele_vector[e_id];

				if (e->GetElementType() != MshElemType::QUAD &&
				    e->GetElementType() != MshElemType::HEXAHEDRON)
					break;

				double ele_v = .0;
				std::map<long, double>::const_iterator itr =
				    map_eleId_val.find(e_id);

				if (itr != map_eleId_val.end())
				{
					ele_v = itr->second;
				}

				ele_values[j] = ele_v;
			}
		}
		break;
	}

	return nodal_val;
}

void convertElementDataToNodalData(
    CRFProcess& pcs,
    const std::vector<double>& vec_ele_data,
    EleToNodeInterpolationMethod::type interpolation_type,
    std::vector<double>& vec_nod_values)
{
	MeshLib::CFEMesh* m_msh = pcs.m_msh;
	const size_t n_nodes = m_msh->nod_vector.size();
	vec_nod_values.resize(n_nodes);

	if (interpolation_type == EleToNodeInterpolationMethod::VOLUME_WEIGHTED)
	{
		for (size_t i = 0; i < n_nodes; i++)
		{
			MeshLib::CNode* m_node = m_msh->nod_vector[i];
			const size_t n_conn_ele = m_node->getConnectedElementIDs().size();

			double nodal_val = 0.0;
			double sum_weights = 0.0;

			// sum the distance weighted data from each connected block
			for (size_t j = 0; j < n_conn_ele; j++)
			{
				const size_t e_id = m_node->getConnectedElementIDs()[j];
				MeshLib::CElem* e = m_msh->ele_vector[e_id];

				double weight = 0.0;

				// calculate representative volume of the considered node in
				// each connected element for weighting
				double volume =
				    e->GetVolume() / ((double)e->GetNodesNumber(false));

				// weight
				weight = 1.0 / volume;
				sum_weights += weight;

				nodal_val += vec_ele_data[e_id] * weight;
			}

			nodal_val /= sum_weights;
			vec_nod_values[i] = nodal_val;
		}
	}
}

#ifdef NEW_EQS
Math_Group::SparseTable* createSparseTable(MeshLib::CFEMesh* a_mesh,
						 bool quadratic,
						 bool symmetry,
						 Math_Group::StorageType storage_type)
{
	Math_Group::SparseTable* st = new Math_Group::SparseTable();
	st->symmetry = symmetry;

	long i = 0, j = 0, ii = 0, jj = 0;
	long lbuff0 = 0, lbuff1 = 0;
	long** larraybuffer;
	larraybuffer = NULL;
	//
	// In sparse table, = number of nodes
	st->rows = a_mesh->GetNodesNumber(quadratic);
	st->size_entry_column = 0;
	st->diag_entry = new long[st->rows];

	if (storage_type == Math_Group::JDS)
	{
		st->row_index_mapping_n2o = new long[st->rows];
		st->row_index_mapping_o2n = new long[st->rows];
	}
	else if (storage_type == Math_Group::CRS)
	{
		st->row_index_mapping_n2o = NULL;
		st->row_index_mapping_o2n = NULL;
	}

	if (symmetry)
	{
		larraybuffer = new long* [st->rows];
		for (i = 0; i < st->rows; i++)
		{
			if (storage_type == Math_Group::JDS) st->row_index_mapping_n2o[i] = i;
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			lbuff1 = (long)a_mesh->nod_vector[i]->getConnectedNodes().size();
			larraybuffer[i] = new long[lbuff1 + 1];
			//
			larraybuffer[i][0] = lbuff1;
			for (j = 0; j < lbuff1; j++)
				larraybuffer[i][j + 1] =
					a_mesh->nod_vector[i]->getConnectedNodes()[j];
			a_mesh->nod_vector[i]->getConnectedNodes().clear();
			for (j = 0; j < lbuff1; j++)
			{
				jj = larraybuffer[i][j + 1];
				if (i <= jj)
					a_mesh->nod_vector[i]->getConnectedNodes().push_back(jj);
			}
		}
	}

	/// CRS storage
	if (storage_type == Math_Group::CRS)
	{
		/// num_column_entries saves vector ptr of CRS
		st->num_column_entries = new long[st->rows + 1];

		std::vector<long> A_index;
		long col_index;

		for (i = 0; i < st->rows; i++)
		{
			st->num_column_entries[i] = (long)A_index.size();

			for (j = 0;
				 j < (long)a_mesh->nod_vector[i]->getConnectedNodes().size();
				 j++)
			{
				col_index = a_mesh->nod_vector[i]->getConnectedNodes()[j];

				/// If linear element is used
				if ((!quadratic) && (col_index >= st->rows)) continue;

				if (i == col_index)st-> diag_entry[i] = (long)A_index.size();
				A_index.push_back(col_index);
			}
		}

		st->size_entry_column = (long)A_index.size();
		st->num_column_entries[st->rows] = st->size_entry_column;

		st->entry_column = new long[st->size_entry_column];
		for (i = 0; i < st->size_entry_column; i++)
			st->entry_column[i] = A_index[i];
	}
	else if (storage_type == Math_Group::JDS)
	{
		//
		//--- Sort, from that has maximum connect nodes to that has minimum
		// connect nodes
		//
		for (i = 0; i < st->rows; i++)
		{
			st->row_index_mapping_n2o[i] = i;
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			st->diag_entry[i] =
				(long)a_mesh->nod_vector[i]->getConnectedNodes().size();
			if (!quadratic)
			{
				lbuff0 = 0;
				for (j = 0; j < st->diag_entry[i]; j++)
					if (a_mesh->nod_vector[i]->getConnectedNodes()[j] <
						static_cast<size_t>(st->rows))
						lbuff0++;
				st->diag_entry[i] = lbuff0;
			}
			st->size_entry_column += st->diag_entry[i];
		}

		//
		for (i = 0; i < st->rows; i++)
		{
			// 'diag_entry' used as a temporary array
			// to store the number of nodes connected to this node
			lbuff0 = st->diag_entry[i];  // Nodes to this row
			lbuff1 = st->row_index_mapping_n2o[i];
			j = i;
			while ((j > 0) && (st->diag_entry[j - 1] < lbuff0))
			{
				st->diag_entry[j] = st->diag_entry[j - 1];
				st->row_index_mapping_n2o[j] = st->row_index_mapping_n2o[j - 1];
				j = j - 1;
			}
			st->diag_entry[j] = lbuff0;
			st->row_index_mapping_n2o[j] = lbuff1;
		}
		// Old index to new one
		for (i = 0; i < st->rows; i++)
			st->row_index_mapping_o2n[st->row_index_mapping_n2o[i]] = i;
		// Maximum number of columns in the sparse table
		st->max_columns = st->diag_entry[0];
		//--- End of sorting
		//
		//--- Create sparse table
		//
		st->num_column_entries = new long[st->max_columns];
		st->entry_column = new long[st->size_entry_column];
		// 1. Count entries in each column in sparse table
		for (i = 0; i < st->max_columns; i++)
			st->num_column_entries[i] = 0;
		for (i = 0; i < st->rows; i++)
			// 'diag_entry' still is used as a temporary array
			// it stores that numbers of nodes connect to this nodes
			for (j = 0; j < st->diag_entry[i]; j++)
				st->num_column_entries[j]++;

		// 2. Fill the sparse table, i.e. store all its entries to
		//    entry_column
		lbuff0 = 0;

		for (i = 0; i < st->max_columns; i++)
			for (j = 0; j < st->num_column_entries[i]; j++)
			{
				// ii is the real row index of this entry in matrix
				ii = st->row_index_mapping_n2o[j];
				// jj is the real column index of this entry in matrix
				jj = a_mesh->nod_vector[ii]->getConnectedNodes()[i];
				st->entry_column[lbuff0] = jj;

				// Till to this stage, 'diag_entry' is really used to store
				// indices of the diagonal entries.
				// Hereby, 'index' refers to the index in entry_column array.
				if (ii == jj) st->diag_entry[ii] = lbuff0;
				//
				lbuff0++;
			}
	}

	// For the case of symmetry matrix
	if (symmetry)
	{
		for (i = 0; i < st->rows; i++)
		{
			lbuff0 = larraybuffer[i][0];
			a_mesh->nod_vector[i]->getConnectedNodes().resize(lbuff0);
			//
			for (j = 0; j < lbuff0; j++)
				a_mesh->nod_vector[i]->getConnectedNodes()[j] =
					larraybuffer[i][j + 1];
		}
		for (i = 0; i < st->rows; i++)
		{
			delete[] larraybuffer[i];
			larraybuffer[i] = 0;
		}
		delete[] larraybuffer;
		larraybuffer = 0;
	}

	return st;
}

void CreateSparseTable(MeshLib::CFEMesh* msh, Math_Group::SparseTable* &sparse_graph, Math_Group::SparseTable* &sparse_graph_H)
{
	Math_Group::StorageType stype;
	stype = Math_Group::JDS;
	for (int i = 0; i < (int)num_vector.size(); i++)
		if (num_vector[i]->ls_storage_method == 100)
		{
			stype = Math_Group::CRS;
			break;
		}

	// Symmetry case is skipped.
	// 1. Sparse_graph_H for high order interpolation. Up to now, deformation
	if (msh->GetNodesNumber(false) != msh->GetNodesNumber(true))
		sparse_graph_H = createSparseTable(msh, true, false, stype);
	// 2. M coupled with other processes with linear element
	if (sparse_graph_H)
	{
		if ((int)pcs_vector.size() > 1)
			sparse_graph = createSparseTable(msh, false, false, stype);
	}
	// 3. For process with linear elements
	else
		sparse_graph = createSparseTable(msh, false, false, stype);

	//  sparse_graph->Write();
	//  sparse_graph_H->Write();
	//
	// ofstream Dum("sparse.txt", ios::out);
	// sparse_graph_H->Write(Dum);
}
#endif  //#ifndef NON_PROCESS

