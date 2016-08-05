/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

// GeoLib.cpp : Definiert den Einsprungpunkt f�r die Konsolenanwendung.
//
/*-------------------------------------------------------------------------
   GeoLib
   Programming:
   09/2003 OK/CC encapsulate Read function
   last modified:
   09/2005 CC GeoLib2
   -------------------------------------------------------------------------*/

#include "geo_lib.h"

#include <stdio.h>
#include <string.h>
// C++
#include <iostream>
#include <string>

// GL
#include "geo_ply.h"
#include "geo_pnt.h"
#include "geo_sfc.h"
#include "geo_vol.h"

using namespace std;

extern void remove_white_space(string* buffer);

enum GEO_TYPES
{
	GEO_POINT = 'P',
	GEO_LINE = 'L',
	GEO_POLYLINE = 'O',
	GEO_SURFACE = 'S',
	GEO_VOLUME = 'V',
	GEO_DOMAIN = 'D'
};

/**************************************************************************
   GeoLib- Funktion: GEOLIB_Read_GeoLib
   Aufgabe: Ansteuerung der GeoLib-Lesefunktionen
   Programmaenderungen:
   09/2003 TK Erste Version
   11/2003 OK Read polylines
   11/2003 WW Read surfaces
   09/2004 TK Path info // please check
   10/2004 OK path_name_slash
   09/2005 CC delete line lesen function
**************************************************************************/
void GEOLIB_Read_GeoLib(const std::string& file_name_path_base)
{
	// Points
	GEORemoveAllPoints();
	GEOReadPoints(file_name_path_base);
	// Polylines
	GEORemoveAllPolylines();
	GEOReadPolylines(file_name_path_base);
	// Surfaces
	GEORemoveAllSurfaces();
	GEOReadSurfaces(file_name_path_base);
	// Volumes
	GEORemoveAllVolumes();
	GEOReadVolumes(file_name_path_base);
	// Determine dependencies between GEO objects
	GEOSurfaceTopology();
	GEOCreateSurfacePointVector();  // OK
}

/*************************************************************************
   GeoLib- Funktion: GEOLIB_Clear_GeoLib_Data

   Aufgabe: L�schen aller GeoLib-Datenstrukturen

   Programmaenderungen:
   09/2003     TK        Erste Version
   01/2005 OK GEORemovePolylines, GEORemoveAllSurfaces, GEORemoveVolumes
 **************************************************************************/
void GEOLIB_Clear_GeoLib_Data()
{
	GEORemoveAllPoints();
	Clear_LineVector();
	GEORemoveAllPolylines();  // OK41
	GEORemoveAllSurfaces();   // OK41 CC change
	GEORemoveAllVolumes();    // OK41
}

