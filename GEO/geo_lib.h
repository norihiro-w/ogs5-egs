/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   GeoLib - Object: Control Level
   Task:
   Programing:
   07/2003 TK/OK Implementation
   11/2003 WW Geo_Read class
   11/2003 OK -> GeoLib class
   09/2005 CC GeoLib2
**************************************************************************/

#ifndef geolib_INC
#define geolib_INC

enum GEO_TYPE
{
	GS_POINT,
	GS_POLYLINE,
	GS_SURFACE,
	GS_VOLUME
};

#include <string>

class GeoLib
{
public:
	GeoLib() {}
	~GeoLib() {}

private:
	friend class CGLPolyline;
	friend class Surface;
};

extern void GEOLIB_Read_GeoLib(const std::string& file_name_path_base);
extern void GEOLIB_Clear_GeoLib_Data();

#define GLI_FILE_EXTENSION ".gli"

#endif
