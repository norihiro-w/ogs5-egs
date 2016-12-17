/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: display.h
 */
/* Aufgabe:
   Enthaelt alle Funktionen fuer Standard Ein- und Ausgabe (Bildschirm,
   Tastatur)
 */
/**************************************************************************/

#ifndef display_INC
#define display_INC

/* Schutz gegen mehrfaches Einfuegen */

/* Andere oeffentlich benutzte Module */
#include <cstdio>
#include <cstring>
//#include <ctype.h>

/*JT: Send output message*/
extern void ScreenMessage(const char* message, ...);
extern void ScreenMessage2(const char* message, ...);

//#ifndef LOG_DEBUG
//#define LOG_DEBUG 0
//#endif
extern int ogs_log_level;
#define ScreenMessaged(fmt, ...)                              \
	do                                                        \
	{                                                         \
		if (ogs_log_level) ScreenMessage(fmt, ##__VA_ARGS__); \
	} while (0)
#define ScreenMessage2d(fmt, ...)                              \
	do                                                         \
	{                                                          \
		if (ogs_log_level) ScreenMessage2(fmt, ##__VA_ARGS__); \
	} while (0)

extern void DisplayMsg(const char* s);
/* Schreibt Zeichenkette ohne Zeilenvorschub auf Standardausgabe */
extern void DisplayMsgLn(const char* s);
/* Schreibt Zeichenkette mit Zeilenvorschub auf Standardausgabe */
/* Schreibt Zeichenkette mit Zeilenruecklauf auf Standardausgabe */
extern void DisplayDouble(double x, int i, int j);
/* Schreibt Double-Wert ohne Zeilenvorschub auf Standardausgabe */
extern void DisplayLong(long x);
/* Schreibt Long-Wert ohne Zeilenvorschub auf Standardausgabe */
extern void DisplayErrorMsg(const char* s);
/* Schreibt Fehlermeldung auf Standardausgabe */
#endif
