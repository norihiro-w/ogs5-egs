/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************/
/* ROCKFLOW - Modul: display.c
 */
/* Aufgabe:
   Enthaelt alle Funktionen fuer Standard Ein- und Ausgabe (Bildschirm,
   Tastatur)
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   10/1999     AH         Warnung entfernt
   01/2002     MK         Umleitung der DisplayX-Funktionen in MSG-Datei
                          Ausnahmen: DisplayStartMsg/DisplayEndMsg */
/**************************************************************************/
#include "display.h"

#include <cstdarg>
#include <cstdlib>
#include <cstdio>
#include "Configure.h"
#include "makros.h"

#ifdef NDEBUG
int ogs_log_level = 0;
#else
int ogs_log_level = 1;
#endif


/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayMsg
 */
/* Aufgabe:
   Schreibt Zeichenkette ohne Zeilenvorschub auf Standardausgabe
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayMsg(const char* s)
{
	printf("%s", s);
}

#ifdef MSVC
int vasprintf(char** strp, const char* format, va_list ap)
{
	va_list args = ap;

	int len = _vscprintf(format, args) +
	          1;  // _vscprintf doesn't count terminating '\0'
	va_end(args);

	if ((*strp = (char*)malloc(len * sizeof(char))) == NULL) return (-1);

	len = vsprintf(*strp, format, ap);  // C4996
	// Note: vsprintf is deprecated; consider using vsprintf_s instead

	return len;
}
#endif

/**************************************************************************
 Task: Output message to screen. Helps to remove so many IFDEFS
 Programming:
  03/2012 JT
**************************************************************************/
void ScreenMessage(const char* format, ...)
{
#if defined(USE_MPI) || defined(USE_PETSC)
	if (myrank > 0) return;
#endif
	va_list arguments;
	va_start(arguments, format);
	char* allocatedBuffer;
	vasprintf(&allocatedBuffer, format, arguments);
	va_end(arguments);
	printf("%s", allocatedBuffer);
	free(allocatedBuffer);
	fflush(stdout);
}

void ScreenMessage2(const char* format, ...)
{
	va_list arguments;
	va_start(arguments, format);
	char* allocatedBuffer = NULL;
	vasprintf(&allocatedBuffer, format, arguments);
	va_end(arguments);
#if defined(USE_MPI) || defined(USE_PETSC)
	printf("rank%d: %s", myrank, allocatedBuffer);
#else
	printf("%s", allocatedBuffer);
#endif
	free(allocatedBuffer);
	fflush(stdout);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayMsgLn
 */
/* Aufgabe:
   Schreibt Zeichenkette mit Zeilenvorschub auf Standardausgabe,
   beginnt immer erst nach 12 Zeichen Einrueckung.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayMsgLn(const char* s)
{
	printf("%s\n            ", s);
}


/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayDouble
 */
/* Aufgabe:
   Schreibt double-Wert ohne Zeilenvorschub formatiert auf
   Standardausgabe. Wird fuer beide Formatangaben 0 angegeben,
   wird im Standardformat geschrieben.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E double x: double-Wert
   E int i: Gesamtstellenzahl
   E int j: Nachkommastellen
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   12/1994     MSR        Erste Version
   12/1995     cb         E-Format
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayDouble(double x, int i, int j)
{
	if ((i == 0) && (j == 0)) /* printf("%f",x); */
		printf("%g", x);
	else
		printf("% *.*g", i, j, x);
}

/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayLong
 */
/* Aufgabe:
   Schreibt long-Wert ohne Zeilenvorschub im Standardformat auf
   Standardausgabe.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E long x: long-Wert
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   12/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayLong(long x)
{
	printf("%ld", x);
}


/**************************************************************************/
/* ROCKFLOW - Funktion: DisplayErrorMsg
 */
/* Aufgabe:
   Schreibt Fehlermeldung auf Standardausgabe.
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E char *s: Zeichenkette (Fehlermeldung)
 */
/* Ergebnis:
   - void -
 */
/* Programmaenderungen:
   03/1994     MSR        Erste Version
   01/2002     MK         Umleitung in MSG-Datei
 */
/**************************************************************************/
void DisplayErrorMsg(const char* s)
{
	printf("\n!!!!!!!!  %s\n\n            ", s);
}
