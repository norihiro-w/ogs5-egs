/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef makros_INC
#define makros_INC

#include <cfloat>
#include <limits>
#include <string>

#include "BuildInfo.h"
#include "Configure.h"

/* Laufzeitausgaben */
#define TESTTIME

/**********************************************************************/
/* Speicher */

//#ifndef NO_ERROR_CONTROL /* Wird ggf. im Makefile gesetzt */
//#define ERROR_CONTROL
/* Fehlertests (Feldgrenzen, Existenz o.ae.), die bei sauberen Netzen und
   einwandfrei funktionierendem Programm nichts bringen und nur Laufzeit
   kosten */
//#endif

#define MEMORY_MANAGEMENT_NOT_ANSI_COMPLIANT
/*  Bei einigen Compilern werden
     malloc(0)
     realloc(NULL,xxx)
     realloc(xxx,0)
     free(NULL)
    nicht ANSI-gerecht gehandhabt. Mit diesem Schalter wird
    ANSI-Verhalten gewaehrleistet. */

#define noMEMORY_ALLOCATION_TEST_SUCCESS
/*  Prueft, ob eine Speicheranforderung erfolgreich absolviert wurde */

#define noMEMORY_TEST_IN_TIME /* fuer Versions-Speichertest */
/* Erstellt waehrend der Laufzeit eine Bilanz des allockierten
   und wieder freigegebenen Speichers. Sehr Zeitintensiv !!! */

#define noMEMORY_STR /* fuer Versions-Speichertest */
                     /* Gibt Informationen bei Memory-Funktionen zur
                        Aufrufstellenlokalisation. Funktioniert nur zusammen mit
                        MEMORY_TEST_IN_TIME. Sehr Speicherintensiv!!! */

#define noMEMORY_SHOW_USAGE
/* Gibt bei MEMORY_TEST_IN_TIME Informationen ueber jeden
   Malloc/Realloc/Free-Vorgang aus. */

#define noMEMORY_REALLOC
/* Ersetzt Realloc durch Malloc und Free und speichert um */

#ifndef MEMORY_TEST_IN_TIME
#ifdef MEMORY_SHOW_USAGE
#undef MEMORY_SHOW_USAGE
#endif
#ifndef MEMORY_ALLOCATION_TEST_SUCCESS
#ifdef MEMORY_STR
#undef MEMORY_STR
#endif
#endif
#endif

#ifdef MEMORY_STR
#define Malloc(a) MAlloc(a, __FILE__, __LINE__)
#define Free(a) FRee(a, __FILE__, __LINE__)
#define Realloc(a, b) REalloc(a, b, __FILE__, __LINE__)
/* Ersetzt Malloc, Free und Realloc in allen *.c-Dateien durch den
   erweiterten Aufruf mit Dateiname und Zeilennummer */
#endif

/**********************************************************************/
/* Protokolle */

/* Definitionen der Dateinamen-Erweiterungen */
#define TEXT_EXTENSION ".rfd"
/* Dateinamen-Erweiterung fuer RF-Output-Dateien */
#define RF_MESSAGE_EXTENSION ".msg"
/* Dateinamen-Erweiterung fuer Text-Eingabedatei (Inverses Modellieren) */
#define CHEM_REACTION_EXTENSION ".pqc"
#define TEC_FILE_EXTENSION ".tec"
#define VTK_FILE_EXTENSION ".vtk"
#define CSV_FILE_EXTENSION ".csv"


/* Eingabeprotokoll ausfuehrlich kommentieren */
#define EXT_RFD_MIN
/* Eingabeprotokoll kommentieren, nur gefundene Schluesselworte */
#ifdef EXT_RFD
#undef EXT_RFD_MIN
#endif

/* Format der Double-Ausgabe ueber FilePrintDouble --> txtinout */
#define FORMAT_DOUBLE
#define FPD_GESAMT 4
#define FPD_NACHKOMMA 14

/**********************************************************************/
#define VEKNORM_BICG MVekNorm2


/**********************************************************************/
/* PCS / C++ */
#define PCS_NUMBER_MAX 30
#define DOF_NUMBER_MAX 6    // JT: max # dof's per process
#define MAX_FLUID_PHASES 2  // JT: max # fluid phases


#ifndef NULL
#define NULL ((void*)0)
#endif
#ifndef TRUE
#define TRUE (0 == 0)
#endif
#ifndef FALSE
#define FALSE (1 == 0)
#endif
#ifndef PI
#define PI 3.14159265358979323846
#endif

/* Feste Zahlen fuer Genauigkeitspruefungen etc. */
#define Mdrittel (1.0 / 3.0)
#define MKleinsteZahl DBL_EPSILON
#define MKleinsteZahlen MKleinsteZahl
#define MFastNull DBL_MIN
#define MSqrt2Over3 sqrt(2.0 / 3.0)


#define MAX_ZEILE 2048
#define MAX_ZEILEN MAX_ZEILE

extern std::string FileName;
extern std::string FilePath;


//---- MPI Parallel --------------
#if defined(USE_MPI) || defined(USE_PETSC)
extern int mysize;
extern int myrank;
#endif
//---- MPI Parallel --------------

extern void ogsAbort(int i);

#endif
