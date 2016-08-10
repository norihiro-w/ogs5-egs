/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#include <cstdlib>
#include <ctime>
#ifdef SUPERCOMPUTER
// kg44 test for buffered output
#include <cstdio>
#include <unistd.h>
#endif

#ifndef WIN32
#include <sys/time.h>
#include <sys/resource.h>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif
#if defined(USE_MPI)
#include <mpi.h>
#endif

#ifdef LIS
#include <lis.h>
#endif
#ifdef USE_PETSC
#include <petscksp.h>
#include <petsctime.h>
#endif

#ifdef USE_PARALUTION
#include <paralution.hpp>
#endif

#include "Configure.h"
#include "BuildInfo.h"

#include "break.h"
#include "display.h"
#include "FileTools.h"
#include "makros.h"
#include "memory.h"
#include "MemWatch.h"
#include "RunTime.h"
#include "timer.h"

#include "files0.h"
#include "problem.h"

#include "ogs_display.h"

#define TEST

#if defined(USE_MPI)
double elapsed_time_mpi;
#endif

#ifndef WIN32
void setmemlimit()
{
	char* env_value = getenv("MAXMEM_GB");
	if(env_value!=NULL)
	{
		long bytes = atol(env_value)*(1024*1024*1024);
		if (bytes > 0)
		{
#ifdef USE_PETSC
			bytes /= mysize;
#endif
			ScreenMessage("-> limit memory size to %g GB\n", atol(env_value));
			struct rlimit memlimit;
			memlimit.rlim_cur = bytes;
			memlimit.rlim_max = bytes;
			setrlimit(RLIMIT_AS, &memlimit);
		}
	}
}
#else
void setmemlimit() {}
#endif

/* Definitionen */

/**************************************************************************/
/* ROCKFLOW - Funktion: main
 */
/* Aufgabe:
   Hauptprogramm
 */
/* Formalparameter: (E: Eingabe; R: Rueckgabe; X: Beides)
   E int argc: Anzahl der Kommandozeilenparameter (incl. Programmname)
   E char *argv[]: Zeiger auf Feld der argc Kommandozeilenparameter
 */
/* Ergebnis:
   Fehlerfreie Bearbeitung: Exit-Code 0
 */
/* Programmaenderungen:
   07/1996     MSR        Erste Version
   08/1999     OK         RF-FEM Applikation
 */
/**************************************************************************/
int main(int argc, char* argv[])
{
	/* parse command line arguments */
	std::string anArg;
	std::string modelRoot;
	bool terminate = false;
	for (int i = 1; i < argc; i++)
	{
		anArg = std::string(argv[i]);
		if (anArg == "--help" || anArg == "-h")
		{
			std::cout << "Usage: ogs [MODEL_ROOT] [OPTIONS]\n"
			          << "Where OPTIONS are:\n"
			          << "  -h [--help]       print this message and exit\n"
			          << "  -b [--build-info] print build info and exit\n"
			          << "  -v [--verbose]    print debug info\n"
			          << "  --version         print ogs version and exit"
			          << "\n";
			terminate = true;
			continue;
		}
		if (anArg == "--build-info" || anArg == "-b")
		{
			std::cout << "ogs version: " << OGS_VERSION << "\n"
			          << "ogs date: " << OGS_DATE << "\n";
#ifdef CMAKE_CMD_ARGS
			std::cout << "cmake command line arguments: " << CMAKE_CMD_ARGS
			          << "\n";
#endif  // CMAKE_CMD_ARGS
#ifdef GIT_COMMIT_INFO
			std::cout << "git commit info: " << GIT_COMMIT_INFO << "\n";
#endif  // GIT_COMMIT_INFO
#ifdef SVN_REVISION
			std::cout << "subversion info: " << SVN_REVISION << "\n";
#endif  // SVN_REVISION
#ifdef BUILD_TIMESTAMP
			std::cout << "build timestamp: " << BUILD_TIMESTAMP << "\n";
#endif  // BUILD_TIMESTAMP
			terminate = true;
			continue;
		}
		if (anArg == "--verbose" || anArg == "-v")
		{
			ogs_log_level = 1;
			continue;
		}
		if (anArg == "--version")
		{
			std::cout << OGS_VERSION << "\n";
			terminate = true;
			continue;
		}
		if (anArg == "--model-root" || anArg == "-m")
		{
			modelRoot = std::string(argv[++i]);
			continue;
		}
		// anything left over must be the model root, unless already found
		if (modelRoot == "") modelRoot = std::string(argv[i]);
	}  // end of parse argc loop

	if (terminate)
	{             // non-interactive mode and no model given
		exit(0);  // e.g. just wanted the build info
	}

	char* dateiname(NULL);
#ifdef SUPERCOMPUTER
	// *********************************************************************
	// buffered output ... important for performance on cray
	// (unbuffered output is limited to 10 bytes per second)
	// georg.kosakowski@psi.ch 11.10.2007

	char buf[1024 * 1024];
	int bsize;

	bsize = 1024 * 1024;  // question: what happens if buffer is full?
	// according to documentation the buffer is flushed when full.
	// If we have a lot of output, increasing buffer is usefull.
	if (bsize > 0)
		//        bufstd = malloc(bsize);
		setvbuf(stdout, buf, _IOFBF, bsize);
//**********************************************************************
#endif
/*---------- MPI Initialization ----------------------------------*/
#if defined(USE_MPI)
	printf("Before MPI_Init\n");
	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);      // 12.09.2007 WW
	elapsed_time_mpi = -MPI_Wtime();  // 12.09.2007 WW
	MPI_Comm_size(MPI_COMM_WORLD, &mysize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	std::cout << "After MPI_Init myrank = " << myrank << '\n';
	time_ele_paral = 0.0;
#endif
/*---------- MPI Initialization ----------------------------------*/

#ifdef USE_PETSC
	PetscLogDouble v1, v2;
	char help[] = "OGS with PETSc \n";
	// PetscInitialize(argc, argv, help);
	PetscInitialize(&argc, &argv, (char*)0, help);
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER >= 3040)
	PetscTime(&v1);
#else
	PetscGetTime(&v1);
#endif
	MPI_Comm_rank(PETSC_COMM_WORLD, &myrank);
	MPI_Comm_size(PETSC_COMM_WORLD, &mysize);
	ScreenMessage("===\nUse PETSc solver (MPI processes = %d)\n", mysize);
#endif

/*---------- LIS solver -----------------------------------------*/
#ifdef LIS
	lis_initialize(&argc, &argv);
#endif
#ifdef USE_PARALUTION
	ScreenMessage("--- PARALUTION INFO BEGIN ----------------------------------\n");
	paralution::init_paralution();
	paralution::info_paralution();
	ScreenMessage("--- PARALUTION INFO END ------------------------------------\n");
#endif

	setmemlimit();

/*========================================================================*/
/* Kommunikation mit Betriebssystem */
/* Ctrl-C ausschalten */
// NoBreak();
/* Timer fuer Gesamtzeit starten */
#ifdef TESTTIME
	TStartTimer(0);
	BaseLib::RunTime runTime;
	runTime.start();
#endif
/* Intro ausgeben */
#if defined(USE_MPI) || defined(USE_PETSC)
	if (myrank == 0)
#endif
		DisplayStartMsg();
	/* Speicherverwaltung initialisieren */
	if (!InitMemoryTest())
	{
		DisplayErrorMsg(
		    "Fehler: Speicherprotokoll kann nicht erstellt werden!");
		DisplayErrorMsg("        Programm vorzeitig beendet!");
		return 1;  // LB changed from 0 to 1 because 0 is indicating success
	}
	if (argc == 1)  // interactive mode
		dateiname = ReadString();
	else  // non-interactive mode
	{
		if (argc == 2)  // a model root was supplied
		{
			dateiname = (char*)Malloc((int)strlen(argv[1]) + 1);
			dateiname = strcpy(dateiname, argv[1]);
		}
		else  // several args supplied
		    if (modelRoot != "")
		{
			dateiname = (char*)Malloc((int)modelRoot.size() + 1);
			dateiname = strcpy(dateiname, modelRoot.c_str());
		}
#if defined(USE_MPI) || defined(USE_PETSC)
		if (myrank == 0)
#endif
			DisplayMsgLn(dateiname);
	}
	// WW  DisplayMsgLn("");
	// WW  DisplayMsgLn("");
	// ----------23.02.2009. WW-----------------

	// LB Check if file exists
	std::string tmpFilename = dateiname;
	tmpFilename.append(".pcs");
	if (!IsFileExisting(tmpFilename))
	{
		std::cout << " Error: Cannot find file " << dateiname << "\n";
		return 1;
	}

	ScreenMessage("\n---------------------------------------------\n");
	ScreenMessage("ogs version    : %s\n", OGS_VERSION);
	ScreenMessage("ogs date       : %s\n", OGS_DATE);
#ifdef CMAKE_CMD_ARGS
	ScreenMessage("cmake command line arguments: %s\n", CMAKE_CMD_ARGS);
#endif  // CMAKE_CMD_ARGS
#ifdef GIT_COMMIT_INFO
	ScreenMessage("git commit info: %s\n", GIT_COMMIT_INFO);
#endif  // GIT_COMMIT_INFO
#ifdef SVN_REVISION
	ScreenMessage("subversion info: %s\n", SVN_REVISION);
#endif  // SVN_REVISION
#ifdef BUILD_TIMESTAMP
	ScreenMessage("build timestamp: %s\n", BUILD_TIMESTAMP);
#endif  // BUILD_TIMESTAMP

	time_t tm =time(NULL );
	struct tm * curtime = localtime ( &tm );
	ScreenMessage("current time   : %s", asctime(curtime));

#ifdef USE_PETSC
	MPI_Barrier(PETSC_COMM_WORLD);
#endif

	FileName = dateiname;
	size_t indexChWin, indexChLinux;
	indexChWin = indexChLinux = 0;
	indexChWin = FileName.find_last_of('\\');
	indexChLinux = FileName.find_last_of('/');
	//
	if (indexChWin != std::string::npos)
		FilePath = FileName.substr(0, indexChWin) + "\\";
	else if (indexChLinux != std::string::npos)
		FilePath = FileName.substr(0, indexChLinux) + "/";
	// ---------------------------WW
	Problem* aproblem = new Problem(dateiname);
#ifndef WIN32
	BaseLib::MemWatch mem_watch;
	ScreenMessage("\tcurrent mem: %d MB\n",
	              mem_watch.getVirtMemUsage() / (1024 * 1024));
#endif
	aproblem->Euler_TimeDiscretize();
#ifdef USE_PETSC
	MPI_Barrier(PETSC_COMM_WORLD);
#endif
	delete aproblem;
	aproblem = NULL;
#ifdef TESTTIME
#if defined(USE_PETSC)
// kg44 quick fix to compile PETSC with version PETSCV3.4
#if (PETSC_VERSION_NUMBER >= 3040)
	PetscTime(&v2);
#else
	PetscGetTime(&v2);
#endif
	PetscPrintf(PETSC_COMM_WORLD, "Elapsed time: %ld s\n", v2 - v1);
#elif defined(USE_MPI)
	elapsed_time_mpi += MPI_Wtime();  // 12.09.2007 WW
	ScreenMessage2("Elapsed time    : %ld s\n", elapsed_time_mpi);
#else
	ScreenMessage2("CPU time    : %ld s\n", TGetTimer(0));
	ScreenMessage2("Elapsed time: %g s\n", runTime.elapsed());
#endif

#endif // TESTTIME
/* Abspann ausgeben */
/* Ctrl-C wieder normal */
// StandardBreak();

	// Finalize libraries
#ifdef LIS
	lis_finalize();
#endif
#if defined(USE_PETSC)
	PetscFinalize();
#elif defined(USE_MPI)
	MPI_Finalize();
#endif
#ifdef USE_PARALUTION
	paralution::stop_paralution();
#endif

	free(dateiname);


	return 0;
}
