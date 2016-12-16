/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "makros.h"

#ifdef USE_PETSC
#include <mpi.h>
#endif

std::string FileName;
std::string FilePath;
#if defined(USE_MPI) || defined(USE_MPI_PARPROC) ||      \
    defined(USE_MPI_REGSOIL) || defined(USE_MPI_GEMS) || \
    defined(USE_MPI_BRNS) || defined(USE_MPI_KRC) || defined(USE_PETSC)
int mysize;
int myrank;
#endif

void ogsAbort(int i)
{
#ifdef USE_PETSC
    MPI_Abort(MPI_COMM_WORLD, i);
#else
    exit(i);
#endif
}
