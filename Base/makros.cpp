/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "makros.h"

std::string FileName;
std::string FilePath;
#if defined(USE_MPI) || defined(USE_PETSC)
int mysize;
int myrank;
#endif
