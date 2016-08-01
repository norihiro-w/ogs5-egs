/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CRFPROCESSTH_H_
#define CRFPROCESSTH_H_

#include <vector>

#include "rf_pcs.h"

#ifdef USE_PETSC
#include "petscsnes.h"
#endif

class CRFProcessTH : public CRFProcess
{
public:
	CRFProcessTH();
	virtual ~CRFProcessTH();

	void Initialization();
	virtual double Execute(int loop_process_number);
#ifdef USE_PETSC
	virtual void setSolver(petsc_group::PETScLinearSolver* petsc_solver);
	void copyVecToNodalValues(Vec x);
	void copyNodalValuesToVec(Vec x);
#endif

protected:
#ifdef USE_PETSC
	double ExecuteNonlinearWithPETsc();
#endif
	void UpdateIterativeStep(const double damp);
private:
	double error_k0;
	std::vector<int> vec_pos;
};

#ifdef USE_PETSC
PetscErrorCode FormFunctionTH(SNES snes, Vec x, Vec f, void* ctx);
#if (PETSC_VERSION_NUMBER >= 3050)
PetscErrorCode FormJacobianTH(SNES snes, Vec x, Mat jac, Mat B, void* ctx);
#else
PetscErrorCode FormJacobianTH(SNES snes, Vec x, Mat* jac, Mat* B,
                              MatStructure* flag, void* ctx);
#endif
#endif

#endif
