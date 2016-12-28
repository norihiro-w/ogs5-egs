/**
 * \copyright
 * Copyright (c) 2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/**************************************************************************
   ROCKFLOW - Object: Process PCS_Deformation
   Task:
   Programing:
   07/2003 WW Implementation
**************************************************************************/
#ifndef pcs_dm_INC
#define pcs_dm_INC

#include <vector>

#include "rf_pcs.h"

namespace FiniteElement
{
class CFiniteElementVec;
}
using FiniteElement::CFiniteElementVec;

namespace process
{
enum InitDataReadWriteType
{
	none,
	read_write,
	read_all_binary,
	write_all_binary,
	read_all_asci,
	write_all_asci,
	read_stress_binary,
	write_stress_binary,
	read_displacement,
	write_displacement,
	read_pressure,
	write_pressure
};

// Elasto-plastic Deformation
class CRFProcessDeformation : public CRFProcess
{
public:
	CRFProcessDeformation();
	virtual ~CRFProcessDeformation();

	void Initialization();

	double Execute(int loop_process_number);
	void ResetTimeStep();
	void Extropolation_GaussValue();

	double const* GetInitialFluidPressure() const { return p0.data(); }

private:
	void InitialMBuffer();
	void InitGauss();

	void AssembleResidual();
	void AssembleJacobian();
	void GlobalAssembly();
	void GlobalAssembly_DM();

	void solveLinear();
	void solveNewton();

	void setDUFromSolution();
	void setPressureFromSolution();
	void incrementNodalDUFromSolution();
	void incrementNodalPressureFromSolution();
	void incrementNodalDisplacement();
	void zeroNodalDU();
	void updateGaussStressStrain();

	double getNormOfCouplingError(int pvar_id_start, int n);

	void StoreLastTimeStepDisplacements();
	void StoreLastCouplingIterationSolution();
	void RecoverLastTimeStepDisplacements();
	void ResetStress();
	void CopyLastTimeStepDisplacementToCurrent();

	std::string GetGaussPointStressFileName();
	void WriteGaussPointStress();
	void ReadGaussPointStress();

private:
	CFiniteElementVec* fem_dm = nullptr;
	std::vector<double> lastTimeStepSolution;
	std::vector<double> lastCouplingSolution;
	std::vector<double> p0;
	bool _isInitialStressNonZero = false;
	double norm_du0_pre_cpl_itr = 0;
	InitDataReadWriteType idata_type = none;
};
}  // end namespace

extern double Tolerance_Local_Newton;
extern int problem_dimension_dm;

#endif
