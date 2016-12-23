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

	// Assemble system equation
	void GlobalAssembly();
	void GlobalAssembly_DM();

	// overloaded
	double Execute(int loop_process_number);

	// Aux. Memory
	double const* GetLastTimeStepSolution() const { return lastTimeStepSolution.data(); }
	double const* GetInitialFluidPressure() const { return p0.data(); }

	void ScalingNodeForce(const double SFactor);
	void InitGauss();
	//
	void SetInitialGuess_EQS_VEC();
    void incrementNodalDUFromSolution();
    void incrementNodalPressureFromSolution();
    void incrementNodalDisplacement();
    double NormOfUpdatedNewton();
    void StoreLastTimeStepDisplacements();
    void StoreLastCouplingIterationSolution();
    void RecoverLastTimeStepDisplacements();
    void CopyLastTimeStepDisplacementToCurrent();

    void zeroNodalDU();
	double NormOfDisp();
#if !defined(USE_PETSC) && \
    !defined(NEW_EQS)  // && defined(other parallel libs)//03~04.3012. WW
	                   //#ifndef NEW_EQS
	double NormOfUnkonwn_orRHS(bool isUnknowns = true);
#endif
    double getNormOfCouplingError(int pvar_id_start, int n);
    // Stress
	// For partitioned HM coupled scheme
    void ResetStress();
	void ResetTimeStep();
	//
	void UpdateStress();
	void UpdateInitialStress(bool ZeroInitialS);
	void Extropolation_GaussValue();

	// Calculate scaling factor for load increment
	double CaclMaxiumLoadRatio();

	// Write stresses
	std::string GetGaussPointStressFileName();
	void WriteGaussPointStress();
	void ReadGaussPointStress();
	void ReadElementStress();

	// Access members
	CFiniteElementVec* GetFEM_Assembler() const { return fem_dm; }

private:
    void solveLinear();
    void solveNewton();
    void setDUFromSolution();
    void setPressureFromSolution();

private:
	CFiniteElementVec* fem_dm = nullptr;
	void InitialMBuffer();
	std::vector<double> lastTimeStepSolution;
	std::vector<double> lastCouplingSolution;
	std::vector<double> p0;
	bool _isInitialStressNonZero = false;

	int counter = 0;
	double InitialNormR0 = 0;
	double InitialNormDU_coupling = 0;
	double InitialNormDU0 = 0;

	InitDataReadWriteType idata_type = none;

	//
	double norm_du0_pre_cpl_itr = 0;

	double getNormOfDisplacements();
};
}  // end namespace

extern double Tolerance_global_Newton;
extern double Tolerance_Local_Newton;
extern int number_of_load_steps;
extern int problem_dimension_dm;
#endif
