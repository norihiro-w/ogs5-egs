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
	std::vector<double> const& GetAuxArray() const { return tempArray; }
	std::vector<double> const& GetInitialFluidPressure() const { return p0; }

	void ScalingNodeForce(const double SFactor);
	void InitGauss();
	//
	void SetInitialGuess_EQS_VEC();
	void UpdateIterativeStep(const double damp, const int u_type);
	void InitializeNewtonSteps(const bool ini_excav = false);
	double NormOfUpdatedNewton();
	void StoreLastSolution(const int ty = 0);
	void RecoverSolution(const int ty = 0);
	double NormOfDisp();
	// Stress
	// For partitioned HM coupled scheme
	void ResetCouplingStep();
	void ResetTimeStep();
	//
	void UpdateStress();
	void UpdateInitialStress(bool ZeroInitialS);
	void Extropolation_GaussValue();

	// Write stresses
	std::string GetGaussPointStressFileName();
	void WriteGaussPointStress();
	void ReadGaussPointStress();
	void ReadElementStress();

	// Access members
	CFiniteElementVec* GetFEM_Assembler() const { return fem_dm; }

private:
	void InitialMBuffer();

	CFiniteElementVec* fem_dm = nullptr;
	std::vector<double> tempArray;
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

extern void CalStressInvariants(const long Node_Inex, double* StressInv);
// For visualization
extern void CalMaxiumStressInvariants(double* StressInv);
extern double LoadFactor;
extern double Tolerance_global_Newton;
extern double Tolerance_Local_Newton;
extern int number_of_load_steps;
extern int problem_dimension_dm;
#endif
