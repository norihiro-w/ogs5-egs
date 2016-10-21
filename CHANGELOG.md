## v5.3.3-egs.1.9
Changes in features
- changed GCC optimization flag from -O3 to -O2 -march=native 
- support Lis configuration via command line arguments 
- print current time for log
- support limiting memory usage via MAXMEM_GB environmental variable (only on Linux)
- support CSV output
- use csvdiff in benchmarking

Changes in development
- now reuires c++11
- support ccache
- fix builds with MSVC 
- add LIS_DIR cmake option
- use trusty in travis
- use Appveyor
- run benchmarks on travis

Bug fix
- fix a parsing error if a comment comes after the number of deactivated subdomains in PCS
- abort if a linear solver fail
- print elasped time instead of CPU time
- bugfix in TECPLOT output for polylines
- bugfix in face integration for applying ST

## v5.3.3-egs.1.8
Changes in features
- New porosity model 13 dependent on volumetric strain and unjacketed volumetric strain (Blöcher 2013) by Antoine Jacquey
- New porosity model 14 dependent on volumetric strain and unjacketed volumetric strain and temperature (case 18 + thermal dependency) by Antoine Jacquey
- New Youngs modulus model 3 for drained model by Antoine Jacquey
- New viscosity model 22 from Ramey et al (1974)
- New viscosity model 30 having exponential dependency (from Fabien Magri)
- New viscosity model 31 from Fabien Magri
- New adaptive time stepping method PID (still experimental)
- Support arbitrary expression in FUNCTION using exprtk (http://www.partow.net/programming/exprtk/)
- Support DIRECT in IC with PETc
- Support PETSc v3.5
- Reactivate Ctr-C break (from Marc)

Changes in the codes
- License headers added
- Removed unused codes for current modeling , e.g. GUI, GEM, PQC

Bug fix
- fixed a bug in TECPLOT output for POINT with PETSc. Don't output for ghost nodes.
- fixed errors in deformation with fracture elements
- fixed a bug in setting 3rd BC with PETSc

## 2014.11.07
Changes
- support DEFORMATION processes with PETSc 
- check line ending of infput files (i.e. Windows or Unix)
- support deactivated domains with Lis and PARDISO
- detect a partitioned mesh file including the partition number

Bug fix
- fix PVTU output that alias of node value names were not used
- fix bugs around quadratic prism elements

## 2014.09.26
Changes
- add -Wno-unused-local-typedefs
- remove  -Woverloaded-virtual

Bug fix
- get PETSc version in cmake

## 2014.04.30
Changes
- support TH_MONOLITHIC process with PETSc nonlinear solve and field split
- delete unused directores: Qt, UTL
- delete unused processes: PTC_FLOW

Bug fix
- avoid having very small time step size due to crtiical time
- correct convergence rate in picard


## 2014.03.11
Changes
- suport BENCHMARK_DIR and BENCHMARK_REF_FIR options in cmake

## 2014.03.10
Changes
- activate compression of LEQS with Lis
- support longlong with MKL

Bug fix
- fix possible memory access error in CElem::MarkingAll

## 2014.02.26
Changes
- merge EoS and stress dependent permeability from trunk
- support $TIME_INTERPOLATION also in BC
- support STEADY in LIQUID_FLOW
- merges works from Fabien's project
  - add $CALC_DIFF_FROM_STRESS0,  $RESET_STRAIN in PCS
  - support INITIAL and ELEMENT distribution types
  - add $STORAGE_DISTRIBUTION
  - improve speed to construct JDS matrix
  - add PETREL data type to output

Bug fix
- fix a bug that STEPS setting in OUT files does not work when multiple #OUTPUT are given
- fix a bug that biot was not multiplied in a strain coupling term

## 2014.01.16
Changes
- support PETSc for liquid flow and heat transport
- support PVTU output 
- change NGPoints for TET from 15 to 5. this saves memory quite lot
- output 13 digits of point coordinates in VTK with PETSc
- support $EXTERNAL_SOLVER_OPTION for PETsc
- add CSV data output type
- add $SOLID_BULK_MODULUS in MSP
- support 3rd BC type:  $TRANSFER_COEFFICIENT
- support $CONSTRAIN in BC
- support coupling iteration based adaptive time stepping
- add time step number to solution file names
- add keyword $LEQS_OUTPUT to PCS
- fix a bug with restart. Initial last_active_time shoud be start_time
- use (non-) symmetric weighted matching for Pardiso
- add OGS_ONLY_TH cmake option to optimize memory for TH simulations
- add new Density models 20,21 (p-T-c dependent) and Viscocity models 20, 21 from FEFLOW
- add MFP output to PVD
- add the new keyword $MAT_ID to BC files

