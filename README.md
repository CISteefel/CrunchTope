k<h1 align='center'>CrunchFlow</h1>

CrunchFlow is reactive transport software developed by Carl I. Steefel, Toshiyuki Bandai, and Sergi Molins at Lawrence Berkeley National Laboratory.

### Copyright and License Agreement

 CrunchFlow, Copyright (c) 2016, The Regents of the University of California, through Lawrence Berkeley National Laboratory
 subject to receipt of any required approvals from the U.S. Dept. of Energy.  All Rights Reserved.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

    (1) Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer

    (2) Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

    (3) Neither the name of the University of California, Lawrence Berkeley National Laboratory, U.S. Department of Energy nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

### CrunchFlow Web Site

[crunchflow.lbl.gov](https://crunch.lbl.gov/)

### News and Updates

Newly added are capabilities for 1D and 2D fully transient unsaturated flow based on the Richards equation (courtesy Toshi Bandai).  Check out Ex13UnsaturatedFlow, Ex14UnsaturatedFlow_evaporation, and Ex15UnsaturatedFlow_transient...

### Installation

Do not use petsc-3.22 with CrunchFlow, since the newly developed PETSc Fortran hooks there do not work.  Instead, use petsc-3.21 (3.21.6). The "git clone" command that follows should put the petsc distribution into ~/petsc (so PETSC_DIR=petsc).

    git clone https://gitlab.com/petsc/petsc.git --branch v3.21.6

Then to download CrunchTope for the first time (i.e., a clone) use

    git clone https://github.com/CISteefel/CrunchTope.git

After the CrunchTope Repo is installed, you can use the various GitHub commands to update (either git Blah_Blah from the Command line, or from GitHub Desktop on Windows)

For those using the pre-built executable (located in the Windows directory in the GitHub repo), you may need to install the Intel oneAPI redistributable libraries:

    https://registrationcenter-download.intel.com/akdlm/IRC_NAS/f6a44238-5cb6-4787-be83-2ef48bc70cba/w_ifort_runtime_p_2024.1.0.968.exe

### Compilation

To run on Windows, first download and install Visual Studio 2022.  Then download and install Intel oneAPI Base Toolkit, choosing the option to integrate with Visual Studio.  The installation script should provide the option to link with Visual Studio 2022, if it doesn't, then something is wrong with the Visual Studio installation.  Then download and install Intel oneAPI HPC toolkit.  Once installed, you should see the option to start a Intel oneAPI command prompt for Intel 64 for Visual Studio.  Starting this will set all of the Environment Variables for the use of oneAPI.  Then run from this same Command Window: 

    C:\cygwin64\bin\mintty.exe -

to start up a Cygwin window (Unix) that inherits the Environment Variables from the oneAPI Command Prompt.  From here set the PETSC_DIR

    export PETSC_DIR=/cygdrive/c/software/petsc

and continue to build PETSC-3.21.6 from the same Cygwin terminal.

A limited set of tests indicate that the "simpler" build of PETSc and Crunch with downloaded Blas-Lapack libraries (no MKL) and no MPI is faster.  So we recommend using this configure.py script to build (first) PETSc, then Crunch:

    export PETSC_DIR=/cygdrive/c/software/petsc

    ./configure PETSC_ARCH=oneAPI-noMPI-ifx-opt \
    --with-cc=/cygdrive/c/software/petsc/lib/petsc/bin/win32fe/win32fe_icx \
    --with-fc=/cygdrive/c/software/petsc/lib/petsc/bin/win32fe/win32fe_ifx \
    --with-cxx=0 \
    --with-mpi=0 \
    --download-fblaslapack \
    --with-debugging=0 \
    --with-shared-libraries=0

And then to use Visual Studio 2022 to build Crunch (use "ifx" Fortran compiler):

    Additional Include Directories
    $(OUTDIR)
    $(PETSC_DIR)
    $(PETSC_DIR)\$(PETSC_ARCH)\lib\petsc\conf
    $(PETSC_DIR)\$(PETSC_ARCH)\include
    $(PETSC_DIR)\include\petsc\finclude
    $(PETSC_DIR)\lib\petsc\conf
    $(PETSC_DIR)\include\
    $(PETSC_DIR)\include\petsc
    $(PETSC_DIR)\$(PETSC_ARCH)\lib
    $(ONEAPI_ROOT)\mkl\mkl\latest\include\mkl\intel64\ilp64
    $(ONEAPI_ROOT)\mpi\latest\include\mpi
    $(ONEAPI_ROOT)\mpi\latest\lib
    ..\

For fully optimized production code, be sure to configure PETSc with "--with-debugging=0" and make sure the CrunchTope Makefile includes the -O3 flag for maximum optimization.

    FFLAGS  = -w -O3 -ffpe-trap=invalid,overflow,zero 

### RunTime

---> Add "H2O" to the end of the list of PRIMARY SPECIES.  
---> Add also to each of the CONDITION blocks  "H2O  55.50843506"

### Pumping Wells
 
---> If there are pumping wells, them "pumpunits" must be set /= 0.0

    ---> Partial Example
    FLOW
    distance_units    centimeters
    time_units        days
    calculate_flow    true
    pumpunits         dm3_min
    pump              0.0006   amendment  1  1  1

    ---> Possible pumpunits:
    cm3_sec
    liter_sec
    dm3_sec
    m3_sec

    cm3_min
    liter_min
    dm3_min
    m3_min

    cm3_hr
    liter_hr
    dm3_hr
    m3_hr

    cm3_day
    liter_day
    dm3_day
    m3_day

    cm3_yr
    liter_yr
    dm3_yr
    m3_yr

### Example of TRANSPORT block

    TRANSPORT
    distance_units         centimeters
    time_units             second
    cementation_exponent   2.00
    fix_diffusion          1.0E-05
    dispersivity           0.0  0.0 
    gas_diffusion          5.E-03
    cementation_exponent   1.00
    MillingtonQuirk        false  
    !!! constant_tortuosity  1.00
    !!! read_tortuosityfile  WhatEver.dat 
    END

### Input File Keyword Blocks

The possible keyword blocks in the input file include
(see markdown files for each Keyword Block [here](https://github.com/CISteefel/CrunchTope/tree/master/docs/input_file)  

[Running CrunchFlow](./docs/input_file/running_crunchflow.md)  
[TITLE](./docs/input_file/title.md)  
[RUNTIME](./docs/input_file/runtime.md)  
[OUTPUT](./docs/input_file/output.md)  
[PRIMARY_SPECIES](./docs/input_file/output.md)  
[SECONDARY_SPECIES](./docs/input_file/secondary_species,md)  
[GASES](./docs/input_file/gases.md)  
[MINERALS](./docs/input_file/minerals.md)  
[AQUEOUS_KINETICS](./docs/input_file/aqueous_kinetics.md)  
[ION_EXCHANGE](./docs/input_file/ion_exchange.md)  
[SURFACE_COMPLEXATION](./docs/input_file/surface_complexation.md)  
[DISCRETIZATION](./docs/input_file/discretization.md)  
[INITIAL_CONDITIONS](./docs/input_file/initial_conditions.md)  
[BOUNDARY_CONDITIONS](./docs/input_file/boundary_conditions.md)  
[INPUT_FILE_FORMATS](./docs/input_file/input_file_formats.md)  
[TRANSPORT](./docs/input_file/transport.md)  
[SATURATED_FLOW](./docs/input_file/saturated_flow.md)  
[RICHARDS_FLOW](./docs/input_file/flow.md)  
[POROSITY](./docs/input_file/porosity.md)  
[TEMPERATURE](./doc/input_file/temperature.md)  
[EROSION](./doc/input_file/erosion.md)  
[PEST](./doc/input_file/pest.md)  
[CONDITION](./doc/input_file/geochemical_conditions.md)  

With the exception of the keyword block CONDITION, the blocks should 
appear only once in the input file. Each keyword block is terminated by 
an END. The keyword block CONDITION is a special case in that it can 
occur multiple times. 

### Shortcourses/Lectures
Check back for announcements of upcoming shortcourses and/or lectures taught by the code developers. In the meantime, check out the shortcourse exercises found in CrunchTope/Exercises and described in "CrunchTope/docs/CrunchCourseDescriptions.pdf" (https://github.com/CISteefel/CrunchTope/blob/master/docs/CrunchCourseExerciseDescriptions.pdf).

### Citations
Steefel, C.I., Appelo, C.A.J., Arora, B., Jacques, D., Kalbacher, T., Kolditz, O., Lagneau, V., Lichtner, P.C., Mayer, K.U., Meeussen, J.C.L. and Molins, S., Moulton, D., Shao, H., Simunek, J., Spycher, N., Yabusaki, S.B., Yeh, G.T. (2015) Reactive transport codes for subsurface environmental simulation. Computational Geosciences, 19, pp.445-478.

Steefel, C.I. (2019) Reactive transport at the crossroads. Reviews in Mineralogy and Geochemistry 85: 1-26.
