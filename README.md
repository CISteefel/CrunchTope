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

To run on Windows, first download and install Visual Studio 2022.  Then download and install Intel oneAPI Base Toolkit, choosing the option to integrate with Visual Studio.  The installation script should provide the option to link with Visual Studio 2022, if it doesn't, then something is wrong with the Visual Studio installation.  Then download and install Intel oneAPI HPC toolkit.  

Once installed, you should see the option to start a Intel oneAPI command prompt for Intel 64 for Visual Studio in the list of Apps under Windows 11 (Intel oneAPI 2025).  Starting this will set all of the Environment Variables for the use of oneAPI.  Then run from this from the same Command Window: 

    C:\cygwin64\bin\mintty.exe -

to start up a Cygwin terminal (Unix) that inherits the Environment Variables from the oneAPI Command Prompt.  From here set the PETSC_DIR

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

And then to use Visual Studio 2022 to build Crunch (use "ifx" Fortran compiler), with maximum optimizations set (unless you want to use the Visual Studio Debugger, in which Optimizations=None).  Attempts to use Interprocedural Optimization doubled the execution time.

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

For fully optimized production code, be sure to configure PETSc with "--with-debugging=0" and make sure the CrunchTope Makefile includes the -O3 flag for maximum optimization (this is for Linux or Mac).

    FFLAGS  = -w -O3 -ffpe-trap=invalid,overflow,zero 

### Tip of the Day

For those wanting to run CrunchFlow directly from Visual Studio, the easiest way to change from one input directory to another is by going to Configuration Properties/Debugging/Working Directory and copying and pasting the directory (AKA folder) from Windows.  This way, Visual Studio will look for your input files there and also write the output files to the same directory.

### RunTime

---> Add "H2O" to the beginning of the list of PRIMARY SPECIES (Species #1).  
---> Add also to each of the CONDITION blocks  "H2O  55.50843506"

### Thermodynamic Database

It turns out that the "Pflotran" thermodynamic database is almost identical to that used by CrunchFlow.  This is because lead Crunch developer Carl Steefel "borrowed" this format from Peter Lichtner way back in 1991 or 1992 when they were both at Universitat Bern in Switzerland.

Now it is possible to easily generate a thermodynamic database for a temperature array and at various pressures using pyGCC (Awolayo and Tutolo, 2022). See the pyGCC website for more details

    https://pygcc.readthedocs.io/en/latest/

Using any Python driver routine like Spyder or Colab, you can run the Python scripts for pyGCC. The following would import pyGCC and indicate to use the thermo.2021 database (that apparently goes with Geochemists Workbench, or GWB). This should successfully create the database with 7 temperature points at a pressure of 230 bars (23 MPa).  This produces a nicely readable thermodynamic database in the GWB format.

    import pygcc
    from pygcc.pygcc_utils import *
    ps = db_reader(sourcedb = 'thermo.2021', sourceformat = 'gwb')

    write_database(T = np.array([0.010, 25, 50, 100, 150, 200, 250]), P = 230, solid_solution = 'Yes', clay_thermo = 'Yes', 
    sourcedb = 'thermo.2021', dataset = 'GWB',print_msg = True)

Then to create the Pflotran-formatted or Crunch-formatted database, you can follow the above (still inside your Python driver) with:

    write_database(T = np.array([0.010, 25, 50, 100, 150, 200, 250]), P = 230, solid_solution = 'Yes', clay_thermo = 'Yes', 
    sourcedb = 'thermo.2021', dataset = 'Pflotran', sourceformat = 'GWB', print_msg = True)

This will create a Pflotran-formatted thermodynamic database for the same temperature range and pressure as the more nicely formatted GWB database. Using Spyder, I find this under my name on Windows and "output/Pflotran". 

Now you just need to add in the Debye-Huckel parameters to the new Crunch database for these temperature and pressures (missing for some reason in Pflotran) by copying and pasting from the GWB-formatted database

    * temperatures (degC)
          0.0100     25.0000     50.0000    100.0000
        150.0000    200.0000    250.0000
    * pressures (bar)
        230.0000    230.0000    230.0000    230.0000
        230.0000    230.0000    230.0000
    * debye huckel a (adh)
          0.4877      0.5052      0.5286      0.5904
          0.6718      0.7772      0.9205
    * debye huckel b (bdh)
          0.3251      0.3285      0.3325      0.3416
          0.3517      0.3626      0.3746
    * bdot
          0.0342      0.0374      0.0409      0.0471
          0.0495      0.0445      0.0288

into the Crunch format:

    'temperature points' 7   0.01  25.  50.  100. 150. 200. 250.
    'Debye-Huckel adh'  0.4877   0.5052  0.5286  0.5904  0.6718  0.7772  0.9205
    'Debye-Huckel bdh'  0.3251   0.3285  0.3325  0.3416  0.3517  0.3626  0.3746
    'Debye-Huckel bdt'  0.0342   0.0374  0.0409  0.0471  0.0495  0.0445  0.0288	 

And then change the end of section reads from 

    'null' 0 0 0
    !:species_name  num (n_i A_i, i=1,num)  log K (1:8)  a0  valence  formula weight [g]

and

    'null' 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
    !:gas_name molar_vol  num (n_i A_i, i=1,num) log K (1:8)  formula weight [g]

and

    'null' 0. 1 1. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.
    !:mineral_name molar_vol  num (n_i A_i, i=1,num) log K (1:8)  formula weight [g]

and

    'null' 0. 1 1. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0.

and so on to:

    'End of primary'   0.0  0.0  0.0
    'End of secondary' 1  0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
    'End of gases'     0.  1  1. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 
    'End of minerals'  0.  1  0. '0' 0. 0. 0. 0. 0. 0. 0. 0.

I deleted the "null" comment at the end of minerals and before oxides altogether (Crunch would consider the oxides as just additional minerals).  So delete this line in the Pflotran-->Crunch database altogether:

    'null' 1 0. '0' 0. 0. 0. 0. 0. 0. 0. 0. 0. 0. 0.
    !:oxide_name molar_vol  num (n_i A_i, i=1,num) log K (1:8)  formula weight [g]

Then append a mineral kinetics section (see Crunch Manual and Shortcourse Exercises for more kinetic options)

    Begin mineral kinetics
    +----------------------------------------------------
    Fo90
      label = default
      type  = tst
      rate(25C) = -6.00
      activation = 15.0  (kcal/mole)
      dependence :
    +----------------------------------------------------
    Lizardite
      label = default
      type  = tst
      rate(25C) = -6.00
      activation = 15.0  (kcal/mole)
      dependence :  H+  -0.13
    +----------------------------------------------------
    Magnetite
      label = default
      type  = tst
      rate(25C) = -6.00
      activation = 15.0  (kcal/mole)
      dependence :
    +----------------------------------------------------
    Brucite
      label = default
      type  = tst
      rate(25C) = -6.00
      activation = 15.0  (kcal/mole)
      dependence :
    +----------------------------------------------------

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
