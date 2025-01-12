<h1 align='center'>CrunchFlow</h1>

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

Do not use petsc-3.22 with CrunchFlow, since the newly developed PETSc Fortran hooks there do not work.  Instead, use petsc-3.21.

    git clone https://gitlab.com/petsc/petsc.git --branch v3.21.6 $PETSC_DIR

### Compilation

For fully optimized production code, be sure to configure PETSc with "--with-debugging=0" and make sure the CrunchTope Makefile includes the -O3 flag for maximum optimization.

    FFLAGS  = -w -O3 -ffpe-trap=invalid,overflow,zero 

### RunTime

---> Add "H2O" to the end of the list of PRIMARY SPECIES.  
---> Add also to each of the CONDITION blocks  "H2O  55.50843506"

### Input File Keyword Blocks

The possible keyword blocks in the input file include
[ see markdown files for each Keyword Block in
(https://github.com/CISteefel/CrunchTope/tree/master/docs/input_file) ].

    TITLE 
    [RUNTIME](./doc/input_file/runtime.md)
    [OUTPUT](./doc/input_file/output.md)
    [PRIMARY_SPECIES](./doc/input_file/output.md)
    [SECONDARY_SPECIES](./doc/input_file/secondary_species,md)
    [GASES](./doc/input_file/gases.md)
    [MINERALS](./doc/input_file/minerals.md)
    [AQUEOUS_KINETICS](./doc/input_file/aqueous_kinetics.md)
    [ION_EXCHANGE](./doc/input_file/ion_exchange.md)
    [SURFACE_COMPLEXATION](./doc/input_file/surface_complexation.md)
    [DISCRETIZATION](./doc/input_file/discretization.md)
    [INITIAL_CONDITIONS](./doc/input_file/initial_conditions.md)
    [BOUNDARY_CONDITIONS](./doc/input_file/boundary_conditions.md)
    [TRANSPORT](./doc/input_file/transport.md)
    [FLOW](./doc/input_file/flow.md)
    [POROSITY](./doc/input_file/porosity.md)
    [TEMPERATURE](./doc/input_file/temperature.md)
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
