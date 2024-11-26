<h1 align='center'>CrunchFlow</h1>

CrunchFlow is reactive transport software developed by Carl I. Steefel, Toshi Bandai, and Sergi Molins at Lawrence Berkeley National Laboratory.

# News and Updates
Newly added are capabilities for 1D and 2D fully transient unsaturated flow based on the Richards equation (courtesy Toshi Bandai).  Check out Ex13UnsaturatedFlow, Ex14UnsaturatedFlow_evaporation, and Ex15UnsaturatedFlow_transient...

# Installation
Do not use petsc-3.22 with CrunchFlow, since the newly developed PETSc Fortran hooks there do not work.  Instead, use petsc-3.21. 
git clone https://gitlab.com/petsc/petsc.git --branch v3.21.6 $PETSC_DIR

# RunTime
---> Add "H2O" to the end of the list of PRIMARY SPECIES.  
---> Add also to each of the CONDITION blocks  "H2O  55.50843506"

# Shortcourses/Lectures
Check back for announcements of upcoming shortcourses and/or lectures taught by the code developers.
In the meantime, check out the shortcourse exercises found in CrunchTope/Exercises and described in docs 
      https://github.com/CISteefel/CrunchTope/blob/master/docs/CrunchCourseExerciseDescriptions.pdf

# Citations
Steefel, C.I., Appelo, C.A.J., Arora, B., Jacques, D., Kalbacher, T., Kolditz, O., Lagneau, V., Lichtner, P.C., Mayer, K.U., Meeussen, J.C.L. and Molins, S., Moulton, D., Shao, H., Simunek, J., Spycher, N., Yabusaki, S.B. (2015) Reactive transport codes for subsurface environmental simulation. Computational Geosciences, 19, pp.445-478.

Steefel, C.I. (2019) Reactive transport at the crossroads. Reviews in Mineralogy and Geochemistry 85: 1-26.
