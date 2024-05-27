# Instructions for Building PETSc and CrunchFlow on Windows
## Install Microsoft Visual Studio with Fortran compilers

You can install Microsoft Visual Studio 2022 from the [link](https://visualstudio.microsoft.com/downloads/). Unfortunately, Microsoft Visual Studio 2019 is not available at the moment.

For Fortran compilers, we need to install [Intel oneAPI Base Kit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html#gs.9yaz3k) and [oneAPI HPC toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html#gs.9yazr2).  It is recommeneded to install Intel oneAPI Base Kit first. The following instructions assume the following versions:
VisualStudio: Visual Studio Community 2022 17.10.0
Intel oneAPI Base Toolkit: 2024.1 (C:\Program Files (x86)\Intel\oneAPI)
Intel HPC Toolkit: 2024.1 (C:\Program Files (x86)\Intel\oneAPI)

Once you stalled them, create a new Visual Studio project to confirm Fortan option is available (you can see "Fortran" in "Empty Project" for example) and run a simple Fortran program:
```
program main      
print *, "Hello, World!"
read(*,*)
end program main
```
 Make sure that Visual Studio can find Fortran compiler. You may need to change the system variable
`IFORT_COMPILER24`
to where the Fortran compiler is. For example, it can be
`C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\bin`

## Building PETSc on Windows

Since the configure scripts for PETSc really work well only with UNIX or Linux type systems, the recommended approach is to Cygwin (see the PETSc documentation on [Microsoft Windows Installation](https://petsc.org/main/install/windows/)).


### Install Cygwin

When installing Cygwin from the [link](https://www.cygwin.com/), make sure you also get the following Cygwin components:
- Python3
- make (Gnu Make within the Devel package)

### Remove Cygwin link.exe

Next, we need to remove Cygwin link.exe to avoid conflict with the Intel ifort compiler. In a Cygwin BASH shell, you type:

`mv /usr/bin/link.exe /usr/bin/link-cygwin.exe`

### Install PETSc

It turns out that much of the difficulty in getting PETSc to build easily is due to the failure to find the right Environmental Variables and compilers in the Cygwin BASH shell.  Open a Windows (or DOS) Command shell (type “cmd” in Windows). Then navigate within the Command window to:

`C:\Program Files (x86)\Intel\oneAPI\`

and give the command:

`setvars.bat intel64`

This will set the various Intel oneAPI flags. After these are set, change directories within the same Command window where you have just set the Intel Environment Variables (no double-clicking) to

`C:\cygwin64\bin`

then login in to a Cygwin (bash) shell with:

`bash.exe -login`

This launches a bash shell in the Cygwin Unix environment, but it has to be done by command line from the same Windows Command window where the environment variables were set.

If everything has worked correctly, the Cygwin bash shell should have inherited the Environment Variable settings from running the “setvars.bat” script.  Test for this by now running within the same Cygwin terminal from a directory other than the one actually containing the files so as to test whether the compilers are in the system search paths:

`which icx`

and

`which ifx`

The location of these compilers should be echoed, if not, the paths have not been set correctly. If not (i.e., you get a message like “No ifort found in …”,), then you will need to add the location of these files manually when configuring PETSc. Cygwin may not know the command `which`. In such a case, follow the direction [here](https://stackoverflow.com/questions/14797194/cygwin-ls-command-not-found).

Then, change directories to where you want to install PETSc, usually something like (note that we use the cygdrive/c address rather than C:\software) when in the Cygwin terminal):

`cd /cygdrive/c/software`

Next, install PETSc on your machine. You can install PETSc from Github repository by following the steps:

`mkdir ~/software`

`$ cd ~/software`

`$ git clone -b release https://gitlab.com/petsc/petsc petsc`
`git pull`

Then, we set a variable `PETSC_DIR` by

`$ cd petsc`
`export PETSC_DIR=$PWD`

This will be using the working directory as PETSC_DIR. Or directly,

`export PETSC_DIR=/cygdrive/c/software/petsc`

### Configure PETSc without MPI

Inside the PETSc folder, run the script below to configure the PETSc. The value set for PETSC_ARCH will override what is set elsewhere (e.g., in Windows Environment Variables, or in .bashrc).  One can create as many PETSC_ARCH as needed, since each configure build will create a separate directory with that name.  The user can then switch between these various PETSC_ARCH options, using either the Windows Environment Variable setting for PETSC_ARCH, or in the user’s .bashrc profile.

For statically linked libraries for optimized version is

`./configure PETSC_ARCH=win64-opt --with-cc=icx --with-fc=ifx --with-cxx=0 --with-mpi=0 --download-fblaslapack --with-debugging=0 --with-shared-libraries=0`

You may need to change the directory for the blas and lapack in the Intel oneAPI.

For MPI,
`./configure PETSC_ARCH=mpi-oneAPI-opt --with-cc=icx --with-fc=ifx --with-cxx=0 --with-mpi-include=/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/include --with-mpi-lib=/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/lib/impi.lib --with-mpiexec=/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/bin/mpiexec.exe --download-fblaslapack --with-debugging=0 --with-shared-libraries=0 FPPFLAGS=-I/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/include/mpi`

### Build PETSc libraries

Follow the direction you got after the build, such as

`make PETSC_DIR=/cygdrive/c/software/petsc-3.21.1 PETSC_ARCH=win64-opt all`

Then, check the build by

`make PETSC_DIR=/cygdrive/c/software/petsc-3.21.1 PETSC_ARCH=win64-opt check`

### Visual Studio Setting

In order to successfully compile codes with PETSc, there are several things to add to the visual studio project setting. By right clicking the Visual Studio Project, and selecting property, you can change the setting.

- Fortran General
Add the following paths to Addtitional include directories:
$(PETSC_DIR)
$(PETSC_DIR)\include
$(PETSC_DIR)\$(PETSC_ARCH)\lib\petsc\conf
$(PETSC_DIR)\$(PETSC_ARCH)\include
$(PETSC_DIR)\$(PETSC_ARCH)\lib
$(PETSC_DIR)\include\petsc\finclude
$(PETSC_DIR)\lib\petsc\conf
$(PETSC_DIR)\include\petsc
$(PETSC_DIR)\include\petsc\mpiuni
$(I_MPI_ONEAPI_ROOT)\include
$(I_MPI_ONEAPI_ROOT)\lib
$(ONEAPI_ROOT)\mkl\latest\lib

- Fortran preprocessor
The key setting here is to turn on preprocessing, since this is necessary for compiling the PETSc instructions included in the dominantly Fortran CrunchFlow.

- Linker General
Put these paths to Addtitioanal library directories:

$(PETSC_DIR)\$(PETSC_ARCH)\lib
$(I_MPI_ONEAPI_ROOT)\lib\release
$(I_MPI_ONEAPI_ROOT)\include
$(ONEAPI_ROOT)\mkl\2024.1.0\lib\intel64

- Linker Input
Add them to Additional dependencies:
Gdi32.lib User32.lib Advapi32.lib Kernel32.lib Ws2_32.lib libpetsc.lib libfblas.lib libflapack.lib

### Test PETSC_DIR

You should be able to compile the program below (this is from the PETSc website).
```
program main
#include <petsc/finclude/petsc.h>
USE petsc
implicit none
!
!  This example demonstrates basic use of the PETSc Fortran interface
!  to vectors.
!
PetscInt  n
PetscErrorCode ierr
PetscBool  flg
PetscScalar      one,two,three,dot
PetscReal        norm,rdot
Vec              x,y,w
PetscOptions     options

n     = 20
one   = 1.0
two   = 2.0
three = 3.0

PetscCallA(PetscInitialize(ierr))
PetscCallA(PetscOptionsCreate(options,ierr))
PetscCallA(PetscOptionsGetInt(options,PETSC_NULL_CHARACTER,'-n',n,flg,ierr))
PetscCallA(PetscOptionsDestroy(options,ierr))

! Create a vector, then duplicate it
PetscCallA(VecCreate(PETSC_COMM_WORLD,x,ierr))
PetscCallA(VecSetSizes(x,PETSC_DECIDE,n,ierr))
PetscCallA(VecSetFromOptions(x,ierr))
PetscCallA(VecDuplicate(x,y,ierr))
PetscCallA(VecDuplicate(x,w,ierr))

PetscCallA(VecSet(x,one,ierr))
PetscCallA(VecSet(y,two,ierr))

PetscCallA(VecDot(x,y,dot,ierr))
rdot = PetscRealPart(dot)
write(6,100) rdot
100  format('Result of inner product ',f10.4)

PetscCallA(VecScale(x,two,ierr))
PetscCallA(VecNorm(x,NORM_2,norm,ierr))
write(6,110) norm
110  format('Result of scaling ',f10.4)

PetscCallA(VecCopy(x,w,ierr))
PetscCallA(VecNorm(w,NORM_2,norm,ierr))
write(6,120) norm
120  format('Result of copy ',f10.4)

PetscCallA(VecAXPY(y,three,x,ierr))
PetscCallA(VecNorm(y,NORM_2,norm,ierr))
write(6,130) norm
130  format('Result of axpy ',f10.4)

PetscCallA(VecDestroy(x,ierr))
PetscCallA(VecDestroy(y,ierr))
PetscCallA(VecDestroy(w,ierr))
PetscCallA(PetscFinalize(ierr))

print *, "Hello, World!"
read(*,*)

end program main
```
### List of Environmental variables
PETSC_ARCH
win64-opt
mpi-oneAPI-opt

PETSC_DIR
C:\software\petsc-3.21.1

I_MPI_ONEAPI_ROOT
C:\Program Files (x86)\Intel\oneAPI\mpi\latest

IFORT_COMPILER24
C:\Program Files (x86)\Intel\oneAPI\compiler\2024.1\bin
