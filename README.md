# Instructions for building PETSc and CrunchFlow on Windows (edited on 2024/10/08)
## Install Microsoft Visual Studio with Fortran compilers

We recommend to use Microsoft Visual Studio Community 2022 to build CrunchFlow on Windows. You can install Microsoft Visual Studio Community 2022 from the [link](https://visualstudio.microsoft.com/downloads/). When installing it, select "Python development" and "Desktop development with C++".

For Fortran compilers, we need to install [Intel oneAPI Base Kit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/base-toolkit.html#gs.9yaz3k) and [oneAPI HPC toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/hpc-toolkit.html#gs.9yazr2).  Install Intel oneAPI Base Kit first and make sure these software are integrated with the Microsoft Visual Studio. If you ever installed the software in the past, we recommend to fully uninstall them because different versions may conflict each other.

The following instructions assume the following versions:

- VisualStudio: Visual Studio Community 2022 17.11.0
- Intel oneAPI Base Toolkit: 2024.2 (C:\Program Files (x86)\Intel\oneAPI)
- Intel HPC Toolkit: 2024.2 (C:\Program Files (x86)\Intel\oneAPI)

### Test your first Fortran Program

Once you install them, create a new Visual Studio project to confirm a Fortran option is available (you can see "Fortran" in "Empty Project" for example) and run a simple Fortran program. You can add a Fortran file by selecting "Add New Item" in "Project". The name of the file can be "test.F90". Copy and paste the codes below to the file:

```
program main      
print *, "Hello, World!"
read(*,*)
end program main
```

Then, build it by selecting "Build Solution" in "Build". You can also test it by pressing 'F5' or selecting 'Start Debugging' in 'Debug' and you can see 'Hello, World!' on the screen. However, if you try to run the executable, which is 'C:\Users\{your_user_name}\source\repos\Console1\Console1\x64\Debug\Console1.exe' for example, your computer may not find "libifcoremdd.dll" file. We will fix this problem after installing PETSc by setting appropriate paths for the environmental variables.

## Building PETSc on Windows

Since the configure scripts for PETSc really work well only with UNIX or Linux type systems, our recommended approach is to Cygwin (follow the section "Native Microsoft/Intel Windows Compilers" in the PETSc documentation on [Microsoft Windows Installation](https://petsc.org/main/install/windows/)).


### Install Cygwin

When installing Cygwin from the [link](https://www.cygwin.com/), make sure you also get the following Cygwin components:
- Python3
- make (Gnu Make within the Devel package)

### Remove Cygwin link.exe

Next, we need to remove Cygwin link.exe to avoid conflict with the Intel ifort compiler. In a Cygwin BASH shell (search on your Windows machine and you can find "Cygwin64 Terminal" and run it as administrator), type the following command:

`mv /usr/bin/link.exe /usr/bin/link-cygwin.exe`

### Install PETSc

It turns out that much of the difficulty in getting PETSc to build easily is due to the failure to find the right Environmental Variables and compilers in the Cygwin BASH shell. The key is to follow the following steps (Note that opening "Intel oneAPI command prompt for Intel 64 for Visual Studio 2022" from Windows as suggested in the PETSc manual did not work for us. In that case, we could not use ifx compiler when buiding PETSc, so we recommend the following steps). Open a Windows (or DOS) Command shell (type “cmd” in Windows). Then navigate within the Command window to:

`C:\Program Files (x86)\Intel\oneAPI\`

and give the command:

`setvars.bat`

This will set the various Intel oneAPI flags. After these are set, run the command within the same Command window where you have just set the Intel Environment Variables (no double-clicking) to

`C:\cygwin64\bin\mintty.exe -`

This launches a bash shell in the Cygwin Unix environment, but it has to be done by command line from the same Windows Command window where the environment variables were set. If everything has worked correctly, the Cygwin bash shell should have inherited the Environment Variable settings from running the “setvars.bat” script. Test for this by now running within the same Cygwin terminal from a directory other than the one actually containing the files so as to test whether the compilers are in the system search paths:

`which icx`

and

`which ifx`

The location of these compilers should be echoed, if not, the paths have not been set correctly. If not (i.e., you get a message like “No ifx found in …”,), then you will need to add the location of these files to the environmental variables of your machine manually. Note that Cygwin may not know the command `which`. In such a case, follow the direction [here](https://stackoverflow.com/questions/14797194/cygwin-ls-command-not-found).

Then, change directories to where you want to install PETSc, usually something like (note that we use the 'cygdrive/c' address rather than 'C:\software') when in the Cygwin terminal):

`cd cygdrive/c/software`

If the directory does not exist, you can create it by

`mkdir software`

after moving to '/cygdrive/c' directory.

Next, install PETSc on your machine. Installing PETSc from Github does not work well. So, download the tar ball from [here link](https://petsc.org/main/install/download/):

`petsc-3.22.0.tar.gz`

and then ideally in C:\software run the following commands:

`gunzip petsc-3.22.0.tar.gz`
`tar xvf petsc-3.22.0.tar`

At this point, you can change the name of the directory from "petsc-3.22.0" to "petsc". Then, we set a variable `PETSC_DIR` by

`cd petsc`
`export PETSC_DIR=$PWD`

This will be using the working directory as PETSC_DIR. Or directly,

`export PETSC_DIR=/cygdrive/c/software/petsc`

### Configure PETSc with MPI (optimized version)

Inside the PETSc folder, run the script below to configure the PETSc. The value set for PETSC_ARCH will override what is set elsewhere (e.g., in Windows Environment Variables, or in .bashrc). One can create as many PETSC_ARCH as needed, since each configure build will create a separate directory with that name. The user can then switch between these various PETSC_ARCH options, using either the Windows Environment Variable setting for PETSC_ARCH, or in the user’s .bashrc profile.

Before configuring PETSc with MPI, run the following command:

`export PETSC_ARCH=oneAPI-MPI-opt`

and run the following command:

```
./configure PETSC_ARCH=oneAPI-MPI-opt \
--with-cc=/cygdrive/c/software/petsc/lib/petsc/bin/win32fe/win32fe_icx \
--with-fc=/cygdrive/c/software/petsc/lib/petsc/bin/win32fe/win32fe_ifx \
--with-cxx=0 \
--with-mpi-include=/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/include \
--with-mpi-lib=/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/lib/impi.lib \
--with-mpiexec=/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/bin/mpiexec.exe  \
--with-blaslapack-dir=/cygdrive/c/PROGRA~2/Intel/oneAPI/2024.2/lib \
--with-debugging=0 \
--with-shared-libraries=0 \
FPPFLAGS=-I/cygdrive/c/PROGRA~2/Intel/oneAPI/mpi/latest/include/mpi
```

You may need to change the version of oneAPI your version from 2024.2 to your version.

### Build PETSc libraries

Follow the direction you got after the build, such as

`make PETSC_DIR=/cygdrive/c/software/petsc PETSC_ARCH=oneAPI-MPI-opt all`

Then, check the build by

`make PETSC_DIR=/cygdrive/c/software/petsc PETSC_ARCH=oneAPI-MPI-opt check`

Once done, you can now use PETSc on Windows! Let's check it with Visual Studio. Restart your comptuer before testing it.

## Adding paths to environmental variables

### Hello World program
We first fix the issue we had when running the executable file of the "Hello, world!". To fix this, you need to add the path to the file "C:\Program Files (x86)\Intel\oneAPI\compiler\2024.2\bin" to the system path by editing the environmental variables. Then, you should be able to run the executable file.

### Buidling CrunchFlow

Next we building CrunchFlow on your computer. If you just create a Visual Studio project yourself and try to build CrunchFlow, it would not work because Visual Studio does not know the path to PETSc! To fix this, we need to change the Visual Studio setting. By right clicking the Visual Studio Project, and selecting property, you can change the setting.

- Fortran General
Add the following paths to Addtitional include directories:
```
$(OUTDIR)
$(PETSC_DIR)
$(PETSC_DIR)\$(PETSC_ARCH)\lib\petsc\conf
$(PETSC_DIR)\$(PETSC_ARCH)\include
$(PETSC_DIR)\include\petsc\flinclude
$(PETSC_DIR)\lib\petsc\conf
$(PETSC_DIR)\include\
$(PETSC_DIR)\include\petsc
$(PETSC_DIR)\$(PETSC_ARCH)\lib
$(ONEAPI_ROOT)\mpi\latest\include
$(ONEAPI_ROOT)\mpi\latest\lib
$(ONEAPI_ROOT)\mpi\latest\include\mpi
..\
```

- Fortran preprocessor
The key setting here is to turn on preprocessing, since this is necessary for compiling the PETSc instructions included in the dominantly Fortran CrunchFlow.

- Linker General
Put these paths to Addtitioanal library directories:

```
$(OUTDIR)
$(PETSC_DIR)
$(PETSC_DIR)\$(PETSC_ARCH)\lib\petsc\conf
$(PETSC_DIR)\$(PETSC_ARCH)\include
$(PETSC_DIR)\include\petsc\flinclude
$(PETSC_DIR)\lib\petsc\conf
$(PETSC_DIR)\include\
$(PETSC_DIR)\include\petsc
$(PETSC_DIR)\include\petsc\mpiuni
$(PETSC_DIR)\$(PETSC_ARCH)\lib
$(ONEAPI_ROOT)\mpi\latest\include
$(ONEAPI_ROOT)\mpi\latest\lib
$(ONEAPI_ROOT)\mpi\latest\include\mpi
..\
```


- Linker Input
Add them to Additional dependencies:
Gdi32.lib User32.lib Advapi32.lib Kernel32.lib Ws2_32.lib libpetsc.lib libfblas.lib libflapack.lib

### Test PETSC_DIR


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

To fix this issue, 