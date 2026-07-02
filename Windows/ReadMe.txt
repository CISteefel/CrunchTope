To run the CrunchFlow.exe executable, place it in a directory/folder where you have full permissions.

Then download the Intel Fortran (ifx) redistributable libraries for Intel oneAPI 2024.1 (2024.1.0)

https://registrationcenter-download.intel.com/akdlm/IRC_NAS/f6a44238-5cb6-4787-be83-2ef48bc70cba/w_ifort_runtime_p_2024.1.0.968.exe

NOTE: Current version of Intel oneAPI is 2024.1 on Carl Steefel's computer (and he built it).  Later versions of the redistributable libraries may work, but no guarantee.  He will upgrade his Intel oneAPI soon, but don't wait for that.

For building CrunchTope from scratch, download first the latest version of Visual Studio (free), and then the free Intel oneAPI Base and HPC packages.  Then build CrunchTope with the instructions provided (somewhat involved insofar as you have to run the Intel environment variable settings followed by launching a Cygwin shell).
