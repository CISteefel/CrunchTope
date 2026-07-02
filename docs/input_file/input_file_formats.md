## Input File Formats

#### Explanation

Various spatially-dependent field variables can be calculated outside of CrunchFlow and read in for use.
These include the temperature, permeability, tortuosity, porosity, liquid saturation, and  erosion/burial rates.
Typically these are specified with a keyword (see below), followed by the file name.
As of August 2006, it is also possible to specify the format of the file.
The basic syntax is:

    keyword     filename    InputFileFormat

where the keyword may be one of the following:

- read_VelocityFile
- read_GasVelocityFile
- read_PermeabilityFile
- read_PorosityFile
- read_SaturationFile
- read_TortuosityFile
- read_BurialFile
- read_TemperatureFile

The possible file formats include

- SingleColumn (the default)
- ContinuousRead
- FullFormat (dummy X and optionally Y and Z coordinates)
- Unformatted (binary)

##### SingleColumn Example

    DO jz = 1,nz
       DO jy = 1,ny
          DO jx = 1,nx
             READ(52,*) tortuosity(jx,jy,jz)
          END DO
       END DO
    END DO

##### ContinuousRead Example

    READ(52,*) (qx(jx,jy,jz),jx=0,nx),jy=1,ny),jz=1,nz)
    READ(53,*) (qy(jx,jy,jz),jx=1,nx),jy=0,ny),jz=1,nz)
    READ(54,*) (qz(jx,jy,jz),jx=1,nx),jy=1,ny),jz=0,nz)

##### FullFormat Example

```
    DO jz = 1,nz
       DO jy = 1,ny
          DO jx = 1,nx
             READ(52,*) xdum, ydum, zdum, tortuosity(jx,jy,jz)
          END DO
       END DO
    END DO
```

##### Unformatted Example (in 1D)

    READ(52) qx

```
