## EROSION

#### Read_BurialFile

Keyword followed by a file name containing erosion/burial values over
the X direction.

Syntax:  &nbsp; read_BurialFile &nbsp; *filename format*

<u> Backwards compatibility </u>: &nbsp; read_burial

<u> Default </u>: &nbsp; *None*

Explanation:  This keyword provides a name for a file
containing burial/erosion rates values in the X coordinate direction.
Other coordinates (Y and Z) are not currently supported. Units are
determined by *space_units* and *time_units* specifications within the
EROSION keyword block. Default units are m/yr.

The format of the file depends on the format specified after the file
name, but the values of the erosion/burial rates are assumed to be
listed with NX varying first, then NY, and then NZ. For those familiar
with FORTRAN and using the *SingleColumn* format (the default), the
actual source code is:

    Do jx = 0,nx
      READ(52,*) FluidBuryX(jx,jy,jz),SolidBuryX(jx,jy,jz)
    End Do