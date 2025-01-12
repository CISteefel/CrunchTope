
## DISCRETIZATION

This keyword block is used to set up the structured grid for a
particular problem. The code assumes (presently) that rectilinear
coordinates are used. The three keyword parameters which are recognized
are:

    xzones
    yzones
    zzones

each being used to set the grid dimensions in one coordinate direction.

#### Xzones

Discretization in the first (X) direction.

Syntax: &nbsp;  xzones &nbsp;*#cells &nbsp; spacing &nbsp; #cells &nbsp; spacing ... &
#cells &nbsp; spacing*

where *#cells* is an integer and *spacing* is a real number.

<u> Default </u>:  &nbsp; *1 &nbsp; 1.0* &nbsp; if **xzones** is not specified.

Explanation:  This keyword parameter specifies the
discretization in the X (first) coordinate direction. In one-dimensional
problems, the X coordinate rather than the Y or Z coordinate should be
used. The discretization is given by specifying pairs of integers and
real numbers which give the number of cells of a given spacing. The
following example would result in 10 cells of 1 meter spacing, 10 cells
of 2 meter spacing and 10 cells of 4 meter spacing:

    xzones 10 1.0 10 2.0 10 4.0

resulting in a total of 30 cells in the X direction. Lines may be
continued on the following line by using an ampersand ("&") at the end
of the line as a continuation marker.

#### Yzones

Discretization in the first (Y) direction.

Syntax: &nbsp;  yzones &nbsp;*#cells &nbsp; spacing &nbsp; #cells &nbsp; spacing ... &
#cells &nbsp; spacing*

where *#cells* is an integer and *spacing* is a real number.

<u> Default </u>:  &nbsp; *1 &nbsp; 1.0* &nbsp; if **yzones** is not specified

#### Zzones

Discretization in the first (Z) direction.

Syntax: &nbsp;  zzones &nbsp;*#cells &nbsp; spacing &nbsp; #cells &nbsp; spacing ... &
#cells &nbsp; spacing*

where *#cells* is an integer and *spacing* is a real number.

<u> Default </u>:  &nbsp; *1 &nbsp; 1.0* &nbsp; if **zzones** is not specified

