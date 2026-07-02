## POROSITY

This keyword block is used to set parameters and options related to the
treatment of the porosity. The main contrast is between simulations in
which the porosity is fixed and those in which it evolves as the mineral
volume fractions evolve.

#### Fix_porosity

Keyword followed by a value giving the default fixed porosity.

    Syntax:  fix_porosity  value

<u> Default:</u> &nbsp; If not provided, the code calculates porosity
based on the mineral volume fractions.

<u>Explanation:</u> &nbsp;  The keyword parameter *fix_porosity* is used
to set a global fixed porosity for a simulation. When provided, it
disables the calculation of porosity based on the sum of the volume
fractions. Accordingly, no porosity update based on these quantities is
possible. When this parameter is set in the POROSITY keyword block,
*set_porosity* keywords within individual geochemical conditions can be
used once the geochemical condition is distributed in space within the
INITIAL_CONDITION keyword block.

#### Minimum_porosity

Keyword giving the minimum porosity value that is allowed as the mineral
volume fractions evolve.

    Syntax:  minimum_porosity  value

<u> Default:</u> &nbsp; $10^{-14}$

<u>Explanation:</u> &nbsp;  This keyword is only used to set a minimum to
which the porosity can evolve as the mineral volume fractions change.
This option is only used where *porosity_update* is true and neither
*fix_porosity* nor *read_porosity* are set.

#### Porosity_update

Keyword followed by a logical (true or false) used to determine whether
the porosity is updated as the mineral volume fractions evolve.

    Syntax:  porosity_update  logical

<u> Default:</u> &nbsp;  False

<u>Explanation:</u> &nbsp; This keyword is only used where neither the
fix_porosity or read_porosity keywords have been set. When true, it
instructs the code to update the porosity as the mineral volume
fractions evolve according to

$$\phi = \sum_{m = 1}^{Nm}{1 - \phi_{m}},$$

where $N_m$ is the number of minerals in the system and $\phi_{m}$ is
the volume fraction of the individual mineral.

#### Read_PorosityFile

Keyword followed by a file name containing permeability values over the
entire spatial domain.

    Syntax:  read_PorosityFile  filename  format

<u>Backwards compatible:</u> &nbsp; read_porosity

<u> Default </u>: &nbsp; None

<u>Explanation: </u> &nbsp;  This keyword provides the name of file
containing porosity values defined at the center of the grid cell for
the entire spatial domain. The format of the file depends on the file
format specified. In none is provided, a single column format consisting
of a single porosity per line, with NX varying first, then NY, and then
NZ is assumed.