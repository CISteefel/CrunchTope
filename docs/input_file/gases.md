### GASES

The **GASES** keyword block specifies the user's choice of gases for a
particular problem. In the case of fully saturated flow, gases are only
used in the initialization procedure at the present time since their
mass is not explicitly accounted for. Future releases of the code will
have an option to specify a kinetic reaction with a gas phase which
could be used, for example, to continually equilibrate a solution with
the atmosphere. Clever users will recognize that the same effect can be
obtained presently by specifying a fictional "mineral" in the
**MINERAL** keyword block whose equilibrium constant corresponds to the
aqueous concentration at which a species would be in equilibrium with a
particular partial pressure of the gas phase (an example is given
below).

During the initialization of the various geochemical conditions, the
option is available for equilibrating the solution with one or more
gases at specified partial pressures. In previous versions of the code,
it was necessary to input the gas O2(g) in order to get most of the
redox reactions to load properly. In the current version, however, the
specification of O2(g) is optional since the code will automatically
load it if the aqueous species O2(aq) is found in the primary or
secondary species list. For fully saturated problems, therefore, gases
need to be input only if they are used in constraining a geochemical
condition.

For unsaturated transport and gas-liquid partitioning, gases must be
explicitly entered. In the case of unsaturated transport and reaction,
an actual mass balance based on the individual gas concentrations is carried
out based on the ideal gas law.

If the *database_sweep* option is specified, then the code will
automatically sweep the thermodynamic database looking for relevant
gases based on the user's choice of primary species. Any entries in the
**GASES** keyword block will be ignored in this case. These will then be
listed in the standard output file "*input filename prefix*.out".

#### Database Format

The format for gases in the database is similar to that of the secondary
species:

    <'Gas'>  <Molar volume> <Number of species in reaction (integer)> <Stoichiometric coefficient> <'SpeciesName'> <Stoichiometric coefficient> <'SpeciesName'> ... <Log K array> <Molecular weight>

#### Input File Entry of Secondary Species

Gases are input as a simple list, with one gas per line of the **GASES**
keyword block. If the database sweep option has been chosen previously,
the list of gases written to the file "*input filename prefix*.out" may
be copied from this file and then pasted directly into the **GASES**
block.
