
## PRIMARY_SPECIES

The **PRIMARY_SPECIES** keyword block is used to specify the primary
component species for a particular problem. Borrowing from the tradition
of EQ3/EQ6, there is a list of "database primary species" given in the
first section of any database file. The list of "database primary
species" serves as the building blocks of reactions in the main part of
the database. From the point of view of CrunchFlow, however, primary
species are any that are so designated in the input file. The user may
list any mathematically valid set of primary species which spans the
concentration space for a particular problem. If aqueous kinetic
reactions are included in the problem, then these must be made up
entirely of the user's choice of primary species (not necessarily those
found in the database). It should be noted, however, that some choices
of primary species will result in non-intuitive components, making it
difficult to specify initial and boundary conditions for a particular
problem.

The primary component species (or "basis species") are the chemical
building blocks for a particular problem. The concentrations of the
primary species are the major unknowns in the chemistry part of a
reactive transport problem. The partial differential equations for the
conservation of chemical mass in the system are written in terms of the
total concentration of these primary component species. This
formulation, which has been widely discussed by a number of workers, is
used instead of equations for the conservation of the various chemical
elements.

#### Database Format

"Database primary species" are listed in the database immediately
following the specification of the temperature field (first line) and 3
lines giving the temperature-dependent Debye-Huckel parameters. The
primary species block in the database is terminated with the line:

    'End of primary ' 0.0 0.0 0.0

The "database primary species" are specified in the following form:

    <'SpeciesName'>  <Debye-Huckel size parameter >  <Charge >  <Molecular weight>

Each "database primary species" is given on a separate line in the
database in free format. The species name is enclosed in single quotes.

#### Input File Entry of Primary Species

The format of the **PRIMARY_SPECIES** keyword block is simple,
consisting of a list of species name to be used. In contrast to other
keyword parameters in the input file, the species names are
case-sensitive. This is to remove the ambiguity which may arise from
various species and compounds have similar names (for example, carbon
monoxide or CO and cobalt or Co).

The list of species in the **PRIMARY_SPECIES** keyword block is
important in that it defines the chemical space for a particular problem
(thus the term "component species"). When geochemical conditions are
given, all of the primary species listed in the **PRIMARY_SPECIES**
keyword block must be specified in each of the conditions. In this
respect, **CrunchFlow** differs from a code like **PHREEQC** (Parkhurst,
1995) where solutions may be mixed which have different sets of non-zero
chemical concentrations specified (**PHREEQc** presumably sets
unmentioned concentrations to 0). Since **CrunchFlow** grew out of a
standard finite-difference approach to solving partial differential
equations rather than out of reaction path approaches, it requires a
consistent set of primary species for each cell within the domain and
for the system boundaries.

The two examples below represent two different basis sets of primary
species which can be used to solve the same problem. They are
mathematically equivalent, therefore, given the correct choice of
constraints. The total concentrations of some of the species will not
change from one representation to the other---for example, total
$HCO_3^-$ will equal total $CO_2(aq)$. However, other total
concentrations will change, including that of $H^+$. Use of the
$Fe^{2+}$/$Fe^{3+}$ redox couple is an alternative way of handling redox than
is given by the use of $Fe^{2+}$ and $O_2$(aq). Note that the definition of
total elemental Fe will change depending upon which formulation is used.

    Example 1:

    PRIMARY_SPECIES
    Na+
    K+
    Ca++
    H+
    Al+++
    Fe++
    SiO2(aq)
    HCO3-
    SO4--
    Cl-
    O2(aq)
    H2O
    END

Example 2:

    PRIMARY_SPECIES
    Na+
    K+
    Ca++
    H+
    Al+++
    Fe++
    SiO2(aq)
    CO2(aq)
    SO4--
    Cl-
    Fe+++
    H2O
    END