### MINERALS

The keyword block **MINERALS** refers to any solid-aqueous phase
heterogeneous reaction, mineral-water reactions being the most common of
these. The keyword parameters for this block are simply the names of the
solid phases involved in the reaction. CrunchFlow currently assumes in
all cases that the reaction involves the dissolution of one mole of the
mineral or solid phase in question (that is, it's stoichiometric
coefficient is --1). Solid phase reactions specified in an input file
involve two separate entries in the database file: 1) a *thermodynamic*
entry which gives the stoichiometry of the reaction, the equilibrium
constants as a function of temperature, the molar volume of the solid
phase, and its molecular weight, and 2) a *kinetic* entry which gives
the rate law and rate coefficients for the reaction. All mineral or
solid-aqueous phase reactions are considered to be kinetic in
CrunchFlow, except in the initialization procedure where equilibrium
with a solid phase may be selected. In practice, equilibrium with
respect to solid phases is attained by choosing rate coefficients that
are large relative to transport rates.

#### Database Formats

<u> Thermodynamic Database Entries </u>

The thermodynamic database for solid-liquid phase reactions begins with
a separate line for each entry, with the mineral or solid name being
enclosed in single quotes. This is followed by the molar volume of the
mineral in units of cm3/mole.

\<'MineralName' \> \<Molar Volume\> \<Number of species in reaction\>
\<Stoichiometric coefficient\> \<'SpeciesName'\> \<Stoichiometric
coefficient\> \<'SpeciesName'\> ...\<Log K array\> \<Molecular Weight\>

In this format, the number of pairs of stoichiometric coefficients and
species names is determined by the preceding value in "Number of species
in reaction" (if there is a mismatch, a read error will result). The
length of the "Log K" array is given by the number following
'temperature points' in the first line of the database file.

An example is given by the entries for quartz and calcite:

> \'Quartz\' 22.6880 1 1.00 \'SiO2(aq)\' -4.6319 -3.9993 -3.4734 -3.0782
> -2.7191 -2.4378 -2.2057 -2.0168 60.0843
>
> \'Calcite\' 36.9340 3 -1.00 \'H+\' 1.00 \'Ca++\' 1.00 \'HCO3-\' 2.2257
> 1.8487 1.3330 0.7743 0.0999 -0.5838 -1.3262 -2.2154 100.0872

<u> Kinetic Database Entries </u>

The kinetic database entries include information on rate formulations
and coefficients to be used in a simulation. Some of these values can be
overwritten in the input file at run time. With the exception of the
case when generic_rates is specified with a value in the **RUNTIME**
block, however, there must be an entry in the database corresponding to
the particular reaction, even if some of the values are not used. The
beginning of the section on solid-liquid phase kinetics in the database
is marked by the string:

> Begin mineral kinetics

which is not enclosed in quotes. Each entry is delimited by a line
beginning with a "+" which speeds up the database searching. This is
followed by a string giving the name of the solid phase. Names given in
the input file (described below) must match these names for a particular
reaction to be loaded. This name must also match an entry in the
thermodynamic part of the database, although this requirement may be
relaxed in the future for the special cases of irreversible reactions
where thermodynamics may be irrelevant. The code also allows for
multiple parallel reactions (multiple reactions running in parallel
which affect the same solid phase) through the database entry on the
line following the solid phase name. This line begins with the string
"label = ", and is to be followed with the label of the rate law. The
following (third) line begins with the string "type = ", which may
include currently up to five different types of rate law (see below).
The fourth line begins with "rate (25C) = ". The fifth line gives the
activation energy for the reaction. While there are different options
after the fifth line depending on the type of reaction involved (see
below), all of the solid-liquid phase kinetic entries have a common
input format for the first five lines following the "+" separator:

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> Calcite\
> label = \<name\>\
> type = \<type of reaction\>\
> rate(25C) = \<rate at 25C in log (base 10) units\>\
> activation = \<activation energy (kcal/mole)\>

The "label" option can be used to construct an overall rate of reaction
that is the sum of a number of parallel reaction pathways. An example
would to use several database entries for "Calcite" with differing
labels and rate laws to capture both the pH independent and pH dependent
rates. After discussing the remainder of the input for the solid-liquid
phase kinetics, an example of multiple parallel reactions involving
calcite will be given.

Several types of reactions may be specified. The five options currently
available are:

type = tst\
type = monod\
type = irreversible\
type = PrecipitationOnly\
type = DissolutionOnly\
type = MonodBiomass

and are described briefly below.

*TST:* These three different reaction types involve different kinds of
input following the "rate(25C)" because of the difference in the rate
formulations. The "tst" rate law is described in detail by Lasaga (1981;
1984) and by Aagaard and Helgeson (1981) and takes the form:

$$\begin{matrix}
R = A_{m}k_{m}\exp\left\lbrack \frac{- E_{a}}{RT} \right\rbrack\prod_{}^{}a_{i}^{n}\left\lbrack 1 - \exp\left( m_{2}g^{m_{3}} \right) \right\rbrack^{m_{1}} \\
g \equiv \frac{\Delta G}{RT} = \ln\left\lbrack \frac{Q}{K_{eq}} \right\rbrack
\end{matrix}$$

where A~m~ is the mineral surface area, k~m~ is the intrinsic rate
constant in units of mol/m^2^/s, E~a~ is the activation energy
(kcal/mole), Q~m~ is the ion activity product for the mineral-water
reaction, K~eq~ is the corresponding equilibrium constant, and
$\prod_{}^{}a_{i}^{n}$is a product representing the inhibition or
catalysis of the reaction by various ions in solution raised to the
power n. Note that this last term operates far from equilibrium, that
is, when Q~m~/K~eq~ is very small (i.e., the mineral is very
undersaturated).

Following the specification of the activation energy is the
specification of the dependences of the far from equilibrium rate on
various species in solution:

dependence : \<species name\> \<exponent\> \<species name\> \<exponent\>
...

The dependence of the far from equilibrium rate is given by pairs of
species names in free format followed by the appropriate exponent.

An example of a database input which provides both pH dependent and pH
independent calcite dissolution rates is given below:

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> Calcite\
> label = default\
> type = tst\
> rate(25C) = -6.19\
> activation = 15.0 (kcal/mole)\
> dependence :\
> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> Calcite\
> label = h+\
> type = tst\
> rate(25C) = -3.00\
> activation = 15.0 (kcal/mole)\
> dependence : H+ 1.0

A general form for a dependence on reaction affinity (or Gibbs free
energy) is given following that given in Burch et al (1993) and Hellmann
and Tisserand (2006). The affinity dependence of the rate is defined
with the parameters *m~1~, m~2~, and m~3~*. A more familiar form would
be:

$$R = \left\lbrack 1 - \left( m3\frac{Q}{K_{eq}} \right)^{m_{2}} \right\rbrack^{m_{1}}$$

The form in the database is

AffinityDependence = m1 m2 m3

As an example, the rate law for albite dissolution proposed by Hellmann
and Tisserand (2006) can be implemented by specifying two parallel rate
laws

> Albite\
> label = HellmannFFE\
> type = tst\
> rate (25C) = -11.9897\
> activation = 0.00 (kcal/mol)\
> dependence :\
> AffinityDependence = 1.00 0.0000798 3.81

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--
>
> Albite\
> label = HellmannCTE\
> type = tst\
> rate (25C) = -13.743\
> activation = 0.00 (kcal/mol)\
> dependence :\
> AffinityDependence = 1.17 1.00 1.00

+\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--

Currently, it is only possible to specify a dependence on reaction
affinity in the database. Only those rate laws in which a dependence on
reaction affinity occurs (TST, PrecipitationOnly, and DissolutionOnly)
will be affected.

<u> Monod</u>:  Monod reactions take a slightly different form than
do TST-type reactions. Following the specification of the activation
energy, the Monod option has a field for inputting various "Monod
terms", that is, the dependence of the reaction rate on electron
acceptors and/or electron donors. The specification of Monod terms takes
the form:

    monod_terms : <electron donor or acceptor> <half-saturation constant> 

The half-saturation constant is assumed to be in concentration units of
moles/kg water. This in turn is followed by a field for the
specification of inhibition terms (also in moles/kgw):

    inhibition_terms : <inhibitor (species name)> <inhibition constant>

The Monod terms take the standard form:

$$R_{m} = k_{\max}\left( \frac{C_{i}}{C_{i} + K_{half}} \right)$$

with the quantities in parentheses representing the "Monod term".
Multiple Monod term may be specified, although the most common approach
is to use a dual Monod form in which a dependence on an electron
acceptor and on an electron donor are provided. An example involving
aerobic respiration is given by:

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> CH2O\
> label = aerobic_respiration\
> type = monod\
> rate(25C) = -9.699\
> activation = 0.0\
> monod_terms : O2(aq) 15.0e-06\
> inhibition :

The inhibition terms take a similar (if inverted) hyperbolic
mathematical form:

$$I_{m} = \frac{K_{in} + C_{i}}{K_{in}}$$

where K~in~ is the inhibition constant and C~i~ refers to the inhibiting
species. An example involving denitrification is given in the following
database input:

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> CH2O\
> label = denitrification\
> type = monod\
> rate(25C) = -10.00 ! Regnier and Steefel (1999)\
> activation = 0.0 (kcal/mole)\
> monod_terms : NO3- 45.0e-06\
> inhibition : O2(aq) 30.0e-06

Note that both of these forms are typical of those used for single as
opposed to dual Monod formulations. In these cases, the organic carbon
is considered to be in excess of any limiting concentration. At this
point, all Monod-type reactions are assumed to be irreversible (that is,
no use is made of the equilibrium constants).

[Irreversible]{.underline}: Irreversible reactions are the simplest kind
in that no use is made of the equilibrium constants. The reactions are
assumed to take the form

$$R_{m} = A_{m}k_{m}\exp\left( \frac{- Ea}{RT} \right)\prod_{}^{}a_{i}^{n}$$

and are therefore similar to the TST form except that a dependence on
the saturation state is missing. Because of the similarity to the TST
reactions, the input is also similar, for example:

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> Iron\
> label = tce\
> type = irreversible\
> rate(25C) = -6.92\
> activation = 15.0 (kcal/mole)\
> dependence : TCE(aq) 1.0

where the rate of dissolution of zero-valent iron is here considered to
have a first-order dependence on the concentration of TCE in solution.

<u>PrecipitationOnly</u>

Invoking this option in the database causes the code to calculate
mineral precipitation only---dissolution is suppressed. An example is
given by:

> KaoliniteYang\
> label = Yang\
> type = PrecipitationOnly\
> rate(25C) = -13.47\
> activation = 0.0 (kcal/mole)\
> dependence :\
> AffinityDependence = 1.00  2.07468 -1.00000
>
> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--

<u>DissolutionOnly</u>

Invoking this option in the database causes the code to calculate
mineral dissolution only---precipitation is suppressed. An example is
given by:

> KaoliniteYang\
> label = Yang\
> type = DissolutionOnly\
> rate(25C) = -12.9393\
> activation = 0.0 (kcal/mole)\
> dependence :\
> AffinityDependence = 1.00 0.5 1.00
>
> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--

This option is useful for specifying different rate laws for
precipitation and dissolution where both PrecipitationOnly and
DissolutionOnly are invoked separately for the same mineral.

<u>MonodBiomass</u>: When a mineral reaction is mediated by
microbes, a monod-type formulation is used with the addition of terms in
the reaction rate expression. One expresses a first order dependence on
biomass concentration (see BIOMASS section below) and the other a
thermodynamic potential factor (Jin and Bethke, 2005):

$$R_{m} = k_{\max}B\frac{C_{i}}{C_{i} + K_{half}}\frac{K_{in} + C_{j}}{K_{in}}\left( 1 - \exp\left( \frac{\Delta G_{redox} + m\Delta G_{p}}{\chi RT} \right) \right)$$

where B is the biomass concentration, $\Delta G_{redox}$ is the Gibbs
free energy of the catabolic pathway and $\Delta G_{p}$ is the
phosphorylation potential (Jin and Bethke, 2005), $m$ is the number of
ATP synthesized and $\chi$ the average stoichiometric number for the
overall reaction.

The kinetic entry is no different than for Monod type reaction. For
example,

> +\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\-\--\
> Ferrihydrite_DIRB_Lc\
> label = DIRB\
  type = MonodBiomass\
> rate(25C) = -11.25\
> activation = 0.0 (kcal/mole)\
> monod_terms : self 0.0005 tot_Lactate 1.0E-3\
> 1 0.0 -1

While the entry in the database file is for the overall stoichiometry of
the reaction including the anabolic (biomass growth) pathway, the
thermodynamic term is calculated exclusively as a function of the
catabolic (or energy-yielding) pathway stoichiometry. As a stop gap
solution, to provide this information, an additional file
(CatabolicPathways.in) must be provided. This file uses the *namelist*
read format available in Fortran 90. For the entry corresponding to the
microbially mediated mineral reaction above, this is required:

> &CatabolicPathway\
> name = Ferrihydrite_DIRB_Lc\
> reaction = -1.833 H+ 1.000 Fe++ 0.250 HCO3- -0.083 Lactate\
> keq = 19.1\
> biomass = \'C5H7O2NFe(s)\'\
> bq = 0.0\
> chi = 1\
> direction = -1

where keq is the equilibrium constant of the catabolic pathway in log10
units, biomass points to the biomass species that mediates the reaction,
bq is $m\Delta G_{p}$, chi is $\chi$ and direction is 1 or -1, being
negative if reactants are written with negative stoichiometric
coefficients. In the future, it is planned that the database entries for
mineral reactions will transition to a format similar to that used for
kinetic aqueous reactions (see AQUEOUS_KINETICS below).

#### Entries in Input File

The simplest possible form for specification of a solid-aqueous phase
reaction in the database is to give the solid or mineral name on a
separate line within the **MINERALS** keyword block without any
accompanying information. If a database sweep has been specified
previously, then the list of minerals may be copied from the file
"'input file name prefix.out" to which they have been written and then
pasted directly into the **MINERALS** keyword block. However, unlike
similar operations for secondary species and gases, the minerals
specified in the input file also require kinetic information which
resides in a separate part of the database from the thermodynamic
entries for minerals or solid phases. Since the kinetic database is not
as complete as the thermodynamic database (letters written to funding
programs within DOE lamenting this situation will be appreciated), some
minerals found during the sweep may not have kinetic database entries
yet. This will be remedied in future releases by 1) augmenting the
kinetic database for CrunchFlow and 2) by sweeping the kinetic database
in addition to the thermodynamic database for solid phase reactions
which can be loaded directly as kinetic reactions. The database sweep
will tend to find a relatively large number of phases which may or may
not actually form at the temperature of interest. It is common, for
example, to find various forms of garnet supersaturated at
25C---geochemists and petrologists know that these phases precipitate
only at high temperature. Thus, editing of the mineral list generated by
a database sweep is essential.

Without any additional information, the code will assume that the
"default" reaction is to be loaded if only the name of the solid phase
is given, that is, it will search for a kinetic entry of that name
labeled as "default".

Following the specification of the reaction, the mineral kinetic options
that can be specified in the input file include:

    -label \<label\>
    -rate \<rate\>
    -activation \<activation energy\>
    -suppress_precipitation \<ε\>
    -local_equilibrium
    -associate \<mineral_name\>

To load solid-liquid phase reactions with labels other than "default",
one uses the form:

> \<Reaction name\> -label \<label\>

where "label" is a string which must match one in the kinetic database.
Similarly, one can overwrite the reaction rate in the database with a
specification in the input file:

> \<Reaction name\> -label \<label\> -rate \<rate\>

where the rate is assumed to be the logarithm of the rate in units of
moles/m^2^/s at 25C. The activation energy may also be overwritten:

> \<Reaction name\> -label \<label\> -activation \<activation energy\>,

where the activation energy is in units of kcal/mol.

Previously, it was possible to specify a threshold for nucleation, but
implemented in the simple fashion in which it was, numerical convergence
was poor due to the discontinuity in rates across the nucleation
threshold. The Newton method for nonlinear equations performs poorly
with functions that are not continuously differentiable. 

For the case of linear kinetics, the code uses an analytical expression
for the approach to equilibrium that is continuously differentiable,
with the rate given by:

$$R_{d} = 0.5 - 0.75\xi + 0.25\frac{\xi^{3}}{\varepsilon^{3}}$$

where the quantity $\xi$is defined by

$$\xi = Log\left( \frac{Q}{K_{eq}} \right) + \varepsilon$$

For nonlinear kinetics, the factor of ε is not used.

Local equilibrium with respect to mineral phases (variously known as the
*local equilibrium* *approximation or fantasy*, depending on one's point
of view) can be specified through a combination of parameters. To
suppress the reduction in reactive surface area as the mineral volume
fraction of 🡪 0, the following option may be specified

    <Reaction name> -local_equilibrium

In other words, the initial reactive surface area remains fixed even
during dissolution as long as the mineral is present at all. This is of
course unrealistic, but the local equilibrium fantasy has been so
drummed into a generation of geochemists and petrologists that it is
impossible to dispense with this option. This alone does not enforce
local equilibrium, but when combined with a sufficiently large rate
constant and reactive surface area, equilibration within about 10^-4^
log(Q/K~eq~) units of equilibrium is possible. The rate constant and/or
reactive surface area should not be set arbitrarily high, however, since
this leads to poor numerical convergence due to the stiffness of the
governing equations. A log rate constant in the range of -3.00 to -5.00
should provide local equilibrium behavior, although this depends in
detail on the magnitude of the transport rates. The following example
should provide local equilibrium with respect to the minerals
*portlandite* and *Ca-oxalate* even under relatively high flow rates

>   MINERALS\
    Portlandite -rate -3.00 -local_equilibrium\
    Ca-Oxalate -rate -3.00 -local_equilibrium\
    END

Another option that has been recently added is given by

    <Reaction name> -associate <mineral_name>

which is used to associate more than one reaction stoichiometry and
equilibrium constant with the same mineral. This is relevant normally
when aqueous kinetics are involved and the reaction for the same mineral
is written as a kinetic reaction involving different sets of primary
species. The most familiar example would be the oxidation of organic
carbon by various electron acceptors, none of which are in equilibrium
with each other (and therefore they are listed as Primary Species---see
below). One example is given by the reductive dissolution of amorphous
iron hydroxide by two different mechanisms: 1) the reductive dissolution
mediated by the bacteria *Geobacter*, and 2) the abiotic reduction by
hydrogen sulfide:

$$\begin{matrix}
FeOOH(s)\  + \ 1.925H^{+} + \ 0.033NH_{4}^{+} + \ 0.208CH_{3}COO^{-}\  \rightarrow \\
Fe^{2 +} + \ 0.033C_{5}H_{7}O_{2}N_{(Fe)}\  + \ 0.25HCO_{3}^{-}\  + \ 1.6H_{2}O
\end{matrix}$$

designated in the input file and database as "Iron_DIRB", while the
reaction

$$2FeOOH + HS^{-} + 5H^{+} \rightarrow 2Fe^{2 +} + S_{0} + 4H_{2}O$$

could be designated as Iron_H2S. The code thinks these are two separate
phases because of the different entries and stoichiometries, unless the
"associate" keyword is used. If one mineral is associated with another,
then only the volume percentage and reactive surface of the "associated"
mineral is updated. So in the example above, the following usage:

> MINERALS\
> Iron_DIRB\
> Iron_H2S -associate Iron_DIRB\
> END

would result in the dissolution according to the Iron_H2S kinetic
pathway affecting the volume fraction and surface area of the Iron_DIRB
phase. To see that this is working, users can check to see which of
these phases is being updated as a function of time in the volume#.out
files.

### BIOMASS

The rate law of microbially-mediated reactions includes a first-order
dependence of the microbial biomass concentration. Currently, only
immobile biomass concentrations are used in this first-order dependence.
As a result, for each microbial biomass considered, one primary species
and one mineral species must be entered in the input file --the former
being the planktonic biomass and the latter the immobile biomass.
Further, a mineral reaction must be included in the database between
this aqueous primary species and the mineral.

#### Database Formats

For example, consider microbial species C5H7O2N mediating one or more
redox reactions.

In the aqueous primary species block, use the format

    'C5H7O2N' 1.0 0.0 113.12

In the mineral species thermodynamic block, use the format

    'C5H7O2N(s)' 1.0 1 1.0 'C5H7O2N' 500.0 -15.0 500.00 500.00 500.00 500.00 500.00 500.00 113.12

where the field used for mineral molar volume is required but ignored.

In the mineral species kinetic block, use the format

> C5H7O2N(s)\
> label = default\
> type = tst\
> rate(25C) = 2.0\
> activation = 15.0 (kcal/mole)\
> dependence :

where a 'tst' rate is used. In this example, the logK for T = 25 C was
set to an arbitrarily small value and the rate constant to an
arbitrarily fast rate constant to ensure that all of the aqueous biomass
produced in microbially mediated reactions is converted instantaneously
to immobile biomass and contributes to the redox reaction rates.

#### Entries in Input File

In the input file, both aqueous and immobile must be included

> PRIMARY_SPECIES\
> [other primary species in the problem]\
> C5H7O2N\
> END
>
> MINERALS\
> [other mineral species in the problem]\
> C5H7O2N(s) -label default -rate 2.0 -biomass\
> END

If a species in the mineral list is associated with immobile biomass, it
must be labeled '-biomass'.

In the geochemical conditions block, use the mineral entry to enter the
initial concentration of biomass in units of mol-biomass/L-H2O, while
use an arbitrarily small initial concentration for the aqueous biomass:

> CONDITION initial\
> units mmol/kg\
> temperature 25.0\
> [other initial conditions]\
> C5H7O2N 1.0E-12\
> C5H7O2N(s) 5.0e-5\
> END

If there is biomass in the aqueous solution used as a boundary
condition, then do use the field for initial aqueous concentrations in
the units specified with the keyword units.