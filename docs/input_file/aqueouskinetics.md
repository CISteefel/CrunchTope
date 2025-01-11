### AQUEOUS_KINETICS

The keyword block and database entries for aqueous phase (homogeneous)
kinetics have a different format than that used for mineral-water (or
solid-liquid) reactions. Aqueous phase reactions must be made up
entirely of primary species---the concept of a "secondary species" in
CrunchFlow is restricted to reactions where an equilibrium relationship
between the species in the reaction applies. In an aqueous phase kinetic
reaction, species which might have been secondary species if equilibrium
applied need to be moved to the list of primary species or the reaction
will not be loaded. Consider the case of the Fe(II)-Fe(III) redox
couple. In the case where the oxidation of Fe(II) by molecular oxygen
occurs rapidly and therefore can be assumed to be at equilibrium, one
could use a list of primary species given by the following to describe
the system:

> PRIMARY_SPECIES\
> H+\
> Fe++\
> O2(aq)\
> END

followed by specification of Fe+++ in the list of secondary species.
This would then load the reaction (neglecting the trailing equilibrium
constants):

    'Fe+++' 4 -0.50 'H2O' 0.25 'O2(aq)' 1.00 'Fe++' 1.00 'H+'

along with its associated equilibrium constants. Fe$^{3+}$ would then
become a secondary species contributing to the total concentrations of
all of the primary species in the reaction (not only Fe$^{2+}$, but also
O$_{2}$(aq) and H$^{+}$). A fuller discussion of this can be found in Steefel
and MacQuarrie (1996) or in Bethke (1996). However, if one wanted to
treat this reaction as a kinetic one, then Fe$^{3+}$ needs to be listed as
a primary species along with the others:

> PRIMARY_SPECIES\
> H+\
> Fe++\
> O2(aq)\
> Fe+++\
> END

If this is done, then the reaction involving Fe$^{3+}$ given above will
not be loaded from the thermodynamic portion of the database. Left like
this, Fe$^{2+}$ and Fe$^{3+}$ would be treated as if they were separate
elements with no reaction relationship linking their concentrations.
This is the so-called "redox disequilibrium" option found in a number of
the reaction path codes. In **CrunchFlow**, however, one can add a
kinetic reaction linking these species using the **AQUEOUS_KINETICS**
keyword block.

Redox reactions are often mediated by microbial species (see **Biomass**
section). The **AQUEOUS_KINETICS** keyword block may also be used for
microbially mediated reactions.

#### Database Formats

While previously the aqueous kinetics inputs were provided in the same
database file as equilibrium and mineral reactions, this option has now
been deprecated. Instead, a new database file specific to aqueous
kinetic reactions must be used. Because this is a separate file, the
name of the file needs to be entered specifically. In the **RUNTIME**
keyword block, use:

>     kinetic_database myAqueousDB.dbs

This file uses the *namelist* read format available in Fortran 90. The
format is rather intuitive and flexible but some simple rules must be
observed. Specifically, it is important to the use of single quotes
('example') to delimit strings such as names of species. Each entry
has a form similar to

> &Name\
> StringVariable = 'example'\
> IntegerVariable = 11\
> /

It is headed by a namelist Name preceded by an ampersand and it is
closed by a forward slash (/). In between, a list of variables are
assigned values via a single equal sign (=).

Two entries must be provided for each aqueous kinetic reaction: one with
the stoichiometry of the reaction and equilibrium constant (&Aqueous
*namelists*), and another one with the kinetic rate law model for that
reaction (&AqueousKinetics *namelists*). Because the namelist format is
flexible, the order in which these are provided is not important.

The &Aqueous *namelists* have the form

> &Aqueous\
> name = \'AceNO3HCO3NO2\'\
> type = \'catabolic\'\
> stoichiometry = 0.125 \'H+\' 0.250 \'HCO3-\' -0.125 \'Acetate\' -0.500
> \'NO3-\' 0.500 \'NO2-\'\
> keq = 12.8969\
> /

where stoichiometry is the reaction stoichiometry and keq is the
logarithm of the equilibrium constant. At this time, the stoichiometry
must be provided as a function only of primary species. Providing the
type of the reaction is optional with two options available: 'anabolic'
and 'catabolic'(see below). The name must match the name of the reaction
or one of the pathways in the input file (see below).

A minimum of one rate law needs to be given in the database along with
the kinetic reaction stoichiometry. . For this purpose, the
&AqueousKinetics namelist is required. In this &AqueousKinetics
namelist, the variable "type" refers to rate laws of the type "tst",
"monod", "MonodBiomass" or "irreversible" just as in the solid-liquid
phase reaction input (although the format of the input clearly differs).
A fuller description of the various reaction "types" is given below
along with an example of the format.

#### Types of Rate Laws

<u>TST</u>: The "tst" rate law is described in detail by Lasaga
(1981; 1984) and by Aagaard and Helgeson (1981) and takes the form:

$$R_{s} = k_{s}\left( 1 - \frac{Q_{s}}{K_{eq}} \right)\prod_{}^{}a_{i}^{n}$$

where $k_s$ is the intrinsic rate constant in units of mol/kg water/yr,
$Q_s$ is the ion activity product for the mineral-water reaction,
$K_{eq}$ is the corresponding equilibrium constant, and
$\prod_{}^{}a_{i}^{n}$is a product representing the inhibition or
catalysis of the reaction by various ions in solution raised to the
power *n*. Unlike the database entries for solid-liquid phase reactions,
there is no dependence on surface area and currently no allowance for a
temperature dependence. As an example, consider the kinetically
controlled oxidation of Fe++ by molecular oxygen which uses two parallel
reactions (one pH dependent and one pH independent) following the rate
law proposed in Singer and Stumm (1970):

>     &Aqueous
>     name = \'FeII_oxidation\'
>     stoichiometry = 1.0 'Fe+++' 0.500 'H2O' -0.25 'O2(aq)' -1.0 'Fe++' -1.0 'H+'
>     keq = 8.4887
>     /
>
>     &AqueousKinetics
>     name = \'FeII_oxidation\'
>     type = \'tst\'
>     rate25C = 1.53e-06
>     dependence = 'tot_Fe++' 1.0 'O2(aq)' 1.0 'H+' -2.0
>     /
>
>     &AqueousKinetics
>     name = \'FeII_oxidation\'
>     type = \'tst\'
>     rate25C = 41.4848
>     dependence = 'tot_Fe++' 1.0 'O2(aq)' 1.0
>     /

As in the case of solid-liquid phase reactions, the reaction rate can be
made to depend on the total concentration of a species (Fe++ in this
case) by prepending the string "tot\_" to the species name. Note that
the name of the reaction must be the same in all entries.

<u>Monod:</u> Monod reactions take a slightly different form
than do TST-type reactions. The Monod option has a field for inputting
various "Monod terms", that is, the dependence of the reaction rate on
electron acceptors and/or electron donors:

$$R_{m} = k_{\max}\left( \frac{C_{i}}{C_{i} + K_{half}} \right)$$

where $C_i$ is the concentration of the electron acceptor or donor,
$K_{half}$ is the half-saturation constant in moles/kg water and $k_{max}$
is the maximum rate of the reaction (mol/kg water/yr) (when the
concentration of the electron donor or acceptor \>\> the half-saturation
constant)

The specification of Monod terms takes the form:

>     &AqueousKinetics
>     name = 'ReactionName'
>     type = 'monod'
>     rate25C = MaxRate
>     monod_terms = 'electron acceptor or donor 1' half-saturation 1 'electron acceptor or donor 2' half-saturation 2 \...
>     inhibition = 'inhibitor' inhibition-constant \...
>     direction = -1
>     /

where the "electron donor or acceptor" is a species name which again can
be prepended with the string "tot\_" to designate a total
concentration.and the half-saturation constant is the concentration of
the electron acceptor or donor at which the rate is ½ of its maximum
value.

Inhibition terms of the following form may also be included:

$$I_{m} = \frac{K_{in} + C_{i}}{K_{in}}$$

where K$_{in}$ is the inhibition constant and C$_i$ is the concentration
of the electron acceptor or donor. At the present time, no parallel
pathways are allowed for monod-type rate laws. Direction is negative if
reactants are written with negative stoichiometric coefficients.

<u>Irreversible:</u> Irreversible reactions are the simplest
kind in that no use is made of the equilibrium constants. The reactions
are assumed to take the form

$$R_{s} = k_{s}\prod_{}^{}a_{i}^{n}$$

and are therefore similar to the TST form except that a dependence on
the saturation state is missing. Such a reaction rate law can be used to
an exponential decay reaction. In combination with a series of other
decay reactions, it can be used to describe an entire decay chain. The
decay chain ^241^Pu 🡪 Am 🡪 Np can be represented by the following
entries in the database:

>     &Aqueous
>     name = 'Pu++++241_decay'
>     stoichiometry = -1.0 'Pu++++241' 1.0 'Am+++'
>     /
>
>     &AqueousKinetics
>     name = 'Pu++++241_decay'
>     type = 'irreversible'
>     rate25C = 4.83E-2
>     dependence = 'tot_Pu++++241' 1.0
>     /
>
>     &Aqueous
>     name = 'Am+++_decay'
>     stoichiometry = -1.0 'Am+++' 1.0 'Np++++'
>     /
>
>     &AqueousKinetics
>     name = 'Am+++_decay'
>     type = 'irreversible'
>     rate25C = 1.60E-3
>     dependence = 'tot_Am+++' 1.0
>     /

Since the reactions are designated as irreversible, no equilibrium
constants are needed and if provided they will be ignored.

<u>Irreversible:</u> Microbially mediated Monod-type reactions
take a form similar to that of Monod reactions but do require additional
parameters associated with the inclusion of biomass and a thermodynamic
potential factor in the rate law (Jin and Bethke, 2005):

$$R_{m} = k_{\max}B\frac{C_{i}}{C_{i} + K_{half}}\frac{K_{in} + C_{j}}{K_{in}}\left( 1 - \exp\left( \frac{\Delta G_{redox} + m\Delta G_{p}}{\chi RT} \right) \right)$$

where B is the biomass concentration, $\Delta G_{redox}$ is the Gibbs
free energy of the catabolic pathway and $\Delta G_{p}$ is the
phosphorylation potential (Jin and Bethke, 2005), $m$ is the number of
ATP synthesized and $\chi$ the average stoichiometric number for the
overall reaction.

For acetate oxidation using nitrate as electron acceptor, one would
write:

>     &Aqueous
>     name = 'AceNO3HCO3NO2'
>     type = 'catabolic'
>     stoichiometry = 0.125 'H+' 0.250 'HCO3-' -0.125 'Acetate' -0.500 'NO3-' 0.500 'NO2-'
>     keq = 12.8969
>     /
>
>     &Aqueous
>     name = 'C5H7O2N_RCH2_Ace_NH4'
>     type = 'anabolic'
>     stoichiometry = -0.075 'H+' -0.125 'Acetate' -0.050 'NH4+' 0.050 'C5H7O2NNO3'
>     /
>
>     &AqueousKinetics
>     name = 'AceNO3HCO3NO2'
>     label = 'default'
>     type = 'MonodBiomass'
>     rate25C = 2000.0
>     monod_terms = 'tot_Acetate' 2.03E-5 'tot_NO3-' 1.06E-5 'tot_NH4+' 1.0e-6
>     biomass = 'C5H7O2NNO3(s)'
>     bq = 4.0
>     chi = 1
>     /

where bq is $m\Delta G_{p}$, chi is $\chi$ and the direction is negative
if the reactants are written with negative stoichiometric coefficients.
The keyword biomass points to the immobile species with the biomass
concentration $B$. Because biomass has units of mol-biomass/bulk-volume,
the units of the rate constant $k_{\max}$ required by the input file
and/or database is in mol/mol-biomass/yr.

#### Entries in Input File

Aqueous kinetic reactions are entered in the input file using the
reaction name in the &Aqueous and &AqueousKinetics namelists. So for the
example for the irreversible reaction above, the **AQUEOUS_KINETICS**
keyword block would look like the following:

>     AQUEOUS_KINETICS
>     Am+++_decay
>     Pu++++241_decay
>     END

If the user desires, the label and the reaction rate constant associated
with any one reaction can be given explicitly in the input file,
overwriting thus the database value of the latter, and pointing to a
specific rate law in case more than one &AqueousKinetics namelist are
available for this reaction:

>     AQUEOUS_KINETICS
>     Am+++_decay -label default -rate 4.83E-2
>     Pu++++241_decay
>     END

The overall stoichiometry of microbially mediated reactions can vary
depending on the how much energy is assigned to the catabolic
(energy-yielding) and anabolic (biomass growth) pathways (Molins et al,
2015). The overall stoichiometry is constructed by adding the catabolic
(energy) and anabolic (biomass growth) pathways weighted by the portion
of electron equivalents used in that pathway (Molins et al, 2015):

$$0 = \sum_{k}^{N_{c}}\nu_{jk}^{e}A_{k}$$

$$0 = \sum_{k}^{N_{c}}\nu_{jk}^{s}A_{k}$$

$$0 = \sum_{k}^{N_{c}}{({f_{e}^{0}\nu}_{jk}^{e}} + f_{s}^{0}\nu_{jk}^{s})A_{k} = \sum_{k}^{N_{c}}{\nu_{jk}^{b}A}_{k}$$

To build the overall stoichiometry of the reaction, the name and portion
of electron-equivalents of each pathway must be given. The names must
correspond each to an &Aqueous namelist, while the name of the reaction
must correspond to an &AqueousKinetics namelist. From the example given
above:

>     AQUEOUS_KINETICS
>     AceNO3HCO3NO2 -pathway AceNO3HCO3NO2 0.60 -pathway C5H7O2N_RCH2_Ace_NH4 0.40 -rate 1261.44
>     END

where -pathway is followed by the name of the &Aqueous namelist
containing the stoichiometry, and the portion of electron equivalents
assign to that pathway. For consistency, reaction stoichiometries in the
database must be derived per each one electron equivalent transferred in
the reaction, and thus the portions assigned to catabolic and anabolic
pathways must add to one. Using the same name for a reaction and a
pathway is permitted (as in the reaction above). As result, if no
anabolic pathway is considered for the reaction, then one can simply
write

>     AQUEOUS_KINETICS
>     AceNO3HCO3NO2 -rate 1261.44
>     END