### SURFACE_COMPLEXATION

Surface complexation follows the approach outlined in Dzombak and Morel
(1990), with either a double layer or non-electrostatic model possible.
Currently, complexation must be on a specific mineral, so a valid
mineral name (listed in the **MINERALS** keyword block) must be given,
for example:

    >FeOH_strong on Fe(OH)3
    >FeOH_weak on Fe(OH)3

To specify a non-electrostatic model, the mineral name should be
followed by the keyword *--no_edl*, for example:

    >FeOH_strong on Fe(OH)3 -no_edl
    >FeOH_weak on Fe(OH)3 -no_edl

The capability for surface complexation on the bulk material will be
added soon.