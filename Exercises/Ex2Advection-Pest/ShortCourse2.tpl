ptf #
TITLE
ShortCourse2.in: Tracer test in one dimension
END
 
RUNTIME
time_units       years
timestep_max     0.010
timestep_init    0.010
time_tolerance   0.1
hindmarsh        true
correction_max   10.0
debye-huckel     false
master           Tracer
database_sweep   no
speciate_only    false
gimrt            false
graphics         kaleidagraph
database         datacom.dbs
screen_output    1000
END
 
OUTPUT
time_units           years
spatial_profile      50.0   
time_series          Tracer1.out 100 1 1
time_series_print    Tracer
time_series_units    Mol/kgw
time_series_interval 1
END

TRANSPORT
distance_units meters
time_units     years
fix_diffusion  0.0
dispersivity   #alpha   #
END

FLOW
time_units       years
distance_units   meters
constant_flow    1.0
END

BOUNDARY_CONDITIONS
x_begin   junk1     flux
x_end     junk2     flux
END

INITIAL_CONDITIONS
junk2   1-100
END

DISCRETIZATION
distance_units  meters
xzones 100 1.00
END

Condition junk1
units            mol/kg
temperature      25.0
Tracer           0.1
END

Condition junk2
units            mol/kg
temperature      25.0
Tracer           0.000001
set_porosity     1.0
END
 
POROSITY
fix_porosity     1.0
END
 
PRIMARY_SPECIES
Tracer
END
 

 
