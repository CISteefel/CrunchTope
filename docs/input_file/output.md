# Output Block


The keyword block OUTPUT controls when and what will be written to disk.
If the OUTPUT block is not given, then no time stepping will occur (see keyword parameter spatial_profile).
The keyword parameters included in the OUTPUT keyword block are:

- time_units
- spatial_profile
- time_series
- time_series_print
- time_series_interval
- MakeMovie

The time_units keyword parameter has the same definition and role as in other keyword blocks where it appears.
It specifies the time units for such keyword parameters as spatial_profile and it determines the units used in the various output files.

Example:
```
OUTPUT
time_units				hours
spatial_profile		250
time_series			50
time_series			100
time_series_print		Cs+  pH
time_series_interval	1
END
```

In this example, a single spatial profile (field variables as a function of space) is written at 250 hours.
The code will time step until 250 hours is reached and then stop.
A total concentration of cesium will be tracked at JX=50 and JX=100 (nodes 50 and 100 in the X direction).
Output to the time series file will be every time step. Since this is a 1D problem, the Y and Z coordinates need not be given.


[Basic keywords](#basic-keywords)
- [spatial_profile](#spatial_profile)
- [time_series](#time_series)
- [time_series_print](#time_series_print)
- [time_series_interval](#time_series_interval)
- [time_series_output](#time_series_output)
- [MakeMovie](#makemovie)

## Basic keywords

### spatial_profile

#### Syntax
```
spatial_profile [time1] [time2  time3 � 
                         & timen]
```
[time1] [time2  time3 � & timen] are real numbers giving the time at which the spatial profile is to be written.
The ampersand (&) can be used to continue listing of output times on multiple lines.
Any one line can be up to 132 characters long.
If not specified, the code assumes that no output file is to be written and therefore no time stepping will occur.

#### Explanation
This keyword parameter instructs the code to output a spatial profile or �snapshot� at the specified time.
The time units are taken from the time_units value specified in the OUTPUT keyword block (if absent, time units of years are assumed).
An ampersand (&) can be used to continue listing of output times on the following line.
This important parameter, in addition to specifying when spatial profiles will be written, controls whether time stepping is carried out at all.
The code will run until the last output time in the list is reached, unless the parameter steady_state in the RUNTIME keyword block is turned on (in this case, the code runs until a steady-state is achieved or until the final spatial profile time is reached, whichever is first).
Field variables written out as part of the spatial profile are given in the Output Variables section.

#### Example

```
spatial_profile 0.001 2e-3 0.1 10 100 1000
```

```
spatial_profile 0.001 2e-3 0.1 10 100 200 300 &
                400 500 600 700 800 900 1000 &
                2000 3000
```


### time_series

#### Syntax
```
time_series [filename] [JX] [optionnal [JY JZ]]
```
for a 1D problem.
```
time_series [filename] [JX] [JY] [optionnal [JZ]]
```
for a 2D problem.
```
time_series [filename] [JX] [JY] [JZ]
```
for a 3D problem.
[filename] is the file name to be used for the time series, and [JX],[JY], and [JZ] are integers giving the node number for which the time series will be tracked.
There are no default values.

#### Explanation
This keyword parameter specifies that a time series should be recorded at the node specified.
Multiple specifications of time_series may occur, each one giving a different filename and node number.
For one-dimensional problems, the Y and Z coordinates are optional.
If in any case a node location is given that is greater than the coordinates specified in the DISCRETIZATION block, execution will terminate and an error statement will be written (e.g., JX > NX, the number of nodes in the X direction).

#### Example

For a 1D problem and three output files at different locations (nodes)
```
time_series MySpatialProfile1.out 1
time_series MySpatialProfile2.out 8 1
time_series MySpatialProfile3.out 11 1 1 
```

For a 2D problem and two output files at different locations (nodes)
```
time_series MySpatialProfile1.out 3 4
time_series MySpatialProfile2.out 8 1 1
```

For a 3D problem and three output files at different locations (nodes)
```
time_series MySpatialProfile10.out 8 4 2
time_series MySpatialProfile40.out 3 20 5
time_series MySpatialProfile05.out 30 2 12
```

### time_series_print

#### Syntax
```
time_series_print [species1] [optionnal [species2 species3 & 
                                      species4 ... speciesn]]
```
or
```
time_series_print all
```

[speciesn] is a character string corresponding to an aqueous species name.
If time_series_print is not specified, no tracking of solutes will be done.
If time_series_print is not followed by a species name, a default of all will be assumed.

#### Explanation
This keyword parameter gives the species to be tracked as a time series.
The time_series_print keyword may be followed with all, in which case all of the total concentrations will be written out in column form in units of mol/kg followed by all of the species concentrations in logarithmic form.
Listing of species names will cause the total concentrations of these species to be written out (not the individual species concentrations).
In addition, pH may be specified to be output.
Species names must match those in the primary species list (case sensitive).

#### Example

```
time_series_print  pH  Ca++  Na+  K+
```

will cause the code to track the time evolution of the solution pH and the total concentrations of Ca, Na, and K.

```
time_series_print  all
```
will cause the code to track the time evolution of all total concentrations, including pH.

### time_series_interval

#### Syntax
```
time_series_interval [interval]
```
[interval] is an integer (Default: 1).

#### Explanation
This keyword parameter gives the time step interval at which the time evolution of the species specified in time_series_print will be written out at the nodal location specified by time_series keyword.

#### Example

```
time_series_interval 10
```

### time_series_output

#### Syntax
```
time_series_output [time1] [time2  time3 � 
                         & timen]
```
[time1] [time2  time3 � & timen] are real numbers giving specific times at which time series is to be written.
The ampersand (&) can be used to continue listing of output times on multiple lines.
Any one line can be up to 132 characters long. Default: none.

#### Explanation
This keyword parameter gives the option of writing out time series at specific times.
This is useful when the objective is to compare discrete data (e.g., from a kinetic experiment) with output from the simulation.

#### Example

```
time_series_output 0.034 2.24 10.1
```

### MakeMovie

#### Syntax
```
MakeMovie [species1] [optionnal [species2 species3 & 
                                      species4 ... speciesn]] 
```
[speciesn] is a character string corresponding to an aqueous primary species name. 
Default: No movie file will be written.

#### Explanation
This keyword parameter provides the option to write out a total concentration for a primary species in 1D or 2D.
The interval at which data is written out is given by time_series_interval.

#### Example

```
MakeMovie  pH  Ca++  Na+  K+
```