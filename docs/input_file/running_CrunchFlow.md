# Explanation

CrunchFlow reads a user-provided input file on startup which provides the necessary physical and chemical parameters needed for a run.
The input file, the name of which is specified by the user, is keyword-based so that the order of appearance does not matter.
In many cases, it is possible to leave various optional parameters out altogether and thus allowing the code to use default parameters.
Keywords are grouped broadly into keyword blocks, which in turn include a variety of keyword parameters.
Keyword blocks may appear anywhere in the input file and keyword parameters may appear anywhere within a particular block, but certain keyword parameters must appear in the appropriate blocks.
The possible keyword blocks include:

```
TITLE
RUNTIME
OUTPUT
PRIMARY_SPECIES
SECONDARY_SPECIES
GASES
MINERALS
AQUEOUS_KINETICS
ION_EXCHANGE
SURFACE_COMPLEXATION
PEST
DISCRETIZATION
INITIAL_CONDITIONS
BOUNDARY_CONDITIONS
TRANSPORT
FLOW
POROSITY
TEMPERATURE
CONDITION
```
 
With the exception of the keyword block CONDITION, the blocks should appear only once in the input file.
Each keyword block is terminated by an END.
The keyword block CONDITION is a special case in that it can occur multiple times.
Each occurrence of CONDITION specifies a separate geochemical condition (these may be boundary or initial conditions or source terms) containing the geochemical input needed to describe a particular problem.
The various keyword blocks will be discussed individually in more detail below.


# Reading Input File Name from a File
Normally, CrunchFlow will prompt the user to enter the name of an input file interactively.
The user can provide the name of the input file, however, by including it in a file names PestControl.ant.
The code will check to see if this file exists in the directly from which CrunchFlow is invoked, and if it does, will attempt to read the file name from it.
Upon successfully reading the input file name, the code will then check to see that this file exists before reading from it.
The use of the PestControl.ant to provide an input file name will then skip the requirement of interactive input from the user, an option that is particularly useful when running PEST.
