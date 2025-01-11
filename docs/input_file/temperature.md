### TEMPERATURE {#temperature-1 .SteefelHeading3}

This keyword block is used to set various temperature parameters. A full
calculation of the temperature will be added soon, but in the meantime,
temperature can be set block by block using individual geochemical
conditions, or by reading temperature values from a file.

#### **Set_temperature** {#set_temperature .SteefelHeading4}

Keyword followed by a value giving the default temperature.

Syntax:  **set_temperature** *value*

<u> Default </u>:  25ºC

Explanation:  The keyword parameter *set_temperature* is
used to set a global value for the temperature. The value provided here
can be overwritten by specifications in individual geochemical
conditions, provided these conditions are distributed in space as
initial conditions.

#### **Temperature_gradient** {#temperature_gradient .SteefelHeading4}

Keyword followed by a value giving the temperature gradient.

Syntax:  **temperature_gradient** *value*

<u> Default </u>:  0.0

Explanation:  The keyword parameter *temperature_gradient*
can be used to set a temperature gradient in units of ºC/m. A positive
value means that the temperature gradient will increase from grid point
1 to grid point NX. The temperature gradient, when input in this way,
operates only in the X direction.

#### **Read_TemperatureFile** {#read_temperaturefile .SteefelHeading4}

Option to specify a file with temperatures to be read.

Syntax:  **read_temperaturefile** *filename format*

*filename* gives the name of the file(up to 132 characters) contining
the liquid saturation distributed over the entire spatial domain.

<u> Default </u>: None.

Explanation:  This provides the name of a file containing
temperatures distributed in space to be read. This option supersedes all
other temperature specifications. The format of the file depends on the
optional format specified (see Section 7.1.2). If no format is
specified, the code assumes a single column (SingleColumn) of values in
the order 1-NX, 1-NY, 1-NZ:

> DO jz = 1,nz\
> DO jy = 1,ny\
> DO jx = 1,nx\
> READ(52,\*) t(jx,jy,jz)\
> END DO\
> END DO\
> END DO