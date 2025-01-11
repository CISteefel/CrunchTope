### PEST {#pest .SteefelHeading3}

#### **CreatePestInstructionFile** {#createpestinstructionfile .SteefelHeading4}

Keyword followed by *true* or *false* which selects whether or not to
create a PEST instruction (*PestExchange.ins*) file.

Syntax:  **CreatePestInstructionFile** *logical*

<u> Default </u>: *False*

Explanation:  This logical determines whether a
*PestExchange.ins* file is created with instructions for PEST to read
output on cation exchange from CrunchFlow. The format specified within
the instruction file will match the format used by CrunchFlow, with the
cation exchange species listed in the order given the *exchange* keyword
within the PEST keyword block. Once created, this option should be
turned off so as to avoid writing over the *PestExchange.ins* file.

#### **CreatePestExchangeFile** {#createpestexchangefile .SteefelHeading4}

Keyword followed by a file name*.*

Syntax:  **CreatePestInstructionFile** *filename*

<u> Default </u>: *PestExchange.out*

Explanation:  This keyword is used to set the filename to
be used for exchange output used in PEST runs. This option is useful
where multiple input files are being run together, thus requiring
different filenames for the output.

#### **Exchange** {#exchange .SteefelHeading4}

Keyword followed by concentration units followed by a list of exchange
species to be written out to a file PestExchange.out.

[Syntax:]{.underline} **Exchange** *units 'exchange species list'*

<u> Default </u>: *None*

[Explanation:]{.underline} This keyword indicates the list of species
participating in exchange reactions (e.g., Na^+^, Ca^2+^, Mg^2+^, ...)
to be written out to the file *PestExchange.out*, which can then be used
as "data" or "observations" in a PEST run using CrunchFlow. If the
keyword *CreatePestInstructionFile* is also set as true, then the format
for the exchange species output will be written to a PEST instruction
file: *PestExchange.ins*. The exchange species list is preceded by a
specification of the units, which may include:

  -----------------------------------------------------------------------
  **Units**  **Alternate Designation**
  ----------------- -----------------------------------------------------
  mol/g  mole/g, moles/g

  mmol/g  mmole/g, mmoles/g

  umol/g  umole/g, umoles/g

  equiv/g  eq/g, equivalents/g

  mequiv/g  milliequiv/g, meq/g, milliequivalents/g

  uequiv/g  microequiv/g, ueq/g, microequivalents/g
  -----------------------------------------------------------------------