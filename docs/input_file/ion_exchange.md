### ION_EXCHANGE

An ion exchange reaction can be described via a mass action expression
with an associated equilibrium constant (Vanselow, 1932; Sposito, 1981;
Appelo and Postma, 1993). The exchange reaction can be written in
generic form as

$$vACl_{u}(aq) + uBX_{\nu}(s) \leftrightarrow uBCl_{\nu}(aq) + vAX_{u}(s)$$

where X refers to the exchange site occupied by the cations $A^{u+}$ and
$B^{v+}$. The equilibrium constant, $K_{eq}$, for this reaction can be
written as (Vanselow, 1932)

$$K_{eq} = \frac{(BCl_{\nu})^{u}(AX_{u})^{\nu}}{(ACl_{u})^{\nu}(BX_{\nu})^{u}}$$

where the parentheses () refer to the thermodynamic activities. Several
activity conventions are in wide use. One possibility is the
Gaines-Thomas activity convention, which assumes a reaction
stoichiometry of the following form (Appelo and Postma,1993), written
here assuming the $Cs^+$ is the relevant cation of interest

$$ Cs^{+} + \frac{1}{m}MX(i)_{m} \leftrightarrow CsX(i) + \frac{1}{m}M^{m +} $$

where *M* is the competing cation ($Na^{+}$, $K^{+}$, $Ca^{2+}$), *m* is its
charge, and *X(i)* refers to the *i*th type of exchange site. In the
Gaines-Thomas convention, each exchange site, *X(i)*, has a charge of
-1. The activities of adsorbed species correspond to the charge equivalent fractions

$$\beta(i)_{M} = \frac{ z_Mq(i)_M } { \sum_M {z_{M}q(i)_{M}} } = [X(i)_{M}]$$  

where $z^{+}$ is the charge of cation *M*, $q(i)_{M}$ is the concentration
of adsorbed cation *M* in exchange site *i* (moles/g), and the square
brackets denote activities. The Gapon activity convention is obtained by
writing the reactions in every case with a single exchanger (Appelo and
Postma, 1993). Alternatively, the Vanselow convention (Vanselow, 1932)
describes the exchanger activity with mole fractions

$$\beta(i)_{M} = \frac{q(i)_{M}}{\sum_{M}^{}{q(i)_{M}}} = [X(i)_{M}]$$

The exchange reactions can then be used to write a mass action equation
for binary Cs-M exchange:

$$K_{M/Cs} = \frac{ \beta(i)_{M}^{1/m} [Cs^{+}] } {\beta(i)_{Cs}[M^{m +}]^{1/m}} = \frac{[X(i)_{M}]^{1/m}[Cs^{+}] }{[X(i)_{Cs}][M^{m +}]^{1/m}}$$

In a single-site ion exchange model, the CEC is equal to the sum of the
charge equivalent concentrations of the adsorbed cations:

$$CEC = \sum_{M}^{}{z_{M}q_{M}}$$

while in a multi-site model, the CEC is the charge summed over all of
the cation exchange sites (Cernik et al., 1996; Voegelin et al., 2000)

$$CEC = \sum_{i}^{}{\sum_{M}^{}{z_{M}q(i)_{M}}}$$

In CrunchFlow, exchangers may be specified with the format

    exchange exchanger_name

where *exchanger_name* is the name of the exchanger given in the
database. Multiple exchangers can be listed. If the exchanger name is
not followed by anything, it is assumed that exchange is on the bulk
material, in which case a CEC is calculated from the combination of the
solid phase density, porosity, and liquid saturation. Alternatively, it is possible to specify exchange on a specific
mineral, in which case the CEC is calculated from the volume fraction of
the mineral, which may change with time, for example:

    exchange Xkaol- on Kaolinite

Also required is the keyword *convention*

    convention activity_convention

where activity_convention must be either *Gaines-Thomas* (of which Gapon
is a variant where the reaction is written with a single exchange site) or *Vanselow*.