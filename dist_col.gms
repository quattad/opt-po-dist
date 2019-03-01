$title Determining Theoretical Trays Required for Propylene Oxide (PO) Separation
$sTitle System --- Multicomponent with PO as light key and Methanol as heavy key at 2bar - 2.6bar.

$onText
This model determines the number of theoretical trays required for PO
separation based on the works of Viswanathan and Grossman.

Variable
         Number of components in feed, c
         Molar flow rate, f
         Feed temperature, tf
         Feed pressure, pf
         Molar specific enthalpy of feed, zfc
         Vector of mole fractions, zf = {zf1, zf2, ..., zfc}

         Molar flow rate for liquid leaving tray i, l(i)
         Mole fraction of liquid stream leaving tray i, x(i)
         Molar specific enthalpy of liquid stream leaving tray i, hl(i)
         Fugacity of component j in liquid stream leaving tray i, fl(i,j)

         Molar flow rate for vapor leaving tray i, v(i)
         Mole fraction of vapor stream leaving tray i, y(i)
         Molar specific enthalpy of vapor stream leaving tray i, hv(i)
         Fugacity of component j in vapor stream leaving tray i, fv(i,j)

         Top Product Rate, P1
         Bottom Product Rate, P2

         Upper bound of liquid and vapor flowrates, fmax

Positive Variable  // used for variables for which negative values are meaningless

Subsets
         Reboiler, R = {1}
         Condenser, C = {N}
         Column Trays, COL = {2, 3, ..., N-1}
Equation
         There are a total of 9 equations that need to be met that mainly
         correspond to the MESH equations.

         1) Phase Equlibrium
         2) Phase Equlibrium Error
         3) Total Material Balance
         4) Component Material Balances
         5) Enthalpy Balances
         6) Reflux single tray constraint
         7) Reboiler single tray constraint
         8) Pressure Profile
         9) Fugacity and enthalpy balances of each component on each tray
$offText

* ----------------------------- STEP 1 : IMPORT THERMODYNAMIC PACKAGE -----------------------------
* Include developed libraries in GAMS

$onecho > NRTLidealID.txt
ID1 1101
ID2 1921
ID3 1442
$offecho

$onecho > NRTLidealaij.txt
a12 610.4
a13 824.8
a21 -48.67
a31 -61.18
a23 1637
a32 1262
$offecho

$onecho > NRTLidealalphaij.txt
alpha12 0.3001
alpha13 0.2986
alpha23 0.4748
$offecho

$funclibin NRTLideal NRTLideal.dll

function f1l     /NRTLideal.f1_liq /;
function f2l     /NRTLideal.f2_liq /;
function f3l     /NRTLideal.f3_liq /;

function f1v     /NRTLideal.f1_vap /;
function f2v     /NRTLideal.f2_vap /;
function f3v     /NRTLideal.f3_vap /;

Set j 'components' / methanol, water, PO, PG /;


* ----------------------------- STEP 2 : DEFINE PHYSICAL PROPERTIES -----------------------------

$ontext
PHYSICAL PROPERTIES CONSTANTS
PG critical temp and pressure not pertinent
$offtext

Table prcon (j, *) 'physical properties'
                 mm      tb      tc      pc
  methanol    32.042  337.85  513.15  79.54
  water       18.015  373.15  647.35 220.89
  PO           58.08  307.38  482.25  49.24
  PG           76.09  460.15                            ;

$ontext
VAPOR PRESSURE CONSTANTS
all constants obtained from Perry's Chemical Engineering handbook except for PO and PG.
Pressures in Pa Formula for Perry: lnP = a + b/T + c*lnT + d*power(T,e)
PO obtained from 1966 paper in the format logP* = 7.658 - 1472/T. Applicable within range of 40-70oC.
$offtext
 
Table vpcon(j, *) "constants for vapor pressure"
                   a        b         c            d        e       tmin    tmax
  methanol    82.718  -6904.5   -8.8622    7.4664e-6        2     175.47   512.5
  water       73.649  -7258.2   -7.3037    4.1653e-6        2     273.16   647.096
  PO           7.658  -1472                                       
  PG           212.8  -15420    -28.109    2.1564e-5        2                       ;

$ontext
ISOBARIC HEAT CAPACITY CONSTANTS
all constants obtained from Perry's Chemical Engineering handbook except for PO.
Specific heat capacity for liquid phase, CpL = a + b*T +c*power(T,2) + d*power(T,3) + e*power(T,4)
For PO in J/mol.K : a + b*T + c*power(T,2)
$offtext

Table cpconl (j, *) "isobaric heat capacity constants for liquid phase"
                    a         b          c              d         e   
  methanol     105800   -362.23     0.9379      
  water        276370   -2090.1     8.125        -0.014116  9.3701e-6
  PO           113.08   -0.15085    6.7e-4
  PG           58080      445.2                                            ;

$ontext
Specific heat capacity for vapor phase, Cp,v = a + b*( (c/t) / (sinh(d/t)) ) + e*( (f / t) / cosh(f / t) )**2

For PO in J/mol.K : a + b*T + c*power(T,2)
Actually do I even need this?
$offtext

Table cpconv (j, *) "isobaric heat capacity constants for vapor phase"
                    a         b          c                 d            e   
  methanol 0.39252e-5   0.879e-5    1.9165e-3     0.53654e-5        896.7
  water    0.33363e-5  0.2679e-5    2.6105e-3     0.08896e-5         1169 
  PO           
  PG                                                       ;

$ontext
TRAY, REBOILER, CONDENSER TEMPERATURES
Define temperatures for the 1 and N-1 trays, reboiler, condenser
Values obtained from HYSYS simulation
$offtext

Scalar
         treb 'guess temperature for reboiler'
         tbot 'guess temperature for bottom tray'
         ttop 'guess temperature for top tray'
         tcon 'guess temperature for condenser';

treb = 359.6;
tbot = 357.5;
ttop = 330.9;
tcon = 326.1;

$ontext
FEED PROPERTIES
Values obtained from HYSYS simulation
$offtext

Scalar
         f       'total number of moles of feed'
         tf      'temperature of feed in K'
         pf      'pressure of feed stream in bar'
         vf      'vapor fraction in feed'
         shf     'specific enthalpy in feed'
         preb    'reboiler pressure'
         pbot    '1st tray pressure'
         ptop    'N-1 tray pressure'
         pcon    'condenser pressure';

f = 6282;
preb = 2.07;
pbot = 1.97;
ptop = 1.90;
pcon = 1.80;
pf = 18.48;
vf = 0;

Parameter xf(j) 'mole fractions in feed stream'
/
         methanol 0.7515,
         water 0.1292,
         PO 0.1073,
         PG 0.0102      / ;

* ----------------------------- STEP 3 : MODELLING DISTLLATION COLUMN -----------------------------
$ontext
Reboiler (1) ~ Condenser (N)
Estimated number of trays is 9 from Gilliand correlation

Sets are as follows.

         Stages, I = {1, 2... N}; N is estimated from Gilliand correlation.
         Possible locations for feed, i_feed = {i_feed + 1, i_feed + 2... N-1}
         Possible locations for reflux, REF = {i_r_min, i_r_min + 1... N -1}
         Possible locations for boil-up, BU = {2, 3..., i_b_max - 1, i_b_max}
         Feed Location = {i_feed}
         Stages above feed, AF = {i|i_feed < i <= (N-1)}; REF subset of AF
         Stages below feed, BF = {i|2<= i < i_feed}; BU subset of BF
         Corresponding index set, J = {1, 2, ... , c}

OPERATOR ord(i)
operator 'ord' returns relative position of member in a set

OPERATOR card(i)
card operator 'card' returns number of elements present in a set

OPERATOR $yes$
yes$ used to define dynamic sets. dynamic sets have elements that are added/removed during runtime

$offText

Set
  i              'stages' / 1*9 /
  reb(i)         'reboiler'
  con(i)         'condenser'
  col(i)         'stages in the column excluding reboiler and condenser'
  above_feed(i)  'stages above the feed stage, excluding feed stage i.e. possible locations for reflux'
  below_feed(i)  'stages below feed stage, including feed stage i.e. possible locations for reboiler'
  iref(i)       'stage at which reflux enters'
  ireb(i)       'stage at which reboiler re-entry stream enters';

Alias(i,i2);
Alias(above_feed, af);
Alias(below_feed, bf);

Parameter
  floc          'feed stage location' / 5 /;

*assigns first element in dynamic set to reboiler
  reb(i) = yes$(ord(i) = 1);

*assigns last element in dynamic set to condenser
  con(i) = yes$(ord(i) = card(i));

*assigns all other elements in dynamic set to the columns
  col(i) = yes - (reb(i) + con(i));

* assigns stages below feed / above feed to dynamic set
  below_feed(i) = yes$( ord(i) le floc ) - reb(i);
  above_feed(i) = yes - below_feed(i) - con(i);

* initial guess for reflux and reboiler stages? should be a set with 1 member each.

$ontext
DEFINE COLUMN PRESSURES
assigns the pressure on a particular tray depending on whether reboiler, condenser or
tray in column
$offtext

Parameter
  p(i) 'pressure prevailing in tray i';

  p(i)$reb(i) = preb;
  p(i)$con(i) = pcon;
  p(i)$col(i) = pbot - ( ((pbot - ptop) / card(i) - 1 - 2) ) * (ord(i) - 2);

Positive Variable
  x(i,j)    'mole fraction of component j in liquid on ith tray'
  y(i,j)    'mole fraction of component j in vapor on ith tray'
  l(i)      'molar flow rate of liquid leaving tray i'
  v(i)      'molar flow rate of vapor leaving tray i'
  t(i)      'temperature of tray i'
  feed(i)   'feed stream entering tray i'
  r         'reflux ratio'
  p1        'top product rate'
  p2        'bottom product rate'
  vref(i)   'amount of reflux entering tray i'
  lbu(i)    'amount of reboiler stream entering tray i'  

Variable
  hl(i)  'molar specific enthalpy of liquid on tray i'
  hv(i)  'molar specific enthalpy of vapor on tray i'

Binary Variable
  zref(i) 'for location of reflux. 1 if i is tray where reflux enters, 0 otherwise.'
  zbu(i) 'for location of reboiler stream. 1 if i is tray where reboiler stream enters, 0 otherwise.'

Equation
  phe1(i),phe2(i), phe3(i)        'phase equilibrium relations'
  errk(i)       'phase equilibrium error function'

  tmbc(i)       'total material balance for entire column'
  tmbf(i)       'total material balance for feed stage'
  tmbref(i)     'total material balance for possible candidates for reflux'
  tmbreb(i)     'total material balance for possible candidates for boil-up re-entry'
  tmb(i)        'total material balance for trays'
  tmbl(i)       'total material balance for reboiler'
  tmbn(i)       'total material balance for condenser'
  defln(i)      'definition of l(n) i.e. l out of condenser'
  defp2(i)      'definition of p2 i.e. l out of reboiler'
  defvn(i)      'definition of v(n) i.e. v out of condenser, which should be zero'
  defref(i)     'definition of reflux re-entering column'

  cmb(i,j)      'component material balance (1 < i < n)'
  cmbf(i,j)     'component material balance for feed stage'
  cmbref(i,j)   'component material balance for possible stages for reflux'
  cmbreb(i,j)   'component material balance for possible stages for boilup reentry'
  cmbl(i,j)     'component material balance on first tray'
  cmbn(i,j)     'component material balance on nth tray'

  defhl(i)      'definition of hl(i)'
  defhv(i)      'definition of hv(i)'

  eb(i)         'enthalpy balance'
  purcon        'purity constraints'
  sumf          'sum of feeds'                                                      ;

* (2) PHASE EQUILIBRIUM RELATIONS fij,v = fij,l
phe1(i) ..    f1l(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) =e= f1v(t(i), p(i), y(i,'methanol'),y(i,'water'),y(i,'PO')); 
phe2(i) ..    f2l(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) =e= f2v(t(i), p(i), y(i,'methanol'),y(i,'water'),y(i,'PO')); 
phe3(i) ..    f3l(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) =e= f3v(t(i), p(i), y(i,'methanol'),y(i,'water'),y(i,'PO')); 

$ontext
Attempted to use VLE equation here, but may not be necessary after including thermodynamic package?
phe(i,j) ..             y(i,j)*p(i)
                        - x(i,j)*exp(
                            vpcon(j, 'a')
                            + (vpcon(j,'b') / t(i))
                            + (vpcon(j,'c')*ln(t(i)))
                            + (vpcon(j,'d')*t(i)*exp(vpcon(j,'e')))
                            ); 
$offtext

* (3) PHASE EQUILIBRIUM ERROR
errk(i) ..               sum(j, x(i,j)) - sum(j, y(i,j)) =e= 0;

* (4) TOTAL MATERIAL BALANCES
* There are a total of 11 material balance equations 1) Entire column
* 2) Reflux 3) After the feed 4) At the feed location 5) Below the feed
* 6) Liquid flow rate from first tray 7) Vapor flow rate from first tray
* 8) Reflux balance 9) Liquid flow rate out of reboiler

* TMB 4a For Entire Column
tmbc(i)$col(i)..            v(i-1) - (sum(i2, vref(i)) + l(i) + p1) =e= 0;

* TMB 4b & 4c for stages above the feed excluding condenser, may / may not include reflux stage
tmbref(i)$above_feed(i)..   l(i) + v(i) - l(i+1)- v(i-1) - vref(i)$iref(i) =e= 0;

* TMB 4d for feed stage balance
tmbf(i) ..                   l(i) + v(i) - l(i+1) -v(i-1) - feed(i)$floc =e= 0;

* TMB 4e,f for stages below feed excluding reboiler, may / may not include boilup re-entrance
tmbreb(i) ..                l(i) + v(i) - l(i+1)- v(i-1) - lbu(i)$ireb(i) =e= 0;

* TMB 4g for constraint on reboiler.
tmbl(i)$reb(i) ..           l(i) + v(i) + sum(bf, lbu(i))- l(i+1) =e= 0;

* TMB 4h for constraint on second product stream = liquid flow rate from reboiler
defp2(i)$reb(i) ..          l(i) - p2 =e= 0;

* TMB 4i for constraint on condenser. no vapor outflow = total condenser
defvn(i)$con(i) ..          v(i) =e= 0;

$ontext
* TMB condenser. no vapor outlet stream. may serve same purpose as 3i, check!
tmbn(i)$con(i) ..           l(i) + p1 - v(i-1) =e= 0;
$offtext

* TMB 4j for constraint on reflux streams back into column. m.b. on reflux = RR * top product
defref(i)$above_feed(i) ..  sum(i2, vref(i)) - (r * p2) =e= 0 ;

* TMB 4k Constraint on liquid flow rate out of condenser to reflux ratio of total product stream
defln(i)$con(i) ..          l(i) - (r * p1) =e= 0;


* (5) COMPONENT MATERIAL BALANCES
* CMB 5a Entire Column
$ontext
Stock Equation
cmb(i,j)$col(i) ..       l(i)*x(i,j) + v(i)*y(i,j) - l(i+1)*x(i+1, j) - v(i-1)*y(i-1,j) - (feed(i)*xf(j))$floc =e= 0;
$offtext

cmb(i,j)$col(i) ..          v(i-1)*y(i-1,j) - (sum(i2,vref(i)) + l(i) + p1) * x(i,j) =e= 0;

* CMB 5b & 5c for stages above the feed excluding condenser, may / may not include reflux stage
* how to settle liq comp for re-entry
cmbref(i,j)$above_feed(i)..   l(i)*x(i,j) + v(i)*y(i,j) - l(i+1)*x(i,j) - v(i-1)*y(i,j) - (vref(i)*x(i,j))$iref(i) =e= 0;

* CMB 5d for feed stage
cmbf(i,j) ..                   l(i)*x(i,j) + v(i)*y(i,j) - l(i+1)*x(i,j) - v(i-1)*y(i,j) - (feed(i)*xf(j))$floc =e= 0;

* CMB 5e,f for stages below feed excluding reboiler, may / may not include boilup re-entrance
cmbreb(i,j) ..                l(i)*x(i,j) + v(i)*y(i,j) - l(i+1)*x(i,j) - v(i-1)*y(i,j) - (lbu(i)*x(i,j))$ireb(i) =e= 0;

* CMB 5g for constraint on reboiler.
cmbl(i,j)$reb(i) ..           l(i)*x(i,j) + v(i)*y(i,j) + (sum(bf, lbu(i)))*y(i,j) - l(i+1)*x(i+1,j) =e= 0;

$ontext
Unnecessary for now
* Component Material for Reboiler
cmbl(i,j)$reb(i) ..         l(i)*x(i,j) + v(i)*y(i,j) - l(i+1)*x(i+1,j) =e= 0;

* Component Material Balance for Condenser
cmbn(i,j)$con(i) ..      l(i)*x(i,j) + p1*x(i,j) - v(i-1)*y(i-1,j) =e= 0;
$offtext

* (6) ENTHALPY BALANCES
* Define equations to calculate enthalpy balances for vapor and liquid on each tray
* Specific heat capacity, CpL = a + b*T +c*power(T,2) + d*power(T,3) + e*power(T,4)
* For PO in J/mol.K : a + b*T + c*power(T,2)

$ontext
defhl(i) ..             hl(i) -( h_liq(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) ) =e= 0;
defhv(i) ..             hv(i) -( h_vap(t(i), p(i), y(i,'methanol'), y(i,'water'), y(i,'PO')) ) =e= 0;
$offtext

$ontext
may not need this, can use thermodynamic librari
defhl(i) ..             hl(i) - (cpconl(j,'a') + cpconl(j,'b')*t(i) + cpconl(j,'c')*(t(i)**2) + cpconl(j,'d')*(t(i)**3) + cpconl(j,'e')*(t(i)**4));                   
defhv(i) ..                
$offtext

* EB 6a,6b for possible candidates for reflux, may / may not include reflux stage

eb(i)$col(i) ..          l(i)*hl(i) + v(i)*hv(i) - l(i+1)*hl(i+1) - v(i-1)*hv(i-1) - (feed(i)*shf)$floc =e= 0;

* EB 6c for feed stage

* EB 6d, 6e for possible candidates for boilup reentry
 
* Reflux entering only on one tray

* Reboiled vapor entering only on one tray

* Pressure profiles

* Purity Constraint on PO
purcon .. x('9', 'PO') =g= 0.85;

* assigns lower, upper bounds and initial conditions to
* reflux ratio, top product rate and bottom product rate

r.l      = 6;
r.up     = 8;
r.lo     = 5;

p1.l     = 2000;
p1.up    = 3000;
p1.lo    = 4000;

p2.l     = 4000;
p2.up    = 6000;
p2.lo    = 7000;

* set upper limits for mole fractions

x.up(i,j) = 1.0;
y.up(i,j) = 1.0;

$onText
* Additional variables / equations for feed tray location problem;
Binary Variable yf(i);
$offText

Variable zf;

Model column 'ideal number of stages' / all /;
Solve column using minlp maximizing zf;
