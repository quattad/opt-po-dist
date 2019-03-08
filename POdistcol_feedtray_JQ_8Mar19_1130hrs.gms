$title Determining Optimal Feed Locations and Feed Flow Rates for 12-Stage PO Distillation Column
$sTitle System --- Multicomponent (Methanol, water, PO ) PO as light key and Methanol as heavy key at 1.8bar - 2.07bar.

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

function h_liq   /NRTLideal.h_liq /;
function h_vap   /NRTLideal.h_vap /;

* ----------------------------- STEP 2 : DEFINE PHYSICAL PROPERTIES -----------------------------
*************************************************************************************************************
$ontext
DECLARE COMPONENT SET
$offtext

Set j 'components' / methanol, water, PO /;

*************************************************************************************************************
$ontext
TRAY, REBOILER, CONDENSER TEMPERATURES
Values obtained from HYSYS simulation
$offtext


Scalar
         treb 'guess temperature for reboiler [K]'
         tbot 'guess temperature for bottom tray [K]'
         ttop 'guess temperature for top tray [K]'
         tcon 'guess temperature for condenser [K]'
         tavg 'avg temp';

treb = 359.6;
tbot = 357.5;
ttop = 330.8;
tcon = 325.5;
tavg = 344.15;
*************************************************************************************************************

*************************************************************************************************************
$ontext
FEED PROPERTIES
Values obtained from HYSYS simulation
$offtext

Scalar
         f          'total number of moles of feed [kmol/h]'
         tf         'temperature of feed [K]'
         pf         'pressure of feed [bar]'
         vf         'vapor fraction in feed [unitless]'
         shf        'specific enthalpy in feed []'
         flow_max   'estimate of maximum flows within system [kmol/h]'
         kmol2mol   'conversion factor for kmol to mol[mol/kmol]'         ;

kmol2mol = 1e3;

f = 6510*kmol2mol;
tf = 388; 
pf = 6.48;
vf = 0;
flow_max = 8000*kmol2mol;

Parameter xf(j) 'mole fractions in feed stream'
/
         methanol 0.75
         water 0.15,
         PO 0.10         / ;
*************************************************************************************************************

*************************************************************************************************************
$ontext
TRAY, REBOILER, CONDENSER PRESSURES
Values obtained from HYSYS simulation
$offtext
Scalar
         preb    'reboiler pressure [bar]'
         pbot    '1st tray pressure [bar]'
         ptop    'N-1 tray pressure [bar]'
         pcon    'condenser pressure [bar]'
         pavg    'avg pressure [bar] between 1st and N-1th tray';

preb = 2.07;
pbot = 1.97;
ptop = 1.90;
pcon = 1.80;
pavg = 2.00;
*************************************************************************************************************


* ----------------------------- STEP 2 : MODELLING DISTLLATION COLUMN -----------------------------

Set
  i              'stages' / 1 * 12 /
  reb(i)         'reboiler'
  con(i)         'condenser'
  col(i)         'stages in the column excluding reboiler and condenser'
  floc(i)        'possible locations for feed stage' / 2 * 11 /
  above_feed(i)  'stages above feed'
  below_feed(i)  'stages below feed'

Parameter
  ifeed          'initial guess for feed stage' / 8 /;

$ontext
DECLARING STAGE NUMBERING
assigns first element in dynamic set to reboiler, last element in dynamic set to condenser and all other elements in dynamic set to the columns
$offtext
  reb(i) = yes$(ord(i) = 1);
  con(i) = yes$(ord(i) = card(i));
  col(i) = yes - (reb(i) + con(i));
  below_feed(i) = yes$( ord(i) < ifeed ) - reb(i);
  above_feed(i) = yes$( ord(i) > ifeed ) - con(i);

$ontext
DECLARING STAGE PRESSURES
$offtext
Parameter
  p(i) 'pressure prevailing in tray i [bar]';

  p(i)$reb(i) = preb;
  p(i)$con(i) = pcon;
  p(i)$col(i) = pavg;

Positive Variable
  x(i,j)        'mole fraction of component j in liquid on ith tray'
  y(i,j)        'mole fraction of component j in vapor on ith tray'
  l(i)          'molar flow rate of liquid leaving tray i [mol/h] after conv. factor'
  v(i)          'molar flow rate of vapor leaving tray i [mol/h] after conv. factor'
  t(i)          'temperature of tray i [K]'
  feed(i)       'feed stream entering tray i [mol/h] after conv. factor'
  r             'reflux ratio'
  p1            'top product rate [mol/h] after conv. factor'
  p2            'bottom product rate [mol/h] after conv factor'

Variable
  hl(i)         'molar specific enthalpy of liquid on tray i [J/mol.K * 1e8; check hscale]'
  hv(i)         'molar specific enthalpy of vapor on tray i [J/mol.K * 1e8; check hscale]'
  hfeed(i)      'molar specific enthalpy of feed [J/mol.K * 1e8; check hscale]'

  qreb          'reboiler duty [J/s or W]'
  qcon          'condenser duty [J/s or W]'

Variable
  objval        'objective function variable to be optimized'       ;

Equation
  phe1(i)       'phase equilibrium relation for methanol'
  phe2(i)       'phase equilibrium relation for water'
  phe3(i)       'phase equilibrium relation for PO'

  errk1(i)      'phase equilibrium check for liquid'
  errk2(i)      'phase equilibrium check for vapour'

  tmb(i)        'total material balance for entire column'
  tmbreb(i)     'total material balance for reboiler'
  tmbcon(i)     'total material balance for condenser'

  cmb(i,j)      'component material balance (1 < i < n)'
  cmbreb(i,j)   'component material balance for reboiler'
  cmbcon(i,j)   'component material balance for condenser'

  defln         'definition of l(n) i.e. l out of condenser'
  defp2(i)      'definition of p2 i.e. l out of reboiler'
  defvn(i)      'definition of v(n) i.e. v out of condenser, which should be zero'

  defhl(i)      'definition of hl(i)'
  defhv(i)      'definition of hv(i)'
  defhfeed      'definition of hfeed(i)'
  
  eb(i)         'enthalpy balance on column'
  ebcon(i)      'enthalpy balance on condenser i.e. condenser duty'
  ebreb(i)      'enthalpy balance on reboiler i.e. reboiler duty'
  
  purcon        'purity constraints'
  sumf          'sum of feeds; NLP approach'

  pofr          'PO flow rate constraint'

  obj           'objective function to minimize'                 ;

Scalar
  hscale        'scaling factor for enthalpy' / 1e8 /;

* EQUATION DECLARATION AND DEFINITIONS
* (2) PHASE EQUILIBRIUM RELATIONS fij,v = fij,l
phe1(i) ..              f1l(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) =e= f1v(t(i), p(i), y(i,'methanol'), y(i,'water'), y(i,'PO'));
phe2(i) ..              f2l(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) =e= f2v(t(i), p(i), y(i,'methanol'), y(i,'water'), y(i,'PO')); 
phe3(i) ..              f3l(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) =e= f3v(t(i), p(i), y(i,'methanol'), y(i,'water'), y(i,'PO'));

* (3) PHASE EQUILIBRIUM NORMALIZATION sum(xij) = sum(yij) = 1
errk1(i) ..             sum(j, x(i,j)) =e= 1;
errk2(i) ..             sum(j, y(i,j)) =e= 1;

* (4) TOTAL MATERIAL BALANCES
defp2(i)$reb(i) ..      l(i) =e= p2;
defvn(i)$con(i) ..      v(i) =e= 0;
defln(i)$con(i) ..      l(i) =e= r*p1;

tmb(i)$col(i)..         l(i+1) + v(i-1) + feed(i)$floc(i) =e= l(i) + v(i);
tmbreb(i)$reb(i)..      l(i+1)                            =e= l(i) + v(i);
tmbcon(i)$con(i) ..     v(i-1)                            =e= l(i)+ p1   ;

* (5) COMPONENT MATERIAL BALANCES
cmb(i,j)$col(i) ..      l(i+1)*x(i+1,j) + v(i-1)*y(i-1,j) + (feed(i)*xf(j))$floc(i) =e= l(i)*x(i,j) + v(i)*y(i,j);
cmbreb(i,j)$reb(i) ..   l(i+1)*x(i+1,j)                                             =e= l(i)*x(i,j) + v(i)*y(i,j);
cmbcon(i,j)$con(i) ..   v(i-1)*y(i-1,j)                                             =e= (l(i)+ p1)*x(i,j)        ;

* (6) ENTHALPY BALANCES
defhl(i)..              hl(i)       =e= h_liq(t(i), p(i), x(i,'methanol'), x(i,'water'), x(i,'PO')) / hscale    ;
defhv(i)..              hv(i)       =e= h_vap(t(i), p(i), y(i,'methanol'), y(i,'water'), y(i,'PO')) / hscale    ;
defhfeed(i)$floc(i)..   hfeed(i)    =e= h_liq(tf, pf, 0.75, 0.15, 0.10) / hscale                                ;

eb(i)$col(i) ..         l(i+1)*hl(i+1) + v(i-1)*hv(i-1) + (feed(i)*hfeed(i))$floc(i) =e= l(i)*hl(i) + v(i)*hv(i);
ebreb(i)$reb(i)..       l(i+1)*hl(i+1) + qreb                                        =e= l(i)*hl(i) + v(i)*hv(i);  
ebcon(i)$con(i)..       v(i-1)*hv(i-1)                                               =e= (l(i)+ p1)*hl(i) + qcon;
 
* (7) FEED CONSTRAINTS
sumf..                  sum(i$floc(i), feed(i)) =e= f;

* (8) PURITY CONSTRAINTS
purcon..                x('12', 'PO') =g= 0.75;

* (9) PO FLOW RATE CONSTRAINT
pofr..                  x('12', 'PO') * p1 =g= 650*kmol2mol;

* DECLARING OBJECTIVE FUNCTION
obj..       objval =e= r;

* PARAMETER DECLARATIONS
* (ord(i) - card(i)) / (1 - card(i)) : 1 --> 0 as i increases.
* Liquid Mole Fractions
x.up(i,j) = 1.0;
x.lo(i,j) = 0;
x.l(i,'methanol') = 0.75;
x.l(i,'water') = 0.15;
x.l(i,'PO') = 0.10;

* Vapor Mole Fractions
y.up(i,j) = 1.0;
y.lo(i,j) = 0;
y.l(i,'methanol') = 0.75;
y.l(i,'water') = 0.15;
y.l(i,'PO') = 0.10;

* Reflux Ratio r
r.l      = 6;
r.up     = 30;
r.lo     = 1;

* Condenser Duty qcon
qcon.l      = 1.468e8;

* Reboiler Duty qreb
qreb.l      = 1.209e8;

* Top Product Molar Flow Rate
p1.l     = 900*kmol2mol;
p1.up    = flow_max;
p1.lo    = 0;

* Bottom Product Molar Flow Rate 
p2.l     = 5500*kmol2mol;
p2.up    = flow_max;
p2.lo    = 0;

* Vapor Flow Rate
v.l(i)    = ((r.l + 1)*p1.l);
v.up(i)   = flow_max;
v.lo(i)   = 0;

* Liquid Flow Rate
l.l(i)$reb(i)    = p2.l;
l.l(i)$(below_feed(i)) = p1.l*r.l + (1 - vf)*f;
l.l(i)$(above_feed(i)) = p1.l*r.l;
l.up(i) = flow_max;
l.lo(i) = 0;    

* Temperature
t.l(i)$reb(i)    = treb;
t.l(i)$con(i)    = tcon;
t.l(i)$col(i)    = tavg;

* Feed
feed.l(i)$floc(i)  = f/10;
feed.up(i)$floc(i) = f;

option nlp = conopt;

Model column 'ideal feed stages' / all /;
Solve column using nlp minimizing objval;