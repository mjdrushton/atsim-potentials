import math

from atsim.potentials import Potential
from atsim.potentials import EAMPotential
from atsim.potentials import writeSetFLFinnisSinclair


#Embedding functions
def AlEmbedFunction(rho):
  if rho == 0.0: return 0.0
  return -math.sqrt(rho)+0.000093283590195398*rho**2-0.0023491751192724*rho*math.log(rho)

def FeEmbedFunction(rho):
  return -math.sqrt(rho) -  0.00067314115586063*rho**2 +  0.000000076514905604792*rho**4

#Density functions
def AlAlDensityFunction(r):
  funcs = [
    (2.5, lambda r: 0.00019850823042883 * (2.5 - r)**4 ),
    (2.6, lambda r: 0.10046665347629    * (2.6 - r)**4 ),
    (2.7, lambda r: 1.0054338881951E-01 * (2.7 - r)**4 ),
    (2.8, lambda r: 0.099104582963213   * (2.8 - r)**4 ),
    (3.0, lambda r: 0.090086286376778   * (3.0 - r)**4 ),
    (3.4, lambda r: 0.0073022698419468  * (3.4 - r)**4 ),
    (4.2, lambda r: 0.014583614223199   * (4.2 - r)**4 ),
    (4.8, lambda r: -0.0010327381407070 * (4.8 - r)**4 ),
    (5.6, lambda r: 0.0073219994475288  * (5.6 - r)**4 ),
    (6.5, lambda r: 0.0095726042919017  * (6.5 - r)**4 )]
  vals = [ func(r) for cutoff, func in funcs if r <= cutoff ]
  return sum(vals)

def FeFeDensityFunction(r):
  funcs = [
    (2.4, lambda r: 11.686859407970*(2.4    - r)**3),
    (3.2, lambda r: -0.014710740098830*(3.2 - r)**3),
    (4.2, lambda r: 0.47193527075943*(4.2   - r)**3)]
  vals = [ func(r) for cutoff, func in funcs if r <= cutoff ]
  return sum(vals)

def FeAlDensityFunction(r):
  funcs = [
    (2.4, lambda r: 0.010015421408039   * (2.4 - r)**4),
    (2.5, lambda r: 0.0098878643929526  * (2.5 - r)**4),
    (2.6, lambda r: 0.0098070326434207  * (2.6 - r)**4),
    (2.8, lambda r: 0.0084594444746494 * (2.8 - r)**4),
    (3.1, lambda r: 0.0038057610928282  * (3.1 - r)**4),
    (5.0, lambda r: -0.0014091094540309 * (5.0 - r)**4),
    (6.2, lambda r: 0.0074410802804324  * (6.2 - r)**4)]
  vals = [ func(r) for cutoff, func in funcs if r <= cutoff ]
  return sum(vals)

def zerowrap(wrapped):
  def f(r):
    if r == 0.0:
      return 0.0
    return wrapped(r)
  return f

#Pair potentials
def ppfuncAlAl(r):
  funcs = [
    ( (0.0, 1.6) , zerowrap(lambda r: (2433.5591473227/r)  * (0.1818 * math.exp(-22.713109144730* r) + 0.5099 * math.exp(-6.6883008584622* r) + 0.2802* math.exp(-2.8597223982536*r) + 0.02817* math.exp(-1.4309258761180*r)) )),
    ( (1.6, 2.25), lambda r: math.exp(6.0801330531321 - 2.3092752322555*r + 0.042696494305190*r**2 - 0.07952189194038*r**3) ),
    ( (2.25, 3.2), lambda r:  17.222548257633*(3.2 - r)**4 - 13.838795389103*(3.2 - r)**5  + 26.724085544227*(3.2 - r)**6 - 4.8730831082596*(3.2 - r)**7 + 0.26111775221382*(3.2 - r)**8),
    ( (2.25, 4.8), lambda r: -1.8864362756631*(4.8 - r)**4 + 2.4323070821980*(4.8 - r)**5 - 4.0022263154653*(4.8 - r)**6 + 1.3937173764119*(4.8 - r)**7 - 0.31993486318965*(4.8 - r)**8),
    ( (2.25, 6.5), lambda r: 0.30601966016455*(6.5 - r)**4 - 0.63945082587403*(6.5 - r)**5 + 0.54057725028875*(6.5 - r)**6 - 0.21210673993915*(6.5 - r)**7 + 0.032014318882870*(6.5 - r)**8) ]

  vals = [ func(r) for ((lowcut, highcut), func) in funcs if lowcut <= r < highcut ]
  return sum(vals)

def ppfuncAlFe(r):
  funcs = [
    ( (0.0, 1.2), zerowrap(lambda r: (4867.1182946454/r) * (0.1818 * math.exp(-25.834107666296*r) + 0.5099*math.exp(-7.6073373918597*r) + 0.2802*math.exp(-3.2526756183596*r) + 0.02817*math.exp(-1.6275487829767*r))) ),
    ( (1.2, 2.2), lambda r: math.exp(6.6167846784367-1.5208197629514*r - 0.73055022396300*r**2 - 0.038792724942647*r**3) ),
    ( (2.2, 3.2), lambda r: -4.1487019439249*(3.2 - r)**4 + 5.6697481153271*(3.2 - r)**5 - 1.7835153896441*(3.2 - r)**6 - 3.3886912738827*(3.2 - r)**7 + 1.9720627768230*(3.2 - r)**8),
    ( (2.2, 6.2), lambda r: 0.094200713038410*(6.2 - r)**4 -0.16163849208165*(6.2 - r)**5 + 0.10154590006100*(6.2 - r)**6 -0.027624717063181*(6.2 - r)**7 + 0.0027505576632627*(6.2 - r)**8) ]
  vals = [ func(r) for ((lowcut, highcut), func) in funcs if lowcut <= r < highcut ]
  return sum(vals)

def ppfuncFeFe(r):
  funcs = [
  ( (0.0, 1.0), zerowrap(lambda r: (9734.2365892908/r)*(0.1818*math.exp(-28.616724320005*r) + 0.5099*math.exp(-8.4267310396064*r) + 0.2802*math.exp(-3.6030244464156*r) + 0.02817*math.exp(-1.8028536321603*r))) ),
  ( (1.0, 2.05), lambda r: math.exp(7.4122709384068 - 0.64180690713367*r - 2.6043547961722*r**2 + 0.62625393931230*r**3) ),
  ( (2.05, 2.2), lambda r: -27.444805994228*(2.2 - r)**3 ),
  ( (2.05, 2.3), lambda r: 15.738054058489*(2.3 - r)**3),
  ( (2.05, 2.4), lambda r: 2.2077118733936*(2.4 - r)**3),
  ( (2.05, 2.5), lambda r: -2.4989799053251*(2.5 - r)**3),
  ( (2.05, 2.6), lambda r: 4.2099676494795*(2.6 - r)**3),
  ( (2.05, 2.7), lambda r: -0.77361294129713*(2.7 - r)**3),
  ( (2.05, 2.8), lambda r: 0.80656414937789*(2.8 - r)**3),
  ( (2.05, 3.0), lambda r: -2.3194358924605*(3.0 - r)**3),
  ( (2.05, 3.3), lambda r: 2.6577406128280*(3.3 - r)**3),
  ( (2.05, 3.7), lambda r: -1.0260416933564*(3.7 - r)**3),
  ( (2.05,4.2), lambda r: 0.35018615891957*(4.2 - r)**3),
  ( (2.05, 4.7), lambda r: -0.058531821042271*(4.7 - r)**3),
  ( (2.05, 5.3), lambda r: -0.0030458824556234*(5.3 - r)**3) ]
  vals = [ func(r) for ((lowcut, highcut), func) in funcs if lowcut <= r < highcut ]
  return sum(vals)

def main():
  # Define list of pair potentials
  pairPotentials = [
    Potential('Al', 'Al', ppfuncAlAl),
    Potential('Al', 'Fe', ppfuncAlFe),
    Potential('Fe', 'Fe', ppfuncFeFe)]

  # Assemble the EAMPotential objects
  eamPotentials = [
    #Al
    EAMPotential('Al', 13, 26.98154, AlEmbedFunction,
      { 'Al' : AlAlDensityFunction,
        'Fe' : FeAlDensityFunction },
      latticeConstant = 4.04527,
      latticeType = 'fcc'),
    #Fe
    EAMPotential('Fe', 26, 55.845, FeEmbedFunction,
      { 'Al': FeAlDensityFunction,
        'Fe' : FeFeDensityFunction},
      latticeConstant = 2.855312,
      latticeType = 'bcc') ]

  # Number of grid points and cut-offs for tabulation.
  nrho = 10000
  drho = 3.00000000000000E-2
  nr   = 10000
  dr   = 6.50000000000000E-4

  with open("Mendelev_Al_Fe.eam.fs", "wb") as outfile:
    writeSetFLFinnisSinclair(
      nrho, drho,
      nr, dr,
      eamPotentials,
      pairPotentials,
      outfile)

if __name__ == '__main__':
    main()
