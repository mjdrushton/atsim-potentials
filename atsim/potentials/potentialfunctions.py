"""Functions for different potential forms.

Most of the potentials in this module are implemented as callable objects. The potential energy
is evaluated by calling one of these objects. By convention the first argument of each is the
atomic separation `r`, with other potential parameters following after. For instance, to
evaluate a Buckingham potential at `r` = 2.0 the following could be called for `A`, `rho` and `C` values 
1000.0, 0.2 and 32.0 respectively:

```#! python

atsim.potentialfunctions.buck(2.0, 1000.0, 0.2, 32.0)

```

The callable objects also have other useful methods. Perhaps most importantly is the `.deriv()` method
this returns the first derivative of the given potential (force). Again using the Buckingham potential as
an example its derivative can be evaluated for `r` = 2.0 as follows:

```#! python

atsim.potentialfunctions.buck.deriv(2.0, 1000.0, 0.2, 32.0)

```

"""

from __future__ import division

import math

class buck(object):
  """Callable object for evaluating the Buckingham potential"""

  def _as_sympy(self):
    import sympy
    r, A, rho, C= sympy.symbols("r A rho C")
    return A * sympy.exp(-r/rho) - C/r**6

  def __call__(self, r, A, rho, C):
    """Buckingham potential form

    .. math ::

      U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right) - \\frac{C}{r_{ij}^6}

    :param r: Atomic separation.
    :param A: Buckingham A parameter
    :param rho: Buckingham rho parameter :math:`\\rho`
    :param C: Buckingham C parameter
    :return: Potential energy."""
    return A * math.exp(-r/rho) - (C/r**6)

  def deriv(self, r, A, rho, C):
    """Return derivative of Buckingham potential at `r`.

    :param r: Atomic separation.
    :param A: Buckingham A parameter
    :param rho: Buckingham rho parameter :math:`\\rho`
    :param C: Buckingham C parameter
    :return: Derivative of Buckingham potential at `r`"""
    return (6.0*C)/r**7 - ((A * math.exp(-r/rho))/rho)

  def deriv2(self, r, A, rho, C):
    """Return 2nd derivative of Buckingham potential at `r`.

    :param r: Atomic separation.
    :param A: Buckingham A parameter
    :param rho: Buckingham rho parameter :math:`\\rho`
    :param C: Buckingham C parameter
    :return: Second derivative of Buckingham potential at `r`"""
    return A*math.exp(-r/rho)/rho**2 - 42.0*C/r**8

buck = buck()

class bornmayer(object):
  """Callable object for evaluating Born-Mayer Potential"""

  def _as_sympy(self):
    import sympy
    r, A, rho = sympy.symbols("r A rho")
    return A * sympy.exp(-r/rho)

  def __call__(self, r, A, rho):
    """Born-Mayer potential form

    .. math ::

      U(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right)

    :param r: Atomic separation.
    :param A: Potential parameter
    :param rho: Potential parameter :math:`\\rho`
    :return: Potential energy"""
    return buck(r, A,rho,0.0)

  def deriv(self, r, A, rho):
    """Return derivative of Born-Mayer potential form at `r`

    :param r: Atomic separation.
    :param A: Potential parameter
    :param rho: Potential parameter :math:`\\rho`
    :return: Derivative at `r`"""
    return buck.deriv(r, A, rho, 0.0)

  def deriv2(self, r, A, rho):
    """Return 2nd derivative of Born-Mayer potential form at `r`

    :param r: Atomic separation.
    :param A: Potential parameter
    :param rho: Potential parameter :math:`\\rho`
    :return: 2nd derivative at `r`"""
    return buck.deriv2(r, A, rho, 0.0)

bornmayer = bornmayer()


class coul(object):
  """Callable representing Coulomb potential (including :math:`4\pi \epsilon_0` term)"""

  def _as_sympy(self):
    import sympy
    r, qi, qj = sympy.symbols("r q_i q_j")
    return (qi * qj)/(4.0*sympy.pi*0.0055264*r)
    

  def __call__(self,r, qi, qj):
    """Coulomb potential (including :math:`4\pi \epsilon_0` term).

    .. math ::

      U(r_{ij}) = \\frac{ q_i q_j }{4 \pi \epsilon_0 r_{ij} }

    .. note:: Constant value appropriate for :math:`r_{ij}` in angstroms and energy in eV.

    :param r: Atomic separation.
    :param qi: Charge on species i
    :param qj: Charge on species j

    :return: Potential energy"""
    return (qi * qj)/(4.0*math.pi*0.0055264*r)

  def deriv(self, r, qi, qj):
    """Return derivative pf Coulomb potential at `r`
    
    :param r: Atomic separation.
    :param qi: Charge on species i
    :param qj: Charge on species j
    :return: Derivative at `r`"""

    return -45.2374059061957*qi*qj/(math.pi*r**2)

  def deriv2(self, r, qi, qj):
    """Return second derivative pf Coulomb potential at `r`
    
    :param r: Atomic separation.
    :param qi: Charge on species i
    :param qj: Charge on species j
    :return: 2nd derivative at `r`"""
    return 90.4748118123914*qi*qj/(math.pi*r**3)

coul = coul()

class constant(object):
  """Callable that returns a constant"""

  def __call__(self,r, constant):
    """Function that simply returns the function that is passed to it.

    :param r: Separation.
    :param constant: The value to be returned by this function.

    :return: The value passed in as `constant`"""
    return constant

  def deriv(self, r, constant):
    """Returns derivative of this potential form. Here this will always be zero

    :param r: Separation.
    :param constant: Constant value

    :return: Derivative (0.0)
    """
    return 0.0

  def deriv2(self, r, constant):
    """Returns 2nd derivative of this potential form. Here this will always be zero

    :param r: Separation.
    :param constant: Constant value

    :return: 2nd derivative (0.0)
    """
    return 0.0

constant = constant()

class hbnd(object):
  """DL_POLY `hbnd` potential type"""
  
  def _as_sympy(self):
    import sympy
    r, A, B = sympy.symbols("r A B")
    return A/r**12 - B/r**10

  def __call__(self, r, A,B):
    """DL_POLY hbnd form:

    .. math::

      U(r_{ij}) = \\frac{A}{r_{ij}^{12}} - \\frac{B}{r_{ij}^{10}}

    :param r: Atomic separation.
    :param A: Potential A parameter
    :param B: Potentials' B parameter
    :return: Potential energy"""
    return (A/r**12) - (B/r**10)

  def deriv(self, r, A,B):
    """Derivative of `hbnd` at `r`

    :param r: Atomic separation.
    :param A: Potential A parameter
    :param B: Potentials' B parameter

    :return: Derivative of `hbnd` at `r`"""

    return (10.0*B)/r**11 - (12.0 * A)/r**13

  def deriv2(self, r, A,B):
    """Second derivative of `hbnd` at `r`

    :param r: Atomic separation.
    :param A: Potential A parameter
    :param B: Potentials' B parameter

    :return: 2nd derivative of `hbnd` at `r`"""
    return 2*(78*A/r**2 - 55*B)/r**12

hbnd = hbnd()

class lj(object):
  """Callable for Lennar-Jones 12-6 potential"""

  def _as_sympy(self):
    import sympy
    r, epsilon, sigma = sympy.symbols("r epsilon sigma")
    return 4.0*epsilon*( (sigma**12/r**12) - (sigma**6/r**6))

  def __call__(self, r, epsilon, sigma):
    """Lennard-Jones 12-6 potential.

    .. math ::

      U(r_{ij}) = 4 \epsilon \left( \\frac{\sigma^{12}}{r_{ij}^{12}} - \\frac{\sigma^6}{r_{ij}^6} \\right)

    :param r: Atomic separation.
    :param epsilon: Epsilon parameter :math:`\epsilon`
    :param sigma: Sigma parameter :math:`\sigma`
    :return: Potential energy"""
    return 4.0*epsilon*( (sigma**12/r**12) - (sigma**6/r**6))

  def deriv(self, r, epsilon, sigma):
    """Derivative of Lennard-Jones 12-6 potential at `r`

    :param r: Atomic separation.
    :param epsilon: Epsilon parameter :math:`\epsilon`
    :param sigma: Sigma parameter :math:`\sigma`
    :return: Derivative of potential at `r`"""
    return (24.0*epsilon*sigma**6)/r**7 - (48.0*epsilon*sigma**12)/r**13

  def deriv2(self, r, epsilon, sigma):
    """Second derivative of Lennard-Jones 12-6 potential at `r`

    :param r: Atomic separation.
    :param epsilon: Epsilon parameter :math:`\epsilon`
    :param sigma: Sigma parameter :math:`\sigma`
    :return: 2nd derivative of potential at `r`"""
    return -168.0*epsilon*sigma**6/r**8 + 624.0*epsilon*sigma**12/r**14

lj = lj()


class morse(object):
  """Callable representing the Morse potential"""

  def _as_sympy(self):
    import sympy
    D, gamma, r_star, r = sympy.symbols("D gamma r_star r")
    return D*(sympy.exp(-2.0*gamma*(r-r_star)) - 2.0*sympy.exp(-gamma*(r-r_star)))

  def __call__(self, r, gamma, r_star, D):
    """Morse potential parametrised with gamma, r_star and D.

    Potential function is (where r is interatomic separation):

    .. math ::

      U(r_{ij}) = D  \left[ \exp \left( -2 \gamma (r_{ij} - r_*) \\right) - 2 \\exp \left( -\gamma (r - r_*) \\right) \\right]

    :param r: Atomic separation.
    :param gamma: Potential parameter :math:`\gamma`
    :param r_star: Potential parameter :math:`r_*`
    :param D: Potential parameter
    :return: Potential energy"""
    return D*(math.exp(-2.0*gamma*(r-r_star)) - 2.0*math.exp(-gamma*(r-r_star)))

  def deriv(self, r, gamma, r_star, D):
    """Evaluate derivative of Morse potential at `r`.
    :param r: Atomic separation.
    :param gamma: Potential parameter :math:`\gamma`
    :param r_star: Potential parameter :math:`r_*`
    :param D: Potential parameter

    :return: Derivative of Morse potential at `r`"""
    return D*(-2.0*gamma*math.exp(-2.0*gamma*(r - r_star)) + 2.0*gamma*math.exp(-gamma*(r - r_star)))

  def deriv2(self, r, gamma, r_star, D):
    """Evaluate derivative of Morse potential at `r`.
    :param r: Atomic separation.
    :param gamma: Potential parameter :math:`\gamma`
    :param r_star: Potential parameter :math:`r_*`
    :param D: Potential parameter

    :return: Derivative of Morse potential at `r`"""
    return D*gamma**2*(4.0*math.exp(-2.0*gamma*(r - r_star)) - 2.0*math.exp(-gamma*(r - r_star)))

morse = morse()


class polynomial(object):
  """Callable for polynomials"""

  def _split_args(self, args):
    return args[0], args[1:]

  def __call__(self, *args):
    """Polynomial of arbitrary order.

    .. math::

      U(r_{ij}) = C_0 + C_1 r_{ij} + C_2 r_{ij}^2 + \dots + C_n r_{ij}^n

    This function accepts a variable number of arguments - the first is :math:`r_{ij}`
    and with the remainder being :math:`C_0, C_1, \dots, C_n` respectively.

    :return: Potential energy for given polynomial"""

    r, coefs = self._split_args(args)
    v = [r**float(i) * c for (i,c) in enumerate(coefs)]
    return sum(v)
  
  def deriv(self, *args):
    """Evaluate polynomial derivative.

    This function accepts a variable number of arguments - the first is :math:`r_{ij}`
    and with the remainder being the polynomial coefficients :math:`C_0, C_1, \dots, C_n` respectively.

    :return: derivative of polynomial at `r` """
    r, coefs = self._split_args(args)
    v = [float(i) * r**float(i-1) * c for (i,c) in enumerate(coefs)][1:]
    return sum([0]+v)

  def deriv2(self, *args):
    """Evaluate polynomial second derivative.

    This function accepts a variable number of arguments - the first is :math:`r_{ij}`
    and with the remainder being the polynomial coefficients :math:`C_0, C_1, \dots, C_n` respectively.

    :return: 2nd derivative of polynomial at `r` """
    r, coefs = self._split_args(args)
    v = [i * float(i-1) * r**float(i-2) * c for (i,c) in enumerate(coefs)][2:]
    return sum([0]+v)

polynomial = polynomial()


class sqrt(object):
  """Callable representing a square root potential form"""

  def _as_sympy(self):
    import sympy
    r, G = sympy.symbols("r, G")
    return G * r**(1/2)

  def __call__(self, r, G):
    """Evaluate square root function at `r`.

    Potential function is:

    .. math ::

      U(r_{ij}) = G\sqrt{r_{ij}}

    :param r: Atomic separation.
    :param r: Variable
    :param G: Parameter :math:`G`
    :return: G multiplied by the square root of r."""
    return G*math.sqrt(r)

  def deriv(self, r, G):
    """Evaluate derivative of square root function at `r`.

    :param r: Atomic separation.
    :param r: Variable
    :param G: Parameter :math:`G`

    :return: Derivative at `r`. """

    return (0.5 * G)/math.sqrt(r)

  def deriv2(self, r, G):
    """Evaluate second derivative of square root function at `r`.

    :param r: Atomic separation.
    :param r: Variable
    :param G: Parameter :math:`G`

    :return: 2nd derivative at `r`. """
    return -G/(4*r**(3/2))
    
sqrt = sqrt()


class zbl(object):
  """Callable representing ZBL potential."""

  Ck1=0.1818
  Ck2=0.5099
  Ck3=0.2802
  Ck4=0.02817
  Bk1=3.2
  Bk2=0.9423
  Bk3=0.4029
  Bk4=0.2016

  def _as_sympy(self):
    import sympy
    Ck1, Ck2, Ck3, Ck4, Bk1, Bk2, Bk3, Bk4 = sympy.symbols("Ck1:5,Bk1:5")
    z1, z2, r = sympy.symbols("z1, z2, r")

    a=(0.8854*0.529)/(z1**0.23 + z2**0.23)
    f = 14.39942 * (z1*z2)/r * (Ck1*sympy.exp((-Bk1*r)/a) + Ck2*sympy.exp((-Bk2*r)/a) + Ck3*sympy.exp((-Bk3*r)/a) + Ck4*sympy.exp((-Bk4*r)/a))
    return f

  def __call__(self, r, z1,z2):
    """ZBL potential.

    Ziegler-Biersack-Littmark screened nuclear repulsion for describing high energy interactions.

    :param r: Atomic separation.
    :param z1: Atomic number of species i
    :param z2: Atomic number of species j
    :return: Potential energy"""
    a=(0.8854*0.529)/(float(z1)**0.23 + float(z2)**0.23)
    return 14.39942 * (z1*z2)/r * (self.Ck1*math.exp((-self.Bk1*r)/a) + self.Ck2*math.exp((-self.Bk2*r)/a) + self.Ck3*math.exp((-self.Bk3*r)/a) + self.Ck4*math.exp((-self.Bk4*r)/a))


  def deriv(self, r, z1, z2):
    """Evaluate derivative of ZBL function at `r`.

    :param r: Atomic separation.
    :param z1: Atomic number of species i
    :param z2: Atomic number of species j
    :return: Derivative of function"""
    v = -14.39942*z1*z2*(self.Ck1*math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        * (self.Bk2 + self.Bk3 + self.Bk4))\
        + self.Ck2*math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk3 + self.Bk4))\
        + self.Ck3*math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk2 + self.Bk4))\
        + self.Ck4*math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk2 + self.Bk3))\
        + 2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1*self.Ck1*math.exp(2.13503407300877\
        *r*(z1**0.23 + z2**0.23)*(self.Bk2 + self.Bk3 + self.Bk4))\
        + self.Bk2*self.Ck2*math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk3 + self.Bk4))\
        + self.Bk3*self.Ck3*math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk2 + self.Bk4)) + self.Bk4*self.Ck4\
        *math.exp(2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk2 + self.Bk3))))\
        *math.exp(-2.13503407300877*r*(z1**0.23 + z2**0.23)\
        *(self.Bk1 + self.Bk2 + self.Bk3 + self.Bk4))/r**2
    return v

  def deriv2(self, r, z1, z2):
    """Evaluate second derivative of ZBL function at `r`.

    :param r: Atomic separation.
    :param z1: Atomic number of species i
    :param z2: Atomic number of species j
    :return: 2nd derivative of function"""
    v = z1*z2*(65.6378912429954*(z1**0.23 + z2**0.23)**2* \
        (self.Bk1**2*self.Ck1*math.exp(-2.13503407300877*self.Bk1*r*(z1**0.23 + z2**0.23)) + \
        self.Bk2**2*self.Ck2*math.exp(-2.13503407300877*self.Bk2*r*(z1**0.23 + z2**0.23)) + \
        self.Bk3**2*self.Ck3*math.exp(-2.13503407300877*self.Bk3*r*(z1**0.23 + z2**0.23)) + \
        self.Bk4**2*self.Ck4*math.exp(-2.13503407300877*self.Bk4*r*(z1**0.23 + z2**0.23))) + \
        28.79884*(2.13503407300877*z1**0.23 + 2.13503407300877*z2**0.23)* \
        (self.Bk1*self.Ck1*math.exp(-2.13503407300877*self.Bk1*r*(z1**0.23 + z2**0.23)) + \
        self.Bk2*self.Ck2*math.exp(-2.13503407300877*self.Bk2*r*(z1**0.23 + z2**0.23)) + \
        self.Bk3*self.Ck3*math.exp(-2.13503407300877*self.Bk3*r*(z1**0.23 + z2**0.23)) + \
        self.Bk4*self.Ck4*math.exp(-2.13503407300877*self.Bk4*r*(z1**0.23 + z2**0.23)))/r + \
        (28.79884*self.Ck1*math.exp(-2.13503407300877*self.Bk1*r*(z1**0.23 + z2**0.23)) + \
        28.79884*self.Ck2*math.exp(-2.13503407300877*self.Bk2*r*(z1**0.23 + z2**0.23)) + \
        28.79884*self.Ck3*math.exp(-2.13503407300877*self.Bk3*r*(z1**0.23 + z2**0.23)) + \
        28.79884*self.Ck4*math.exp(-2.13503407300877*self.Bk4*r*(z1**0.23 + z2**0.23)))/r**2)/r
    return v

zbl = zbl()


class zero(object):
  """Callable that returns 0.0 for any given input"""

  def __call__(self, r):
    """Returns 0.0 for any value of r.

    :param r: Separation.
    
    :return: Returns 0.0"""
    return 0.0

  def deriv(self, r):
    """Derivative function - always returns 0.0"""
    return 0.0

  def deriv2(self, r):
    """Second derivative function - always returns 0.0"""
    return 0.0

zero = zero()

class exp_spline(object):
  """Callable representing exponential spline"""

  def _as_sympy(self):
    import sympy
    r, B0, B1, B2, B3, B4, B5, C = sympy.symbols("r, B0, B1, B2, B3, B4, B5, C")
    return sympy.exp(B0 + B1*r + B2*r**2 + B3*r**3 + B4*r**4 + B5*r**5) + C

  def __call__(self, r, B0, B1, B2, B3, B4, B5, C):
    """Exponential spline function (as used in splining routines).

      .. math::

              U(r_{ij}) = \exp \left( B_0 + B_1 r_{ij} + B_2 r_{ij}^2 + B_3 r_{ij}^3 + B_4 r_{ij}^4 + B_5 r_{ij}^5 \\right) + C

    :param r: Atomic separation
    :param B0: Spline coefficient
    :param B1: Spline coefficient
    :param B2: Spline coefficient
    :param B3: Spline coefficient
    :param B4: Spline coefficient
    :param B5: Spline coefficient
    :param C: C parameter

    :return: Splined values"""
    return math.exp(B0 + B1*r + B2*r**2 + B3*r**3 + B4*r**4 + B5*r**5) + C

  def deriv(self, r, B0, B1, B2, B3, B4, B5, C):
    """Derivative of exponential spline. 

      :param r: Atomic separation
      :param B0: Spline coefficient
      :param B1: Spline coefficient
      :param B2: Spline coefficient
      :param B3: Spline coefficient
      :param B4: Spline coefficient
      :param B5: Spline coefficient
      :param C: C parameter

      :return: Derivative of spline at `r`"""

    return (B1 + r*(2*B2 + r*(3*B3 + 4*B4*r + 5*B5*r**2))) * math.exp(B0 + r*(B1 + r*(B2 + r*(B3 + r*(B4 + B5*r)))))

  def deriv2(self, r, B0, B1, B2, B3, B4, B5, C):
    """Second derivative of exponential spline. 

      :param r: Atomic separation
      :param B0: Spline coefficient
      :param B1: Spline coefficient
      :param B2: Spline coefficient
      :param B3: Spline coefficient
      :param B4: Spline coefficient
      :param B5: Spline coefficient
      :param C: C parameter

      :return: 2nd derivative of spline at `r`"""

    return (2*B2 + 6*B3*r + 12*B4*r**2 + 20*B5*r**3 \
      + (B1 + 2*B2*r + 3*B3*r**2 + 4*B4*r**3 + 5*B5*r**4)**2) \
      * math.exp(B0 + B1*r + B2*r**2 + B3*r**3 + B4*r**4 + B5*r**5)

exp_spline = exp_spline()

