"""Functions for different potential forms.

Most of the potentials in this module are implemented as callable _Potential_Function_Bases. The potential energy
is evaluated by calling one of these objects. By convention the first argument of each is the
atomic separation `r`, with other potential parameters following after. For instance, to
evaluate a Buckingham potential at `r` = 2.0 the following could be called for `A`, `rho` and `C` values 
1000.0, 0.2 and 32.0 respectively:

.. code-block:: python

  atsim.potentialfunctions.buck(2.0, 1000.0, 0.2, 32.0)


The callable objects also have other useful methods. Perhaps most importantly is the `.deriv()` method
this returns the first derivative of the given potential (force). Again using the Buckingham potential as
an example its derivative can be evaluated for `r` = 2.0 as follows:

.. code-block:: python

  atsim.potentialfunctions.buck.deriv(2.0, 1000.0, 0.2, 32.0)


See :ref:`list-of-potential-forms` for descriptions of these potential forms.

"""

import math

class _Potential_Function_Base(object):
  is_potential = True


class _buck(_Potential_Function_Base):
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


buck = _buck()

class _bornmayer(_Potential_Function_Base):
  """Callable object for evaluating Born-Mayer Potential"""

  def _as_sympy(self):
    import sympy
    r, A, rho = sympy.symbols("r A rho")
    return A * sympy.exp(-r/rho)

  def __call__(self, r, A, rho):
    """Born-Mayer potential form

    .. math ::

      V(r_{ij}) = A \exp \left( \\frac{-r_{ij}}{\\rho} \\right)

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

bornmayer = _bornmayer()


class _coul(_Potential_Function_Base):
  """Callable representing Coulomb potential (including :math:`4\pi \epsilon_0` term)"""

  def _as_sympy(self):
    import sympy
    r, qi, qj = sympy.symbols("r q_i q_j")
    return (qi * qj)/(4.0*sympy.pi*0.0055264*r)
    

  def __call__(self,r, qi, qj):
    """Coulomb potential (including :math:`4\pi \epsilon_0` term).

    .. math ::

      V(r_{ij}) = \\frac{ q_i q_j }{4 \pi \epsilon_0 r_{ij} }

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

coul = _coul()

class _constant(_Potential_Function_Base):
  """Callable that returns a constant"""

  def __call__(self,r, constant):
    """Function that simply returns the value that is passed to it.

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

constant = _constant()

class _exponential(_Potential_Function_Base):
  """`exponential` potential type A*r^n"""
  
  def _as_sympy(self):
    import sympy
    r, A, n = sympy.symbols("r A n")
    return A*r**n

  def __call__(self, r, A,n):
    """General exponential form:

    .. math::

      V(r_{ij}) = A r^n

    :param r: Atomic separation.
    :param A: Potential A parameter
    :param n: n
    :return: Potential energy"""
    return A*r**n

  def deriv(self, r, A, n):
    """Derivative of `exponential` at `r`

    :param r: Atomic separation.
    :param A: Potential A parameter
    :param n: Potentials' B parameter

    :return: Derivative of `exponential` at `r`"""
    return A*n*r**(n-1)

  def deriv2(self, r, A,n):
    """Second derivative of `exponential` at `r`

    :param r: Atomic separation.
    :param A: Potential A parameter
    :param n: Potentials' B parameter

    :return: 2nd derivative of `exponential` at `r`"""
    return A*n*(n-1)*r**(n-2) 

exponential = _exponential()

class _hbnd(_Potential_Function_Base):
  """DL_POLY `hbnd` potential type"""
  
  def _as_sympy(self):
    import sympy
    r, A, B = sympy.symbols("r A B")
    return A/r**12 - B/r**10

  def __call__(self, r, A,B):
    """DL_POLY hbnd form:

    .. math::

      V(r_{ij}) = \\frac{A}{r_{ij}^{12}} - \\frac{B}{r_{ij}^{10}}

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

hbnd = _hbnd()

class _lj(_Potential_Function_Base):
  """Callable for Lennard-Jones 12-6 potential"""

  def _as_sympy(self):
    import sympy
    r, epsilon, sigma = sympy.symbols("r epsilon sigma")
    return 4.0*epsilon*( (sigma**12/r**12) - (sigma**6/r**6))

  def __call__(self, r, epsilon, sigma):
    """Lennard-Jones 12-6 potential.

    .. math ::

      V(r_{ij}) = 4 \epsilon \left( \\frac{\sigma^{12}}{r_{ij}^{12}} - \\frac{\sigma^6}{r_{ij}^6} \\right)

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

lj = _lj()


class _morse(_Potential_Function_Base):
  """Callable representing the Morse potential"""

  def _as_sympy(self):
    import sympy
    D, gamma, r_star, r = sympy.symbols("D gamma r_star r")
    return D*(sympy.exp(-2.0*gamma*(r-r_star)) - 2.0*sympy.exp(-gamma*(r-r_star)))

  def __call__(self, r, gamma, r_star, D):
    """Morse potential parametrised with gamma, r_star and D.

    Potential function is (where r is interatomic separation):

    .. math ::

      V(r_{ij}) = D  \left[ \exp \left( -2 \gamma (r_{ij} - r_*) \\right) - 2 \\exp \left( -\gamma (r - r_*) \\right) \\right]

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

morse = _morse()


class _polynomial(_Potential_Function_Base):
  """Callable for polynomials"""

  def _split_args(self, args):
    return args[0], args[1:]

  def __call__(self, *args):
    """Polynomial of arbitrary order.

    .. math::

      V(r_{ij}) = C_0 + C_1 r_{ij} + C_2 r_{ij}^2 + \dots + C_n r_{ij}^n

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

polynomial = _polynomial()


class _sqrt(_Potential_Function_Base):
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
    
sqrt = _sqrt()

class _tang_toennies(_Potential_Function_Base):
  """Callable for Tang-Toennies potential form"""

  def _f2n(self, x, n):
    import sympy 
    v = None
    for k in range(2*n+1):
      j = x**k / sympy.factorial(k)
      if v is None:
        v = j
      else:
        v = v+j
    return 1.0 - sympy.exp(-x) * v

  def _sp_second_term(self, r, b, C_6, C_8, C_10):
    n_3 = self._f2n(b*r, 3) * C_6/r**6
    n_4 = self._f2n(b*r, 4) * C_8/r**8
    n_5 = self._f2n(b*r, 5) * C_10/r**10
    return n_3 + n_4 + n_5

  def _as_sympy(self):
    import sympy

    A,b,C_6,C_8,C_10,r = sympy.symbols("A,b,C_6,C_8,C_10,r")

    r = r/0.5292

    f = A*sympy.exp(-b*r)
    s = self._sp_second_term(r,b,C_6, C_8, C_10)
    v = (f-s)*27.211

    return v

  def __call__(self, r, A,b,C_6,C_8,C_10):
    """Evaluate the Tang-Toennies potential form at a separation `r`.

    This potential form describes the Van der Waal's interactions between the noble gases (He to Rn) and is described in:

      * K.T. Tang and J.P. Toennies, “The van der Waals potentials between all the rare gas atoms from He to Rn”, *J. Chem. Phys.*\ , **118** (2003) 4976.

    The potential has the following form:

    .. math::

      V(r) = A \\exp(-br) - \\sum_{n=3}^N f_{2N} (bR) \\frac{C_{2N}}{R^{2N}}

    Where:

    .. math::

      f_{2N}(x) = 1- \\exp(-x) \\sum_{k=0}^{2n} \\frac{x^k}{k!}


    :param r: Separation.
    :param A: Potential parameter.
    :param b: Potential parameter.
    :param C_6: Potential parameter.
    :param C_8: Potential parameter.
    :param C_10: Potential parameter.

    :return: Potential evaluated at `r`."""


    return 27.210999999999998522*A*math.exp(-1.8896447467876*b*r) -\
     0.046875170184330419709*C_10*(-(0.00015997002473914140102*b**10*r**10 +\
     0.00084656137091953641977*b**9*r**9 + 0.0040320024974155677447*b**8*r**8 +\
     0.017069885773058547651*b**7*r**7 + 0.063233684857718089334*b**6*r**6 +\
     0.2007795961602264756*b**5*r**5 + 0.53126281143995923717*b**4*r**4 +\
     1.1245771192561058172*b**3*r**3 + 1.7853786345309936578*b**2*r**2 +\
     1.8896447467876038573*b*r + 1.0)*math.exp(-1.8896447467876*b*r) + 1.0)/r**10 -\
     0.59767283277249427798*C_6*(-(0.063233684857718089334*b**6*r**6 +\
     0.2007795961602264756*b**5*r**5 + 0.53126281143995923717*b**4*r**4 +\
     1.1245771192561058172*b**3*r**3 + 1.7853786345309936578*b**2*r**2 +\
     1.8896447467876038573*b*r + 1.0)*math.exp(-1.8896447467876*b*r) + 1.0)/r**6 -\
     0.16737985467421556685*C_8*(-(0.0040320024974155677447*b**8*r**8 +\
     0.017069885773058547651*b**7*r**7 + 0.063233684857718089334*b**6*r**6 +\
     0.2007795961602264756*b**5*r**5 + 0.53126281143995923717*b**4*r**4 +\
     1.1245771192561058172*b**3*r**3 + 1.7853786345309936578*b**2*r**2 +\
     1.8896447467876038573*b*r + 1.0)*math.exp(-1.8896447467876*b*r) + 1.0)/r**8

  def deriv(self, r, A,b,C_6,C_8,C_10):
    """Evaluate the first derivative of Tang-Toennies potential form at a separation `r`.

    :param r: Separation.
    :param A: Potential parameter.
    :param b: Potential parameter.
    :param C_6: Potential parameter.
    :param C_8: Potential parameter.
    :param C_10: Potential parameter.

    :return: Derivative evaluated at `r`."""


    return (-51.419123204837482888*A*b*r**11 - 0.000014169731923731671203*C_10*b**11*r**11 - 0.000074986221340387997341*C_10*b**10*r**10 -\
      0.00039682708333333333072*C_10*b**9*r**9 - 0.0018900080324999999661*C_10*b**8*r**8 - 0.0080015380063920005238*C_10*b**7*r**7 -\
      0.029640897390878523376*C_10*b**6*r**6 - 0.094115777395517505322*C_10*b**5*r**5 - 0.24903034698853929174*C_10*b**4*r**4 -\
      0.5271474385053400713*C_10*b**3*r**3 - 0.83689927337107783423*C_10*b**2*r**2 - 0.88577419093594889077*C_10*b*r +\
      0.46875170184330416934*C_10*math.exp(1.8896447467876*b*r) - 0.46875170184330416934*C_10 - 0.071415448895607608337*C_6*b**7*r**11 -\
      0.22675833333333325625*C_6*b**6*r**10 - 0.720003059999999806*C_6*b**5*r**9 - 1.9051280967599995009*C_6*b**4*r**8 -\
      4.0327751552215671538*C_6*b**3*r**7 - 6.4024338364297603832*C_6*b**2*r**6 - 6.7763359724772591619*C_6*b*r**5 +\
      3.5860369966349656679*C_6*r**4*math.exp(1.8896447467876*b*r) - 3.5860369966349656679*C_6*r**4 - 0.0012752758731358502728*C_8*b**9*r**11 -\
      0.0053990079365079353749*C_8*b**8*r**10 - 0.022857239999999997421*C_8*b**7*r**9 - 0.084672359855999981826*C_8*b**6*r**8 -\
      0.2688516770147711954*C_8*b**5*r**7 - 0.71138153738108456103*C_8*b**4*r**6 - 1.5058524383282798631*C_8*b**3*r**5 -\
      2.3906913310899771119*C_8*b**2*r**4 - 2.5303077048256321646*C_8*b*r**3 + 1.3390388373937245348*C_8*r**2*math.exp(1.8896447467876*b*r) -\
      1.3390388373937245348*C_8*r**2)*math.exp(-1.8896447467876*b*r)/r**11

  def deriv2(self, r, A,b,C_6,C_8,C_10):
    """Evaluate the second derivative of Tang-Toennies potential form at a separation `r`.

    :param r: Separation.
    :param A: Potential parameter.
    :param b: Potential parameter.
    :param C_6: Potential parameter.
    :param C_8: Potential parameter.
    :param C_10: Potential parameter.

    :return: Second derivative evaluated at `r`."""

    return 1.0*(97.163876048445729339*A*b**2*r**12 +\
     0.000026775759493068161437*C_10*b**12*r**12 + 0.00014169731923731671203*C_10*b**11*r**11 +\
     0.0008248484347442676997*C_10*b**10*r**10 + 0.0043650979166666662584*C_10*b**9*r**9 +\
     0.02079008835750000006*C_10*b**8*r**8 + 0.088016918070311991884*C_10*b**7*r**7 +\
     0.32604987129966372938*C_10*b**6*r**6 + 1.0352735513506925447*C_10*b**5*r**5 + 2.7393338168739322924*C_10*b**4*r**4 +\
     5.7986218235587401182*C_10*b**3*r**3 + 9.2058920070818572867*C_10*b**2*r**2 + 9.7435161002954373544*C_10*b*r -\
     5.1562687202763459737*C_10*math.exp(1.8896447467876*b*r) + 5.1562687202763459737*C_10 + 0.13494982784506348583*C_6*b**8*r**12 +\
     0.428492693373645539*C_6*b**7*r**11 + 1.587308333333332433*C_6*b**6*r**10 + 5.040021419999998642*C_6*b**5*r**9 +\
     13.335896677319995618*C_6*b**4*r**8 + 28.229426086550969188*C_6*b**3*r**7 + 44.817036855008325347*C_6*b**2*r**6 +\
     47.434351807340810581*C_6*b*r**5 - 25.102258976444758787*C_6*r**4*math.exp(1.8896447467876*b*r) + 25.102258976444758787*C_6*r**4 +\
     0.0024098183543761341092*C_8*b**10*r**12 + 0.010202206985086802182*C_8*b**9*r**11 + 0.048591071428571420976*C_8*b**8*r**10 +\
     0.2057151599999999525*C_8*b**7*r**9 + 0.76205123870399982255*C_8*b**6*r**8 + 2.4196650931329406475*C_8*b**5*r**7 +\
     6.4024338364297603832*C_8*b**4*r**6 + 13.552671944954518324*C_8*b**3*r**5 + 21.516221979809795783*C_8*b**2*r**4 +\
     22.772769343430688593*C_8*b*r**3 - 12.051349536543520813*C_8*r**2*math.exp(1.8896447467876*b*r) + 12.051349536543520813*C_8*r**2)*math.exp(-1.8896447467876*b*r)/r**12

tang_toennies = _tang_toennies()

class _zbl(_Potential_Function_Base):
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

zbl = _zbl()


class _zero(_Potential_Function_Base):
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

zero = _zero()

class _exp_spline(_Potential_Function_Base):
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

exp_spline = _exp_spline()

