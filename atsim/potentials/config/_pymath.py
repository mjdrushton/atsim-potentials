import math

import sys

def ceil(x):
  return math.ceil(x)

def copysign(a,b):
  return math.copysign(a,b)

def fabs(x):
  return math.fabs(x)

def floor(x):
  return math.floor(x)

def fmod(a,b):
  return math.fmod(a,b)

def fsum(*args):
  return math.fsum(args)

def ldexp(a,b):
  return math.ldexp(a,int(b))

def trunc(x):
  return math.trunc(x)

def exp(x):
  return math.exp(x)

def log(*args):
  return math.log(*args)

def log1p(x):
  return math.log1p(x)

if hasattr(math, "log2"):
  def log2(x):
    return math.log2(x)
else:
  def log2(x):
    return math.log(x,2)

def log10(x):
  return math.log10(x)

def pow(x,a):
  return math.pow(x,a)

def sqrt(x):
  return math.sqrt(x)

def acos(x):
  return math.acos(x)

def atan(x):
  return math.atan(x)

def atan2(x,y):
  return math.atan2(x,y)

def cos(x):
  return math.cos(x)

def hypot(x,y):
  return math.hypot(x,y)

def sin(x):
  return math.sin(x)

def tan(x):
  return math.tan(x)

def radians(x):
  return math.radians(x)

def degrees(x):
  return math.degrees(x)

def acosh(x):
  return math.acosh(x)

def asinh(x):
  return math.asinh(x)

def atanh(x):
  return math.atanh(x)

def cosh(x):
  return math.cosh(x)

def sinh(x):
  return math.sinh(x)

def tanh(x):
  return math.tanh(x)

def factorial(x):
  return math.factorial(int(x))

if hasattr(math, "gcd"):
  _gcd = math.gcd
else:
  import fractions
  _gcd = fractions.gcd

def gcd(a,b):
  return _gcd(int(a),int(b))


