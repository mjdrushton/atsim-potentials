import math

import types

"""Module containing functions and classes that are shared between different tests"""

def assertFloatWithinPercentage(testCase, expect, actual, percenttolerance = 0.5, places = 5):
  """Assert that actual is within percenttolerance of expect. Percenttolerance is specified as a percentage of
  expect. If expect is == 0.0 then revert to using testCase.assertAlmostEqual() method and check
  that actual is the same as expect for given number of places"""

  if expect == 0.0:
    testCase.assertAlmostEqual(expect, actual, places = places)
  else:
    percentDifference = math.fabs( ((actual-expect)/float(expect)) * 100.0 )
    msg = "Actual %f != Expect %f within tolerance of %f%% (difference = %f%%)" % (actual, expect, percenttolerance, percentDifference)
    testCase.assertTrue(percentDifference <= percenttolerance, msg)

def checkVector(tc, expect, actual, msg=None, tolerance = 0.0025):
  """Assert that two vectors are within tolerance of each other

  @param tc unittest.TestCase object
  @param expect Expected vector
  @param actual Actual vector
  @param msg Test fail message
  @param tolerance Acceptable distance between expected and actual vectors"""
  tc.assertEqual(len(expect), len(actual), msg=msg)
  diff = (expect[0] - actual[0], expect[1] - actual[1], expect[2]-actual[2])
  dist = math.sqrt( diff[0]**2 + diff[1]**2 + diff[2]**2)
  if msg == None:
    msg = "%s != %s" % (expect, actual)
  tc.assertTrue(dist <= dist, msg = msg)

def _compareCollection(path, testCase, expect, actual, places, percenttolerance):
  expectType = type(expect)
  if expectType == list or expectType == tuple:
    #Compare lists
    try:
      testCase.assertEqual(len(expect), len(actual))
    except AssertionError as e:
      raise AssertionError("%s at '%s'" % (str(e), path))

    for i,(e,a) in enumerate(zip(expect, actual)):
      _compareCollection(path+'[%d]'% i, testCase, e,a, places, percenttolerance)
  elif expectType == dict:
    #Compare dictionaries
    ekeys = list(expect.keys())
    akeys = list(actual.keys())
    ekeys.sort()
    akeys.sort()
    testCase.assertEqual(ekeys, akeys)
    for k,v in expect.items():
      _compareCollection(path+'[%s]'% (k,), testCase, v, actual[k], places, percenttolerance)
  elif expectType == float:
    #Compare float type in a fuzzy manner
    try:
      if math.isnan(expect):
        testCase.assertTrue(math.isnan(actual))
      elif percenttolerance != None:
        assertFloatWithinPercentage(testCase, expect, actual, percenttolerance = percenttolerance, places = places)
      else:
        testCase.assertAlmostEqual(expect, actual, places = places)
    except AssertionError as e:
      raise AssertionError("%s at '%s'" % (str(e), path))
  else:
    #Compare anything else
    try:
      testCase.assertEqual(expect,actual)
    except AssertionError as e:
      raise AssertionError("%s at '%s'" % (str(e), path))

def compareCollection(testCase, expect, actual, places = 5, percenttolerance = None):
  """Check two collections are the same"""
  path = "collection"
  _compareCollection(path, testCase, expect, actual, places, percenttolerance)

