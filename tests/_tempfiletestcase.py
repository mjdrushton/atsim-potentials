import unittest
import shutil

class TempfileTestCase(unittest.TestCase):
  """TestCase base class that provides a temporary directory that is
  torn-down between tests.

  The location of the temporary directory is accessible through the
  tempdir attribute"""

  def setUp(self):
    import tempfile
    self.tempdir = tempfile.mkdtemp()

  def tearDown(self):
    shutil.rmtree(self.tempdir, ignore_errors = True)
