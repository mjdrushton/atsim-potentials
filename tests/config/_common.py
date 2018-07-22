import py

def _get_resource_dir(subdir):
  rd = py.path.local(__file__).dirpath()
  rd = rd.parts()[-2].join(subdir)
  return rd

def _get_lammps_resource_dir():
  return _get_resource_dir('lammps_resources')

def _get_dlpoly_resource_dir():
  return _get_resource_dir('dl_poly_resources')
