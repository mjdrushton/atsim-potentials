import pathlib

def _get_resource_dir(subdir) -> pathlib.Path:
  rd = pathlib.Path(__file__).parent.parent / subdir
  return rd

def _get_lammps_resource_dir() -> pathlib.Path:
  return _get_resource_dir('lammps_resources')

def _get_dlpoly_resource_dir() -> pathlib.Path:
  return _get_resource_dir('dl_poly_resources')
