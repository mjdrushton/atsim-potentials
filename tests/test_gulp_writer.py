import pytest

import io

from atsim.potentials import potentialforms, pair_tabulation, Potential, Multi_Range_Defn, create_Multi_Range_Potential_Form, writePotentials

from ._rungulp import needsGULP, runGULP, extractGULPEnergy

@needsGULP
def test_single_pair(tmpdir):
  gulp_input = u"""single

cell
50.0 50.0 50.0 90.0 90.0 90.0

cartesian
Au 25.0 25.0 25.0
B  {} 25.0 25.0

species
Au 0.0
B 0.0

include potentials.lib"""

  buck = potentialforms.buck(1000.0, 0.23, 11.6)

  mrdefn = Multi_Range_Defn(">", 0, buck)
  mrpot = create_Multi_Range_Potential_Form(mrdefn)

  pot = Potential(u"Au", u"B", mrpot)

  tabulation = pair_tabulation.GULP_PairTabulation([pot], 10.0, 500)

  with tmpdir.join("potentials.lib").open("w") as potfile:
    tabulation.write(potfile)

  # Reproduce a buckingham potential at sepns of 0.1 1.1 2.1 and 3.1
  for x in [0.6, 1.1, 2.1, 3.1]:
    gulp_infile = io.StringIO(gulp_input.format(x+25.0))
    gulp_infile.seek(0)

    expect  = pytest.approx(buck(x), rel = 1e-4)

    gulp_outfile = io.StringIO()
    runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

    gulp_outfile.seek(0)
    actual = extractGULPEnergy(gulp_outfile)

    assert expect == actual


@needsGULP
def test_structure(tmpdir):
  gulp_input = u"""single

cell
5.468 5.468 5.468 90.0 90.0 90.0

frac
U 0 0 0
U 1/2 1/2 0
U 1/2 0 1/2
U 0 1/2 1/2

O 1/4 1/4 1/4
O 1/4 3/4 1/4
O 3/4 3/4 1/4
O 3/4 1/4 1/4

O 1/4 1/4 3/4
O 1/4 3/4 3/4
O 3/4 3/4 3/4
O 3/4 1/4 3/4


species
U 4.0
O -2.0

include potentials.lib"""

  # First calculate the expected energy using GULP's built-in analytical potentials
  with tmpdir.join("potentials.lib").open("w") as potfile:
    potfile.write("buck\n")
    potfile.write("O O 9547.96 0.2192 32.0 15.0\n")
    potfile.write("O U 1761.775 0.35642 0.0 15.0\n")

  gulp_infile = io.StringIO(gulp_input)
  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  expect = extractGULPEnergy(gulp_outfile)

  tmpdir.join("potentials.lib").remove()
  assert not tmpdir.join("potentials.lib").exists()

  # Now build a potential model and tabulate it - then re-run the calculation and check the energies match.

  pot_OO = Potential("O", "O", create_Multi_Range_Potential_Form(
    Multi_Range_Defn(">", 0, potentialforms.buck(9547.96, 0.2192, 32.0))))

  pot_UO = Potential("U", "O", create_Multi_Range_Potential_Form(
    Multi_Range_Defn(">", 0, potentialforms.bornmayer(1761.775, 0.35642))))

  potlist = [pot_OO, pot_UO]

  tabulation = pair_tabulation.GULP_PairTabulation(potlist, 15.0, 500)

  with tmpdir.join("potentials.lib").open("w") as potfile:
    tabulation.write(potfile)

  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  actual = extractGULPEnergy(gulp_outfile)
  assert pytest.approx(expect) == actual

  tmpdir.join("potentials.lib").remove()
  assert not tmpdir.join("potentials.lib").exists()

  # Now do the same again but use the procedural interface
  with tmpdir.join("potentials.lib").open("w") as potfile:
    writePotentials("GULP", potlist, 15.0, 500, out = potfile)

  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  actual = extractGULPEnergy(gulp_outfile)
  assert pytest.approx(expect) == actual

@needsGULP
def test_badly_conditioned_spline(tmpdir):
  """Test for GULP 'badly_conditioned_spline' error encountered when writing Basak potential"""

  gulp_input = u"""single

cell
5.468 5.468 5.468 90.0 90.0 90.0

frac
U 0 0 0
U 1/2 1/2 0
U 1/2 0 1/2
U 0 1/2 1/2

O 1/4 1/4 1/4
O 1/4 3/4 1/4
O 3/4 3/4 1/4
O 3/4 1/4 1/4

O 1/4 1/4 3/4
O 1/4 3/4 3/4
O 3/4 3/4 3/4
O 3/4 1/4 3/4


species
U 2.4
O -1.2

include potentials.lib"""

  # First calculate the expected energy using GULP's built-in analytical potentials
  with tmpdir.join("potentials.lib").open("w") as potfile:
    potfile.write("buck\n")
    potfile.write("O O 1633.01 0.3270196735 3.94879 10.0\n")
    potfile.write("U U 294.640906285709 0.327022 0.0 10.0\n")
    potfile.write("O U 693.650933805978 0.327022 0.0 10.0\n")
    potfile.write("\n")
    potfile.write("morse\n")
    potfile.write("O U 0.577189831995 1.65 2.369 10.0\n")


  gulp_infile = io.StringIO(gulp_input)
  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  expect = extractGULPEnergy(gulp_outfile)

  tmpdir.join("potentials.lib").remove()
  assert not tmpdir.join("potentials.lib").exists()

  # Now build a potential model and tabulate it - then re-run the calculation and check the energies match.
  aspot = io.StringIO(u"""
[Tabulation]
target :  GULP
cutoff : 10.0
nr : 1000

[Pair]
O-O = as.buck 1633.010242995040 0.327022 3.948787
U-U = as.buck 294.640906285709 0.327022 0.0
O-U = sum(as.buck 693.650933805978 0.327022 0.0, 
		  as.morse 1.65 2.369 0.577189831995)
"""
  )

  aspot.seek(0)

  from atsim.potentials.config import Configuration
  tabulation =  Configuration().read(aspot)

  with tmpdir.join("potentials.lib").open("w") as potfile:
    tabulation.write(potfile)

  gulp_infile.seek(0)

  gulp_outfile = io.StringIO()
  runGULP(gulp_infile, gulp_outfile, cwd = tmpdir.strpath)

  gulp_outfile.seek(0)
  actual = extractGULPEnergy(gulp_outfile)
  assert pytest.approx(expect, rel=1e-4) == actual

  tmpdir.join("potentials.lib").remove()
  assert not tmpdir.join("potentials.lib").exists()
