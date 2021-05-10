import io

from pytest import approx, raises, fail

from ._runlammps import extractLAMMPSEnergy, needsLAMMPS, runLAMMPS

from atsim.potentials.config import Configuration
from atsim.potentials.config._common import ConfigurationException

from .config._common import _get_lammps_resource_dir
from atsim.potentials.eam_tabulation import ADP_EAMTabulation


def test_adp_tabulation(tmpdir):
    cfg_file = _get_lammps_resource_dir().join("Al_Cu_adp.aspot")
    configuration = Configuration()
    tabulation = configuration.read(cfg_file.open('r'))

    assert type(tabulation) is ADP_EAMTabulation

    # assert "eam_adp" == tabulation.tabulation_type
    assert 10000 == tabulation.nrho
    assert tabulation.drho == approx(2.2770502180000001e-03)
    assert 10000 == tabulation.nr
    assert tabulation.dr == approx(6.2872099999999995e-04)

    assert 2 == len(tabulation.eam_potentials)

    eam_dict = dict([(ep.species, ep) for ep in tabulation.eam_potentials])

    assert ['Al', 'Cu'] == sorted(eam_dict.keys())
    assert eam_dict['Al'].embeddingValue(0.1001902096E+01) == approx(-0.2418739157E+01)
    assert eam_dict['Cu'].embeddingValue(0.1009492263E+01) == approx(-0.9390266026E+00)

    assert eam_dict['Al'].electronDensity(0.1609791553E+01) == approx(0.7617792157E-01)
    assert eam_dict['Cu'].electronDensity(0.1512859370E+01) == approx(0.1824604072E+00)

    assert 3 == len(tabulation.potentials)
    assert [('Al', 'Al'), ('Al','Cu'), ('Cu', 'Cu')] == sorted([tuple(sorted((p.speciesA, p.speciesB))) for p in tabulation.potentials])

    pair_dict = dict([(tuple(sorted((p.speciesA, p.speciesB))), p) for p in tabulation.potentials])
    assert pair_dict[('Al', 'Al')].energy(0.1123368233E+01) == approx(0.5178467078E+01)
    assert pair_dict[('Al', 'Cu')].energy(0.1500522547E+01) == approx(0.2913519771E+01)
    assert pair_dict[('Cu', 'Cu')].energy(0.1817755147E+01) == approx(0.2103976024E+01)

    assert 1 == len(tabulation.dipole_potentials)
    pair_dict = dict([(tuple(sorted((p.speciesA, p.speciesB))), p) for p in tabulation.dipole_potentials])
    assert pair_dict[('Al', 'Cu')].energy(0.1634465200E+01) == approx(-0.5117209620E-01)

    assert 1 == len(tabulation.quadrupole_potentials)
    pair_dict = dict([(tuple(sorted((p.speciesA, p.speciesB))), p) for p in tabulation.quadrupole_potentials])
    assert pair_dict[('Al', 'Cu')].energy(0.5044715650E+01) == approx(-0.4386367758E-02)


def test_throw_configuration_exception():
    cfg = u"""[Tabulation]
target : eam_adp
drho : 2.2770502180000001e-03
nrho : 10000

nr : 10000
dr : 6.2872099999999995e-04

[EAM-Embed]
Al : as.zero
Cu : as.zero

[EAM-Density]
Al : as.zero
Cu : as.zero

[Pair]
Al-Al : as.zero
Cu-Al : as.zero
Cu-Cu : as.zero
"""

    configuration = Configuration()

    with raises(ConfigurationException):
        configuration.read(io.StringIO(cfg))

    cfg = u"""[Tabulation]
target : eam_adp
drho : 2.2770502180000001e-03
nrho : 10000

nr : 10000
dr : 6.2872099999999995e-04

[EAM-Embed]
Al : as.zero
Cu : as.zero

[EAM-Density]
Al : as.zero
Cu : as.zero

[Pair]
Al-Al : as.zero
Cu-Al : as.zero
Cu-Cu : as.zero

[EAM-ADP-Dipole]
Al-Al : as.zero
Al-Cu : as.zero
Cu-Cu : as.zero

"""

    with raises(ConfigurationException):
        configuration.read(io.StringIO(cfg))

    cfg = u"""[Tabulation]
target : eam_adp
drho : 2.2770502180000001e-03
nrho : 10000

nr : 10000
dr : 6.2872099999999995e-04

[EAM-Embed]
Al : as.zero
Cu : as.zero

[EAM-Density]
Al : as.zero
Cu : as.zero

[Pair]
Al-Al : as.zero
Cu-Al : as.zero
Cu-Cu : as.zero

[EAM-ADP-Quadrupole]
Al-Al : as.zero
Al-Cu : as.zero
Cu-Cu : as.zero

"""

    with raises(ConfigurationException):
        configuration.read(io.StringIO(cfg))


@needsLAMMPS
def test_adp_in_lammps(tmpdir):
    expect_dir = tmpdir.join("expect")
    expect_dir.ensure(dir=True)
    
    lmpin =  _get_lammps_resource_dir().join("AlCu3.lmpstruct")
    lmpin.copy(expect_dir.join("structure.lmpstruct"))

    lmpin =  _get_lammps_resource_dir().join("calc_energy.lmpin")
    lmpin.copy(expect_dir.join("calc_energy.lmpin"))

    # Copy existing table file
    _get_lammps_resource_dir().join("AlCu.adp").copy(expect_dir.join("AlCu.adp"))

    with expect_dir.join("potentials.lmpinc").open("w") as potentials:
        potentials.write("pair_style adp\n")
        potentials.write("pair_coeff * * AlCu.adp Al Cu\n")

    runLAMMPS(cwd=expect_dir.strpath)
    expect_energy = extractLAMMPSEnergy(cwd=expect_dir.strpath)

    actual_dir = tmpdir.join("actual")
    actual_dir.ensure(dir=True)

    lmpin =  _get_lammps_resource_dir().join("AlCu3.lmpstruct")
    lmpin.copy(actual_dir.join("structure.lmpstruct"))

    lmpin =  _get_lammps_resource_dir().join("calc_energy.lmpin")
    lmpin.copy(actual_dir.join("calc_energy.lmpin"))

    with actual_dir.join("potentials.lmpinc").open("w") as potentials:
        potentials.write("pair_style adp\n")
        potentials.write("pair_coeff * * AlCu.adp Al Cu\n")

    # Tabulate potential
    cfg_file = _get_lammps_resource_dir().join("Al_Cu_adp.aspot")
    configuration = Configuration()
    tabulation = configuration.read(cfg_file.open('r'))
    with actual_dir.join("AlCu.adp").open("w") as outfile:
        tabulation.write(outfile)

    runLAMMPS(cwd=actual_dir.strpath)
    actual_energy = extractLAMMPSEnergy(cwd=actual_dir.strpath)

    assert expect_energy == approx(actual_energy)
