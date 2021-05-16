import pytest

from atsim.potentials.config._config_parser import _RawConfigParser

import io

def test_variables():
    tst_cfg = """[Variables]
A = 1
B = 2
C = 3

[Species]
Al.charge = 1.2

[Pair]
Al-Cu = ${A} + ${B} + ${C}
A-B = ${Species:Al.charge}

"""

    parser = _RawConfigParser()
    parser.read_file(io.StringIO(tst_cfg))

    assert parser["Pair"]["Al-Cu"] == "1 + 2 + 3"
    assert parser["Pair"]["A-B"] == "1.2"


