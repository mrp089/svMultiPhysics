import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "gr"

# Fields to test
fields = [
    "Displacement",
    "Velocity",
    "Jacobian",
    "Stress",
    "Strain",
    "Cauchy_stress",
    "Def_grad",
    "VonMises_stress",
    "GR",
]


@pytest.mark.parametrize("n_proc", [1])
def test_gr_equilibrated(n_proc):
    test_folder = "equilibrated"
    t_max = 4
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max)
