import os
import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "fsi"

# Fields to test
fields = ["Displacement", "Pressure", "Velocity"]


@pytest.mark.parametrize("name_inp", ["svFSI.xml", "svFSI_petsc.xml"])
def test_pipe_3d(name_inp, n_proc):
    folder = os.path.join(base_folder, "pipe_3d")
    fields = ["Displacement", "Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 5, name_inp=name_inp)
