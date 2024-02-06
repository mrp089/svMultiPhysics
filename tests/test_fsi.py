<<<<<<< HEAD
import os
=======
>>>>>>> equilibrated_gr_49
import pytest

from .conftest import run_with_reference

# Common folder for all tests in this file
base_folder = "fsi"

# Fields to test
fields = ["Displacement", "Pressure", "Velocity"]


@pytest.mark.parametrize("name_inp", ["svFSI.xml", "svFSI_petsc.xml"])
def test_pipe_3d(name_inp, n_proc):
<<<<<<< HEAD
    folder = os.path.join(base_folder, "pipe_3d")
    fields = ["Displacement", "Pressure", "Velocity"]
    run_with_reference(folder, fields, n_proc, 5, name_inp=name_inp)
=======
    test_folder = "pipe_3d"
    t_max = 5
    run_with_reference(base_folder, test_folder, fields, n_proc, t_max, name_inp=name_inp)
>>>>>>> equilibrated_gr_49
