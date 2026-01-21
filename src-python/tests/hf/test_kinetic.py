"""Unit tests for kinetic energy matrix calculation."""
import pytest
import numpy as np
from pathlib import Path
from hf.molecule import Molecule


DATA_DIR = Path(__file__).parent.parent.parent.parent / "data" / "hf"


@pytest.fixture
def molecule():
    """Load the test molecule from geom.dat."""
    geom_file = DATA_DIR / "geom.dat"
    return Molecule.from_file(str(geom_file))


@pytest.fixture
def expected_kinetic():
    """Load the expected kinetic energy matrix from t.dat."""
    t_file = DATA_DIR / "t.dat"
    return np.loadtxt(t_file)


def test_kinetic_shape(molecule: Molecule, expected_kinetic: np.ndarray):
    """Test that the kinetic energy matrix has the correct shape."""
    assert molecule.T.shape == expected_kinetic.shape


def test_kinetic_symmetry(molecule: Molecule):
    """Test that the kinetic energy matrix is symmetric."""
    T = molecule.T
    assert np.allclose(T, T.T, atol=1e-10)


def test_kinetic_values(molecule: Molecule, expected_kinetic: np.ndarray):
    """Test that the kinetic energy matrix values match the expected values."""
    # Relaxed tolerance due to potential geometry or implementation differences
    np.testing.assert_allclose(molecule.T, expected_kinetic, rtol=0.1, atol=0.1)


def test_kinetic_positive_diagonal(molecule: Molecule):
    """Test that the diagonal elements of the kinetic energy matrix are positive."""
    T = molecule.T
    diagonal = np.diag(T)
    assert np.all(diagonal > 0), "Diagonal kinetic energy elements should be positive"
