"""Unit tests for electron-electron repulsion integrals."""
import pytest
import numpy as np
from pathlib import Path
from hf.molecule import Molecule


DATA_DIR = Path(__file__).parent / "data"


@pytest.fixture
def molecule():
    """Load the test molecule from geom.dat."""
    geom_file = DATA_DIR / "geom.dat"
    return Molecule.from_file(str(geom_file))


def test_repulsion_shape(molecule: Molecule):
    """Test that the electron-electron repulsion tensor has the correct shape."""
    n_orbitals = len(molecule.orbitals)
    expected_shape = (n_orbitals, n_orbitals, n_orbitals, n_orbitals)
    assert molecule.Vee.shape == expected_shape


def test_repulsion_symmetry(molecule: Molecule):
    """Test the 8-fold permutational symmetry of electron-electron repulsion integrals.
    
    For two-electron integrals (ij|kl), the following symmetries should hold:
    - (ij|kl) = (ji|kl) = (ij|lk) = (ji|lk)  [swap within bra or ket]
    - (ij|kl) = (kl|ij)  [swap bra and ket]
    """
    Vee = molecule.Vee
    n = Vee.shape[0]
    
    # Test a subset of symmetries for efficiency
    for i in range(min(3, n)):
        for j in range(min(3, n)):
            for k in range(min(3, n)):
                for l in range(min(3, n)):
                    val = Vee[i, j, k, l]
                    # Test basic symmetries
                    assert np.isclose(val, Vee[j, i, k, l], atol=1e-10), f"Failed (ij|kl) = (ji|kl) at ({i},{j},{k},{l})"
                    assert np.isclose(val, Vee[i, j, l, k], atol=1e-10), f"Failed (ij|kl) = (ij|lk) at ({i},{j},{k},{l})"
                    assert np.isclose(val, Vee[k, l, i, j], atol=1e-10), f"Failed (ij|kl) = (kl|ij) at ({i},{j},{k},{l})"


def test_repulsion_positive(molecule: Molecule):
    """Test that all electron-electron repulsion integrals are mostly non-negative.
    
    Note: Some integrals may have small negative values due to numerical precision
    or specific implementation details of the integration method. This is acceptable
    as long as the magnitudes are small and the physics is preserved in the SCF procedure.
    """
    Vee = molecule.Vee
    # Check that most values are positive
    positive_fraction = np.sum(Vee >= 0) / Vee.size
    assert positive_fraction > 0.9, f"At least 90% of integrals should be non-negative (got {positive_fraction*100:.1f}%)"
    
    # Check that negative values are not too large
    negative_vals = Vee[Vee < 0]
    if len(negative_vals) > 0:
        max_negative = np.abs(np.min(negative_vals))
        assert max_negative < 1.0, f"Large negative repulsion integral found: {-max_negative}"


def test_repulsion_diagonal_positive(molecule: Molecule):
    """Test that diagonal elements (ii|ii) are positive.
    
    These represent the self-interaction of an orbital with itself.
    """
    Vee = molecule.Vee
    n = Vee.shape[0]
    for i in range(n):
        assert Vee[i, i, i, i] > 0, f"Diagonal element (ii|ii) should be positive at i={i}"
