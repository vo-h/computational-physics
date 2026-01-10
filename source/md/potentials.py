import math


##### NON-BONDING POTENTIALS #####
def lennard_jones(r, epsilon=1.0, sigma=1.0):
    """
    Calculate the Lennard-Jones potential between two particles.

    Parameters:
    r (float): The distance between the two particles.
    epsilon (float): The depth of the potential well (default is 1.0).
    sigma (float): The finite distance at which the inter-particle potential is zero (default is 1.0).

    Returns:
    float: The Lennard-Jones potential energy.
    """
    if r <= 0:
        raise ValueError("Distance r must be greater than zero.")
    
    term1 = (sigma / r) ** 12
    term2 = (sigma / r) ** 6
    potential = 4 * epsilon * (term1 - term2)
    
    return potential


##### BONDING POTENTIALS #####
def morse(r, D_e=1.0, r_e=1.0, a=1.0):
    """
    Calculate the Morse potential between two atoms.

    Parameters:
    r (float): The distance between the two atoms.
    D_e (float): The depth of the potential well (default is 1.0).
    r_e (float): The equilibrium bond length (default is 1.0).
    a (float): The width of the potential well (default is 1.0).

    Returns:
    float: The Morse potential energy.
    """
    assert r > 0, f"Distance r must be greater than zero but got {r}."
    return D_e * (1 - math.exp(-a * (r - r_e))) ** 2


def harmonic(r, r0, k):
    """
    Calculate the harmonic bond potential between two atoms.

    Parameters:
    r (float): The current distance between the two atoms.
    r0 (float): The equilibrium bond length.
    k (float): The force constant of the bond.

    Returns:
    float: The harmonic bond potential energy.
    """
    assert r > 0, f"Distance r must be greater than zero but got {r}."
    return 0.5 * k * (r - r0) ** 2