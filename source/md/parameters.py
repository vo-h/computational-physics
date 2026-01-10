import math
from scipy import constants

def force_constant_from_IR(wave_number: float, m1: float, m2: float) -> float:
    """Compute force constant from IR spectroscopy
    wave_number: The wave number of the IR absorption peak in cm^-1
    m1: Mass of the first atom in atomic mass units (amu)
    m2: Mass of the second atom in atomic mass units (amu)
    """
    reduced_mass = (m1 * m2) / (m1 + m2)
    return (2 * math.pi * constants.c * wave_number * 100) ** 2 * reduced_mass * constants.atomic_mass