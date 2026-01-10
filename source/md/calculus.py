from typing import Callable

def central_diff_quotient(func: Callable[[float], float], x: float, h: float = 1e-11) -> float:
    """Calculate the numerical gradient of a function at a given point using central difference."""
    return (func(x + h) - func(x - h)) / (2 * h)