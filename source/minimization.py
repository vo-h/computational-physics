from typing import Callable
from source.calculus import central_diff_quotient

def minimize_by_steepest_descent(func: Callable[[float], float], guess: float, lr: float = 1e-6, tolerance: float = 1e-8):
    """Minimize a function using the steepest descent method."""
    # curr = guess
    # for _ in range(int(max_iteration)):
    #     gradient = central_diff_quotient(func, curr)
    #     new_guess = curr - lr * gradient
    #     if abs(func(new_guess) - func(curr)) < tolerance:
    #         return new_guess
    #     curr = new_guess
    # return curr

    gradient = central_diff_quotient(func, guess)
    new_guess = guess - lr * gradient
    if abs(func(new_guess) - func(guess)) < tolerance:
        return new_guess
    return minimize_by_steepest_descent(func, new_guess, lr, tolerance)
