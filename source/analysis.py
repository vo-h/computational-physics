import numpy as np
import plotly.express as px
import pandas as pd
from typing import Callable
from plotly.graph_objs import Figure


def get_coords(plane: tuple[float, float], start: int = 0, stop: int = 10, step: float = 0.1, index=2) -> np.ndarray:
    """Generate a list of coordinates along a plane"""
    points = np.linspace(start, stop, int((stop - start) / step) + 1) 
    planes = np.vstack([np.array(plane)] * len(points))
    return np.insert(planes, index, points, axis=1)

def get_scatter(coords: np.ndarray, func1: tuple[str, Callable[[np.ndarray], np.ndarray]], func2: tuple[str, Callable[[np.ndarray], np.ndarray]] = None, index: int = 2, orbital: int = 2) -> Figure:
    """Generate a DataFrame for comparing the results of two functions."""
    
    func1_results: np.ndarray = func1[1](coords=coords)
    func1_results = pd.DataFrame({'x': coords[:,index], 'y': func1_results[:,orbital], 'method': func1[0]})

    if func2 is None:
        return px.scatter(func1_results, x='x', y='y', color='method', title='AO Density along axis')
    
    func2_results: np.ndarray = func2[1](coords=coords)
    func2_results = pd.DataFrame({'x': coords[:,index], 'y': func2_results[:,orbital], 'method': func2[0]})
    return px.scatter(pd.concat([func1_results, func2_results], ignore_index=True), x='x', y='y', color='method', title='Comparison of AO Density along axis')