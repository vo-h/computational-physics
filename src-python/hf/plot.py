"""
Molecular orbital visualization utilities for 3D plotting using Plotly.

This module provides functions to visualize molecular orbitals after Hartree-Fock
calculations, including isosurface plots and electron density visualizations.
"""

from typing import Tuple, Optional, Dict
import numpy as np
import plotly.graph_objects as go
import plotly.colors
from numpy.typing import NDArray
from hf.molecule import Molecule
from mendeleev import get_all_elements

def create_3d_grid(
    molecule: Molecule,
    grid_size: int = 50,
    padding: float = 3.0
) -> Tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """
    Create a 3D grid around the molecule for orbital evaluation.
    
    Parameters
    ----------
    molecule : Molecule
        Molecule object containing atomic coordinates
    grid_size : int, optional
        Number of points along each axis (default: 50)
    padding : float, optional
        Extra space around the molecule in Bohr (default: 3.0)
    
    Returns
    -------
    X : NDArray[np.float64]
        3D meshgrid X coordinates
    Y : NDArray[np.float64]
        3D meshgrid Y coordinates
    Z : NDArray[np.float64]
        3D meshgrid Z coordinates
    """
    # Get molecule bounds
    coords = np.array([atom.coords for atom in molecule.atoms])
    x_min, y_min, z_min = coords.min(axis=0) - padding
    x_max, y_max, z_max = coords.max(axis=0) + padding
    
    # Create coordinate ranges
    x_range = np.linspace(x_min, x_max, grid_size)
    y_range = np.linspace(y_min, y_max, grid_size)
    z_range = np.linspace(z_min, z_max, grid_size)
    
    # Create 3D meshgrid
    X, Y, Z = np.meshgrid(x_range, y_range, z_range, indexing='ij')
    
    return X, Y, Z


def compute_mo_on_grid(
    molecule: Molecule,
    mo: int,
    X: NDArray[np.float64],
    Y: NDArray[np.float64],
    Z: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Compute molecular orbital values on a 3D grid.
    
    Parameters
    ----------
    molecule : Molecule
        Molecule object with orbital coefficients
    mo : int
        Index of the molecular orbital to compute
    X : NDArray[np.float64]
        3D meshgrid X coordinates
    Y : NDArray[np.float64]
        3D meshgrid Y coordinates
    Z : NDArray[np.float64]
        3D meshgrid Z coordinates
    
    Returns
    -------
    NDArray[np.float64]
        3D array of molecular orbital values at each grid point
    """
    values = np.zeros_like(X)
    
    # Iterate over grid points
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            for k in range(X.shape[2]):
                values[i, j, k] = molecule(X[i, j, k], Y[i, j, k], Z[i, j, k], mo)
    return values


def get_atom_color_map() -> Dict[str, str]:
    """
    Get color mapping for atomic elements using Plotly's color palette.
    
    Returns
    -------
    Dict[str, str]:
        Dictionary mapping element symbols to Plotly colors
    """
    elements = [elem.symbol for elem in get_all_elements()]
    colors = plotly.colors.qualitative.Plotly
    return {elem: colors[i % len(colors)] for i, elem in enumerate(elements)}


def plot_molecular_orbital(
    molecule: Molecule,
    mo: int,
    grid_size: int = 40,
    padding: float = 3.0,
    isovalue: float = 0.1,
    show_atoms: bool = True,
    title: Optional[str] = None
) -> go.Figure:
    """
    Plot a molecular orbital in 3D using Plotly isosurfaces.
    
    Creates an interactive 3D visualization showing positive (blue) and negative (red)
    phases of the molecular orbital using isosurface plots. Optionally displays atomic
    positions.
    
    Parameters
    ----------
    molecule : Molecule
        Molecule object with SCF-calculated orbital coefficients
    mo : int
        Index of the molecular orbital to plot (0-indexed)
    grid_size : int, optional
        Number of grid points along each axis. Lower values are faster but less
        smooth (default: 40)
    padding : float, optional
        Extra space around molecule in Bohr (default: 3.0)
    isovalue : float, optional
        Isosurface threshold value. Smaller values show more diffuse orbitals,
        larger values show more compact orbitals (default: 0.02)
    show_atoms : bool, optional
        Whether to display atomic positions as spheres (default: True)
    title : str, optional
        Custom title for the plot. If None, generates default title
    
    Returns
    -------
    go.Figure
        Plotly figure object that can be displayed with fig.show() or saved
        
    Examples
    --------
    >>> from hf.molecule import Molecule
    >>> from hf.atom import Atom
    >>> 
    >>> # Create a water molecule
    >>> atoms = [
    ...     Atom(atom="O", coords=(0.0, 0.0, 0.0), basis='STO-3G'),
    ...     Atom(atom="H", coords=(0.757, 0.586, 0.0), basis='STO-3G'),
    ...     Atom(atom="H", coords=(-0.757, 0.586, 0.0), basis='STO-3G'),
    ... ]
    >>> mol = Molecule(atoms=atoms)
    >>> 
    >>> # Plot HOMO (highest occupied molecular orbital)
    >>> n_electrons = sum(atom.Z for atom in mol.atoms)
    >>> homo_index = n_electrons // 2 - 1
    >>> fig = plot_molecular_orbital(mol, homo_index, grid_size=30)
    >>> fig.show()
    """
    print(f"Computing molecular orbital {mo} for {len(molecule.atoms)}-atom molecule...")
    
    # Create 3D grid
    X, Y, Z = create_3d_grid(molecule, grid_size, padding)
    
    # Compute MO values on grid
    print("Evaluating orbital on grid...")
    mo_values = compute_mo_on_grid(molecule, mo, X, Y, Z)
    
    # Create figure
    fig = go.Figure()
    
    # Add positive isosurface (blue)
    if mo_values.max() >= isovalue:
        fig.add_trace(go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=mo_values.flatten(),
        isomin=isovalue,
        isomax=mo_values.max(),
        surface_count=1,
        colorscale='Blues',
        showscale=False,
        caps=dict(x_show=False, y_show=False, z_show=False),
        name='Positive'
        ))
    
    # Add negative isosurface (red)
    if mo_values.min() <= -isovalue:
        fig.add_trace(go.Isosurface(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=mo_values.flatten(),
        isomin=mo_values.min(),
        isomax=-isovalue,
        surface_count=1,
        colorscale='Reds',
        showscale=False,
        caps=dict(x_show=False, y_show=False, z_show=False),
        name='Negative'
        ))
    
    # Add atoms as spheres
    if show_atoms:
        atom_coords = np.array([atom.coords for atom in molecule.atoms])
        atom_symbols = [atom.atom for atom in molecule.atoms]
        
        # Get colors for atoms
        color_map = get_atom_color_map()
        colors = [color_map.get(symbol, 'gray') for symbol in atom_symbols]
        
        fig.add_trace(go.Scatter3d(
            x=atom_coords[:, 0],
            y=atom_coords[:, 1],
            z=atom_coords[:, 2],
            mode='markers+text',
            marker=dict(size=10, color=colors, line=dict(color='black', width=2)),
            text=atom_symbols,
            textposition='top center',
            name='Atoms'
        ))
    
    # Determine orbital type for title
    n_electrons = sum(atom.Z for atom in molecule.atoms)
    n_occupied = n_electrons // 2
    
    if title is None:
        orbital_type = ""
        if mo == n_occupied - 1:
            orbital_type = " (HOMO)"
        elif mo == n_occupied:
            orbital_type = " (LUMO)"
        elif mo < n_occupied:
            orbital_type = " (Occupied)"
        else:
            orbital_type = " (Virtual)"
        
        title = f'Molecular Orbital {mo}{orbital_type} (Isovalue: ±{isovalue})'
    
    # Update layout
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X (Bohr)',
            yaxis_title='Y (Bohr)',
            zaxis_title='Z (Bohr)',
            aspectmode='data'
        ),
        width=800,
        height=800
    )
    
    print("Done!")
    return fig