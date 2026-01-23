import click
from hf.molecule import Molecule
import numpy as np
import sys

@click.command()
@click.option('--input-file', type=click.Path(exists=True), help='Input file for the molecule')
@click.option('--method', type=click.Choice(['1e-mat', '2e-mat', 'guesses', 'scf']), default='scf', help='Method to use')
def hf(input_file: str, method: str):
    """Perform SCF on a molecule. Have options for input file, and method used like compute overlap, hamiltonian and final SCF."""
    
    # Red input file and create molecule object
    molecule = Molecule.from_file(input_file)
    with np.printoptions(linewidth=sys.maxsize):
        if method == '1e-mat':
            S = molecule.S
            T = molecule.T
            V = molecule.Vne
            print(f"Overlap matrix:\n{S}\n")
            print(f"Kinetic energy matrix:\n{T}\n")
            print(f"Nuclear attraction matrix:\n{V}")
        elif method == '2e-mat':
            ERI = molecule.Vee
            with np.printoptions(linewidth=sys.maxsize):
                for i in range(ERI.shape[0]):
                    for j in range(ERI.shape[1]):
                        print(f"Two-electron repulsion integrals (ERI) tensor slice [{i},{j},:,:]:\n{ERI[i,j,:,:]}\n")

        elif method == 'guesses':
            F0 = molecule.F0
            C0 = molecule.C0
            D0 = molecule.D0
            E0 = molecule.E0
            print(f"Initial Fock matrix:\n{F0}\n")
            print(f"Initial coefficients:\n{C0}\n")
            print(f"Initial density matrix:\n{D0}\n")
            print(f"Initial electronic energy: {E0}")
        elif method == 'scf':
            molecule.CHF()
            print("SCF converged. Final energy: {:.6f} Hartree".format(molecule.E))
