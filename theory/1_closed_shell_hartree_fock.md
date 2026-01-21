# Hartree-Fock Theory for Closed-Shell Systems
Electronic structure theory, at its core, tries to solve

$$
\hat{H}\Psi = E\Psi
$$

where $\hat{H}$ is the Hamiltonian operator, $\Psi$ is the wavefunction of the multi-electron system, and $E$ is the energy eigenvalue. The Hartree-Fock method guesses $\Psi$ takes the below form (where N is the number of electrons; $1/\sqrt{N!}$ is a normalization constant; $\chi_i(j)$ is the molecular spin-orbital i evaluated at electron j):

$$
\Psi = \frac{1}{\sqrt{N!}}
\begin{vmatrix}
\chi_1(1) & \chi_2(1) & \cdots & \chi_N(1) \\
\chi_1(2) & \chi_2(2) & \cdots & \chi_N(2) \\
\vdots & \vdots & \ddots & \vdots \\
\chi_1(N) & \chi_2(N) & \cdots & \chi_N(N)
\end{vmatrix}
$$

for several reasons:
1. The determinant ensures that the wavefunction is antisymmetric with respect to the exchange of any two electrons, satisfying the Pauli exclusion principle, because exchanging any two electrons is equivalent to swapping two rows of the determinant, which changes its sign.
2. Writing the wavefunction as a sum of products simplifies the mathematical treatment of many-electron systems (at the cost of ignoring electron correlation effects).

## Born-Oppenheimer Approximation
In atomic units, we can write the Hamiltonian for the whole system as:

$$
\hat{H} = -\frac{1}{2}\sum_{i}\nabla_i^2 - \frac{1}{2} \sum_{A}\nabla_A^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i<j} \frac{1}{r_{ij}} + \sum_{A<B} \frac{Z_A Z_B}{R_{AB}}
$$  

where capital letters (A, B) denote the nuclei and lowercase letters (i, j) denote the electrons. The first two terms represent the kinetic energy of the electrons and nuclei, respectively. The third term represents the electron-nucleus attraction, the fourth term represents the electron-electron repulsion, and the last term represents the nucleus-nucleus repulsion.

The Born-Oppenheimer approximation simplifies the problem by assuming that the nuclei are stationary relative to the electrons due to their much larger mass. This allows us to neglect the nuclear kinetic energy term, leading to a simplified Hamiltonian:

$$
\hat{H} = -\frac{1}{2}\sum_{i}\nabla_i^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i<j} \frac{1}{r_{ij}} + \sum_{A<B} \frac{Z_A Z_B}{R_{AB}} = -\frac{1}{2}\sum_{i}\nabla_i^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i<j} \frac{1}{r_{ij}} + V_{NN}
$$

Ignoring the nuclear repulsion term $V_{NN}$ for now, we can compute the Hartree-Fock version of the electronic energy as the expectation value of the Hamiltonian with respect to the wavefunction:

$$
E_{HF} 
= \langle \Psi | \hat{H} | \Psi \rangle 
= \langle \Psi | -\frac{1}{2}\sum_{i}\nabla_i^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i<j} \frac{1}{r_{ij}} | \Psi \rangle 
= \sum_{i=1}^{N} \langle \chi_i | \hat{h}(i) | \chi_i \rangle + \sum_{i<j}^{N} \left[ \langle \chi_i \chi_j | \hat{J}_j | \chi_i \chi_j \rangle - \langle \chi_i \chi_j | \hat{K}_j | \chi_i \chi_j \rangle \right]
$$

where 

$$ \hat{h}(i) = -\frac{1}{2}\nabla_i^2 + \sum_{A} \frac{Z_A}{r_{iA}} $$

$$
\langle \chi_i(\vec{x}_1) \chi_j (\vec{x}_2) | \hat{J}_j | \chi_i (\vec{x}_1) \chi_j (\vec{x}_2) \rangle 
= \langle \chi_i (\vec{x}_1) \chi_j (\vec{x}_2) | \frac{1}{r_{12}} | \chi_i (\vec{x}_1) \chi_j (\vec{x}_2) \rangle 
$$

$$
\langle \chi_i(\vec{x}_1) \chi_j (\vec{x}_2) | \hat{K}_j | \chi_i (\vec{x}_1) \chi_j (\vec{x}_2) \rangle 
= \langle \chi_i (\vec{x}_1) \chi_j (\vec{x}_2) | \frac{1}{r_{12}} | \chi_j (\vec{x}_1) \chi_i (\vec{x}_2) \rangle 
$$

## Constrained Optimization
In order to find the ground state of the system, we need to minimize $E_{HF}$ with respect to the molecular orbitals $\chi_i$, subject to the constraint that the orbitals remain orthonormal $\langle \chi_i | \chi_j \rangle = \delta_{ij}$, which is what Lagrange multipliers are for (spoiler alert—the multipliers will turn out to be the orbital energies $\epsilon_i$):

$$
\mathcal{L}[\{\chi_i\}] = E_{HF}[\{\chi_i\}]  - \sum_{i,j} \epsilon_{ij} \left( \langle \chi_i | \chi_j \rangle - \delta_{ij} \right)
$$

After substituting $E_{HF}$ into the Lagrangian and taking the functional derivative with respect to $\chi_i$, we arrive at the non-canonical (i.e., non-eigenvalue) Hartree-Fock equations:

$$
\hat{F} |\chi_i \rangle = \sum_j^N \epsilon_{ij} |\chi_i \rangle
$$

## Unitary Transformation to Canonical Form
To get to the canonical form, we can define a unitary transformation matrix U that transforms the molecular orbitals $\chi_i$ into a new set of orbitals $\chi_j'$ while preserving their orthonormality:

$$
|\chi_j' \rangle = \sum_j U_{lj} |\chi_j \rangle \rightarrow 
\langle  \chi_j'| \hat{F} | \chi_i' \rangle 

= \sum_k^N \sum_i^N U_{ki}^* \epsilon_{kl} U_{lj} 
= (U^{\dagger} \epsilon U)_{ij} = \epsilon_{ij}'

$$

By choosing U such that it diagonalizes the matrix of Lagrange multipliers $\epsilon$, we can rewrite the Hartree-Fock equations in their canonical form:

$$
\hat{F} |\chi_i' \rangle = \epsilon_i' |\chi_i' \rangle
$$

## Basis Set Expansion and Roothaan-Hall Equations
In practice, we express the molecular orbitals as linear combinations of a finite set of known atomic basis functions $\{\phi_\mu\}$:

$$
|\chi_i \rangle = \sum_{\mu} c_{\mu i} |\phi_\mu \rangle 
\rightarrow \hat{F} |\sum_{\nu} c_{\nu i} \phi_\nu \rangle = \epsilon_i |\sum_{\nu} c_{\nu i} \phi_\nu \rangle
$$

Multiplying both sides by $\langle \phi_\mu |$ and rearranging, we obtain the Roothaan-Hall equations:

$$
\sum_{\nu} F_{\mu \nu} c_{\nu i} = \epsilon_i \sum_{\nu} S_{\mu \nu} c_{\nu i}
\rightarrow \mathbf{F} \mathbf{c} = \mathbf{S} \mathbf{c} \mathbf{\epsilon}
$$

where $\mathbf{F}$ is the Fock matrix with elements $F_{\mu \nu} = \langle \phi_\mu | \hat{F} | \phi_\nu \rangle$, $\mathbf{S}$ is the overlap matrix with elements $S_{\mu \nu} = \langle \phi_\mu | \phi_\nu \rangle$, $\mathbf{c}$ is the coefficient matrix with elements $c_{\nu i}$, and $\mathbf{\epsilon}$ is the diagonal matrix of orbital energies $\epsilon_i$. To once again obtain an eigenvalue problem, we can apply a transformation using the symmetric orthogonalization matrix $\mathbf{S}^{-1/2}$:

$$
\mathbf{F}' = \mathbf{S}^{-1/2} \mathbf{F} \mathbf{S}^{-1/2}
$$

$$
\mathbf{c}' = (\mathbf{S}^{-1/2})^{-1} \mathbf{c} 
$$

$$
\mathbf{c} = \mathbf{S}^{-1/2} \mathbf{c}'
$$

Applying these transformations to the Roothaan-Hall equations, we get:

$$
\mathbf{F}' \mathbf{c}' = \mathbf{c}' \mathbf{\epsilon}
$$



# Resource(s)
* [In-depth math behind Hartree-Fock method](https://azchemie.com/Hartree_Fock.pdf)
* [Cheat Sheet on Hartree-Fock Algorithm](https://github.com/CrawfordGroup/ProgrammingProjects/tree/master/Project%2303)
* Hartree-Fock:
    * [Computed Basis Sets For Most Elements](https://www.basissetexchange.org/)
    * Example Calculations of Water
        * [Using Gauss09](https://chemistry.montana.edu/callis/courses/chmy564/460hf.pdf)
        * [Analytical Evaluation of Overlap Integrals](https://content.wolfram.com/sites/19/2012/02/Ho.pdf)
        * [Analytical Evaluation of Kinetic Energy Integrals](https://content.wolfram.com/sites/19/2013/01/Ho_Kinetic.pdf)
        * [Analytical Evaluation of Nuclear Attraction Integrals](https://content.wolfram.com/sites/19/2014/12/Ho_Nuclear.pdf)
        * [Survey on Integral Evaluation Methods](https://rsc.anu.edu.au/~pgill/papers/045Review.pdf)
    * [PySCF](https://pyscf.org/) & [Associated Tutorials](https://pyscf.org/pyscf_tutorials/index.html) for comparing results against.
* https://kthpanor.github.io/echem/docs/mol_struct/hessians.html
* [Project Ideas in Computational Physics/Chemistry](https://github.com/CrawfordGroup/ProgrammingProjects)

