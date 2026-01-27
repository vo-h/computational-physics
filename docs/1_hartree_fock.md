# Hartree-Fock Theory

*   [Theory](#Theory)
    *   [The Slater Determinant $`\Psi`$ ](#the-slater-determinant-psi)
    *   [The Born-Oppenheimer Approximation of $`\hat{H}`$ ](#the-born-oppenheimer-approximation-of-hath)
    *   [Constrained Optimization To Obtain Ground State Energy: $`E_{HF}`$ ](#constrained-optimization-to-obtain-ground-state-energy-e_hf)
    *   [Basis Set Expansion and Roothaan-Hall Equations](#Basis-Set-Expansion-and-Roothaan-Hall-Equations)
*   [Computational Implementation](#Computational-Implementation)
    *   [Contracted Gaussian Basis Sets: STO-nG](#Contracted-Gaussian-Basis-Sets-STO-nG)
    *   [Overlap, Kinetic, and Nuclear Attraction Integrals](#Overlap,-Kinetic,-and-Nuclear-Attraction-Integrals)
    *   [Two-Electron Repulsion Integrals](#Two-Electron-Repulsion-Integrals)
    *   [Self-Consistent Field (SCF) Procedure](#Self-Consistent-Field-(SCF)-Procedure)

## Theory
Electronic structure theory, at its core, tries to solve (where $\hat{H}$ is the Hamiltonian operator, $\Psi$ is the wavefunction of the multi-electron system, and $E$ is the energy eigenvalue):

$$
\hat{H}\Psi = E\Psi
$$

### The Slater Determinant $\Psi$
The Hartree-Fock method guesses $\Psi$ takes the below form (where N is the number of electrons; $1/\sqrt{N!}$ is a normalization constant; $\chi_i(j)$ is the molecular spin-orbital i evaluated at electron j):

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
1. The determinant ensures that the wavefunction is antisymmetric with respect to the exchange of any two electrons, satisfying the Pauli exclusion principle (exchanging any two electrons is equivalent to swapping two rows of the determinant, which changes its sign).
2. Writing the wavefunction as a sum of products simplifies the mathematical treatment of many-electron systems (at the cost of ignoring electron correlation effects).


### The Born-Oppenheimer Approximation of $\hat{H}$

In atomic units, we can write the Hamiltonian for the whole system as:

$$
\hat{H} = -\frac{1}{2}\sum_{i}\nabla_i^2 - \frac{1}{2} \sum_{A}\nabla_A^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i\lt j} \frac{1}{r_{ij}} + \sum_{A\lt B} \frac{Z_A Z_B}{R_{AB}}
$$  

where capital letters (A, B) denote the nuclei and lowercase letters (i, j) denote the electrons. The first two terms represent the kinetic energy of the electrons and nuclei, respectively. The third term represents the electron-nucleus attraction, the fourth term represents the electron-electron repulsion, and the last term represents the nucleus-nucleus repulsion.

The Born-Oppenheimer approximation simplifies the problem by assuming that the nuclei are stationary relative to the electrons due to their much larger mass. This allows us to neglect the nuclear kinetic energy term, leading to a simplified Hamiltonian:

$$
\hat{H} = -\frac{1}{2}\sum_{i}\nabla_i^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i\lt j} \frac{1}{r_{ij}} + \sum_{A\lt B} \frac{Z_A Z_B}{R_{AB}} = -\frac{1}{2}\sum_{i}\nabla_i^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i\lt j} \frac{1}{r_{ij}} + V_{NN}
$$

### Constrained Optimization To Obtain Ground State Energy $E_{HF}$

Ignoring the trivial nuclear repulsion term $V_{NN}$ for now, we can the Hartree-Fock electronic energy expression is:

$$
E_{HF} 
= \langle \Psi | \hat{H} | \Psi \rangle 
= \langle \Psi | -\frac{1}{2}\sum_{i}\nabla_i^2 + \sum_{i,A} \frac{Z_A}{r_{iA}} + \sum_{i\lt j} \frac{1}{r_{ij}} | \Psi \rangle 
= \sum_{i=1}^{N} \langle \chi_i | \hat{h}(i) | \chi_i \rangle + \sum_{i\lt j}^{N} \left[ \langle \chi_i| \hat{J}_j | \chi_i \rangle - \langle \chi_i| \hat{K}_j | \chi_i \rangle \right]
$$

where 

$$ \hat{h}(i) = -\frac{1}{2}\nabla_i^2 + \sum_{A} \frac{Z_A}{r_{iA}} $$

$$
\hat{J}_j | \chi_i (\vec{x}_1) \rangle 
= \left[ \int \chi_j^* (\vec{x}_2) \frac{1}{r_{12}} \chi_j (\vec{x}_2) \right] \chi_i (\vec{x}_1) d\vec{x}_2
$$

$$
\hat{K}_j | \chi_i (\vec{x}_1) \rangle 
= \left[ \int \chi_j^* (\vec{x}_2) \frac{1}{r_{12}} \chi_i (\vec{x}_2) \right] \chi_j (\vec{x}_1) d\vec{x}_2
$$

In order to find the ground state of the system, we need to minimize $E_{HF}$ with respect to the molecular orbitals $\chi_i$, subject to orthonormality $\langle \chi_i | \chi_j \rangle = \delta_{ij}$, which is what Lagrange multipliers are for (spoiler alert—the multipliers will turn out to be the orbital energies $\epsilon_i$):

$$
\mathcal{L}[\lbrace\chi_i\rbrace] = E_{HF}[\lbrace\chi_i\rbrace]  - \sum_{i,j} \epsilon_{ij} \left( \langle \chi_i | \chi_j \rangle - \delta_{ij} \right)
$$

After substituting $E_{HF}$ into the Lagrangian and taking the functional derivative with respect to $\chi_i$, we arrive at the non-canonical (i.e., non-eigenvalue) Hartree-Fock equations:

$$
\left[ \hat{h} + \sum_{j}^{N} \left( \hat{J}_j - \hat{K}_j \right) \right] |\chi_i \rangle = \hat{F} |\chi_i \rangle  =\sum_j^N \epsilon_{ij} |\chi_i \rangle
$$

It can be shown that we can always find a unitary matrix that transforms the above equations into an eigenvalue problem by forming new orbitals from a linear combination of the original orbitals, but I'm not showing that here. The end result is the canonical Hartree-Fock equations:

$$
\hat{F} |\chi_i \rangle = \epsilon_i |\chi_i \rangle
$$


### Basis Set Expansion and Roothaan-Hall Equations
In practice, we express the unknown molecular orbitals as linear combinations of known atomic orbitals $\lbrace\phi_\nu\rbrace$:

$$
|\chi_i \rangle = \sum_{\nu} c_{\nu i} |\phi_\nu \rangle 
\rightarrow \hat{F} |\sum_{\nu} c_{\nu i} \phi_\nu \rangle = \epsilon_i |\sum_{\nu} c_{\nu i} \phi_\nu \rangle
$$

Multiplying both sides by $\langle \phi_\mu |$ and rearranging, we obtain the Roothaan-Hall equations:

$$
\sum_{\nu} F_{\mu \nu} c_{\nu i} = \epsilon_i \sum_{\nu} S_{\mu \nu} c_{\nu i}
\rightarrow \mathbf{F} \mathbf{c} = \mathbf{S} \mathbf{c} \mathbf{\epsilon}
$$

where $\mathbf{F}$ is the Fock matrix with elements $F_{\mu \nu} = \langle \phi_\mu | \hat{F} | \phi_\nu \rangle$, $\mathbf{S}$ is the overlap matrix with elements $S_{\mu \nu} = \langle \phi_\mu | \phi_\nu \rangle$, $\mathbf{c}$ is the coefficient matrix with elements $c_{\nu i}$, and $\mathbf{\epsilon}$ is the diagonal matrix of orbital energies $\epsilon_i$. Note that since we're working with atomic orbitals, which are generally not orthogonal (to be orthogonal, the electrons in two orbitals need to not interact with one another and physically speaking, you kind of need interaction to form bonds), $\mathbf{S} \neq \delta_{ij}$ .

To once again get an eigenvalue problem, we can apply a transformation using the symmetric orthogonalization matrix $\mathbf{S}^{-1/2}$ such that:

$$
\mathbf{F}' = \mathbf{S}^{-1/2} \mathbf{F} \mathbf{S}^{-1/2} \text{ ; } \mathbf{c}' = (\mathbf{S}^{-1/2})^{-1} \mathbf{c} \text{ ; } \mathbf{c} = \mathbf{S}^{-1/2} \mathbf{c}'
$$

where ($\mathbf{\Lambda}^{-1/2}$ is the diagonal matrix of the inverse square roots of the eigenvalues of $\mathbf{S}$ and $\mathbf{L}$ is the matrix of eigenvectors of $\mathbf{S}$):

$$
\mathbf{S}^{-1/2} = \mathbf{L} \mathbf{\Lambda}^{-1/2} \mathbf{L}^T
$$

Applying these transformations to the Roothaan-Hall equations, we get:

$$
\mathbf{F} \mathbf{c} = \mathbf{S} \mathbf{c} \mathbf{\epsilon} \\

\mathbf{F} (\mathbf{S}^{-1/2} \mathbf{c}') = \mathbf{S} (\mathbf{S}^{-1/2} \mathbf{c}') \mathbf{\epsilon} \\

(\mathbf{S}^{-1/2} \mathbf{F} \mathbf{S}^{-1/2}) \mathbf{c}' = (\mathbf{S}^{-1/2} \mathbf{S} \mathbf{S}^{-1/2}) \mathbf{c}' \mathbf{\epsilon} \\
\mathbf{F}' \mathbf{c}' = \mathbf{I} \mathbf{c}' \mathbf{\epsilon} \\
\mathbf{F}' \mathbf{c}' = \mathbf{c}' \mathbf{\epsilon}
$$

Given these ingredients, we can roughly see how the Hartree-Fock method can be implemented: guess $\mathbf{F}$, transform it to $\mathbf{F}'$, diagonalize to get $\mathbf{c}'$ and $\mathbf{\epsilon}$, back-transform to get $\mathbf{c}$, use $\mathbf{c}$ to compute a new $\mathbf{F}$ somehow, and repeat until convergence.


### The Density Matrix

Starting from the Hartree-Fock energy for a closed-shell system where each spatial orbital contains 2 electrons:

$$
E_{\text{elec}} = \sum_{i=1}^{N/2} \langle \chi_i | 2\hat{h} | \chi_i \rangle + \sum_{i=1}^{N/2} \sum_{j=1}^{N/2} \left[ 2\langle \chi_i \chi_i | \chi_j \chi_j \rangle - \langle \chi_i \chi_j | \chi_j \chi_i \rangle \right]
$$

**Expanding using LCAO:** Substitute the LCAO expansion $|\chi_i \rangle = \sum_{\nu} c_{\nu i} |\phi_\nu \rangle$ into the one-electron term:

$$
\langle \chi_i | \hat{h} | \chi_i \rangle = \left\langle \sum_{\mu} c_{\mu i} \phi_\mu \middle| \hat{h} \middle| \sum_{\nu} c_{\nu i} \phi_\nu \right\rangle = \sum_{\mu \nu} c_{\mu i}^* c_{\nu i} \langle \phi_\mu | \hat{h} | \phi_\nu \rangle = \sum_{\mu \nu} c_{\mu i} c_{\nu i} H^{\text{core}}_{\mu \nu}
$$

(assuming real coefficients). Summing over occupied orbitals:

$$
\sum_{i=1}^{N/2} \langle \chi_i | 2\hat{h} | \chi_i \rangle = 2\sum_{i=1}^{N/2} \sum_{\mu \nu} c_{\mu i} c_{\nu i} H^{\text{core}}_{\mu \nu} = 2\sum_{\mu \nu} H^{\text{core}}_{\mu \nu} \sum_{i=1}^{N/2} c_{\mu i} c_{\nu i} = 2\sum_{\mu \nu} H^{\text{core}}_{\mu \nu} D_{\mu \nu}
$$

where we identify the **density matrix** $D_{\mu \nu} = \sum_{i=1}^{N/2} c_{\mu i} c_{\nu i}$.

**For the two-electron terms:** Similarly, expand the two-electron integrals:

$$
\langle \chi_i \chi_i | \chi_j \chi_j \rangle = \sum_{\mu \nu \lambda \sigma} c_{\mu i} c_{\nu i} c_{\lambda j} c_{\sigma j} (\mu \nu | \lambda \sigma)
$$

Summing over occupied orbitals and using the density matrix:

$$
\sum_{i,j=1}^{N/2} \langle \chi_i \chi_i | \chi_j \chi_j \rangle = \sum_{\mu \nu \lambda \sigma} (\mu \nu | \lambda \sigma) \sum_{i=1}^{N/2} c_{\mu i} c_{\nu i} \sum_{j=1}^{N/2} c_{\lambda j} c_{\sigma j} = \sum_{\mu \nu \lambda \sigma} (\mu \nu | \lambda \sigma) D_{\mu \nu} D_{\lambda \sigma}
$$

The exchange term follows analogously. Combining all terms with proper factors:

$$
E_{\text{elec}} = \sum_{\mu \nu} D_{\mu \nu} H^{\text{core}}_{\mu \nu} + \frac{1}{2} \sum_{\mu \nu \lambda \sigma} D_{\mu \nu} D_{\lambda \sigma} \left[ 2(\mu \nu | \lambda \sigma) - (\mu \lambda | \nu \sigma) \right]
$$

Recognizing that the Fock matrix is $F_{\mu \nu} = H^{\text{core}}_{\mu \nu} + \sum_{\lambda \sigma} D_{\lambda \sigma} \left[ 2(\mu \nu | \lambda \sigma) - (\mu \lambda | \nu \sigma) \right]$, we can simplify:

$$
E_{\text{elec}} = \sum_{\mu \nu} D_{\mu \nu} \left[ H^{\text{core}}_{\mu \nu} + F_{\mu \nu} \right]
$$

This shows explicitly how the LCAO expansion coefficients $c_{\mu i}$ lead to the density matrix $D_{\mu \nu}$, which in turn enables a compact expression for the electronic energy in terms of basis function integrals


### The Self-Consistent Field (SCF) Procedure: Closed Shell Variant

The following procedure is specific to a closed-shell system, where each spatial orbital is doubly occupied (i.e., two electrons with opposite spins occupy the same spatial orbital):

**Step 1: Compute Basis Integrals**

Calculate all one-electron and two-electron integrals over the atomic orbital basis $\lbrace\phi_\mu\rbrace$:

$$
T_{\mu \nu} = \langle \phi_\mu | -\frac{1}{2}\nabla^2 | \phi_\nu \rangle \\
V_{\text{nuc}, \mu \nu} = \langle \phi_\mu | \sum_A \frac{Z_A}{r_{iA}} | \phi_\nu \rangle \\
S_{\mu \nu} = \langle \phi_\mu | \phi_\nu \rangle \\
(\mu \nu | \lambda \sigma) = \iint \phi_\mu^*(\vec{r}_1) \phi_\nu(\vec{r}_1) \frac{1}{r_{12}} \phi_\lambda^*(\vec{r}_2) \phi_\sigma(\vec{r}_2) d\vec{r}_1 d\vec{r}_2
$$

Form the core Hamiltonian matrix: $\mathbf{H}^{\text{core}} = \mathbf{T} + \mathbf{V}_{\text{nuc}}$.

**Step 2: Compute Orthogonalization Matrix**

Diagonalize the overlap matrix $\mathbf{S}$ to obtain eigenvalues $\mathbf{\Lambda}$ and eigenvectors $\mathbf{L}$. Then compute:

$$
\mathbf{S}^{-1/2} = \mathbf{L} \mathbf{\Lambda}^{-1/2} \mathbf{L}^T
$$

**Step 3: Initial Guess for Fock Matrix**

For the first iteration, use the core Hamiltonian as the initial guess:

$$
\mathbf{F}^{(0)} = \mathbf{H}^{\text{core}}
$$

**Step 4: Transform to Orthogonal Basis**

Apply the orthogonalization transformation:

$$
\mathbf{F}'^{(n)} = \mathbf{S}^{-1/2} \mathbf{F}^{(n)} \mathbf{S}^{-1/2}
$$

**Step 5: Solve Eigenvalue Problem**

Diagonalize $\mathbf{F}'^{(n)}$ to obtain the coefficient matrix $\mathbf{c}'^{(n)}$ and orbital energies $\mathbf{\epsilon}^{(n)}$:

$$
\mathbf{F}'^{(n)} \mathbf{c}'^{(n)} = \mathbf{c}'^{(n)} \mathbf{\epsilon}^{(n)}
$$

**Step 6: Back-Transform Coefficients**

Transform coefficients back to the original non-orthogonal basis:

$$
\mathbf{c}^{(n)} = \mathbf{S}^{-1/2} \mathbf{c}'^{(n)}
$$

**Step 7: Build Density Matrix**

Construct the density matrix from the coefficient matrix, summing over the occupied orbitals (N/2 spatial orbitals for a closed-shell system):

$$
D_{\mu \nu}^{(n)} = \sum_{i=1}^{N/2} c_{\mu i}^{(n)} c_{\nu i}^{(n)}
$$

**Step 8: Build New Fock Matrix**

Using the density matrix, construct the next Fock matrix:

$$
F_{\mu \nu}^{(n+1)} = H^{\text{core}}_{\mu \nu} + \sum_{\lambda \sigma} D_{\lambda \sigma}^{(n)} \left[ 2(\mu \nu | \lambda \sigma) - (\mu \lambda | \nu \sigma) \right]
$$

**Step 9: Check Convergence**

Calculate the electronic energy:

$$
E_{\text{elec}}^{(n)} = \sum_{\mu \nu} D_{\mu \nu}^{(n)} \left[ H^{\text{core}}_{\mu \nu} + F_{\mu \nu}^{(n)} \right]
$$

If $|E_{\text{elec}}^{(n)} - E_{\text{elec}}^{(n-1)}| < \text{threshold}$, convergence is achieved. Otherwise, return to Step 4 with the new Fock matrix $\mathbf{F}^{(n+1)}$.

**Step 10: Compute Total Energy**

Add the nuclear repulsion energy:

$$
E_{\text{total}} = E_{\text{elec}} + V_{NN} = E_{\text{elec}} + \sum_{A < B} \frac{Z_A Z_B}{R_{AB}}
$$

## Computational Implementation
In practice, we once again expand the atomic basis functions $\lbrace\phi_\mu\rbrace$ with another basis set for computational efficiency. 

### Contracted Gaussian Basis Sets: STO-nG
For the STO-nG basis set (very popular and simple), each Gaussian function has the following form:

$$
g(\alpha, \vec{r}, \vec{R}) = N (x-X)^{l} (y-Y)^{m} (z-Z)^{n} e^{-\alpha |\vec{r} - \vec{R}|^2}
$$

where `(x,y,z)` are the Cartesian coordinates of the electron, `(X,Y,Z)` are the coordinates of the nucleus, $\alpha$ is the exponent that determines the width of the Gaussian, and N is a normalization constant. The integers l, m, n determine the angular momentum of the orbital. For an s-orbital, `l = m = n = 0`. For a p-orbital, `(n,l,m) = (0,0,1), (0,1,0), or (1,0,0)`, correpsonding to px,py,pz. For a d-orbital, `(n,l,m) = (0,2,0), (2,0,0), (0,0,2), (1,1,0), or (1,0,1)`, corresponding to dxy, dxz, dyz, dx2-y2, dz2.

$\phi_\mu$ can then be expressed as:

$$
\phi_\mu(\vec{r}) = \sum_{k=1}^{n} c_k g(\alpha_k, \vec{r}, \vec{R})
$$

where $c_k$ are the contraction coefficients, $\alpha_k$ are the exponent coefficients, and n is the number of primitive Gaussians in the contraction (e.g., n=3 for STO-3G). All these coefficients are tabulated for each element and basis set (see [here](https://www.basissetexchange.org/)).

As for the normalization constant N, it can be computed as:

$$
N = \left( \frac{2\alpha}{\pi} \right)^{3/4} 
\left[ \frac{(8\alpha)^{n+l+m}n!l!m!}{(2n)!(2l)!(2m)!} \right]^{1/2}
$$



### Overlap, Kinetic, and Nuclear Attraction Integrals

[]

### Two-Electron Repulsion Integrals
[]

### Self-Consistent Field (SCF) Procedure
[]

## Resource(s)
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

