
# Core Plan
1. &#x2714; Write Python code to compute the overlap matrix
2. Write C code for numerical integration & differentiation to speed up overall calculations
3. Write C & Python code to compute the kinetic matrix
4. Write C & Python code to compute the nuclear attraction matrix
5. Write C & Python code to compute the two-electron repulsion integrals
6. Write C & Python code to perform the SCF procedure and compute the total energy of the system
7. Write C & Python code to compute the forces on the nuclei and perform geometry optimization

# Stretch Goals
1. Write C & Python code to compute the vibrational frequencies of the molecule
2. Write C & Python code to compute the electronic excited states of the molecule using time-dependent Hartree-Fock (TDHF) or configuration interaction singles (CIS) methods
3. Write C & Python code to compute the properties of the molecule, such as dipole moment, polarizability, etc.
4. Write C & Python code to perform molecular dynamics simulations using the computed forces on the nuclei.

# Resource(s)

* [Computed Basis Sets For Most Elements](https://www.basissetexchange.org/)
* Example Calculations of H2O: [Lecture Notes With Results from Gaussian for Each Matrix](https://chemistry.montana.edu/callis/courses/chmy564/460water.pdf) & [Hand Calculations For the Overlap Matrix Using GTO-3G](https://content.wolfram.com/sites/19/2012/02/Ho.pdf)
* [PySCF](https://pyscf.org/) & [Associated Tutorials](https://pyscf.org/pyscf_tutorials/index.html) for comparing results against.
* [Integral Formulas](https://booksite.elsevier.com/9780444594365/downloads/16755_10036.pdf) for computing the integrals analytically.
* Series of Integral Evaluations:
    * [Overlap Integral](https://content.wolfram.com/sites/19/2012/02/Ho.pdf)
    * [Kinetic Energy Integral](https://content.wolfram.com/sites/19/2013/01/Ho_Kinetic.pdf)
    * [Nuclear Attraction Integral](https://content.wolfram.com/sites/19/2014/12/Ho_Nuclear.pdf)
    * [Two-Electron Repulsion Integral](https://www.youtube.com/watch?v=9n2s8Xo7l3k&t=1h30m)

## Additional Resources
* [PyQuint](https://ifilot.github.io/pyqint/electronic_structure_calculations.html) & [Associated Exercises](https://github.com/ifilot/hfhsl2021/)
* Some More Notebooks & Scripts: [link1](https://joshuagoings.com/2017/04/28/integrals/) & [link2](https://cwagen.substack.com/p/teach-yourself-quantum-chemistry)

# Integral Evaluation

## Overlap Integral
The overlap integral between two atomic orbitals $\phi_i$ and $\phi_j$ is given by
$$S_{ij} = \int \phi_i(\mathbf{r}) \phi_j(\mathbf{r}) d\mathbf{r}$$
where 

$$\phi_i = \sum_\mu N_\mu c_\mu g_\mu$$

* $N_\mu$ is the normalization constant for the primitive Gaussian function $g_\mu$.
* $c_\mu$ is the contraction coefficient for the primitive Gaussian function $g_\mu$.
* $g_\mu$ is the primitive Gaussian function defined as (where $(X, Y, Z)$ is the center of the orbital):
$$
g_\mu(\mathbf{r}) = (x - X)^{n_x} (y - Y)^{n_y} (z - Z)^{n_z} e^{-\alpha [(x - X)^2 + (y - Y)^2 + (z - Z)^2]}
$$

$N_\mu$ can be computed using the formula
$$N_\mu = \left( \frac{2\alpha}{\pi} \right)^{3/4} \sqrt{\frac{(8\alpha)^{n_x + n_y + n_z} n_x! n_y! n_z!}{(2n_x)! (2n_y)! (2n_z)!}}$$
where
* For the s orbital, $n_x = n_y = n_z = 0$.
* For the p orbital, ($n_x$, $n_y$, $n_z$) is one of (1, 0, 0), (0, 1, 0), or (0, 0, 1) depending on whether it's a $p_x$, $p_y$, or $p_z$ orbital.


### $S_{ij}$ Computation Deep Dive
For each pair of orbitals, $S_{ij}$ involves summing over all pairs of primitive Gaussian functions from the two orbitals.

$$

S_{ij} = \sum_{\mu=1}^{N_i} \sum_{\nu=1}^{N_j} N_\mu c_\mu N_\nu c_\nu \int g_\mu(\mathbf{r_i}) g_\nu(\mathbf{r_j}) d\mathbf{r}$$

### S-S Overlap
For s-s overlap, the integral simplifies to
$$

\int_{s-s}^{overlap} g_\mu(\mathbf{r_i}) g_\nu(\mathbf{r_j}) d\mathbf{r}
\\\\

=
\int_{-\infty}^\infty \int_{-\infty}^\infty \int_{-\infty}^\infty
e^{-\alpha_\mu [(x - X_i)^2 + (y - Y_i)^2 + (z - Z_i)^2]} *
e^{-\alpha_\nu [(x - X_j)^2 + (y - Y_j)^2 + (z - Z_j)^2]} \\\\

=

\int_{-\infty}^\infty e^{-\alpha_\mu (x - X_i)^2 - \alpha_\nu (x - X_j)^2} dx
\int_{-\infty}^\infty e^{-\alpha_\mu (y - Y_i)^2 - \alpha_\nu (y - Y_j)^2} dy
\int_{-\infty}^\infty e^{-\alpha_\mu (z - Z_i)^2 - \alpha_\nu (z - Z_j)^2} dz

\\\\

= 
\int_{-\infty}^\infty e^{-(\alpha_\mu + \alpha_\nu) x^2 + (2\alpha_\mu X_i + 2\alpha_\nu X_j) x - (\alpha_\mu X_i^2 + \alpha_\nu X_j^2)} dx
\int_{-\infty}^\infty e^{-(\alpha_\mu + \alpha_\nu) y^2 + (2\alpha_\mu Y_i + 2\alpha_\nu Y_j) y - (\alpha_\mu Y_i^2 + \alpha_\nu Y_j^2)} dy
\int_{-\infty}^\infty e^{-(\alpha_\mu + \alpha_\nu) z^2 + (2\alpha_\mu Z_i + 2\alpha_\nu Z_j) z - (\alpha_\mu Z_i^2 + \alpha_\nu Z_j^2)} dz \\\\

=
\boxed{
\left( \sqrt{\frac{\pi}{\alpha_\mu + \alpha_\nu}} * e^{\frac{(\alpha_\mu X_i + \alpha_\nu X_j)^2}{4*(\alpha_\mu + \alpha_\nu)}+(\alpha_\mu*X_i^2- \alpha_\nu*X_j^2)} \right)
\left( \sqrt{\frac{\pi}{\alpha_\mu + \alpha_\nu}} * e^{\frac{(\alpha_\mu Y_i + \alpha_\nu Y_j)^2}{4*(\alpha_\mu + \alpha_\nu)}+(\alpha_\mu*Y_i^2- \alpha_\nu*Y_j^2)} \right)
\left( \sqrt{\frac{\pi}{\alpha_\mu + \alpha_\nu}} * e^{\frac{(\alpha_\mu Z_i + \alpha_\nu Z_j)^2}{4*(\alpha_\mu + \alpha_\nu)}+(\alpha_\mu*Z_i^2- \alpha_\nu*Z_j^2)} \right)
}
$$

### S-P Overlap
For s-p overlap, the integral is more complex due to the presence of the Cartesian component (x, y, or z) in the p orbital. For example, for a $p_x$ orbital centered at $X_j$, the integral separates into three 1D integrals:
$$

\int_{s-p}^{overlap} g_\mu(\mathbf{r_i}) g_\nu(\mathbf{r_j}) d\mathbf{r}
= 
\int_{-\infty}^\infty (x - X_j) e^{-\alpha_\mu (x - X_i)^2 - \alpha_\nu (x - X_j)^2} dx \cdot I_y^{ss} \cdot I_z^{ss}
$$
where 


$$
\boxed{
I_y^{ss} = \int_{-\infty}^\infty e^{-\alpha_\mu (y - Y_i)^2 - \alpha_\nu (y - Y_j)^2} dy = \sqrt{\frac{\pi}{\alpha_\mu + \alpha_\nu}} * e^{\frac{(\alpha_\mu Y_i + \alpha_\nu Y_j)^2}{4*(\alpha_\mu + \alpha_\nu)}+(\alpha_\mu*Y_i^2- \alpha_\nu*Y_j^2)}
}
$$

and

$$
\boxed{
I_z^{ss} = \int_{-\infty}^\infty e^{-\alpha_\mu (z - Z_i)^2 - \alpha_\nu (z - Z_j)^2} dz = \sqrt{\frac{\pi}{\alpha_\mu + \alpha_\nu}} * e^{\frac{(\alpha_\mu Z_i + \alpha_\nu Z_j)^2}{4*(\alpha_\mu + \alpha_\nu)}+(\alpha_\mu*Z_i^2- \alpha_\nu*Z_j^2)}
}
$$

Using the weighted center $P_x = \frac{\alpha_\mu X_i + \alpha_\nu X_j}{\alpha_\mu + \alpha_\nu}$ and the substitution $(x - X_j) = (x - P_x) + (P_x - X_j)$, where the $(x - P_x)$ term vanishes due to symmetry, we get:

$$
\boxed{
\int_{s-p}^{overlap} = \frac{\alpha_\mu (X_i - X_j)}{\alpha_\mu + \alpha_\nu} \cdot \int_{s-s}^{overlap}
}
$$

### P-P Overlap
For p-p overlap (e.g., $p_x$ on center i with $p_x$ on center j), we have $(x - X_i)(x - X_j)$ in the integrand. Expanding around the weighted center and evaluating the Gaussian integrals yields:

$$
\boxed{
\int_{p_x-p_x}^{overlap} = \left[\frac{1}{2(\alpha_\mu + \alpha_\nu)} - \frac{\alpha_\mu \alpha_\nu (X_i - X_j)^2}{(\alpha_\mu + \alpha_\nu)^2}\right] \cdot \int_{s-s}^{overlap}
}
$$

For orthogonal p orbitals (e.g., $p_x$ with $p_y$):

$$
\boxed{
\int_{p_x-p_y}^{overlap} = \frac{\alpha_\mu(X_i - X_j)}{\alpha_\mu + \alpha_\nu} \cdot \frac{\alpha_\mu(Y_i - Y_j)}{\alpha_\mu + \alpha_\nu} \cdot \int_{s-s}^{overlap}
}
$$

## $T_{ij}$ Computation Deep Dive
The kinetic energy integral between two atomic orbitals $\phi_i$ and $\phi_j$ is given by
$$T
_{ij} = -\frac{1}{2} \int \phi_i(\mathbf{r}) \nabla^2 \phi_j(\mathbf{r})
$$

With the same expansion of $\phi_i$ and $\phi_j$ in terms of primitive Gaussian functions, we have
$$
T_{ij} = -\frac{1}{2} \sum_{\mu=1}^{N_i} \sum_{\nu=1}^{N_j} N_\mu c_\mu N_\nu c_\nu \int g_\mu(\mathbf{r_i}) \nabla^2 g_\nu(\mathbf{r_j}) d\mathbf{r}
$$

### S-S Overlap

$$
\nabla^2 g_\nu(\mathbf{r_j}) = \frac{\partial^2 g_\nu}{\partial x^2} + \frac{\partial^2 g_\nu}{\partial y^2} + \frac{\partial^2 g_\nu}{\partial z^2} \\\\
$$

Looking only at the double derivative in x and focusing on the relevant terms, we get

$$
\frac{\partial^2 g_\nu}{\partial x^2} =
e^{-\alpha [(y - Y)^2 + (z - Z)^2]}
\frac{\partial}{\partial x^2} \left[ e^{-\alpha (x - X)^2}
\right]
\\\\
= e^{-\alpha [(y - Y)^2 + (z - Z)^2]} * \left[ -2\alpha e^{-\alpha (x - X)^2} + 4\alpha^2 (x - X)^2 e^{-\alpha (x - X)^2} \right] \\\\
= e^{-\alpha [(y - Y)^2 + (z - Z)^2 + (x - X)^2]} * \left[ -2\alpha + 4\alpha^2 (x - X)^2 \right]
$$

Now, applying the integral, we get
$$
\int_{s-s}^{kinetic,x} = 
I_y^{ss} \cdot I_z^{ss} \cdot \int_{-\infty}^\infty e^{-\alpha_\mu (x - X_i)^2 - \alpha_\nu (x - X_j)^2} * \left[ -2\alpha_\nu + 4\alpha_\nu^2 (x - X_j)^2 \right] dx

$$

The integral in x can be computed then refactored in terms of $I_x^{ss}$ to give:

$$
I_x^{ss} * \frac{2\alpha_\nu \alpha_\mu}{\alpha_\nu + \alpha_\mu} * \left[ \frac{2\alpha_\nu \alpha_\mu (X_i - X_j)^2}{(\alpha_\nu + \alpha_\mu)^2} - 1\right]
$$

Combining all these expressions gives

$$
\int_{s-s}^{kinetic,x} = \frac{2\alpha_\nu \alpha_\mu}{\alpha_\nu + \alpha_\mu} * \left[ \frac{2\alpha_\nu \alpha_\mu (X_i - X_j)^2}{(\alpha_\nu + \alpha_\mu)^2} - 1\right] * \int_{s-s}^{overlap}
$$

The same expression applies for the y and z components, with the appropriate substitutions for the coordinates. Summing over all three components gives the final expression for the kinetic energy integral between two s orbitals:

$$\boxed{
\int_{s-s}^{kinetic} =  \left[ \left[ \frac{2\alpha_\nu \alpha_\mu (X_i - X_j)^2}{(\alpha_\nu + \alpha_\mu)^2} - 1\right] +  \left[ \frac{2\alpha_\nu \alpha_\mu (Y_i - Y_j)^2}{(\alpha_\nu + \alpha_\mu)^2} - 1\right] + \left[ \frac{2\alpha_\nu \alpha_\mu (Z_i - Z_j)^2}{(\alpha_\nu + \alpha_\mu)^2} - 1\right] \right] * \frac{2\alpha_\nu \alpha_\mu}{\alpha_\nu + \alpha_\mu} * \int_{s-s}^{overlap}
}
$$



$$
2*\alpha*e^{-\alpha (x - X)^2} 
*
(-1 + 2x\alpha (x-X) + 2\alpha X * (x - X))

$$