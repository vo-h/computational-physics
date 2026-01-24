use crate::hf::atom::Atom;
use crate::hf::basis_stog::stog_orbital::STOGOrbital;
use crate::hf::basis_stog::stog_integrator::compute_Sij;
use crate::hf::basis_stog::stog_integrator::compute_Tij;
use crate::hf::basis_stog::stog_integrator::compute_VijR;
use crate::hf::basis_stog::stog_integrator::compute_Vijkl;
use nalgebra::DMatrix;

#[derive(Clone, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub orbitals: Vec<STOGOrbital>,
}

struct EigenPairs {
    pub values: Vec<f64>,
    pub vectors: DMatrix<f64>,
}

fn sort_by_eigenvalues(eigvals: &Vec<f64>, eigvecs: &DMatrix<f64>) -> EigenPairs {
    // Sort eigenvectors by eigenvalues ascending
    let mut eig_pairs: Vec<(f64, Vec<f64>)> = eigvals.iter().cloned()
        .zip(eigvecs.column_iter().map(|col| col.iter().cloned().collect()))
        .collect();
    eig_pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());   
    
    return EigenPairs {
        values: eig_pairs.iter().map(|pair| pair.0).collect(),
        vectors: DMatrix::<f64>::from_columns(
            &eig_pairs.iter().map(|pair| DMatrix::<f64>::from_column_slice(eigvals.len(), 1, &pair.1).column(0).clone_owned()).collect::<Vec<_>>()
        ),
    };
}

impl Molecule {
    pub fn new(atoms: Vec<Atom>) -> Self {
        let mut orbitals = Vec::new();
        for atom in &atoms {
            for orb in &atom.orbitals {
                orbitals.push((*orb).clone());
            }
        }
        Molecule {
            atoms: atoms,
            orbitals: orbitals,
        }
    }

    pub fn from_string(lines: &str) -> Self {
        let mut atoms: Vec<Atom> = Vec::new();
        let items = lines.split('\n').collect::<Vec<&str>>();
        for line in items.iter() {
            let parts = line.trim().split_whitespace().collect::<Vec<&str>>();
            if parts.len() >= 5 {
                let symbol = parts[0];
                let x: f64 = parts[1].parse().unwrap();
                let y: f64 = parts[2].parse().unwrap();
                let z: f64 = parts[3].parse().unwrap();
                let Z: u8 = parts[4].parse().unwrap();
                // For simplicity, we use a fixed basis here; in practice, this could be parameterized
                let basis = parts[5].to_string();
                let atom = Atom::new(symbol, (x, y, z), Z, basis).unwrap();
                atoms.push(atom);
            }
        }
        Molecule::new(atoms)
    }

    pub fn compute_S(&self) -> DMatrix<f64> {
        let n = self.orbitals.len();
        let mut S = DMatrix::<f64>::zeros(n, n);

        for i in 0..n {
            for j in i..n {
                S[(i, j)] = compute_Sij(&self.orbitals[i], &self.orbitals[j]);
                S[(j, i)] = S[(i, j)];
            }
        }
        return S;
    }

    pub fn compute_T(&self) -> DMatrix<f64> {
        let n = self.orbitals.len();
        let mut T = DMatrix::<f64>::zeros(n, n);

        for i in 0..n {
            for j in i..n {
                T[(i, j)] = compute_Tij(&self.orbitals[i], &self.orbitals[j]);
                T[(j, i)] = T[(i, j)];
            }
        }
        return T;
    }

    pub fn compute_Vne(&self) -> DMatrix<f64> {
        let n = self.orbitals.len();
        let mut Vne = DMatrix::<f64>::zeros(n, n);
        for i in 0..n {
            for j in i..n {
                let mut val = 0.0;
                for atom in &self.atoms {
                    val += compute_VijR(&self.orbitals[i], &self.orbitals[j], atom.coords) * (-(atom.Z as f64));
                }
                Vne[(i, j)] = val;
                Vne[(j, i)] = Vne[(i, j)];
            }
        }
        return Vne;
    }

    pub fn compute_Vee(&self) -> Vec<Vec<Vec<Vec<f64>>>> {
        let n = self.orbitals.len();
        let mut Vee = vec![vec![vec![vec![0.0; n]; n]; n]; n];

        for i in 0..n {
            for j in 0..n {
                for k in 0..n {
                    for l in 0..n {
                        if (i*n + j) >= (k*n + l) {
                            Vee[i][j][k][l] = compute_Vijkl(&self.orbitals[i], &self.orbitals[j], &self.orbitals[k], &self.orbitals[l]);
                            Vee[j][i][k][l] = Vee[i][j][k][l];
                            Vee[i][j][l][k] = Vee[i][j][k][l];
                            Vee[j][i][l][k] = Vee[i][j][k][l];
                            Vee[k][l][i][j] = Vee[i][j][k][l];
                            Vee[k][l][j][i] = Vee[i][j][k][l];
                            Vee[l][k][i][j] = Vee[i][j][k][l];
                            Vee[l][k][j][i] = Vee[i][j][k][l];
                        }
                    }
                }
            }
        }
        return Vee;
    }

    pub fn compute_H(&self) -> DMatrix<f64> {
        let T = self.compute_T();
        let Vne = self.compute_Vne();
        let H = T + Vne;
        return H;
    }

    pub fn compute_Vnn(&self) -> f64 {
        let mut Vnn = 0.0;
        for i in 0..self.atoms.len() {
            for j in (i+1)..self.atoms.len() {
                let Zi = self.atoms[i].Z as f64;
                let Zj = self.atoms[j].Z as f64;
                let dx = self.atoms[i].coords.0 - self.atoms[j].coords.0;
                let dy = self.atoms[i].coords.1 - self.atoms[j].coords.1;
                let dz = self.atoms[i].coords.2 - self.atoms[j].coords.2;
                let r_ij = (dx*dx + dy*dy + dz*dz).sqrt();
                Vnn += Zi * Zj / r_ij;
            }
        }
        return Vnn;
    }

    pub fn compute_S12(&self) -> DMatrix<f64> {
        let S = self.compute_S();
        let eig = S.symmetric_eigen();
        let mut eigvals = eig.eigenvalues;
        let eigvecs = eig.eigenvectors;
        eigvals = eigvals.map(|x| 1.0 / x.sqrt());
        let eigvecst = eigvecs.transpose();
        let res = eigvecs * eigvals * eigvecst;
        return res;
    }

    pub fn compute_F0(&self) -> DMatrix<f64> {
        let S12 = self.compute_S12();
        let S12T = S12.transpose();
        let H = self.compute_H();
        return S12T * H * S12;
    }

    pub fn compute_C0(&self) -> DMatrix<f64> {
        let F0 = self.compute_F0();
        let S12 = self.compute_S12();
        let eig = F0.symmetric_eigen();
        let eigvals = eig.eigenvalues;
        let eigvecs = eig.eigenvectors;
        
        let sorted_eig =  sort_by_eigenvalues(&eigvals.data.as_vec(), &eigvecs);
        return S12 * sorted_eig.vectors;
    }

    pub fn compute_D0(&self) -> DMatrix<f64> {
        let occ =  self.orbitals.len() / 2;
        let mut C0 = self.compute_C0();
        C0 = C0.view_mut((0, 0), (occ, occ)).into_owned();
        let C0T = C0.transpose();
        let D0 = C0 * C0T;
        return D0;
    }

    pub fn compute_E0(&self) -> f64 {
        let D0 = self.compute_D0();
        let H = self.compute_H();
        return 2.0 * (D0.component_mul(&H)).sum();
    }

    pub fn compute_CHF(&self) -> f64 {
        // Placeholder implementation
        0.0
    }
}
