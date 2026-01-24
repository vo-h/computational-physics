use crate::hf::atom::Atom;
use crate::hf::basis_stog::stog_orbital::STOGOrbital;
use crate::hf::basis_stog::stog_integrator::compute_Sij;
use crate::hf::basis_stog::stog_integrator::compute_Tij;
use crate::hf::basis_stog::stog_integrator::compute_VijR;
use crate::hf::basis_stog::stog_integrator::compute_Vijkl;

#[derive(Clone, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub orbitals: Vec<STOGOrbital>,
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
                let basis = "3-21g".to_string();
                let atom = Atom::new(symbol, (x, y, z), Z, basis).unwrap();
                atoms.push(atom);
            }
        }
        Molecule::new(atoms)
    }

    pub fn compute_S(&self) -> Vec<Vec<f64>> {
        let n = self.orbitals.len();
        let mut S = vec![vec![0.0; n]; n];

        for i in 0..n {
            for j in i..n {
                S[i][j] = compute_Sij(&self.orbitals[i], &self.orbitals[j]);
                S[j][i] = S[i][j];
            }
        }
        return S;
    }

    pub fn compute_T(&self) -> Vec<Vec<f64>> {
        let n = self.orbitals.len();
        let mut T = vec![vec![0.0; n]; n];

        for i in 0..n {
            for j in i..n {
                T[i][j] = compute_Tij(&self.orbitals[i], &self.orbitals[j]);
                T[j][i] = T[i][j];
            }
        }
        return T;
    }

    pub fn compute_Vne(&self) -> Vec<Vec<f64>> {
        let n = self.orbitals.len();
        let mut Vne = vec![vec![0.0; n]; n];
        for i in 0..n {
            for j in i..n {
                let mut val = 0.0;
                for atom in &self.atoms {
                    val += compute_VijR(&self.orbitals[i], &self.orbitals[j], atom.coords) * (atom.Z as f64);
                }
                Vne[i][j] = val;
                Vne[j][i] = Vne[i][j];
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

    pub fn compute_S12(&self) -> Vec<Vec<f64>> {
        vec![vec![0.0]]
    }

    pub fn compute_F0(&self) -> Vec<Vec<f64>> {
        // Placeholder implementation
        vec![vec![0.0]]
    }

    pub fn compute_C0(&self) -> Vec<Vec<f64>> {
        // Placeholder implementation
        vec![vec![0.0]]
    }


    pub fn compute_D0(&self) -> Vec<Vec<f64>> {
        // Placeholder implementation
        vec![vec![0.0]]
    }

    pub fn compute_E0(&self) -> f64 {
        // Placeholder implementation
        0.0
    }

    pub fn compute_CHF(&self) -> f64 {
        // Placeholder implementation
        0.0
    }
}
