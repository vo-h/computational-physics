use std::f64::consts::PI;

pub mod stog_primitive {
    pub struct StoGPrimitive {
        pub coords: (f64, f64, f64),
        pub alpha: f64,
        pub cc: f64,
        pub nx: u32,
        pub ny: u32,
        pub nz: u32,
        pub l: u32,
        pub N: f64,
    }

    impl StoGPrimitive {
        pub fn new(
            coords: (f64, f64, f64),
            alpha: f64,
            cc: f64,
            nx: u32,
            ny: u32,
            nz: u32,
        ) -> Self {

            fn factorial(num: u32) -> f64 {
                (1..=num).map(|x| x as f64).product()
            }
            let prefactor = (2.0 * alpha / PI).powf(0.75);
            let numerator = (8.0 * alpha).powi((nx + ny + nz) as i32) * factorial(nx) * factorial(ny) * factorial(nz);
            let denominator = factorial(2 * nx) * factorial(2 * ny) * factorial(2 * nz);
            let N = prefactor * (numerator / denominator).sqrt();
            StoGPrimitive {
                coords,
                alpha,
                cc,
                nx,
                ny,
                nz,
                nx + ny + nz,
                N,
            }
        }
    }
}

pub mod stog_orbital {
    use super::stog_primitive::StoGPrimitive;

    pub struct StoGOrbital {
        pub atom: String,
        pub shell: String,
        pub coords: (f64, f64, f64),
        pub gtos: Vec<StoGPrimitive>,
    }

    impl StoGOrbital {
        pub fn new(atom: String, shell: String, coords: (f64, f64, f64), cc: Vec<f64>, alpha: Vec<f64>, nx: Vec<u32>, ny: Vec<u32>, nz: Vec<u32>) -> Self {
            let mut gtos = Vec::new();
            for i in 0..cc.len() {
                gtos.push(StoGPrimitive::new(coords, alpha[i], cc[i], nx[i], ny[i], nz[i]));
            }
            StoGOrbital { atom, shell, coords, gtos }
        }
    }
}

pub mod stog_integrator {
    use super::stog_orbital::StoGOrbital;
    
    pub fn compute_Sij(orb1: &StoGOrbital, orb2: &StoGOrbital) -> f64 {
        let mut sij = 0.0;
        for u in 0..orb1.gtos.len() {
            for v in 0..orb2.gtos.len() {
                let prim1 = &orb1.gtos[u];
                let prim2 = &orb2.gtos[v];
                let sx = compute_sx(prim1.coords.0, prim2.coords.0, prim1.alpha, prim2.alpha, prim1.nx as i32, prim2.nx as i32);
                let sy = compute_sx(prim1.coords.1, prim2.coords.1, prim1.alpha, prim2.alpha, prim1.ny as i32, prim2.ny as i32);
                let sz = compute_sx(prim1.coords.2, prim2.coords.2, prim1.alpha, prim2.alpha, prim1.nz as i32, prim2.nz as i32);
                let coeff = orb1.gtos[u].cc * orb2.gtos[v].cc * orb1.gtos[u].N * orb2.gtos[v].N;
                let prefactor = prim1.N * prim2.N * prim1.cc * prim2.cc;  
                let E_AB = (-prim1.alpha * prim2.alpha / (prim1.alpha + prim2.alpha) * (
                    (prim1.coords.0 - prim2.coords.0).powi(2) +
                    (prim1.coords.1 - prim2.coords.1).powi(2) +
                    (prim1.coords.2 - prim2.coords.2).powi(2)
                )).exp();
                sij += coeff * E_AB * prefactor * sx * sy * sz;
            }
        }


    }

    pub fn compute_Tij(orb1: &StoGOrbital, orb2: &StoGOrbital) -> f64 {
        // Placeholder for kinetic energy integral computation between STO-G orbitals
    }

    pub fn compute_VijR(orb1: &StoGOrbital, orb2: &StoGOrbital, Z: f64) -> f64 {
        // Placeholder for nuclear attraction integral computation between STO-G orbitals
    }

    pub fn compute_Vijkl(orb1: &StoGOrbital, orb2: &StoGOrbital, orb3: &StoGOrbital, orb4: &StoGOrbital) -> f64 {
        // Placeholder for electron repulsion integral computation between STO-G orbitals
    }

    fn compute_sx(&self, A: f64, B: f64, alpha: f64, beta: f64, ai: i32, bi: i32) -> f64 {
        let P = (alpha * A + beta * B)/(alpha + beta);
        if ai < 0 || bi < 0 {return 0.0;}
        if (ai, bi) == (0, 0) {return 1.0;}
        if (ai, bi) == (1, 0) {return -(P-A);}
        if ai > 1 && bi == 0 {
            let term1 = -(A-P)*self.compute_sx(A, B, alpha, beta, ai-1, 0);
            let term2 = (ai-1) as f64/(2.0*(alpha+beta))*self.compute_sx(A, B, alpha, beta, ai-2, 0);
            return term1 + term2;
        }
        let term1 = self.compute_sx(A, B, alpha, beta, ai+1, bi-1);
        let term2 = self.compute_sx(A, B, alpha, beta, ai, bi-1);
        return term1 + (A-B) * term2;
    }

    fn compute_tx(&self, A: f64, B: f64, alpha: f64, beta: f64, ai: i32, bi: i32) -> f64 {
        let P = (alpha * A + beta * B)/(alpha + beta);
        if ai < 0 || bi < 0 {return 0.0;}
        if (ai, bi) == (0, 0) {return 2.0 *alpha * beta / self.compute_sx(A, B, alpha, beta, 1, 1);}
        if bi == 0 {
            let term1 = -ai * beta * self.compute_sx(A, B, alpha, beta, ai-1, 1);
            let term2 = 2.0 * alpha * beta * self.compute_sx(A, B, alpha, beta, ai+1, 1);
            return term1 + term2;
        }
        if ai == 0 {
            let term1 = -bi * alpha * self.compute_sx(A, B, alpha, beta, 1, bi-1);
            let term2 = 2.0 * alpha * beta * self.compute_sx(A, B, alpha, beta, 1, bi+1);
            return term1 + term2;
        }
        let term1 = ai * bi * self.compute_sx(A, B, alpha, beta, ai-1, bi-1);
        let term2 = 2 * ai * beta * self.compute_sx(A, B, alpha, beta, ai-1, bi+1);
        let term3 = 2 * alpha * bi * self.compute_sx(A, B, alpha, beta, ai+1, bi-1);
        let term4 = 4 * alpha * beta * self.compute_sx(A, B, alpha, beta, ai+1, bi+1);
        return 0.5 * (term1 - term2 - term3 + term4);
    }
    

}