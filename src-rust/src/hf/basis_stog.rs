use std::f64::consts::PI;
use xsf::hyp1f1;

pub mod stog_primitive {
    pub struct STOGPrimitive {
        pub coords: (f64, f64, f64),
        pub alpha: f64,
        pub cc: f64,
        pub nx: i32,
        pub ny: i32,
        pub nz: i32,
        pub l: i32,
        pub N: f64,
    }

    impl STOGPrimitive {
        pub fn new(
            coords: (f64, f64, f64),
            alpha: f64,
            cc: f64,
            nx: i32,
            ny: i32,
            nz: i32,
        ) -> Self {

            fn factorial(num: u32) -> f64 {
                (1..=num).map(|x| x as f64).product()
            }
            let prefactor = (2.0 * alpha / PI).powf(0.75);
            let numerator = (8.0 * alpha).powi((nx + ny + nz) as i32) * factorial(nx as u32) * factorial(ny as u32) * factorial(nz as u32);
            let denominator = factorial(2 * nx as u32) * factorial(2 * ny as u32) * factorial(2 * nz as u32);
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
    use super::stog_primitive::STOGPrimitive;

    pub struct STOGOrbital {
        pub atom: String,
        pub shell: String,
        pub coords: (f64, f64, f64),
        pub gtos: Vec<STOGPrimitive>,
    }

    impl STOGOrbital {
        pub fn new(atom: String, shell: String, coords: (f64, f64, f64), cc: Vec<f64>, alpha: Vec<f64>, nx: Vec<u32>, ny: Vec<u32>, nz: Vec<u32>) -> Self {
            let mut gtos = Vec::new();
            for i in 0..cc.len() {
                gtos.push(STOGPrimitive::new(coords, alpha[i], cc[i], nx[i], ny[i], nz[i]));
            }
            STOGOrbital { atom, shell, coords, gtos }
        }
    }
}

pub mod stog_integrator {
    use super::stog_primitive::STOGPrimitive;
    use super::stog_orbital::STGOrbital;
    
    pub fn compute_Sij(orb1: &STOGOrbital, orb2: &STOGOrbital) -> f64 {

        fn compute_Sab(gto1: &STOGPrimitive, gto2: &STOGPrimitive) -> f64 {
            // Overlap integral between two STO-G primitives
            let prefactor = (PI / (gto1.alpha + gto2.alpha)).powf(1.5);
            let s_x = compute_E(gto1.coords.0, gto2.coords.0, gto1.alpha, gto2.alpha, gto1.nx as i32, gto2.nx as i32, 0);
            let s_y = compute_E(gto1.coords.1, gto2.coords.1, gto1.alpha, gto2.alpha, gto1.ny as i32, gto2.ny as i32, 0);
            let s_z = compute_E(gto1.coords.2, gto2.coords.2, gto1.alpha, gto2.alpha, gto1.nz as i32, gto2.nz as i32, 0);
            return prefactor * s_x * s_y * s_z
        }

        let mut sij = 0.0;
        for u in 0..orb1.gtos.len() {
            for v in 0..orb2.gtos.len() {
                let gto1 = &orb1.gtos[u];
                let gto2 = &orb2.gtos[v];
                sij += gto1.N * gto2.N * gto1.cc * gto2.cc * compute_Sab(gto1, gto2);
            }
        }
        return sij;
    }

    pub fn compute_Tij(orb1: &STOGOrbital, orb2: &STOGOrbital) -> f64 {
        let tij = 0.0;

        fn compute_D(A: f64, B: f64, alpha: f64, beta: f64, ai: i32, bi: i32, t: i32) -> f64 {
            let p = alpha + beta;
            if t == 0 {return compute_E(A, B, alpha, beta, ai, bi, 0) * (PI/p).powf(0.5);}
            let term1 = bi as f64 * compute_D(A, B, alpha, beta, ai, bi-1, t-1);
            let term2 = 2.0 * beta * compute_D(A, B, alpha, beta, ai, bi+1, t-1);
            return term1 - term2;
        }

        fn compue_Tab(orb1: &STOGPrimitive, orb2: &STOGPrimitive) -> f64 {
            let term1 = compute_D(orb1.coords.0, orb2.coords.0, orb1.alpha, orb2.alpha, orb1.nx, orb2.nx, 2);
            let term2 = compute_D(orb1.coords.1, orb2.coords.1, orb1.alpha, orb2.alpha, orb1.ny, orb2.ny, 2);
            let term3 = compute_D(orb1.coords.2, orb2.coords.2, orb1.alpha, orb2.alpha, orb1.nz, orb2.nz, 2);
            let term4 = compute_D(orb1.coords.0, orb2.coords.0, orb1.alpha, orb2.alpha, orb1.nx, orb2.nx, 0);
            let term5 = compute_D(orb1.coords.1, orb2.coords.1, orb1.alpha, orb2.alpha, orb1.ny, orb2.ny, 0);
            let term6 = compute_D(orb1.coords.2, orb2.coords.2, orb1.alpha, orb2.alpha, orb1.nz, orb2.nz, 0);
            return -0.5 * (term1*term5*term6 + term4*term2*term6 + term4*term5*term3);
        }
        
        for u in 0..orb1.gtos.len() {
            for v in 0..orb2.gtos.len() {
                let gto1 = &orb1.gtos[u];
                let gto2 = &orb2.gtos[v];
                tij += gto1.N * gto2.N * gto1.cc * gto2.cc * compue_Tab(gto1, gto2);
            }
        }
        return tij;
    }

    pub fn compute_VijR(orb1: &STOGOrbital, orb2: &STOGOrbital, R: (f64, f64, f64)) -> f64 {
        // Placeholder for nuclear attraction integral computation between STO-G orbitals
        fn compute_Vab(gto1: &STOGPrimitive, gto2: &STOGPrimitive, R: (f64, f64, f64)) -> f64 {
            let p = gto1.alpha + gto2.alpha;
            let P = (
                (gto1.alpha * gto1.coords.0 + gto2.alpha * gto2.coords.0) / p,
                (gto1.alpha * gto1.coords.1 + gto2.alpha * gto2.coords.1) / p,
                (gto1.alpha * gto1.coords.2 + gto2.alpha * gto2.coords.2) / p,
            );

            let mut val = 0.0;
            for t in 0..=(gto1.nx + gto2.nx) {
                for u in 0..=(gto1.ny + gto2.ny) {
                    for v in 0..=(gto1.nz + gto2.nz) {
                        let E_x = compute_E(gto1.coords.0, gto2.coords.0, gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, t);
                        let E_y = compute_E(gto1.coords.1, gto2.coords.1, gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, u);
                        let E_z = compute_E(gto1.coords.2, gto2.coords.2, gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, v);
                        let R_tuv = compute_R(t as i32, u as i32, v as i32, 0, p, P, R);
                        val += E_x * E_y * E_z * R_tuv;
                    }
                }
            }
            return (2.0 * PI / p) * val;
        }

        let mut vijr = 0.0;
        for u in 0..orb1.gtos.len() {
            for v in 0..orb2.gtos.len() {
                let gto1 = &orb1.gtos[u];
                let gto2 = &orb2.gtos[v];
                vijr += gto1.N * gto2.N * gto1.cc * gto2.cc * compute_Vab(gto1, gto2, R);
            }
        }
        return vijr;
    }

    pub fn compute_Vijkl(orb1: &STOGOrbital, orb2: &STOGOrbital, orb3: &STOGOrbital, orb4: &STOGOrbital) -> f64 {
        fn compute_Vabcd(gto1: &STOGPrimitive, gto2: &STOGPrimitive, gto3: &STOGPrimitive, gto4: &STOGPrimitive) -> f64 {
            p = gto1.alpha + gto2.alpha;
            q = gto3.alpha + gto4.alpha;
            P = (
                (gto1.alpha * gto1.coords.0 + gto2.alpha * gto2.coords.0) / p,
                (gto1.alpha * gto1.coords.1 + gto2.alpha * gto2.coords.1) / p,
                (gto1.alpha * gto1.coords.2 + gto2.alpha * gto2.coords.2) / p,
            );
            Q = (
                (gto3.alpha * gto3.coords.0 + gto4.alpha * gto4.coords.0) / q,
                (gto3.alpha * gto3.coords.1 + gto4.alpha * gto4.coords.1) / q,
                (gto3.alpha * gto3.coords.2 + gto4.alpha * gto4.coords.2) / q,
            );

            let mut val = 0.0;
            for t in 0..=(gto1.nx + gto2.nx) {
                for u in 0..=(gto1.ny + gto2.ny) {
                    for v in 0..=(gto1.nz + gto2.nz) {
                        for tau in 0..=(gto3.nx + gto4.nx) {
                            for nu in 0..=(gto3.ny + gto4.ny) {
                                for phi in 0..=(gto3.nz + gto4.nz) {
                                    let term1 = compute_E(gto1.coords.0, gto2.coords.0, gto1.alpha, gto2.alpha, gto1.nx, gto2.nx, t);
                                    let term2 = compute_E(gto1.coords.1, gto2.coords.1, gto1.alpha, gto2.alpha, gto1.ny, gto2.ny, u);
                                    let term3 = compute_E(gto1.coords.2, gto2.coords.2, gto1.alpha, gto2.alpha, gto1.nz, gto2.nz, v);
                                    let term4 = compute_E(gto3.coords.0, gto4.coords.0, gto3.alpha, gto4.alpha, gto3.nx, gto4.nx, tau);
                                    let term5 = compute_E(gto3.coords.1, gto4.coords.1, gto3.alpha, gto4.alpha, gto3.ny, gto4.ny, nu);
                                    let term6 = compute_E(gto3.coords.2, gto4.coords.2, gto3.alpha, gto4.alpha, gto3.nz, gto4.nz, phi);
                                    let term7 = (-1.0 as f64).powi((tau + nu + phi) as i32);
                                    let term8 = compute_R((t+tau) as i32, (u + nu) as i32, (v + phi) as i32, 0, p * q / (p + q), P, Q);
                                    val += term1 * term2 * term3 * term4 * term5 * term6 * term7 * term8;
                                }
                            }
                        }
                    }
                }
            }
            return (2.0 * PI.powf(2.5) / (p * q * (p + q).sqrt())) * val;
        }

        let mut vijkl = 0.0;
        for a in 0..orb1.gtos.len() {
            for b in 0..orb2.gtos.len() {
                for c in 0..orb3.gtos.len() {
                    for d in 0..orb4.gtos.len() {
                        let gto1 = &orb1.gtos[a];
                        let gto2 = &orb2.gtos[b];
                        let gto3 = &orb3.gtos[c];
                        let gto4 = &orb4.gtos[d];
                        vijkl += gto1.N * gto2.N * gto3.N * gto4.N * gto1.cc * gto2.cc * gto3.cc * gto4.cc * compute_Vabcd(gto1, gto2, gto3, gto4);
                    }
                }
            }
        }
        return vijkl;
    }

    fn compute_E(A: f64, B: f64, alpha: f64, beta: f64, ai: i32, bi: i32, t: i32) -> f64 {
        let p = alpha + beta;
        let q = alpha * beta / p;
        let Qx = A - B;
        if t < 0 || t > ai + bi {return 0.0;}
        if ai == 0 && bi == 0 && t == 0 {return (-q * Qx * Qx).exp();}
        if bi == 0 {
            let term1 = (1.0/(2.0*p)) * compute_E(A, B, alpha, beta, ai-1, 0, t-1);
            let term2 = - (q * Qx / alpha) * compute_E(A, B, alpha, beta, ai-1, 0, t);
            let term3 = (t + 1) * compute_E(A, B, alpha, beta, ai-1, 0, t+1);
            return term1 + term2 + term3;   
        }
        let term1 = (1.0/(2.0*p)) * compute_E(A, B, alpha, beta, ai, bi-1, t-1);
        let term2 = (q * Qx / beta) * compute_E(A, B, alpha, beta, ai, bi-1, t);
        let term3 = (t + 1) * compute_E(A, B, alpha, beta, ai, bi-1, t+1);
        return term1 + term2 + term3;   
    }

    fn compute_R(t: i32, u: i32, v: i32, n: i32, p: f64, P: (f64, f64, f64), C: (f64, f64, f64)) -> f64 {

        fn boys(n: u32, T: f64) -> f64 {
            return hyp1f1(n as f64 + 0.5, n as f64 + 1.5, -T) / (2.0 * n as f64 + 1.0);
        }
        
        let RPC = ((P.0 - C.0).powi(2) + (P.1 - C.1).powi(2) + (P.2 - C.2).powi(2)).sqrt();
        let mut val = 0.0;

        if t < 0 || u < 0 || v < 0 {
            return 0.0;
        }
        if t == 0 && u == 0 && v == 0{
            val += (-2.0 * p).powi(n + 1) * boys(n,p * RPC * RPC);
        }
        else if t == 0 && u == 0 {
            val += (v-1) as f64 * compute_R(0, 0, v-2, n+1, p, P, C);
            val += (P.2 - C.2) * compute_R(0, 0, v-1, n+1, p, P, C);
        }
        else if t == 0 {
            val += (u-1) as f64 * compute_R(0, u-2, v, n+1, p, P, C);
            val += (P.1 - C.1) * compute_R(0, u-1, v, n+1, p, P, C);
        }
        else {
            val += (t-1) as f64 * compute_R(t-2, u, v, n+1, p, P, C);
            val += (P.0 - C.0) * compute_R(t-1, u, v, n+1, p, P, C);
        }
        return val;
    }

}