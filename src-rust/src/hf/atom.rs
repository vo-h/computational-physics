pub mod atom {
    use crate::hf::basis_stog::stog_orbital::STOGOrbital;

    struct Atom {
        symbol: String,
        Z: u8,
        coords: (f64, f64, f64),
        coords_units: String,
        basis: String,
        orbitals: Vec<STOGOrbital>,
    }

    impl Atom {
        pub fn new(symbol: String, Z: u8, coords: (f64, f64, f64), coords_units: String, basis: String) -> Self {
            Atom {
                symbol,
                Z,
                coords,
                coords_units,
                basis,
                orbitals: Vec::new(),
            }
        }
    } 
}