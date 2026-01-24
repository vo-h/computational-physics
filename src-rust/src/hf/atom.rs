pub mod atom {
    struct Atom {
        symbol: String,
        Z: u8,
        coords: (f64, f64, f64),
        coords_units: String,
        basis: String,
    }

    impl Atom {
        pub fn new(symbol: String, Z: u8, coords: (f64, f64, f64), coords_units: String, basis: String) -> Self {
            Atom {
                symbol,
                Z,
                coords,
                coords_units,
                basis,
            }
        }
    } 
}