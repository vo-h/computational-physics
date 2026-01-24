use crate::hf::basis_stog::stog_orbital::STOGOrbital;
use reqwest::blocking::get;
use reqwest::Error;

#[derive(Clone, Debug)]
pub struct Atom {
    pub symbol: String,
    pub Z: u8,
    pub coords: (f64, f64, f64),
    pub coords_units: String,
    pub basis: String,
    pub orbitals: Vec<STOGOrbital>,
}

impl Atom {

    pub fn new(symbol: &str, coords: (f64, f64, f64), Z: u8, basis: String) -> Result<(Self), Error>  {
        let url = format!("https://www.basissetexchange.org/basis/{}/format/gaussian94/?version=1&elements={}&uncontract_spdf=true", basis, symbol);
        let res = get(url)?.text()?;
        let sections = res.split("\n\n").collect::<Vec<&str>>();
        let n = sections.len();

        let info = sections[n-3].to_string().trim().to_string();
        let lines = info.split('\n').collect::<Vec<&str>>();
        let mut orbitals = Vec::<STOGOrbital>::new();
        
        let mut curr: String = String::new();
        for line in lines[1..].iter() {
            // Case: new orbital section
            if curr.len() == 0 && ! line.starts_with(' ') {
                curr += &line.to_string();

            // Case: continuation of current orbital type
            } else if curr.len() > 0 && line.starts_with(' ') {
                curr += &format!("\n{}", line);
            
            // Case: new orbital type, but we have already collected one
            } else if curr.len() > 0 && ! line.starts_with(' ') {
                let mut new_orbitals = STOGOrbital::from_string("O".to_string(), coords, &curr);
                orbitals.append(&mut new_orbitals);
                curr = line.to_string();
            }
        }

        if curr.len() > 0 && !curr.contains("****") {
            let mut new_orbitals = STOGOrbital::from_string("O".to_string(), coords, &curr);
            orbitals.append(&mut new_orbitals);
        }

        Ok(Atom {
            symbol: symbol.to_string(),
            Z: Z,
            coords: coords,
            coords_units: "bohr".to_string(),
            basis: format!("sto-{}g", basis),
            orbitals: orbitals,
        })
    }
} 
