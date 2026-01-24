#![allow(non_snake_case)]
#[allow(dead_code)]

mod hf;
use reqwest::blocking::get;
use reqwest::Error;
use clap::Parser;
use crate::hf::atom::Atom;
use crate::hf::basis_stog::stog_orbital::STOGOrbital;
use crate::hf::basis_stog::stog_primitive::STOGPrimitive;

/// Simple program to greet a person
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Location of input file containing molecule specification
    #[arg(short, long)]
    input: std::path::PathBuf,

    /// Method to compute
    #[arg(short, long)]
    method: String,
}

fn main() -> () {
    let args = Args::parse();
    let input = std::fs::read_to_string(args.input).expect("Failed to read input file");
    
    let molecule = hf::molecule::Molecule::from_string(&input);
    if args.method == "overlap" {
        let S = molecule.compute_S();
        println!("Overlap matrix S:\n {:#?}", S);
        return;
    }
}
