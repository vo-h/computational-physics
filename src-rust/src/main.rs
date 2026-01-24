#![allow(non_snake_case)]
#[allow(dead_code)]

mod hf;
use clap::Parser;

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
    if args.method == "1e-mat" {
        let S = molecule.compute_S();
        println!("Overlap matrix S:\n {a:.5}", a=S);

        let T = molecule.compute_T();
        println!("Kinetic energy matrix T:\n {a:.5}", a=T);

        let V = molecule.compute_Vne();
        println!("Nuclear attraction matrix Vne:\n {a:.5}", a=V);
        return;
    } else if args.method == "2e-mat" {
        let Vee = molecule.compute_Vee();
        for i in 0..Vee.len() {
            for j in 0..Vee.len() {
                // Print slice Vee[i][j][:][:]
                println!("Vee[{}][{}]:", i, j);
                for k in 0..Vee.len() {
                    for l in 0..Vee.len() {
                        print!("{:.5} ", Vee[i][j][k][l]);
                    }
                    println!();
                }
            }
        }
    } else if args.method == "guesses" {
        let F0 = molecule.compute_F0();
        println!("Initial Fock matrix F0:\n {a:.5}", a=F0);
        let C0 = molecule.compute_C0();
        println!("Initial coefficient matrix C0:\n {a:.5}", a=C0);
        let D0 = molecule.compute_D0();
        println!("Initial density matrix D0:\n {a:.5}", a=D0);
        let E0 = molecule.compute_E0();
        println!("Initial electronic energy E0: {:.5}", E0);
    } else if args.method == "scf" {
        molecule.compute_CHF(1e-6, 100);
    } else {
        println!("Unknown method: {}", args.method);
    }
}
