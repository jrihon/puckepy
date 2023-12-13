use std::{ffi::OsString, fs::File, io::Write};

use pyo3::pyfunction;
use pyo3::PyErr;

pub mod sixring;
pub mod fivering;








/// https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html : PDB format
///
/// https://doc.rust-lang.org/std/fmt/index.html#syntax : Formatting syntax in Rust
#[pyfunction]
pub fn write_pdb(fname: OsString,  array: [[f32;3];6], resname: String) -> Result<(), PyErr> {
    
    if resname.len() > 3 {
        panic!("Residue name cannot be larger than three characters")

    };
    let mut atomnames: Vec<&str> = vec!["O4'", "C1'", "C2'", "C3'", "C4'"];
    let mut atomnumbs: Vec<&str> = vec!["1", "2", "3", "4", "5"];

    if array.len() == 6 {
        atomnames.push("C5'");
        atomnumbs.push("6");
    };


    let mut buffer = File::create(fname).expect("Cannot open file!");

    for (i, aname) in atomnames.iter().enumerate() { // Range type

        let coordinates = &array[i];
        let content = format!(
            "ATOM   {:>4} {:<4} {:>3} A   1    {:width$.precision$}{:width$.precision$}{:width$.precision$}  {:>22}\n",
            atomnumbs[i], aname, resname, coordinates[0], coordinates[1], coordinates[2], aname.chars().nth(0).unwrap(), width=8, precision=3 
            // Atom number, Atom name, residue name, x coord, y, coord, z coord, element symbol
            );
        buffer.write_all(content.as_bytes()).expect("Cannot convert String to bytes");
    }
    
    Ok(())
}
