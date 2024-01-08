use std::{ffi::OsString, fs::File, io::Write};

use pyo3::pyfunction;
use pyo3::PyErr;


/// https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html : PDB format
///
/// https://doc.rust-lang.org/std/fmt/index.html#syntax : Formatting syntax in Rust
#[pyfunction]
pub fn write_to_pdb(fname: OsString,  array: Vec<[f64;3]>, resname: String) -> Result<(), PyErr> {
    
    let mut fname: String = fname.to_str().expect("Passed argument `fname` contains invalid UTF-8").to_owned();
    if !fname.ends_with(".pdb"){ 
        fname.push_str(".pdb")
    };

    // Residue Name limitations of PBB format
    if resname.len() > 3 {
        panic!("Residue name cannot be larger than three characters")
    };

    let mut atomnames: Vec<&str> = vec!["O4'", "C1'", "C2'", "C3'", "C4'"];
    let mut atomnumbs: Vec<&str> = vec!["1", "2", "3", "4", "5"];

    // If array is a sixring system
    if array.len() == 6 {
        atomnames.push("C5'");
        atomnumbs.push("6");
    };

    let mut buffer = File::create(fname).expect("Cannot open file!");

    // Iterate over array of Coordinates and format the pdb file correctly
    for (i, aname) in atomnames.iter().enumerate() {

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

#[pyfunction]
pub fn write_to_xyz(fname: OsString, array: Vec<[f64;3]>) -> Result<(), PyErr> {

    let mut fname: String = fname.to_str().expect("Passed argument `fname` contains invalid UTF-8").to_owned();
    if !fname.ends_with(".xyz"){ 
        fname.push_str(".xyz")
    };
    let mut buffer = File::create(fname).expect("Cannot open file!");

    let mut elements: Vec<&str> = vec!["O", "C", "C", "C", "C"];
    if array.len() == 6 {
        elements.push("C");
    };

    buffer.write_all(format!("{}\n", array.len()).as_bytes()).expect("Cannot convert &str to bytes");
    buffer.write_all("Coordinates generated by pucke.py\n".as_bytes()).expect("Cannot convert &str to bytes");

    // Iterate over array of Coordinates and format the xyz file correctly
    for i in 0..array.len() {

        let coordinates = &array[i];
        let content = format!(
            "{:>2} {:width$.precision$}   {:width$.precision$}   {:width$.precision$}\n",
            elements[i], coordinates[0], coordinates[1], coordinates[2], width=19, precision=14 
            // Element symbol, x coord, y, coord, z coord
            );
        buffer.write_all(content.as_bytes()).expect("Cannot convert String to bytes");
    }
    
    Ok(())
}
