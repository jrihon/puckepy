use std::fs;
use pyo3::{pyfunction, pyclass,PyErr};

/// The only thing we need from the pdb is 
/// Atom names Vec<String>
/// Coordinates, best to do as Vec<[f32;3]>
#[pyclass(get_all)]
pub struct Pdb {
    atom_names: Vec<String>,
    coordinates: Vec<[f32;3]>
}

/// The only thing we need from the xyz is 
/// Coordinates, best to do as Vec<[f32;3]>
#[pyclass(get_all)]
pub struct Xyz {
    coordinates: Vec<[f32;3]>
}

/// Parses an xyz-file format
/// This means a format that looks like this
/// ```
/// ATOM      1  O6'  MA    41      24.802  52.534  40.016  1.00  0.00           O  
/// ATOM      2  C6'  MA    41      24.803  51.735  41.199  1.00  0.00           C  
/// ATOM      3 H6'1  MA    41      25.476  50.878  41.168  1.00  0.00           H  
/// ATOM      4 H6'2  MA    41      23.806  51.294  41.182  1.00  0.00           H  
/// ATOM      5  C5'  MA    41      25.097  52.567  42.397  1.00  0.00           C  
///
/// ```
#[pyfunction]
pub fn from_pdb(fname: String) -> Result<Pdb, PyErr> {

    if !fname.ends_with(".pdb") {
        panic!("The {} is not a valid `.pdb` file format ", &fname)
    };

    let file = match fs::read_to_string(&fname) {
        Ok(a) => a,
        Err(e)=> panic!("{}", e),
    };

    if file.is_empty() { 
        panic!("The file {} is empty.", &fname)
    }

    let pdblines = file.lines().map(|s| s.into()).collect::<Vec<String>>();

    let mut atom_names: Vec<String> = vec![];
    let mut coordinates: Vec<[f32;3]> = vec![];

    for l in pdblines.iter() {
        if l.starts_with("ATOM") { 
            atom_names.push(l[12..16].trim().into());

            let x = match l[31..39].trim().parse::<f32>() {
                Ok(a) => a,
                Err(e) => panic!("{}", e)
            };
            let y = match l[39..47].trim().parse::<f32>() {
                Ok(a) => a,
                Err(e) => panic!("{}", e)
            };
            let z = match l[47..55].trim().parse::<f32>() {
                Ok(a) => a,
                Err(e) => panic!("{}", e)
            };
            coordinates.push([x,y,z]);

        }
    }
    
    Ok(Pdb {
        atom_names,
        coordinates
    })
}




/// Parses an xyz-file format
/// This means a format that looks like this
/// ```
/// 31
///Coordinates from ORCA-job Conformation_X
///  H   4.01196826057662      2.03352821967286      2.01847309650732
///  O   3.76770440038636      1.71999235396699      1.14581624607411
///  C   2.53548022010070      2.32709191442346      0.78140278302649
///  H   2.69801965937301      3.28480341404723      0.28455391459758
///
/// ```
#[pyfunction]
pub fn from_xyz(fname: String) -> Result<Xyz, PyErr> {

    if !fname.ends_with(".xyz") {
        panic!("The {} is not a valid `.xyz` file format ", &fname)
    };

    let file = match fs::read_to_string(&fname) {
        Ok(a) => a,
        Err(e)=> panic!("{}", e),
    };

    if file.is_empty() { 
        panic!("The file {} is empty.", &fname)
    }

    let mut coordinates: Vec<[f32;3]> = vec![];
    let xyz_lines = file.lines().map(|s| s.into()).collect::<Vec<String>>();
    let mut xyz_iter = xyz_lines.iter();

    // Two next calls, because xyz files always start with two lines of header data
    // We just discard this
    let _ = &xyz_iter.next();
    let _ = &xyz_iter.next();

    for l in xyz_iter {
        let spl: Vec<&str> = l.split_whitespace().collect();

        if spl.len() != 4 {
            continue
        };

        let x = match spl[1].parse::<f32>() {
            Ok(a) => a,
            Err(e) => panic!("{}", e)
        };
        let y = match spl[2].trim().parse::<f32>() {
            Ok(a) => a,
            Err(e) => panic!("{}", e)
        };
        let z = match spl[3].trim().parse::<f32>() {
            Ok(a) => a,
            Err(e) => panic!("{}", e)
        };
        coordinates.push([x,y,z]);

    }

    Ok(Xyz { coordinates})
}
