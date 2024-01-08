use std::fs;
use pyo3::{pyclass, pymethods, PyErr};

/// The only thing we need from the pdb is 
/// Atom names Vec<String>
/// Coordinates, best to do as Vec<[f64;3]>
#[pyclass(get_all)]
pub struct Pdb {
    pub fname : String,
    pub atomnames: Vec<String>,
    pub coordinates: Vec<[f64;3]>
}
/// Parses an pdb-file format
/// This means a format that looks like this
/// ```
/// ATOM      1  O6'  MA    41      24.802  52.534  40.016  1.00  0.00           O  
/// ATOM      2  C6'  MA    41      24.803  51.735  41.199  1.00  0.00           C  
/// ATOM      3 H6'1  MA    41      25.476  50.878  41.168  1.00  0.00           H  
/// ATOM      4 H6'2  MA    41      23.806  51.294  41.182  1.00  0.00           H  
/// ATOM      5  C5'  MA    41      25.097  52.567  42.397  1.00  0.00           C  
/// ```
#[pymethods]
impl Pdb {

    // Result<Pdb,PyErr>  
    // This is a Result type because the user might mistype the name of the file,
    // causing the function to appropriately crash
    #[new]
    fn new(fname: String) -> Result<Pdb, PyErr> {

        if !fname.ends_with(".pdb") {
            panic!("The {} is not a valid `.pdb` file format ", &fname)
        };

        let file = match fs::read_to_string(&fname) {
            Ok(a) => a,
            Err(e)=> panic!("{}", e),
        };

        if file.is_empty() { 
            panic!("The file {} is empty.", &fname)
        }; 

        Ok(Pdb {
            fname,
            atomnames: vec![],
            coordinates: vec![],
        })
    }

    fn parse(&self) -> Pdb {

        let pdblines = self.fname.lines().map(|s| s.into()).collect::<Vec<String>>();

        let mut atomnames: Vec<String> = vec![];
        let mut coordinates: Vec<[f64;3]> = vec![];

        for lines in pdblines.iter() {
            if lines.starts_with("ATOM") { 
                atomnames.push(lines[12..16].trim().into());

                let x = match lines[31..39].trim().parse::<f64>() {
                    Ok(a) => a,
                    Err(e) => panic!("{}", e)
                };
                let y = match lines[39..47].trim().parse::<f64>() {
                    Ok(a) => a,
                    Err(e) => panic!("{}", e)
                };
                let z = match lines[47..55].trim().parse::<f64>() {
                    Ok(a) => a,
                    Err(e) => panic!("{}", e)
                };
                coordinates.push([x,y,z]);

            }
        }
        
        Pdb {
            atomnames,
            coordinates,
            fname: self.fname.to_owned()
        }
    }

    // TODO: 
    // Go over the molecular structure and parse by the change of residue numbers. 
    // For iterate if there are multiple residue numbers to begin with, 
    // then store a Vec of Pdb structs and return this
    fn parse_by_monomers(&self) -> Vec<Pdb> {
        todo!()
    }
}











/// The only thing we need from the xyz is 
/// Coordinates, best to do as Vec<[f64;3]>
#[pyclass(get_all)]
pub struct Xyz {
    pub coordinates: Vec<[f64;3]>
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
/// ```
#[pymethods]
impl Xyz {
    
    // Result<Pdb,PyErr>  
    // This is a Result type because the user might mistype the name of the file,
    // causing the function to appropriately crash
    #[new]
    fn new(fname: String) -> Result<Xyz, PyErr> {

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

        let mut coordinates: Vec<[f64;3]> = vec![];
        let xyz_lines = file.lines().map(|s| s.into()).collect::<Vec<String>>();
        let mut xyz_iter = xyz_lines.iter();

        // Two next calls, because xyz files always start with two lines of header data
        // We just discard this
        let _ = &xyz_iter.next();
        let _ = &xyz_iter.next();

        for l in xyz_iter {
            let splits: Vec<&str> = l.split_whitespace().collect();

            if splits.len() != 4 {
                continue
            };

            let x = match splits[1].parse::<f64>() {
                Ok(a) => a,
                Err(e) => panic!("{}", e)
            };
            let y = match splits[2].trim().parse::<f64>() {
                Ok(a) => a,
                Err(e) => panic!("{}", e)
            };
            let z = match splits[3].trim().parse::<f64>() {
                Ok(a) => a,
                Err(e) => panic!("{}", e)
            };
            coordinates.push([x,y,z]);

        }

        Ok(Xyz { coordinates})
    }
}
