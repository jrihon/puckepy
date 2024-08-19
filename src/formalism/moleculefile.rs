use core::panic;
use std::fs::read_to_string;
use pyo3::{pyclass, pymethods, PyErr, pyfunction};

use std::{ffi::OsString, fs::File, io::Write};



// Read in queried file and see if it is valid
fn validate_contents(fname: &String, fileformat: &str) -> String {

    if !fname.ends_with(fileformat) {
        panic!("The {} is not a valid `{}` file format ", &fname, fileformat)
    };

    // Read contents once
    let filecontents = match read_to_string(&fname) {
        Ok(contents) => contents,
        Err(e) => panic!("{}", e)
    };

    // Check if contents exists
    if filecontents.is_empty() { panic!("The {} is empty!", &filecontents)};

    filecontents
    
}

// anames : atomnames Vector
// coords : coordinates Vector
fn populate_pdbfields(lines: &String, anames_container: &mut Vec<String>, coords_container: &mut Vec<[f64;3]>) {

    anames_container.push(lines[12..16].trim().into());

    let x = match lines[31..39].trim().parse::<f64>() {
        Ok(a) => a,
        Err(e) => panic!("Cannot parse x-coordinate : {}. At line\n{}", e, &lines)
    };
    let y = match lines[39..47].trim().parse::<f64>() {
        Ok(a) => a,
        Err(e) => panic!("Cannot parse y-coordinate : {}. At line\n{}", e, &lines)
    };
    let z = match lines[47..55].trim().parse::<f64>() {
        Ok(a) => a,
        Err(e) => panic!("Cannot parse z-coordinate : {}. At line\n{}", e, &lines)
    };
    coords_container.push([x,y,z]);



}
//pub struct FileContents {
//    fc: String // filecontents
//}
//
/// The only thing we need from the pdb is 
/// Atom names Vec<String>
/// Coordinates, best to do as Vec<[f64;3]>
#[pyclass(get_all)]
pub struct Pdb {
    pub data : String,
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
    fn new(filename: String) -> Result<Pdb, PyErr> {

        let filecontents = validate_contents(&filename, ".pdb");

        Ok(Pdb {
            data: filecontents,
            atomnames: vec![],
            coordinates: vec![],
        })
    }

    fn parse(&self) -> Pdb {

        if self.atomnames.len() > 0 || self.coordinates.len() > 0 {
            panic!("This Pdb object has already been populated. Will not parse again.")
        };

        let pdblines = &self.data.lines()
                                 .map(|s| s.into())
                                 .collect::<Vec<String>>();

        let mut atomnames: Vec<String> = vec![];
        let mut coordinates: Vec<[f64;3]> = vec![];

        for lines in pdblines.iter() {
            if lines.starts_with("ATOM") || lines.starts_with("HETATM") { 
                populate_pdbfields(lines, &mut atomnames, &mut coordinates);
            }
        };
        
        Pdb { 
            data: self.data.to_string(),
            coordinates,
            atomnames,
        }
    }


    // Go over the molecular structure and parse by the change of residue numbers. 
    // Iterate if there are multiple residue numbers to begin with, 
    // then store a Vec of Pdb structs and return this
    fn parse_by_monomers(&self) -> Vec<Pdb> {

        let pdblines = &self.data.lines()
                                 .map(|s| s.into())
                                 .collect::<Vec<String>>();

        let mut pdbs: Vec<Pdb> = vec![];
        let mut resnumber: u16 = 42069; // residue names can only go to 9999, so this is safe :^)

        let mut atomnames_container: Vec<String> = vec![];
        let mut coordinates_container: Vec<[f64;3]> = vec![];


        let mut iter_lines = pdblines.iter();
        while let Some(lines) = iter_lines.next() {

            if lines.starts_with("ATOM") || lines.starts_with("HETATM") { 

                // Check residue number first
                let parsed_resname: u16 = match lines[22..26].trim().parse() {
                    Ok(a) => a,
                    Err(_) => panic!("Residue number cannot be parsed as an integer, at\n{}", &lines)
                };
                // If there is a valid u16, set it as the current residue number, for the first one
                // parsed from the file
                if resnumber == 42069 { 
                    resnumber = parsed_resname 
                };

                if resnumber != parsed_resname {
                    // Drain the atomnames and coordinates Vecs into a Pdb 
                    // Push the Pdb onto the Vec<Pdb>
                    pdbs.push( Pdb {
                                   data: "monomer_".to_string() + &resnumber.to_string(),
                                   atomnames: atomnames_container.drain(..).collect(),
                                   coordinates: coordinates_container.drain(..).collect(),
                               }
                    );

                    resnumber = parsed_resname; // reset the parsed residuename to the residue name

                    // Start pushing to the cleared Vecs at the current line for a new Pdb struct
                    populate_pdbfields(lines, &mut atomnames_container, &mut coordinates_container);

                } else {
                    // Populate the current containers as normal
                    populate_pdbfields(lines, &mut atomnames_container, &mut coordinates_container);
                }
            }
        }
        // Drain the final atomnames and coordinates Vecs into the last Pdb 
        pdbs.push( Pdb {
                       data: "monomer_".to_string() + &resnumber.to_string(),
                       atomnames: atomnames_container.drain(..).collect(),
                       coordinates: coordinates_container.drain(..).collect(),
                   }
        );

        pdbs // return Vec<Pdb>
    }
}











/// The only thing we need from the xyz is 
/// Coordinates, best to do as Vec<[f64;3]>
#[pyclass]
pub struct Xyz {
    filecontents: String,
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
    fn new(filename: String) -> Result<Xyz, PyErr> {

        let filecontents = validate_contents(&filename, ".xyz");

        Ok(Xyz {
            filecontents,
        })

    }

    // Parses filecontents and returns an array of coordinates
    fn parse(&self) -> Vec<[f64;3]> {

        let mut coordinates: Vec<[f64;3]> = vec![];
        let xyz_lines = &self.filecontents.lines()
                                  .map(|s| s.into())
                                  .collect::<Vec<String>>();
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

        coordinates
    }
}






/// https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html : PDB format
///
/// https://doc.rust-lang.org/std/fmt/index.html#syntax : Formatting syntax in Rust
#[pyfunction]
pub fn write_to_pdb(filename: OsString,  coordinates: Vec<[f64;3]>, residuename: String) -> Result<(), PyErr> {
    
    let mut filename: String = filename.to_str().expect("Passed argument `filename` contains invalid UTF-8").to_owned();
    if !filename.ends_with(".pdb"){ 
        filename.push_str(".pdb")
    };

    // Residue Name limitations of PBB format
    if residuename.len() > 3 {
        panic!("Residue name cannot be larger than three characters")
    };

    let mut atomnames: Vec<&str> = vec!["O4'", "C1'", "C2'", "C3'", "C4'"];
    let mut atomnumbs: Vec<&str> = vec!["1", "2", "3", "4", "5"];

    // If coordinates is a sixring system
    if coordinates.len() == 6 {
        atomnames.push("C5'");
        atomnumbs.push("6");
        atomnames[0] = "O5'"; 
    };

    let mut buffer = File::create(filename).expect("Cannot open file!");

    // Iterate over coordinates of Coordinates and format the pdb file correctly
    for (i, aname) in atomnames.iter().enumerate() {

        let coordinate = &coordinates[i];
        let content = format!(
            "ATOM   {:>4} {:<4} {:>3} A   1    {:width$.precision$}{:width$.precision$}{:width$.precision$}  {:>22}\n",
            atomnumbs[i], aname, residuename, coordinate[0], coordinate[1], coordinate[2], aname.chars().nth(0).unwrap(), width=8, precision=3 
            // Atom number, Atom name, residue name, x coord, y, coord, z coord, element symbol
            );
        buffer.write_all(content.as_bytes()).expect("Cannot convert String to bytes");
    }
    
    Ok(())
}

#[pyfunction]
pub fn write_to_xyz(filename: OsString, coordinates: Vec<[f64;3]>) -> Result<(), PyErr> {

    let mut filename: String = filename.to_str().expect("Passed argument `filename` contains invalid UTF-8").to_owned();
    if !filename.ends_with(".xyz"){ 
        filename.push_str(".xyz")
    };
    let mut buffer = File::create(filename).expect("Cannot open file!");

    let mut elements: Vec<&str> = vec!["O", "C", "C", "C", "C"];
    if coordinates.len() == 6 {
        elements.push("C");
    };

    buffer.write_all(format!("{}\n", coordinates.len()).as_bytes()).expect("Cannot convert &str to bytes");
    buffer.write_all("Coordinates generated by pucke.py\n".as_bytes()).expect("Cannot convert &str to bytes");

    // Iterate over coordinates of Coordinates and format the xyz file correctly
    for i in 0..coordinates.len() {

        let coordinate = &coordinates[i];
        let content = format!(
            "{:>2} {:width$.precision$}   {:width$.precision$}   {:width$.precision$}\n",
            elements[i], coordinate[0], coordinate[1], coordinate[2], width=19, precision=14 
            // Element symbol, x coord, y, coord, z coord
            );
        buffer.write_all(content.as_bytes()).expect("Cannot convert String to bytes");
    }
    
    Ok(())
}
