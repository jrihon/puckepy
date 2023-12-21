use pyo3::{pymethods,pyclass, PyErr};

use std::fs;
//use crate::molfile::Xyz;

use self::cremerpople::{
    geometric_center_of_molecule,
    molecular_axis,
    local_elevation,
    return_cp_coordinates
};
//
pub mod cremerpople;


/// The CP tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct CP {
    pub amplitude: f32,
    pub phase_angle: f32,
}

#[pymethods]
impl CP {

    #[new]
    fn new(amplitude: f32, phase_angle: f32) -> CP {
        if amplitude > 1. {
            panic!("amplitude value is larger than 1.")
        }

        if !(0.0..=360.0).contains(&phase_angle) {
            panic!("phase_angle value should be within the range of 0 -> 360")
        }

        CP { amplitude, phase_angle }
    }

    fn to_as(&self) -> AS {
        // let mut phase_angle = self.1 - 90.;
        // if phase_angle < 0. { 
        //     phase_angle += 360.
        // }; => Original code

        // If the value is smaller than 0 after decreasing 90, it is already smaller than 90
        // This means that we will do two operations, a -90 and then +360
        // This cuts out an operation or two down the line
        let new_angle = if self.phase_angle < 90. { self.phase_angle + 270. } else { self.phase_angle - 90. };

        AS { amplitude: self.amplitude, phase_angle : new_angle }

    }
    
}


/// The AS tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct AS {
    pub amplitude: f32,
    pub phase_angle: f32,
}

#[pymethods]
impl AS {

    #[new]
    fn new(amplitude: f32, phase_angle: f32) -> AS {
        if amplitude > 1. {
            panic!("amplitude value is larger than 1.")
        }

        if !(0.0..=360.0).contains(&phase_angle) {
            panic!("phase_angle value should be within the range of 0 -> 360")
        }
        AS { amplitude, phase_angle }
    }
    
    fn to_cp(&self) -> CP {
        // let mut phase_angle = self.1 + 90.;
        // if phase_angle > 360. { 
        //     phase_angle -= 360.
        // }; => Original code

        // If the value is larger than 360 after adding 90, it is already larger than 270
        // This means that we will do two operations, a +90 and then -360
        // This cuts out an operation or two down the line
        let new_angle = if self.phase_angle > 270. { self.phase_angle - 270. } else { self.phase_angle + 90. };

        CP {amplitude: self.amplitude, phase_angle: new_angle }

    }
}





/// The only thing we need from the pdb is 
/// Atom names Vec<String>
/// Coordinates, best to do as Vec<[f32;3]>
#[pyclass(get_all)]
pub struct Pdb {
    pub atom_names: Vec<String>,
    pub coordinates: Vec<[f32;3]>
}
/// Parses an pdb-file format
/// This means a format that looks like this
/// ```
/// ATOM      1  O6'  MA    41      24.802  52.534  40.016  1.00  0.00           O  
/// ATOM      2  C6'  MA    41      24.803  51.735  41.199  1.00  0.00           C  
/// ATOM      3 H6'1  MA    41      25.476  50.878  41.168  1.00  0.00           H  
/// ATOM      4 H6'2  MA    41      23.806  51.294  41.182  1.00  0.00           H  
/// ATOM      5  C5'  MA    41      25.097  52.567  42.397  1.00  0.00           C  
///
/// ```



#[pymethods]
impl Pdb {

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


    fn cp_from_indices(&self, indices: Vec<usize>) -> CP {
        
        let mut molarray = vec![];

        for idx in indices {
            molarray.push(self.coordinates[idx])
        };

        geometric_center_of_molecule(&mut molarray);
        let mol_axis = molecular_axis(&molarray);
        let zj = local_elevation(&molarray, mol_axis);
        let (amplitude, phase_angle) = return_cp_coordinates(zj);

        CP { amplitude, phase_angle }
    }
    
    fn cp_from_atomnames(&self, query_names: Vec<String>) -> CP {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        // Search for the indices of the atom names
        for name in query_names.iter() {
            match self.atom_names.at(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find {} in the queried pdb.", name)
            }
        }

        // Call cp_from_indices
        self.cp_from_indices(indices)
    }
}

trait FindString {
    fn at(&self, pattern: &str) -> Result<usize, ()> ;
}

impl FindString for Vec<String> {

    fn at(&self, pattern: &str) -> Result<usize, ()> {

        let mut c = 0;

        for name in self {
            if name != pattern {
                c += 1
            } else { 
                return Ok(c as usize) 
            }
        };

        Err(())

//        let a = self.iter().try_fold( 0_u32, |mut prev: u32, atomname| 
//            {
//                if atomname == pattern {
//                    ControlFlow::Continue(prev += 1_u32)
//                } else {
//                    ControlFlow::Break(prev)
//                }
//            });
//
//        Ok(4)
    }
    
}


/// The only thing we need from the xyz is 
/// Coordinates, best to do as Vec<[f32;3]>
#[pyclass(get_all)]
pub struct Xyz {
    pub coordinates: Vec<[f32;3]>
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
#[pymethods]
impl Xyz {
    
    #[new]
    fn from_xyz(fname: String) -> Result<Xyz, PyErr> {

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

    fn cp_from_indices(&self, indices: Vec<usize>) -> CP {
        
        let mut molarray = vec![];

        for idx in indices {
            molarray.push(self.coordinates[idx])
        };

        geometric_center_of_molecule(&mut molarray);
        let mol_axis = molecular_axis(&molarray);
        let zj = local_elevation(&molarray, mol_axis);
        let (amplitude, phase_angle) = return_cp_coordinates(zj);

        CP { amplitude, phase_angle }
    }
    
}
