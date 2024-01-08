use pyo3::{pymethods,pyclass, PyErr};
use std::fs;

pub mod cremerpople;
pub mod altonasund;

mod search_atomname;
use self::search_atomname::FindString; // Find name of &str 'String'.at_position()

const PIS_IN_180: f64 = 57.2957795130823208767981548141051703_f64;              // taken from <f64>.to_degrees() for i in range(NUM):




// Enum to control the which type of n-membered ring system is being produced and 
// returns the correct one to the user.
// Acts as an addition safety measure whenever users prompt incorrect amount of values in function
// calls
pub enum MemberedRing {
    Five(CP5),
    Six(CP6)
}

/// The CP tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct CP5 {
    pub amplitude: f64,
    pub phase_angle: f64,
}

#[pyclass(get_all)]
pub struct CP6 {
    pub amplitude: f64,
    pub phase_angle: f64,
    pub theta: f64,
}

#[pymethods]
impl CP5 {

    #[new]
    fn new(amplitude: f64, phase_angle: f64) -> CP5 {
        if amplitude > 1. {
            panic!("amplitude value is larger than 1.")
        }

        if !(0.0..=360.0).contains(&phase_angle) {
            panic!("phase_angle value should be within the range of 0 -> 360")
        }

        CP5 { amplitude, phase_angle }
    }

    // Calculate Cremer-Pople formalism by prompted indices
    fn cp_from_indices(&self, coordinates : Vec<[f64; 3]>, indices: Vec<usize>) -> CP5 {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        for idx in indices {
            molarray.push(coordinates[idx])
        }

       match cremerpople::cremer_pople(&mut molarray) {
           MemberedRing::Five(a) => a,
           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
       }
    }
    
    // Find indices of atomnames and pass them to self.cp_from_indices()
    fn cp_from_atomnames(&self, pdb : &Pdb, query_names: Vec<String>) -> CP5 {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        // Search for the indices of the atom names
        for name in query_names.iter() {
            match pdb.atomnames.at_position(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
            }
        }

        self.cp_from_indices(pdb.coordinates.clone(), indices)
    }

    fn to_as_angle(&self) -> f64 {
        // let mut phase_angle = self.1 - 90.;
        // if phase_angle < 0. { 
        //     phase_angle += 360.
        // }; => Original code

        // If the value is smaller than 0 after decreasing 90, it is already smaller than 90
        // This means that we will do two operations, a -90 and then +360
        // This cuts out an operation or two down the line
        if self.phase_angle < 90. { self.phase_angle + 270. } else { self.phase_angle - 90. }

    }
    
}


#[pymethods]
impl CP6 {

    #[new]
    fn new(amplitude: f64, phase_angle: f64, theta: f64) -> CP6 {
        if amplitude > 1. {
            panic!("amplitude value is larger than 1.")
        }

        if !(0.0..=360.0).contains(&phase_angle) {
            panic!("phase_angle value should be within the range of 0 -> 360")
        }

        if !(0.0..=180.0).contains(&theta) {
            panic!("theta value should be within the range of 0 -> 180")
        }

        CP6 { amplitude, phase_angle, theta }
    }

    // Calculate Cremer-Pople formalism by prompted indices
    fn cp_from_indices(&self, coordinates : Vec<[f64; 3]>, indices: Vec<usize>) -> CP6 {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        for idx in indices {
            molarray.push(coordinates[idx])
        }

       match cremerpople::cremer_pople(&mut molarray) {
           MemberedRing::Six(a) => a,
           _ => panic!("An amount, not equal to 6, has been queried. Expected 6 elements.")
       }
    }
    
    // Find indices of atomnames and pass them to self.cp_from_indices()
    fn cp_from_atomnames(&self, pdb : &Pdb, query_names: Vec<String>) -> CP6 {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        // Search for the indices of the atom names
        for name in query_names.iter() {
            match pdb.atomnames.at_position(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
            }
        }

        self.cp_from_indices(pdb.coordinates.clone(), indices)
    }

}


/// The AS tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct AS {
    pub amplitude: f64,
    pub phase_angle: f64,
}

#[pymethods]
impl AS {

    #[new]
    fn new(amplitude: f64, phase_angle: f64) -> AS {
        if amplitude > 1. {
            panic!("amplitude value is larger than 1.")
        }

        if !(0.0..=360.0).contains(&phase_angle) {
            panic!("phase_angle value should be within the range of 0 -> 360")
        }
        AS { amplitude, phase_angle }
    }
    
    // Calculate Altona Sundaralingam formalism by the indices
    fn as_from_indices(&self, pdb: &Pdb, indices: Vec<usize>) -> AS {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        for idx in indices {
            molarray.push(pdb.coordinates[idx])
        }

        altonasund::altona_sundaralingam(&mut molarray)

//       match cremerpople::cremer_pople(&mut molarray) {
//           MemberedRing::Five(a) => a.to_as(),
//           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
//       }
    }
    
    // Find the indices of the atomnames and pass them to self.as_from_indices()
    fn as_from_atomnames(&self, pdb: &Pdb, query_names: Vec<String>) -> AS {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        // Search for the indices of the atom names
        for name in query_names.iter() {
            match pdb.atomnames.at_position(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
            }
        }

        // Call cp_from_indices
        self.as_from_indices(pdb, indices)
    }

    fn to_cp(&self) -> f64 {
        // let mut phase_angle = self.1 + 90.;
        // if phase_angle > 360. { 
        //     phase_angle -= 360.
        // }; => Original code

        // If the value is larger than 360 after adding 90, it is already larger than 270
        // This means that we will do two operations, a +90 and then -360
        // This cuts out an instruction or two down the line
        if self.phase_angle > 270. { self.phase_angle - 270. } else { self.phase_angle + 90. }

    }
}





/// The only thing we need from the pdb is 
/// Atom names Vec<String>
/// Coordinates, best to do as Vec<[f64;3]>
#[pyclass(get_all)]
pub struct Pdb {
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
        
        Ok(Pdb {
            atomnames,
            coordinates
        })
    }


//    // Calculate Cremer-Pople formalism by prompted indices
//    fn cp_from_indices(&self, indices: Vec<usize>) -> MemberedRing {
//        
//        let mut molarray: Vec<[f64; 3]> = vec![];
//
//        for idx in indices {
//            molarray.push(self.coordinates[idx])
//        }
//
//       cremerpople::cremer_pople(&mut molarray)
//    }
//    
//    // Find indices of atomnames and pass them to self.cp_from_indices()
//    fn cp_from_atomnames(&self, query_names: Vec<String>) -> MemberedRing {
//
//        // Make empty vec :
//        let mut indices: Vec<usize> = Vec::with_capacity(6);
//
//        // Search for the indices of the atom names
//        for name in query_names.iter() {
//            match self.atomnames.at_position(name) {
//                Ok(a) => indices.push(a),
//                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
//            }
//        }
//
//        self.cp_from_indices(indices)
//    }
//
//    // Calculate Altona Sundaralingam formalism by the indices
//    fn as_from_indices(&self, indices: Vec<usize>) -> AS {
//        
//        let mut molarray: Vec<[f64; 3]> = vec![];
//
//        for idx in indices {
//            molarray.push(self.coordinates[idx])
//        }
//
//        altonasund::altona_sundaralingam(&molarray)
//    }
//    
//    // Find the indices of the atomnames and pass them to self.as_from_indices()
//    fn as_from_atomnames(&self, query_names: Vec<String>) -> AS {
//
//        // Make empty vec :
//        let mut indices: Vec<usize> = Vec::with_capacity(6);
//
//        // Search for the indices of the atom names
//        for name in query_names.iter() {
//            match self.atomnames.at_position(name) {
//                Ok(a) => indices.push(a),
//                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
//            }
//        }
//
//        // Call cp_from_indices
//        self.as_from_indices(indices)
//    }
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

//    // Cremer Pople for indices
//    fn cp_from_indices(&self, indices: Vec<usize>) -> MemberedRing {
//        
//        let mut molarray = vec![];
//
//        for idx in indices {
//            molarray.push(self.coordinates[idx])
//        };
//
//        cremerpople::cremer_pople(&mut molarray)
//    }
//
//    // Altona Sundaralingam for indices
//    fn as_from_indices(&self, indices: Vec<usize>) -> AS {
//        
//        let mut molarray = vec![];
//
//        for idx in indices {
//            molarray.push(self.coordinates[idx])
//        };
//
//        altonasund::altona_sundaralingam(&mut molarray)
//    }
    
}
