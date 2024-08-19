use std::f64::consts::PI;
use pyo3::{pyclass, pymethods};

use crate::conf_sampling::sixring::TWOPI;
use crate::formalism::{
    moleculefile::Pdb,
    PIS_IN_180,
    search_atomname::FindString
};

use crate::geometry::fundamental_ops::{normalise_vector, cross_product, dot_product};
use crate::formalism::inversion;


// Enum to control the which type of n-membered ring system is being produced and 
// returns the correct one to the user.
// Acts as an addition safety measure whenever users prompt incorrect amount of values in function
// calls too
// Only for internal usage in the backend
enum MemberedRing {
    Five(CP5),
    Six(CP6)
}

/// The CP tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct CP5 {
    amplitude: f64,
    phase_angle: f64,
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

    // Find indices of atomnames and pass them to self.cp_from_indices()
    fn from_atomnames(&self, pdb : &Pdb, query_names: Vec<String>) -> (f64, f64) {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        let _ = match query_names.len() {
            5 => 5,
           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
        };

        // Search for the indices of the atom names
        for name in query_names.iter() {
            match pdb.atomnames.at_position(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
            }
        }

        self.from_indices(pdb.coordinates.clone(), indices) // I have to make this a clone()
                                                               // because it does not satisfy the
                                                               // PyObject trait bound for a reason
                                                               // unbeknownst to me
    }

    // Calculate Cremer-Pople formalism by prompted indices
    fn from_indices(&self, coordinates : Vec<[f64; 3]>, indices: Vec<usize>) -> (f64, f64) {
        
        let mut molarray: Vec<[f64; 3]>  = indices.iter().map(|i| coordinates[*i]).collect();

        match cremer_pople(&mut molarray) {
           MemberedRing::Five(cp) => (cp.amplitude, cp.phase_angle),
           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
        }
    }
    
    fn invert(&self) -> [[f64;3]; 5] {
        inversion::fivering::invert_fivering(self.amplitude, self.phase_angle)
    }

}


#[pyclass]
pub struct CP6 {
    amplitude: f64,
    phase_angle: f64,
    theta: f64,
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
    fn from_indices(&self, coordinates : Vec<[f64; 3]>, indices: Vec<usize>) -> (f64, f64, f64) {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        for idx in indices {
            molarray.push(coordinates[idx])
        }

       match cremer_pople(&mut molarray) {
           MemberedRing::Six(cp) => (cp.amplitude, cp.phase_angle, cp.theta),
           _ => panic!("An amount, not equal to 6, has been queried. Expected 6 elements.")
       }
    }
    
    // Find indices of atomnames and pass them to self.cp_from_indices()
    fn from_atomnames(&self, pdb : &Pdb, query_names: Vec<String>) -> (f64, f64, f64) {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        let _ = match query_names.len() {
            6 => 6,
           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
        };


        // Search for the indices of the atom names
        for name in query_names.iter() {
            match pdb.atomnames.at_position(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
            }
        }

        self.from_indices(pdb.coordinates.clone(), indices)
    }

    fn invert(&self) -> [[f64;3]; 6] {
        inversion::sixring::invert_sixring(self.amplitude, self.phase_angle, self.theta)
    }

}



// The Cremer-Pople algorithm; the main function
fn cremer_pople(molarray: &mut Vec<[f64; 3]>) -> MemberedRing {
    
    geometric_center_of_molecule(molarray);
    let mol_axis = molecular_axis(&molarray);
    let zj = local_elevation(&molarray, mol_axis);
    
    return_cp_coordinates(zj)

}

// Copied the coordinates over to make the array a mutable reference
// works for all any-membered ring systems
fn geometric_center_of_molecule(molarray : &mut Vec<[f64;3]>) { 

    let (x, y, z) = calculate_average_per_dimension(&molarray);

    // move molecule array to geometric center
    for coord in molarray { // already a mutable reference.
                            // molarray.iter_mut().for_each() was being annoying
            coord[0] -= x;
            coord[1] -= y;
            coord[2] -= z  // doing it like this is more readable than mapping
    };

}

// works for all any-membered ring systems
fn molecular_axis(molarray : &Vec<[f64;3]>) -> [f64;3] { 

    let (cos_uv, sin_uv) = unit_vector(molarray.len());

    // Calculate R prime
    let rp = molarray.iter().zip(cos_uv.iter()).map({ |(arr, ci)|
                    arr.map(|x| x * ci)
                }).collect::<Vec<[f64;3]>>();

    // Calculate R prime prime
    let rpp = molarray.iter().zip(sin_uv.iter()).map({ |(arr, si)|
                    arr.map(|x| x * si)
                }).collect::<Vec<[f64;3]>>();

    let (x0, y0, z0) = calculate_average_per_dimension(&rp);
    let (x1, y1, z1) = calculate_average_per_dimension(&rpp);

    // return molecular axis
    cross_product(
        normalise_vector([x0, y0, z0]),
        normalise_vector([x1, y1, z1])
    )    

}

// Calculate local elevation by taking the dot product of the 
// centered molecule's array and doing a dot(a, b) every 
// coordinates and the molecular_axis
// works for all any-membered ring systems
fn local_elevation(molarray : &Vec<[f64;3]>, mol_axis: [f64;3]) -> Vec<f64> {

    // iterate over the array and get the local elevation for every coordinate
    molarray.iter()
        .map(|coord| dot_product(*coord, mol_axis) )
        .collect()
}

// Calculate the Cremer Pople Coordinates based on the local elevation
fn return_cp_coordinates(zj : Vec<f64>) -> MemberedRing { 

    // constant values for the calculations 
    let size = zj.len();
    let cos_uv2: Vec<f64> = (0..size).map(|i| ((4. * PI * i as f64) / size as f64).cos() ).collect();     // cos(2pi * m * i / 5) (Eq. 12)
    let sin_uv2: Vec<f64> = (0..size).map(|i| ((4. * PI * i as f64) / size as f64).sin() ).collect();     // sin(2pi * m * i / 5) (Eq. 12)

    // We are not using multiplying by sqrt_cst value (sqrt(2/N)), because the factor cancels out when
    // calculating the phase_angle -> saves an operation here and there ...
    let sum1 = zj.iter().zip(cos_uv2.iter()).fold(0., |acc, (x, c)| acc + (x * c)); // q_2 * cos(phi_2) = sqrt_cst * sum1 (Eq. 12)
    let sum2 = zj.iter().zip(sin_uv2.iter()).fold(0., |acc, (y, s)| acc - (y * s)); // q_2 * sin(phi_2) = sqrt_cst * sum2 (Eq. 13)

    // By summing all zj^2 values and sqrting the result
    let amplitude = zj.iter().map(|i| i * i).sum::<f64>().sqrt();

    // (sum2/sum1) = sin(phase_angle) / cos(phase_angle) -> atan2(sum2/sum1) = phase_angle
    let mut phase_angle = sum2.atan2(sum1); 

    // Some mirroring and subtractions are needed to make everything come out right
//    if sum1 <= 0.0 { phase_angle = PI - phase_angle }; 
//    if sum1 < 0.0 { phase_angle = TWOPI - phase_angle }; 
    if sum1 <= 0.0 { phase_angle = PI + phase_angle }; 
    if sum1 > 0.0 { phase_angle -= PI }

    if phase_angle < 0.0 { phase_angle += TWOPI }; // radians range
    phase_angle *= PIS_IN_180; // <f64>.to_degrees() takes a self, not &mut self

    match size {
        5 => {
            MemberedRing::Five(CP5::new(amplitude, phase_angle))
        },
        6 => {
            let q3: f64 = zj.iter().zip([1., -1., 1., -1., 1., -1.])
                                    .map( |(z, factor)| z * factor).sum::<f64>() / (size as f64)
                                    .sqrt();

            // For some reason, it is necessary to mirror the value over PI
            let theta = (PI - (q3/amplitude).acos()) * PIS_IN_180; // acos -> to_degrees()

            MemberedRing::Six(CP6::new(amplitude, phase_angle, theta))
        },
        _ => panic!("Ringsystem prompted is not FIVE-membered or SIX-membered.")
    }
}


// Returns (cosined array, sined array)
// works for up to six-membered ring systems
fn unit_vector(size: usize) -> (Vec<f64>, Vec<f64>) {

        let cos_uv = (0..size).map(|x| ((2. * PI * x as f64) / size as f64 ).cos() ).collect();
        let sin_uv = (0..size).map(|x| ((2. * PI * x as f64) / size as f64 ).sin() ).collect();

        (cos_uv, sin_uv)
}

// works for up to six-membered ring systems
fn calculate_average_per_dimension(molarray: &Vec<[f64;3]>) -> (f64, f64, f64) {

    let size = molarray.len() as f64;

    // with_capacity allows pushing onto the heap for six elements without reallocation
    let mut xvec: Vec<f64> = Vec::with_capacity(6);
    let mut yvec: Vec<f64> = Vec::with_capacity(6);
    let mut zvec: Vec<f64> = Vec::with_capacity(6);

    for i in molarray.iter() {
        xvec.push(i[0]); yvec.push(i[1]); zvec.push(i[2]);
    }

    // Calculate averages of coordinate to define geometric center
    let x = xvec.iter().fold(0., |acx, xi| acx + xi) / size;
    let y = yvec.iter().fold(0., |acy, yi| acy + yi) / size;
    let z = zvec.iter().fold(0., |acz, zi| acz + zi) / size;

    (x, y, z)


}
