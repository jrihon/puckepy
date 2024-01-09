use std::f64::consts::PI;
use pyo3::{pymethods, pyclass};

use crate::conf_sampling::sixring::TWOPI;
use crate::geometry::molecule_ops::dihedral;
use crate::formalism::{
    moleculefile::Pdb,
    PIS_IN_180,
    search_atomname::FindString
};

const PIOVERFIVE: f64 = 0.628318530718;

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
    
    
    // Find the indices of the atomnames and pass them to self.as_from_indices()
    fn from_atomnames(&self, pdb: &Pdb, query_names: Vec<String>) -> (f64, f64) {

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
        self.from_indices(pdb.coordinates.clone(), indices)
    }

    // Calculate Altona Sundaralingam formalism by the indices
    fn from_indices(&self, coord_array: Vec<[f64;3]>, indices: Vec<usize>) -> (f64, f64) {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        for idx in indices {
            molarray.push(coord_array[idx])
        }

        altona_sundaralingam(&mut molarray)

//       match cremerpople::cremer_pople(&mut molarray) {
//           MemberedRing::Five(a) => a.to_as(),
//           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
//       }
    }

//    fn to_cp(&self) -> f64 {
//        // let mut phase_angle = self.1 + 90.;
//        // if phase_angle > 360. { 
//        //     phase_angle -= 360.
//        // }; => Original code
//
//        // If the value is larger than 360 after adding 90, it is already larger than 270
//        // This means that we will do two operations, a +90 and then -360
//        // This cuts out an instruction or two down the line
//        if self.phase_angle > 270. { self.phase_angle - 270. } else { self.phase_angle + 90. }
//
//    }
}





/// (b) For abbreviated nomenclature see M. Sundaralingam, J. A". Chem. SOC.,93, 6644 (1971). and references therein.

// tan(P) = (theta2 + theta4 - theta1 - theta3 ) / (2 * theta0 * (sin36° + sin72°))
// theta_M = theta0 / cos(P)
// altona sundaralingam
// Where P => phase angle (expressed in degrees);  tau_M => amplitude of pucker (expressed in degrees)
// Who the fuck thought that amplitude should be expressed in degrees and not in radians?
// Note that : 
// theta2 = nu0
// theta3 = nu1
// theta4 = nu2
// theta0 = nu3
// theta1 = nu4
// Here, we will assume that we start from O4' -> C1' -> C2' -> C3' -> C4', like CP
//     While AS assumes C2' -> C3' -> C4' -> O4' -> C1'  
// 
// Function courtesy of Cpptraj Github : https://github.com/Amber-MD/cpptraj/blob/master/src/TorsionRoutines.cpp
fn altona_sundaralingam(coordinates: &Vec<[f64;3]>) -> (f64, f64) {
    
    let theta4 = dihedral(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
    let theta0 = dihedral(coordinates[1], coordinates[2], coordinates[3], coordinates[4]);
    let theta1 = dihedral(coordinates[2], coordinates[3], coordinates[4], coordinates[0]);
    let theta2 = dihedral(coordinates[3], coordinates[4], coordinates[0], coordinates[1]);
    let theta3 = dihedral(coordinates[4], coordinates[0], coordinates[1], coordinates[2]);
//    println!("{}, {}, {}, {}, {}", theta0, theta1, theta2, theta3, theta4);

    let a = (theta0 + // theta0 + cos(0.) == theta0
             (theta1 * (4. * PIOVERFIVE).cos()) +
             (theta2 * (8. * PIOVERFIVE).cos()) +
             (theta3 * (12. * PIOVERFIVE).cos()) +
             (theta4 * (16. * PIOVERFIVE).cos())) * 0.4;

    let b = ((theta0 + 1.) + // theta0 + sin(0.) == theta0 + 1
             (theta1 * (4. * PIOVERFIVE).sin()) +
             (theta2 * (8. * PIOVERFIVE).sin()) +
             (theta3 * (12. * PIOVERFIVE).sin()) +
             (theta4 * (16. * PIOVERFIVE).sin())) * -0.4;

    let amplitude = ((a * a) + (b * b)).sqrt().to_radians() * 0.65_f64 ;
    let mut phase_angle = ((theta2 + theta4) - (theta1 + theta3)).atan2(2. * theta0 * ((36_f64.to_radians().sin() + 72_f64.to_radians().sin())));

    // if the amplitude is roughly equal to 0.0 , then that means that the conformer has all 
    // five atoms in the same plane. This makes the phase_angle undefined and therefore we 
    // equate it to 0.0
//    let mut phase_angle = 0.;
//    if amplitude > 0.0 || amplitude < 0.0 {
//        phase_angle = b.atan2(a)
//    }

    phase_angle += PI; 
//
    if phase_angle < 0.0 {
        phase_angle += TWOPI; 
    }

    phase_angle *= PIS_IN_180; // to_degrees()

    (amplitude, phase_angle)
}
