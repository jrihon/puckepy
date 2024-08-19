use pyo3::{pymethods, pyclass};

use crate::geometry::molecule_ops::dihedral;
use crate::formalism::{
    moleculefile::Pdb,
    search_atomname::FindString,
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

        // Call cp_from_indices
        self.from_indices(pdb.coordinates.clone(), indices)
    }

    // Calculate Altona Sundaralingam formalism by the indices
    fn from_indices(&self, coordinates: Vec<[f64;3]>, indices: Vec<usize>) -> (f64, f64) {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        let _ = match indices.len() {
            5 => 5,
           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
        };

        for idx in indices {
            molarray.push(coordinates[idx])
        }

        altona_sundaralingam(&mut molarray)

//       match cremerpople::cremer_pople(&mut molarray) {
//           MemberedRing::Five(a) => a.to_as(),
//           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
//       }
    }

}





/// (b) For abbreviated nomenclature see M. Sundaralingam, J. A". Chem. SOC.,93, 6644 (1971). and references therein.

// tan(P) = (theta2 + theta4 - theta1 - theta3 ) / (2 * theta0 * (sin36° + sin72°))
// tan(P) = (nu4 + nu1 - nu3 - nu0 ) / (2 * nu2 * (sin36° + sin72°))
// theta_M = theta0 / cos(P)
// altona sundaralingam
// The AS formalism is both the dumbest and the mathmetically most inconsistent formula.
//  The torsion angles are in degrees and [-180 -> 180] ... who calculates with that????
//
// Also, who the frick thought that amplitude should be expressed in degrees and not in radians?
//
// Note that : 
// theta2 = nu0
// theta3 = nu1
// theta4 = nu2
// theta0 = nu3
// theta1 = nu4
// Here, we will assume that we start from O4' -> C1' -> C2' -> C3' -> C4', like CP
//     Instead of AS's assumption of C2' -> C3' -> C4' -> O4' -> C1'  
// 
// Function courtesy of Cpptraj Github : https://github.com/Amber-MD/cpptraj/blob/master/src/TorsionRoutines.cpp
fn altona_sundaralingam(coordinates: &Vec<[f64;3]>) -> (f64, f64) {
    
    //  we follow the order of O4' - C1' - C2' - C3' - C4' when the atoms are being passed to the function
    //  NB: cpptraj follows  C1' - C2' - C3' - C4' - O4' when the atoms are being passed to the function
    let nu0 = dihedral(coordinates[4], coordinates[0], coordinates[1], coordinates[2]); // nu0 -> theta3 -> v4
    let nu1 = dihedral(coordinates[0], coordinates[1], coordinates[2], coordinates[3]); // nu1 -> theta4 -> v5
    let nu2 = dihedral(coordinates[1], coordinates[2], coordinates[3], coordinates[4]); // nu2 -> theta0 -> v1
    let nu3 = dihedral(coordinates[2], coordinates[3], coordinates[4], coordinates[0]); // nu3 -> theta1 -> v2
    let nu4 = dihedral(coordinates[3], coordinates[4], coordinates[0], coordinates[1]); // nu4 -> theta2 -> v3

    // this part is courtesy of the AS_Pucker() function in https://github.com/Amber-MD/cpptraj/blob/master/src/TorsionRoutines.cpp 
    let a = ((nu2 + 1. ) +// nu0 + cos(1.) == nu0 + 1
             (nu3 * (4. * PIOVERFIVE).cos()) +
             (nu4 * (8. * PIOVERFIVE).cos()) +
             (nu0 * (12. * PIOVERFIVE).cos()) +
             (nu1 * (16. * PIOVERFIVE).cos())) * 0.4;

    let b =  (nu2 + // nu0 + sin(0.) == nu0 + 0
             (nu3 * (4. * PIOVERFIVE).sin()) +
             (nu4 * (8. * PIOVERFIVE).sin()) +
             (nu0 * (12. * PIOVERFIVE).sin()) +
             (nu1 * (16. * PIOVERFIVE).sin())) * -0.4;

    let amplitude = ((a * a) + (b * b)).sqrt().to_radians() ;

    // swap thetas for nus
    let mut phase_angle = (nu4 + nu1 - nu3 - nu0 ).atan2(2. * nu2 * (PIOVERFIVE.sin() + (PIOVERFIVE * 2.).sin())).to_degrees();

    if phase_angle < 0. { phase_angle += 360. };

    (amplitude, phase_angle)
}
