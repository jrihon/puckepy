#![allow(unused)]
use crate::conf_sampling::sixring::TWOPI;
use crate::geometry::dihedral;
use crate::formalism::AS;
use std::f64::consts::PI;


use crate::formalism::PIS_IN_180;


const PIOVERFIVE: f64 = 0.628318530718;


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
pub fn altona_sundaralingam(coordinates: &Vec<[f64;3]>) -> AS {
    
//    let theta1 = dihedral(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
//    let theta2 = dihedral(coordinates[1], coordinates[2], coordinates[3], coordinates[4]);
//    let theta3 = dihedral(coordinates[2], coordinates[3], coordinates[4], coordinates[0]);
//    let theta4 = dihedral(coordinates[3], coordinates[4], coordinates[0], coordinates[1]);
//    let theta0 = dihedral(coordinates[4], coordinates[0], coordinates[1], coordinates[2]);
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
                               //
//  amp = sqrt(a*a + b*b); // => this is how you calculate tau_m
//
//  if (amp != 0.0)
//    pucker = atan2(b,a);
//  if (pucker < 0) pucker += Constants::TWOPI;
//

    AS::new(amplitude, phase_angle)
}
