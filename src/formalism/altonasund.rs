use crate::conf_sampling::sixring::TWOPI;
use crate::geometry::dihedral;
use crate::formalism::AS;



// Redo these, because it is not making sense. Sin(72°) is most definitely not 0.95
const PIOVERFIVE: f32 = 0.628318530718;


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
pub fn altona_sundaralingam(coordinates: &Vec<[f32;3]>) -> AS {
    
    let theta2 = dihedral(coordinates[0], coordinates[1], coordinates[2], coordinates[3]);
    let theta3 = dihedral(coordinates[1], coordinates[2], coordinates[3], coordinates[4]);
    let theta4 = dihedral(coordinates[2], coordinates[3], coordinates[4], coordinates[0]);
    let theta0 = dihedral(coordinates[3], coordinates[4], coordinates[0], coordinates[1]);
    let theta1 = dihedral(coordinates[4], coordinates[0], coordinates[1], coordinates[2]);

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

    let amplitude = (a.powi(2) + b.powi(2)).sqrt();
    let mut phase_angle = 0.0;

    if amplitude > 0.0 || amplitude < 0.0 {
        phase_angle = b.atan2(a)
    }

    if phase_angle < 0.0 {
        phase_angle += TWOPI; 
    }
//  amp = sqrt(a*a + b*b); // => this is how you calculate tau_m
//
//  if (amp != 0.0)
//    pucker = atan2(b,a);
//  if (pucker < 0) pucker += Constants::TWOPI;
//

    AS::new(amplitude, phase_angle)
}
