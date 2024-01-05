use std::f64::consts::PI;

use crate::conf_sampling::sixring::TWOPI;

use crate::formalism::{CP5, CP6, MemberedRing};

use crate::geometry::{normalise_vector, cross_product, dot_product};

// The Cremer-Pople algorithm; the main function
pub fn cremer_pople(molarray: &mut Vec<[f64; 3]>) -> MemberedRing {
    
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
    const PIS_IN_180: f64 = 57.2957795130823208767981548141051703_f64;              // taken from
                                                                                    // <f64>.to_degrees()
                                                                                    //     for i in range(NUM):

    // We are not using multiplying by sqrt_cst value (sqrt(2/N)), because the factor cancels out when
    // calculating the phase_angle -> saves an operation here and there ...
    let sum1 = zj.iter().zip(cos_uv2.iter()).fold(0., |acc, (x, c)| acc + (x * c)); // q_2 * cos(phi_2) = sqrt_cst * sum1 (Eq. 12)
    let sum2 = zj.iter().zip(sin_uv2.iter()).fold(0., |acc, (y, s)| acc - (y * s)); // q_2 * sin(phi_2) = sqrt_cst * sum2 (Eq. 13)

    // By summing all zj^2 values and sqrting the result
    let amplitude = zj.iter().map(|i| i * i).sum::<f64>().sqrt();

    // (sum2/sum1) = sin(phase_angle) / cos(phase_angle) -> atan2(sum2/sum1) = phase_angle
    let mut phase_angle = sum2.atan2(sum1); 

    // Some mirroring and subtractions are needed to make everything come out right
    if sum1 <= 0.0 { phase_angle = PI - phase_angle }; 
    if sum1 < 0.0 { phase_angle = TWOPI - phase_angle }; 
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
