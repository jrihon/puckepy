use std::f32::consts::PI;
use crate::conf_sampling::sixring::TWOPI;

use crate::CP;

use crate::geometry::{normalise_vector, cross_product, dot_product};

// The Cremer-Pople algorithm
pub fn cremer_pople(molarray: &mut Vec<[f32; 3]>) -> CP {
    
    geometric_center_of_molecule(molarray);
    let mol_axis = molecular_axis(&molarray);
    let zj = local_elevation(&molarray, mol_axis);
    let (amplitude, phase_angle) = return_cp_coordinates(zj);

    CP::new(amplitude, phase_angle)
}

// Copied the coordinates over to make the array a mutable reference
fn geometric_center_of_molecule(molarray : &mut Vec<[f32;3]>) { 

    let (x, y, z) = calculate_average_per_dimension(&molarray);

    // move molecule array to geometric center
    for coord in molarray { // already a mutable reference.
                            // molarray.iter_mut().for_each() was being annoying
            coord[0] -= x;
            coord[1] -= y;
            coord[2] -= z
    };

}

fn molecular_axis(molarray : &Vec<[f32;3]>) -> [f32;3] { 

    let (cos_uv, sin_uv) = unit_vector(molarray.len() as f32);

    // Calculate R prime
    let rp = molarray.iter().zip(cos_uv.iter()).map({ |(arr, ci)|
                    arr.map(|x| x * ci)
                }).collect::<Vec<[f32;3]>>();

    // Calculate R prime prime
    let rpp = molarray.iter().zip(sin_uv.iter()).map({ |(arr, si)|
                    arr.map(|x| x * si)
                }).collect::<Vec<[f32;3]>>();

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
fn local_elevation(molarray : &Vec<[f32;3]>, mol_axis: [f32;3]) -> Vec<f32> {

    // iterate over the array and get the local elevation for every coordinate
    molarray.iter().map(|coord| dot_product(*coord, mol_axis) ).collect()
}

// Calculate the Cremer Pople Coordinates based on the local elevation
fn return_cp_coordinates(zj : Vec<f32>) -> (f32, f32) { 

    // constant values for the calculations 
    let size = zj.len() as f32;
    let cos_uv2 = [0., 1., 2., 3., 4.].map(|i| ((4. * PI * i) / size ).cos() );     // cos(2pi * m * i / 5) (Eq. 12)
    let sin_uv2 = [0., 1., 2., 3., 4.].map(|i| ((4. * PI * i) / size ).sin() );     // sin(2pi * m * i / 5) (Eq. 13)
    const PIS_IN_180: f32 = 57.2957795130823208767981548141051703_f32;              // taken from
                                                                                    // <f32>.to_degrees()
                                                                                    //     for i in range(NUM):

    let sqrt_cst = (2_f32/5_f32).sqrt();                                            // sqrt(2/5)

    let sum1 = zj.iter().zip(cos_uv2.iter()).fold(0., |acc, (x, c)| acc + (x * c)); // q_2 * cos(phi_2) = sqrt_cst * sum1 (Eq. 12)
    let sum2 = zj.iter().zip(sin_uv2.iter()).fold(0., |acc, (y, s)| acc - (y * s)); // q_2 * sin(phi_2) = sqrt_cst * sum2 (Eq. 13)

    // The norm here is the same as the amplitude, if we calculate it this way
    // By summing all zj^2 values and sqrting the result
    let amplitude = zj.iter().map(|i| i * i).sum::<f32>().sqrt();

    
    let mut phase_angle = ((sum2 * sqrt_cst) / amplitude).asin() + PI;
    
    if sum1 < 0.0 {
        phase_angle = PI - phase_angle
    }; 

    if phase_angle < 0.0 {
        phase_angle += TWOPI // make phase_angle in range 0. -> 2PI
    };

    phase_angle *= PIS_IN_180; // <f32>.to_degrees() takes a self, not &mut self

    (amplitude, phase_angle)
}


// Returns (cosined array, sined array)
fn unit_vector(size: f32) -> ([f32;5], [f32;5]) {
//fn unit_vector(size: f32) -> ([f32;6], [f32;6]) {

        // Fiverings
        let cos_uv = [0., 1., 2., 3., 4.].map(|x| ((2. * PI * x) / size ).cos() );
        let sin_uv = [0., 1., 2., 3., 4.].map(|x| ((2. * PI * x) / size ).sin() );

        (cos_uv, sin_uv)

//        // Sixrings
//        let cos_uv = [0., 1., 2., 3., 4., 5.].map(|x| ((2. * PI * x) / size ).cos() );
//        let sin_uv = [0., 1., 2., 3., 4., 5.].map(|x| ((2. * PI * x) / size ).sin() );
//
//        (cos_uv, sin_uv)



}

fn calculate_average_per_dimension(molarray: &Vec<[f32;3]>) -> (f32, f32, f32) {

    let size = molarray.len() as f32;

    let mut xvec: Vec<f32> = Vec::with_capacity(6);
    let mut yvec: Vec<f32> = Vec::with_capacity(6);
    let mut zvec: Vec<f32> = Vec::with_capacity(6);

    for i in molarray.iter() {
        xvec.push(i[0]); yvec.push(i[1]); zvec.push(i[2]);
    }

    // Calculate averages of coordinate to define geometric center
    let x = xvec.iter().fold(0., |acx, xi| acx + xi) / size;
    let y = yvec.iter().fold(0., |acy, yi| acy + yi) / size;
    let z = zvec.iter().fold(0., |acz, zi| acz + zi) / size;

    (x, y, z)


}
