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
#[allow(unused_variables)]
fn return_cp_coordinates(zj : Vec<f64>) -> MemberedRing { 

    // constant values for the calculations 
    let size = zj.len();
    let cos_uv2: Vec<f64> = (0..size).map(|i| ((4. * PI * i as f64) / size as f64).cos() ).collect();     // cos(2pi * m * i / 5) (Eq. 12)
    let sin_uv2: Vec<f64> = (0..size).map(|i| ((4. * PI * i as f64) / size as f64).sin() ).collect();     // sin(2pi * m * i / 5) (Eq. 12)
//    let cos_uv2 = [0., 1., 2., 3., 4.].map(|i| ((4. * PI * i) / size ).cos() );     // cos(2pi * m * i / 5) (Eq. 12)
//    let sin_uv2 = [0., 1., 2., 3., 4.].map(|i| ((4. * PI * i) / size ).sin() );     // sin(2pi * m * i / 5) (Eq. 13)
    const PIS_IN_180: f64 = 57.2957795130823208767981548141051703_f64;              // taken from
                                                                                    // <f64>.to_degrees()
                                                                                    //     for i in range(NUM):

    let sqrt_cst = (2_f64/ size as f64).sqrt();                                            // sqrt(2/5) -- sqrt(2/6)

    let sum1 = zj.iter().zip(cos_uv2.iter()).fold(0., |acc, (x, c)| acc + (x * c)); // q_2 * cos(phi_2) = sqrt_cst * sum1 (Eq. 12)
    let sum2 = zj.iter().zip(sin_uv2.iter()).fold(0., |acc, (y, s)| acc - (y * s)); // q_2 * sin(phi_2) = sqrt_cst * sum2 (Eq. 13)

    // The norm here is the same as the amplitude, if we calculate it this way
    // By summing all zj^2 values and sqrting the result
    let amplitude = zj.iter().map(|i| i * i).sum::<f64>().sqrt();

    let mut phase_angle = ((sum1 * sqrt_cst) / amplitude).acos() + PI; // need to offset by ' + PI', this works for fiverings
    
//    println!("{}", sum1);
//    println!("{}", sum2);
//    if sum1 < 0.0 {
//        phase_angle = PI - phase_angle
//    }; 
//
//    if phase_angle < 0.0 {
//        phase_angle += TWOPI // make phase_angle in range 0. -> 2PI
//    };
//
//    phase_angle *= PIS_IN_180; // <f64>.to_degrees() takes a self, not &mut self

    // Nadenken en testen of de waarde van amplitude nog geupdate moet worden of niet 
    // In Cpptraj wordt de amplitude (Q) eerst bepaald adhv summation over sum1 en sum2, 
    // en niet als een sum van de zj's . 
    // Op het einde wordt amplitude, bij sixrings, nog eens gecorrigeerd met de voor N = 6 ? 
    //     Maar moeten wij dit nog doen dan?

    match size {
        5 => {
            if sum1 < 0.0 { phase_angle = PI - phase_angle }; 
            if phase_angle < 0.0 { phase_angle += TWOPI }; // radians range
            phase_angle *= PIS_IN_180; // <f64>.to_degrees() takes a self, not &mut self

            MemberedRing::Five(CP5::new(amplitude, phase_angle))
        },
        6 => {
            let q3: f64 = zj.iter().zip([1., -1., 1., -1., 1., -1.]).map( |(z, factor)| z * factor).sum::<f64>() / (size as f64).sqrt();

            if sum1 < 0.0 || sum2 < 0.0 { phase_angle = PI - phase_angle }; 
            if phase_angle < 0.0 { phase_angle += TWOPI }; // radians range
            phase_angle *= PIS_IN_180; // <f64>.to_degrees() takes a self, not &mut self
            // For some reason, it is necessary to prepend with 'PI - ' to make theta come out right
            let theta = (PI - (q3/amplitude).acos()) * PIS_IN_180; // acos -> to_degrees()

            MemberedRing::Six(CP6::new(amplitude, phase_angle, theta))
        },
        _ => panic!("Ringsystem prompted is not FIVE-membered or SIX-membered.")
    }

//    (amplitude, phase_angle)
//      if (N == 6) {
//    double q3 = 0.0;
//    double mult = 1.0;
//    for (int i = 0; i < N; i++) {
//      q3 += mult * Zn[i]; // mult ~ pow( -1.0, i )
//      mult = -mult;
//    }
//    q3 /= sqrt( dN );
//    // Calculate theta
//    theta = atan2( amplitude, q3 );

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
