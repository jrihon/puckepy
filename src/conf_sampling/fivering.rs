use pyo3::{pyclass, pymethods};
use ndarray::Array1; // not public, useful for the linspace function

const FOURPIOVERFIVE : f32 = (4. * PI) / 5.;
use std::f32::consts::PI;




/// Struct to keep the nu_1 and nu_3 generated dihedrals
#[pyclass(get_all)]
pub struct Fivering {
    pub nu1: Vec<f32>,
    pub nu3: Vec<f32>,
}

#[pymethods]
impl Fivering {
    
    // implement the __init__() function for the wrapper
    #[new]
    fn new(num: u16) -> Self {
        
        // Derive torsion angles from the given axes
        let polars = FiveringAxes::new(num as usize);

        // Setup variable
        let amount : u16 = num * num;
        let num_f32 : f32 = num as f32;

        // Initialise equation-specific constants
        let denominator_x : f32 = FOURPIOVERFIVE.cos();
        let denominator_y : f32 = FOURPIOVERFIVE.sin();

        // Instance Fivering struct
        let mut nu1: Vec<f32> = Vec::with_capacity(amount as usize);
        let mut nu3: Vec<f32> = Vec::with_capacity(amount as usize);

        let mut x : f32;
        let mut y : f32;

        for i in 0..amount as usize {
            // Calculate indexes for the array axises
            x = (i as f32 / num_f32).floor(); // X axis, returns with floor
            y = i as f32 % num_f32; // Y axis, return with modulo

            // fill out the array
            nu1.push((polars.zx[x as usize] * denominator_x ) + ( polars.zy[y as usize] * denominator_y));
            nu3.push((polars.zx[x as usize] * denominator_x ) - ( polars.zy[y as usize] * denominator_y));
        }

        // Make values ORCA-ready
        Self {
            nu1 : nu1.iter().map(|x| if x < &0. { x + 360.} else {*x}).collect(),
            nu3 : nu3.iter().map(|x| if x < &0. { x + 360.} else {*x}).collect()
        }

    }
}

#[pyclass(get_all)]
pub struct FiveringAxes {
    zx : Vec<f32>,
    zy : Vec<f32>,

}

// Do not make the methods available to the user 
// We only return the FiveringAxes to the user to allow them to use the attributes it hold.
// The module has no purpose in allowing the user to generate their own Fivering Class
impl FiveringAxes {
    /// Initialise the struct with a near-empty array
    fn new(num: usize) -> FiveringAxes {
        FiveringAxes {
            zx: Array1::linspace(-60., 60., num).into_raw_vec(),
            zy: Array1::linspace(-60., 60., num).into_raw_vec(),
        }
        
    }
    
}
