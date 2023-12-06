use std::f32::consts::PI;

use ndarray::Array1;
use pyo3::{pyclass, pymethods};

const FOURPIOVERFIVE : f32 = (4. * PI) / 5.;

#[pyclass(get_all)]
pub struct Fivering {
    nu1: Vec<f32>,
    nu3: Vec<f32>,
}

#[pymethods]
impl Fivering {
    
    #[new]
    fn new(num: u16) -> Self {
        
        // Derive torsion angles from the given axes
        let polars = FuranoseAxes::new(num as usize);

        // Setup variable
        let amount : u16 = num * num;
        let num_f32 : f32 = num as f32;

        // Initialise equation-specific constants
        let denominator_x : f32 = FOURPIOVERFIVE.cos();
        let denominator_y : f32 = FOURPIOVERFIVE.sin();

        // Instance Furanose struct
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
//            f.nu1[i] = ( polars.zx[x as usize] * denominator_x ) + ( polars.zy[y as usize] * denominator_y);
//            f.nu3[i] = ( polars.zx[x as usize] * denominator_x ) - ( polars.zy[y as usize] * denominator_y);
        }

        // Make values ORCA-ready
        Self {
            nu1 : nu1.iter().map(|x| if x < &0. { x + 360.} else {*x}).collect(),
            nu3 : nu3.iter().map(|x| if x < &0. { x + 360.} else {*x}).collect()
        }

    }
}

struct FuranoseAxes {
    zx : Array1<f32>,
    zy : Array1<f32>,

}

impl FuranoseAxes {
    /// Initialise the struct with a near-empty array
    fn new(num: usize) -> FuranoseAxes {
        FuranoseAxes {
            zx: Array1::linspace(-60., 60., num),
            zy: Array1::linspace(-60., 60., num),
        }
        
    }
    
}
