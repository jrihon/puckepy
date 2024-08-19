use pyo3::{pyclass, pymethods};
use ndarray::Array1; // not public, useful for the linspace function

const FOURPIOVERFIVE : f64 = (4. * PI) / 5.;
use std::f64::consts::PI;




/// Struct to keep the nu_1 and nu_3 generated dihedrals
#[pyclass(get_all)]
pub struct Fivering {
    pub nu1: Vec<f64>,
    pub nu3: Vec<f64>,
}

#[pymethods]
impl Fivering {
    
    #[new]
    fn new(interval: u16) -> Self {
        
        // Derive torsion angles from the given axes
        let polars = FAxes::new(interval as usize);

        // Setup variable
        let amount : u16 = interval * interval;
        let interval_f64 : f64 = interval as f64;

        // Initialise equation-specific constants
        let denominator_x : f64 = FOURPIOVERFIVE.cos();
        let denominator_y : f64 = FOURPIOVERFIVE.sin();

        // Instance Fivering struct
        let mut nu1: Vec<f64> = Vec::with_capacity(amount as usize);
        let mut nu3: Vec<f64> = Vec::with_capacity(amount as usize);

        let mut x : f64;
        let mut y : f64;

        for i in 0..amount as usize {
            // Calculate indexes for the array axises
            x = (i as f64 / interval_f64).floor(); // X axis, returns with floor
            y = i as f64 % interval_f64; // Y axis, return with modulo

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
    zx : Vec<f64>,
    zy : Vec<f64>,
}

#[pymethods]
impl FiveringAxes {

    #[new]
    fn new(interval: u16) -> Self {

        // Derive torsion angles from the given axes
        let polars = FAxes::new(interval as usize);

        // Setup variable
        let amount : u16 = interval * interval;
        let interval_f64 : f64 = interval as f64;

        // Instance Fivering struct
        let mut zx: Vec<f64> = Vec::with_capacity(amount as usize);
        let mut zy: Vec<f64> = Vec::with_capacity(amount as usize);
        let mut x : f64;
        let mut y : f64;

        for i in 0..amount as usize {
            // Calculate indexes for the array axises
            x = (i as f64 / interval_f64).floor(); // X axis, returns with floor
            y = i as f64 % interval_f64; // Y axis, return with modulo

            // fill out the array
            zx.push(polars.zx[x as usize]);
            zy.push(polars.zy[y as usize]);
        }

        // Make values ORCA-ready
        Self {
            zx,
            zy
        }
        
    }
    
}



struct FAxes {
    zx : Vec<f64>,
    zy : Vec<f64>,
}

// Do not make the methods available to the user 
// We only return the FiveringAxes to the user to allow them to use the attributes it hold.
// The module has no purpose in allowing the user to generate their own Fivering Class
impl FAxes {

    fn new(interval: usize) -> FAxes {
        FAxes {
            zx: Array1::linspace(-60., 60., interval).into_raw_vec(),
            zy: Array1::linspace(-60., 60., interval).into_raw_vec(),
        }
        
    }
    
}
