use ndarray::Array1;
use pyo3::{pyclass, pymethods};


/// the `phi-psi` dihedrals, which are the peptide backbone dihedrals in proteins
/// public `phi` field : Vec<f64>
/// public `psi` field : Vec<f64>
#[pyclass(get_all)]
pub struct Peptide {
    phi : Vec<f64>,
    psi : Vec<f64>,
}

#[pymethods]
impl Peptide {

    #[new]
    fn new(interval: u16) -> Self {

        let amount = (interval * interval) as usize;

        let axes = PAxes::new(interval as usize);

        let mut phi = Vec::with_capacity(amount);
        let mut psi = Vec::with_capacity(amount);
        
        let mut xi : f64;
        let mut yi : f64;
        for i in 0..amount as usize {

            // For every x value, return all y values
            xi = (i as f64 / interval as f64).floor(); // floor, to return x axis value
            yi = i as f64 % interval as f64; // return with modulo, to return y axis value

            // fill out the array
            phi.push(axes.x[xi as usize]); 
            psi.push(axes.y[yi as usize]); 
        }

        Self {
            phi,
            psi,
        }
    }
}




#[pyclass(get_all)]
pub struct PeptideAxes {
    x : Vec<f64>,
    y : Vec<f64>,
}

#[pymethods]
impl PeptideAxes {

    #[new]
    fn new(interval: u16) -> Self {
        
        let amount = (interval * interval) as usize;

        let axes = PAxes::new(interval as usize);

        let mut phi = Vec::with_capacity(amount);
        let mut psi = Vec::with_capacity(amount);
        
        let mut xi : f64;
        let mut yi : f64;
        for i in 0..amount as usize {

            // For every x value, return all y values
            xi = (i as f64 / interval as f64).floor(); // floor, to return x axis value
            yi = i as f64 % interval as f64; // return with modulo, to return y axis value

            // fill out the array
            phi.push(axes.x[xi as usize]); 
            psi.push(axes.y[yi as usize]); 
        }

        Self {
            x: phi,
            y: psi,
        }
    }
}
//
/// The axes to iterate over for peptide-like molecules : 
/// Its extent is : [0 , 2pi] (rad)
/// Its extent is : [0 , 360] (degrees)
/// public `x` field : Vec<f64>
/// public `y` field : Vec<f64>
/// Can remain a private struct, as this only is required to build the Peptide struct
struct PAxes {
    x : Vec<f64>,
    y : Vec<f64>,
}

impl PAxes {
    /// Initialise the struct with an array of zeroes
    fn new(interval: usize) -> PAxes {
        PAxes {
            x: Array1::linspace(0., 360., interval).into_raw_vec(),
            y: Array1::linspace(0., 360., interval).into_raw_vec(),
        }
    }
}

