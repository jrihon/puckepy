use ndarray::Array1;
use pyo3::{pyclass, pymethods};


/// the `phi-psi` dihedrals, which are the peptide backbone dihedrals in proteins
/// public `phi` field : Vec<f32>
/// public `psi` field : Vec<f32>
#[pyclass(get_all)]
pub struct Peptide {
    phi : Vec<f32>,
    psi : Vec<f32>,
}

#[pymethods]
impl Peptide {

    #[new]
    fn new(num: u16) -> Peptide {

        let amount = (num * num) as usize;

        let axes = PeptideAxes::new(num as usize);

        let mut phi = Vec::with_capacity(amount);
        let mut psi = Vec::with_capacity(amount);
        
        let mut xi : f32;
        let mut yi : f32;
        for i in 0..amount as usize {

            // For every x value, return all y values
            xi = (i as f32 / num as f32).floor(); // floor, to return x axis value
            yi = i as f32 % num as f32; // return with modulo, to return y axis value

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
//
/// The axes to iterate over for peptide-like molecules : 
/// Its extent is : [0 , 2pi] (rad)
/// Its extent is : [0 , 360] (degrees)
/// public `x` field : Vec<f32>
/// public `y` field : Vec<f32>
/// Can remain a private struct, as this only is required to build the Peptide struct
//#[pyclass(get_all)]
pub struct PeptideAxes {
    x : Vec<f32>,
    y : Vec<f32>,

}

impl PeptideAxes {
    /// Initialise the struct with an array of zeroes
    fn new(num: usize) -> PeptideAxes {
        PeptideAxes {
            x: Array1::linspace(0., 360., num).into_raw_vec(),
            y: Array1::linspace(0., 360., num).into_raw_vec(),
        }
    }
}

