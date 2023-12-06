use ndarray::Array1;
use pyo3::{pyclass, pymethods};


/// the `phi-psi` dihedrals, which are the peptide backbone dihedrals in proteins
/// public `phi` field : Array1<f32>
/// public `psi` field : Array1<f32>
/// Needs to be a public struct is it gets returned to the user
#[pyclass]
pub struct Peptide {
    phi : Vec<f32>,
    psi : Vec<f32>,
//    phi : Array1<f32>,
//    psi : Array1<f32>,
}

#[pymethods]
impl Peptide {

    #[new]
    fn new(num: u16) -> Peptide {

        let amount = (num * num) as usize;

        let bb = PeptideAxes::new(num as usize);

        let mut phi = Vec::with_capacity(amount);
        let mut psi = Vec::with_capacity(amount);
//        let mut phi = Array1::zeros(amount);
//        let mut psi = Array1::zeros(amount);
        
        let mut x : f32;
        let mut y : f32;
        for i in 0..amount as usize {

            // For every x value, return all y values
            x = (i as f32 / num as f32).floor(); // floor, to return x axis value
            y = i as f32 % num as f32; // return with modulo, to return y axis value

            // fill out the array
            phi.push(bb.x[x as usize]); 
            psi.push(bb.y[y as usize]); 
//            phi[i as usize] = bb.x[x as usize]; 
//            psi[i as usize] = bb.y[y as usize]; 
        }

        Self {
            phi,
            psi,
        }
    }

    fn get_phi(&self) -> Vec<f32> {
        self.phi.clone()
    }

    fn get_psi(&self) -> Vec<f32> {
        self.psi.clone()
    }
}
//
/// The axes to iterate over for peptide-like molecules : 
/// Its extent is : [0 , 2pi] (rad)
/// Its extent is : [0 , 360] (degrees)
/// public `x` field : Array1<f32>
/// public `y` field : Array1<f32>
/// Can remain a private struct, as this only is required to build the Peptide struct
struct PeptideAxes {
    x : Array1<f32>,
    y : Array1<f32>,

}

impl PeptideAxes {
    /// Initialise the struct with an array of zeroes
    fn new(num: usize) -> PeptideAxes {
        PeptideAxes {
            x: Array1::linspace(0., 360., num),
            y: Array1::linspace(0., 360., num),
        }
    }
}

