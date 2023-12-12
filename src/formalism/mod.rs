use pyo3::{pymethods,pyclass};

use crate::molfile::{Xyz, Pdb};

use self::cremerpople::{
    geometric_center_of_molecule,
    molecular_axis,
    local_elevation,
    return_cp_coordinates
};

mod cremerpople;


/// The CP tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct CP {
    amplitude: f32,
    phase_angle: f32,
}

#[pymethods]
impl CP {

    #[new]
    fn new(amplitude: f32, phase_angle: f32) -> CP {
        CP { amplitude, phase_angle }
    }

    fn to_as(&self) -> AS {
        // let mut phase_angle = self.1 - 90.;
        // if phase_angle < 0. { 
        //     phase_angle += 360.
        // }; => Original code

        // If the value is smaller than 0 after decreasing 90, it is already smaller than 90
        // This means that we will do two operations, a -90 and then +360
        // This cuts out an operation or two down the line
        let new_angle = if self.phase_angle < 90. { self.phase_angle + 270. } else { self.phase_angle - 90. };

        AS { amplitude: self.amplitude, phase_angle : new_angle }

    }
    
}


/// The AS tuple-struct holds the (amplitude, phase_angle) parameters
#[pyclass(get_all)]
pub struct AS {
    amplitude: f32,
    phase_angle: f32,
}

#[pymethods]
impl AS {

    #[new]
    fn new(amplitude: f32, phase_angle: f32) -> AS {
        AS { amplitude, phase_angle }
    }
    
    fn to_cp(&self) -> CP {
        // let mut phase_angle = self.1 + 90.;
        // if phase_angle > 360. { 
        //     phase_angle -= 360.
        // }; => Original code

        // If the value is larger than 360 after adding 90, it is already larger than 270
        // This means that we will do two operations, a +90 and then -360
        // This cuts out an operation or two down the line
        let new_angle = if self.phase_angle > 270. { self.phase_angle - 270. } else { self.phase_angle + 90. };

        CP {amplitude: self.amplitude, phase_angle: new_angle }

    }
}


#[pymethods]
impl Xyz {

    fn cp_from_indices(&self, indices: Vec<usize>) -> CP {
        
        let mut molarray = vec![];

        for idx in indices {
            molarray.push(self.coordinates[idx])
        };

        geometric_center_of_molecule(&mut molarray);
        let mol_axis = molecular_axis(&molarray);
        let zj = local_elevation(&molarray, mol_axis);
        let (amplitude, phase_angle) = return_cp_coordinates(zj);

        CP { amplitude, phase_angle }
    }
    
}


#[pymethods]
impl Pdb {

    fn cp_from_indices(&self, indices: Vec<usize>) -> CP {
        
        let mut molarray = vec![];

        for idx in indices {
            molarray.push(self.coordinates[idx])
        };

        geometric_center_of_molecule(&mut molarray);
        let mol_axis = molecular_axis(&molarray);
        let zj = local_elevation(&molarray, mol_axis);
        let (amplitude, phase_angle) = return_cp_coordinates(zj);

        CP { amplitude, phase_angle }
    }
    
    fn cp_from_atomnames(&self, query_names: Vec<String>) -> CP {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        // Search for the indices of the atom names
        for name in query_names.iter() {
            match self.atom_names.at(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find {} in the queried pdb.", name)
            }
        }

        // Call cp_from_indices
        self.cp_from_indices(indices)
    }
}

trait FindString {
    fn at(&self, pattern: &str) -> Result<usize, ()> ;
}

impl FindString for Vec<String> {

    fn at(&self, pattern: &str) -> Result<usize, ()> {

        let mut c = 0;

        for name in self {
            if name != pattern {
                c += 1
            } else { 
                return Ok(c as usize) 
            }
        };

        Err(())

//        let a = self.iter().try_fold( 0_u32, |mut prev: u32, atomname| 
//            {
//                if atomname == pattern {
//                    ControlFlow::Continue(prev += 1_u32)
//                } else {
//                    ControlFlow::Break(prev)
//                }
//            });
//
//        Ok(4)
    }
    
}
