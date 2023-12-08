use pyo3::pymethods;

use crate::molfile::Xyz;

use self::cremerpople::{
    geometric_center_of_molecule,
    molecular_axis,
    local_elevation,
    return_cp_coordinates
};

mod cremerpople;



#[pymethods]
impl Xyz {

    fn return_crempop(&self, indices: Vec<usize>) -> (f32, f32) {
        
        let mut molarray = vec![];

        for idx in indices {
            molarray.push(self.coordinates[idx])
        };

        geometric_center_of_molecule(&mut molarray);
        let mol_axis = molecular_axis(&molarray);
        let zj = local_elevation(&molarray, mol_axis);
        let (amp, phase_angle) = return_cp_coordinates(zj);

        (amp, phase_angle)
    }
    
}
