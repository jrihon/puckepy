use pyo3::{pyclass,pymethods};

use crate::geometry::molecule_ops::{dihedral, bondangle};
use crate::formalism::{
    moleculefile::Pdb,
    search_atomname::FindString
};

#[pyclass(get_all)]
pub struct SP { }

#[pymethods]
impl SP {

    #[new]
    fn new() -> SP {
        SP { }

    }

    // Calculate Cremer-Pople formalism by prompted indices
    fn from_indices(&self, coordinates : Vec<[f64; 3]>, indices: Vec<usize>) -> ([f64;3], [f64;3]) {
        
        let mut molarray: Vec<[f64; 3]> = vec![];

        for idx in indices {
            molarray.push(coordinates[idx])
        }

        strauss_pickett(molarray)
    }
    
    // Find indices of atomnames and pass them to self.cp_from_indices()
    fn from_atomnames(&self, pdb : &Pdb, query_names: Vec<String>) -> ([f64;3], [f64;3])  {

        // Make empty vec :
        let mut indices: Vec<usize> = Vec::with_capacity(6);

        let _ = match query_names.len() {
            6 => 6,
           _ => panic!("An amount, not equal to 5, has been queried. Expected 5 elements.")
        };


        // Search for the indices of the atom names
        for name in query_names.iter() {
            match pdb.atomnames.at_position(name) {
                Ok(a) => indices.push(a),
                Err(()) => panic!("Could not find \"{}\" atomname in the queried pdb.", name)
            }
        }

        self.from_indices(pdb.coordinates.clone(), indices)
    }

}







fn strauss_pickett(molarray: Vec<[f64;3]>) -> ([f64;3], [f64;3]) {


    ([
        dihedral(molarray[4], molarray[0], molarray[2], molarray[1]), //Alpha_1
        dihedral(molarray[0], molarray[2], molarray[4], molarray[3]), //Alpha_2
        dihedral(molarray[2], molarray[4], molarray[0], molarray[5]), //Alpha_3
        ],
        [
        bondangle(molarray[0], molarray[1], molarray[2]), //Beta_1
        bondangle(molarray[2], molarray[3], molarray[4]), //Beta_2
        bondangle(molarray[4], molarray[5], molarray[0]), //Beta_3
    ])
}
