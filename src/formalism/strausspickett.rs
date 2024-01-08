
use crate::geometry::{dihedral, angle};


pub fn strauss_pickett(molarray: Vec<[f64;3]>) -> ([f64;3], [f64;3]) {


    ([
        dihedral(molarray[4], molarray[0], molarray[2], molarray[1]), //Alpha_1
        dihedral(molarray[0], molarray[2], molarray[4], molarray[3]), //Alpha_2
        dihedral(molarray[2], molarray[4], molarray[0], molarray[5]), //Alpha_3
        ],
        [
        angle(molarray[0], molarray[1], molarray[2]), //Beta_1
        angle(molarray[2], molarray[3], molarray[4]), //Beta_2
        angle(molarray[4], molarray[5], molarray[0]), //Beta_3
    ])
}
