//! PUBLIC FUNCTIONS FOR THE API : DIHEDRAL, BONDANGLE, BONDLENGTH
//!
//!
//!
//! CREMER-POPLE Calculations : dot_product, cross_product, normalise_vector, subtract_arr, RotationMatrix
//!
//! General Linear Algebra and math stuff I do not want to import from different libraries
//! so I write it myself and this way I can implement them on primitives

use std::f64::consts::PI;
use pyo3::pyfunction;

use crate::geometry::fundamental_ops::LinAlg;
use crate::formalism::PIS_IN_180;

// 3D coordinates-types and -matrices from primitives
pub type Coordinate = [f64; 3];

/// Calculate the dihedral between four Coordinate points
/// A dihedral is an angle between four points.
/// Returns a value in degrees [-180. -> 180.]
// Essentially : 
//     get three vector from the four points; b0, b1 and b2
//     from cross(b0, b1) and cross(b1, b2) we get two direction axes
//     -> The dot product between those to direction axes results in the dihedral angle
//
//     Here we use the praxeolitic formula, which involves 1 sqrt and 1 cross product
//     This does not use the description above, but it is more performant than this description
//     See : https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
//

#[pyfunction]
/// Returns the angle of the dihedral in `degrees`
/// Calculate the dihedral between four Coordinate points
/// A dihedral is an angle between four points.
/// Returns a value in degrees [-180. -> 180.]
pub fn dihedral(p0 : Coordinate, p1 : Coordinate, p2 : Coordinate, p3 : Coordinate) -> f64 {
    let b0 = p0.subtract_arr(&p1);
    let b1 = p2.subtract_arr(&p1).normalise_vector(); // do not let magnitude affect subsequent operations
    let b2 = p3.subtract_arr(&p2);
    
    // Vector rejections (as opposed to projections)
    // the b0/b2 will receive a rhs-subtraction from the b1-vector, which has been scaled
    let v = b0.subtract_arr( 
                &b1.scale_vector( 
                    b0.dot_product(&b1)
                )
            );
    let w = b2.subtract_arr( 
                &b1.scale_vector( 
                    b2.dot_product(&b1)
                )
            );
    
    let x = v.dot_product(&w);
    let y = b1.cross_product(&v).dot_product(&w);
    
    // return in degrees : 
    y.atan2(x) * (180. / PI)
    
    
}

#[pyfunction]
pub fn bondangle(p0 : Coordinate, p1 : Coordinate, p2 : Coordinate) -> f64 {

    let a = p0.subtract_arr(&p1).normalise_vector();
    let b = p2.subtract_arr(&p1).normalise_vector();
    b.dot_product(&a).acos() * PIS_IN_180
}


/// -> sqrt(x² + y² + z²)
/// This gives the length of the vector
#[pyfunction]
pub fn bondlength(p0 : Coordinate, p1: Coordinate) -> f64 {
    let c = p0.subtract_arr(&p1);
    c.map(|x: f64| x.powi(2))
        .into_iter()
        .sum::<f64>()
        .sqrt()
}
