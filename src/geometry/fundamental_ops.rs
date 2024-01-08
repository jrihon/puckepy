//! PUBLIC FUNCTIONS FOR THE API : DIHEDRAL, BONDANGLE, BONDLENGTH
//!
//!
//!
//! CREMER-POPLE Calculations : dot_product, cross_product, normalise_vector, subtract_arr, RotationMatrix
//!
//! General Linear Algebra and math stuff I do not want to import from different libraries
//! so I write it myself and this way I can implement them on primitives


// 3D coordinates-types and -matrices from primitives
pub type Coordinate = [f64; 3];
pub type DirectionAxis = [f64; 3];
pub type RotationMatrix = [[f64; 3]; 3];

//
// Implement standard linear algebra on the Coordinate types
// This part is not public to the user
pub trait LinAlg {
    fn dot_product(&self, rhs : &Coordinate) -> f64;
    fn cross_product(&self, rhs : &Coordinate) -> DirectionAxis;
    fn subtract_arr(&self, rhs : &Coordinate) -> Coordinate;
    fn add_arr(&self, rhs : &Coordinate) -> Coordinate;
    fn normalise_vector(&self) -> Coordinate;
    fn norm(&self) -> f64;
    fn scale_vector(&self, factor: f64) -> Coordinate;

}


impl LinAlg for Coordinate {

    /// Calculate the scalar product between two vectors
    fn dot_product(&self, rhs : &Coordinate) -> f64 {
      (self[0] * rhs[0]) + (self[1] * rhs[1]) + (self[2] * rhs[2])
    }

    /// Calculate the cross product between two vectors
    fn cross_product(&self, rhs : &Coordinate) -> DirectionAxis {

        [
            (self[1] * rhs[2]) - (self[2]) * (rhs[1]),
           (-self[0] * rhs[2]) + (self[2]) * (rhs[0]),
            (self[0] * rhs[1]) - (self[1]) * (rhs[0]),
        ]

    }

    /// Subtract one Coordinate from another
    fn subtract_arr(&self, rhs : &Coordinate) -> Coordinate {
        [ self[0] - rhs[0], self[1] - rhs[1], self[2] - rhs[2] ]
    }
    ///
    /// Add one Coordinate from another
    fn add_arr(&self, rhs : &Coordinate) -> Coordinate {
        [ self[0] + rhs[0], self[1] + rhs[1], self[2] + rhs[2] ]
    }

    /// -> sqrt(x² + y² + z²)
    /// This gives the length of the vector
    fn norm(&self) -> f64 {
        self.map(|x: f64| x.powi(2))
            .into_iter()
            .sum::<f64>()
            .sqrt()
    }

    /// Normalise the size of the Coordinate
    fn normalise_vector(&self) -> Coordinate {
       let d = 1. / self.norm(); 

       self.map(|x: f64| d * x) // apply the factor `d` to all elements of the coordinate
    }

    /// Scale the vector by a certain factor
    fn scale_vector(&self, factor: f64) -> Coordinate {
        self.map(|x| x * factor)
    }


}





//
// The following functions are not public to the user


// Dot product of two vectors
pub fn dot_product(a : Coordinate, b : Coordinate) -> f64 {
    (a[0] * b[0]) + (a[1] * b[1]) + (a[2] * b[2])
}

// Cross product of two vectors
pub fn cross_product(a : Coordinate, b : Coordinate) -> DirectionAxis {

    [
        (a[1] * b[2]) - (a[2]) * (b[1]),
       (-a[0] * b[2]) + (a[2]) * (b[0]),
        (a[0] * b[1]) - (a[1]) * (b[0]),
    ]
}

pub fn norm(c : [f64;3]) -> f64 {
    c.map(|x: f64| x.powi(2))
        .into_iter()
        .sum::<f64>()
        .sqrt()
}

// Normalise the vector
pub fn normalise_vector(c : Coordinate) -> Coordinate {
    let d = 1. / norm(c); 

    c.map(|x: f64| d * x) // apply the factor `d` to all elements of the coordinate
}

// Subtract two vectors
pub fn subtract_arr(a : Coordinate, b : Coordinate) -> Coordinate {
    [ a[0] - b[0], a[1] - b[1], a[2] - b[2] ]
}



// The following part is not public to the user

/// Only used for Cremer-Pople calculations
/// Custom trait to extend primitive type :
/// Make extension trait on the primitive type `RotationMatrix`
pub trait RotMatrix {
    fn new(phi: f64) -> RotationMatrix; 
    fn apply_rotation(&self, p : Coordinate) -> Coordinate; 
    fn apply_rotation_around_g(&self, p : Coordinate, idx: usize) -> f64; 
}

impl RotMatrix for RotationMatrix {

    /// make a RotationMatrix out of an angle `phi`
    fn new(phi: f64) -> RotationMatrix {

        [ [ phi.cos(), -phi.sin(), 0., ],   // rotation around i-hat (x-axis)
          [ phi.sin(),  phi.cos(), 0., ],   // rotation around j-hat (y-axis)
          [        0.,         0., 1.  ] ]  // rotation around k-hat (z-axis)
        
    }

    /// apply a Rotation to a Coordinate
    fn apply_rotation(&self, p : Coordinate) -> Coordinate {
        let mut c = [0.,0.,0.];

        for (i, arr) in self.iter().enumerate() { 
            c[i] = (arr[0] * p[0]) + (arr[1] * p[1]) + (arr[2] * p[2])
        }

        c
    }

    /// apply a Rotation to a Coordinate, with a specified axis of the RotationMatrix
    fn apply_rotation_around_g(&self, p : Coordinate, idx: usize) -> f64 {
            (self[idx][0] * p[0]) + (self[idx][1] * p[1]) + (self[idx][2] * p[2])
        }
}


#[cfg(test)]
mod test_linalg {

    use assert_float_eq::*;
    use super::*;
    use crate::geometry::molecule_ops::dihedral;

    #[test]
    pub fn planarity() {
        //N1   6.105   8.289   4.633    
        //C2   7.360   7.768   4.827    
        //N3   7.390   6.603   5.551    
        //C4   6.301   5.942   6.079    

        let planar = dihedral(
                        [6.105, 8.289, 4.633],
                        [7.360, 7.768, 4.827],
                        [7.390, 6.603, 5.551],
                        [6.301, 5.942, 6.079]
                        );
        // assert up until the 3rd decimal
        assert_float_absolute_eq!(planar, -0.07, 0.001)
    }

    #[test]
    pub fn chi_angle() {
        //O4'   5.157  10.381   4.681  
        //C1'   5.981   9.551   3.863  
        //N1    6.105   8.289   4.633  
        //C2    7.360   7.768   4.827  

        let chi = dihedral(
                        [5.157, 10.381, 4.681],
                        [5.981, 9.551, 3.863],
                        [6.105, 8.289, 4.633],
                        [7.360, 7.768, 4.827]
                        );
        // assert up until the 3rd decimal
        assert_float_absolute_eq!(chi, -130.214, 0.001)
    }

    #[test]
    pub fn subtract_points() {

        let a = [1., 2., 3.];
        let b = [4., 5., 6.];

        let c = subtract_arr(b, a);
        assert_float_absolute_eq!(c.iter().sum::<f64>(), [3., 3., 3.,].iter().sum::<f64>(), 0.001)

    }

}

