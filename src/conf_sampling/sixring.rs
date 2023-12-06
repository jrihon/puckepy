use crate::conf_sampling::equidistance_sphere;
use crate::conf_sampling::local_elevation;
use crate::conf_sampling::reconstruct_ring;
use crate::conf_sampling::ring_partition::RingPartition;
use crate::geometry::dihedral;

use ndarray::Array1;
use pyo3::{pyclass, pymethods};

use std::f32::consts::PI;
// two pi; constant

pub const TWOPI : f32 = 2. * PI; 
pub const Z_SIZE: usize = 6;

/// the `alpha` dihedrals according to the Strauss-Piccket (SP) pyranose puckering formalism
/// public `alpha1` field : Vec<f32>
/// public `alpha2` field : Vec<f32>
/// public `alpha3` field : Vec<f32>
#[pyclass(get_all)]
pub struct Sixring {
    alpha1 : Vec<f32>,
    alpha2 : Vec<f32>,
    alpha3 : Vec<f32>,
}

#[pymethods]
impl Sixring {
    /// Initialise the struct with an array of zeroes
    #[new]
    pub fn new(num : usize) -> Self {
        let sphere_size = equidistance_sphere::equidistance_sphere(num as u16);

        let zj = local_elevation::cremerpople_evelation(&sphere_size);

        let projection = zj.projection_and_partition(sphere_size.amount);

        let mut a1 = Vec::with_capacity(sphere_size.amount);
        let mut a2 = Vec::with_capacity(sphere_size.amount);
        let mut a3 = Vec::with_capacity(sphere_size.amount);

        let vec_of_pyranoses = reconstruct_ring::reconstruct_coordinates(
                                &projection,
                                sphere_size.amount,
                                zj, 
                                );

        for pyr in vec_of_pyranoses.iter() {
            a1.push(dihedral(pyr.p5, pyr.p1, pyr.p3, pyr.p2));
            a2.push(dihedral(pyr.p1, pyr.p3, pyr.p5, pyr.p4));
            a3.push(dihedral(pyr.p3, pyr.p5, pyr.p1, pyr.p6));
        };

        Self { 
            alpha1: a1.iter().map(|x| if x < &0. {x + 360.} else {*x}).collect(),
            alpha2: a2.iter().map(|x| if x < &0. {x + 360.} else {*x}).collect(),
            alpha3: a3.iter().map(|x| if x < &0. {x + 360.} else {*x}).collect(),
        }

    }
}
/// The axes to iterate over for sixring molecules : 
/// public `rho` field : f64 . Standard value of 0.67
/// public `theta` field : Array1<f64>. [0, pi] or [0, 180]
/// public `phi` field : Array1<f64>. [0, 2pi] or [0, 360]
/// public `amount` field : Array1<f64>. The corrected amount of points to sample
pub struct SphericalAxes {
    pub rho : f32,
    pub theta : Array1<f32>,
    pub phi : Array1<f32>,
    pub amount : usize,
}


impl SphericalAxes {
    pub fn new(amount: usize, m_theta : usize, rho: f32) -> SphericalAxes {
        SphericalAxes {
            rho, // shorthand initialisation
            theta : Array1::<f32>::zeros(m_theta),
            phi : Array1::<f32>::zeros(amount),
            amount // shorthand initialisation
        }
    }
}

