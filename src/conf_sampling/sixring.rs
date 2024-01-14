use pyo3::{pyclass, pymethods};
use ndarray::{Array1, Array2, ArrayBase, DataOwned, Ix2};

// Crate imports
use crate::geometry::fundamental_ops::{subtract_arr, RotationMatrix, Coordinate, RotMatrix, LinAlg};
use crate::geometry::molecule_ops::dihedral;

use std::f64::consts::PI;
use crate::formalism::inversion::{RIJ,RIJSQ,COSBIJK};
pub const TWOPI : f64 = 2. * PI; 
pub const Z_SIZE: usize = 6;
pub const RHO : f64 = 0.67; // radius of the sphere; constant

/// the `alpha` dihedrals according to the Strauss-Piccket (SP) pyranose puckering formalism
/// public `alpha1` field : Vec<f64>
/// public `alpha2` field : Vec<f64>
/// public `alpha3` field : Vec<f64>
#[pyclass(get_all)]
pub struct Sixring {
    pub alpha1 : Vec<f64>,
    pub alpha2 : Vec<f64>,
    pub alpha3 : Vec<f64>,
}

#[pymethods]
impl Sixring {

    #[new]
    pub fn new(amount : usize) -> Self {
        let (globe_points, corrected_amount) = equidistance_sphere(amount as u16);

        let zj = cremerpople_evelation(&globe_points, corrected_amount);

        let projection = zj.projection_and_partition(corrected_amount);

        let mut a1 = Vec::with_capacity(corrected_amount);
        let mut a2 = Vec::with_capacity(corrected_amount);
        let mut a3 = Vec::with_capacity(corrected_amount);

        let vec_of_pyranoses = reconstruct_coordinates(
                                projection,
                                corrected_amount,
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

#[pyclass(get_all)]
pub struct SixringAxes {
    pub rho : f64,
    pub theta : Vec<f64>,
    pub phi : Vec<f64>,
}

#[pymethods]
impl SixringAxes {

    #[new]
    fn new(amount : usize) -> Self {

        let (globe_points, corrected_amount) = equidistance_sphere(amount as u16);
        let mut theta_vec: Vec<f64> = Vec::with_capacity(corrected_amount);
        let mut phi_vec: Vec<f64> = Vec::with_capacity(corrected_amount);

        let mut idx_theta: usize = 0;
        for i in 0..corrected_amount { 
            if (globe_points.phi[i] == 0.0) && i != 0 {
                idx_theta += 1 
            };

            theta_vec.push(globe_points.theta[idx_theta]);
            phi_vec.push(globe_points.phi[i]);
        }
        
        Self {
            rho: globe_points.rho,
            theta: theta_vec,
            phi: phi_vec,
        }
    }
}



/// The axes to iterate over for sixring molecules : 
/// public `rho` field : f64 . Standard value of 0.67
/// public `theta` field : Array1<f64>. [0, pi] or [0, 180]
/// public `phi` field : Array1<f64>. [0, 2pi] or [0, 360]
struct SAxes {
    pub rho : f64,
    pub theta : Array1<f64>,
    pub phi : Array1<f64>,
}


impl SAxes {
    pub fn new(m_theta : usize, rho: f64, amount: usize) -> Self {
        Self {
            rho, 
            theta : Array1::<f64>::zeros(m_theta),
            phi : Array1::<f64>::zeros(amount),
        }
    }
}





/// Generate a sphere with equidistant points on its surface. This represents the 
/// Cremer-Pople sphere, where every point represents a specific conformation
///
///
fn equidistance_sphere(amount : u16 ) -> (SAxes, usize) {
    // Set a value as surface area / points
    let corrected_amount: f64 = corrected_amount_by_sphere_radius(amount as f64);
    let a: f64 = ( 4. * PI * RHO.powi(2)) / corrected_amount;

    let mut idx : u32 = 0; // indexing the arrays

    // Set d as the square root of a
    let d: f64 = a.sqrt();

    // Round of the ratio between PI and the value of d
    let m_theta: f64 = (PI / d).round();

    // Set d_theta and d_phi
    let d_theta: f64 = PI / m_theta;
    let d_phi: f64 = a / d_theta;

    let amount_sizeof: usize = corrected_amount_to_size_up_arrays(m_theta, d_phi);
    let mut globe = SAxes::new(m_theta as usize, RHO, amount_sizeof as usize);

    for m in 0..m_theta as u32 {
        globe.theta[m as usize] = (PI * (m as f64 + 0.5)) / m_theta;
        let m_phi: f64 = (TWOPI * globe.theta[m as usize].sin() / d_phi).round();

        for n in 0..m_phi as u32 {
            globe.phi[idx as usize] = (TWOPI * n as f64) / m_phi;
            idx += 1;
            
        }
    }

    (globe, amount_sizeof) // return the struct containing all the points on the surface of our globe
}


fn corrected_amount_to_size_up_arrays(m_theta : f64, d_phi : f64) -> usize {
    // Counting the amount of points that are actually generated
    let mut size_array: u32 = 0;

    for m in 0..m_theta as u32 {
        let theta: f64 = (PI * (m as f64 + 0.5)) / m_theta;
        let m_phi: f64 = (TWOPI * theta.sin() / d_phi).round();
        size_array += m_phi as u32;

    };

    size_array as usize // return exact amount of points that will be sampled over
}

/// Markus Deserno's mathematics only works out if we commit to a radius = 1 unit
/// 
/// Since the Rho value (the radius of the sphere) is set to be 0.67 for our purposes (see
/// Cremer-Pople standard puckering values), we need to correct the amount of prompted points.
/// 
/// What we need is the ratio of the surface are at rho(0.67) and rho(1.00)
/// --> (0.67^2 * PI * 4) / (1.00^2 * PI * 4) => 0.67^2
fn corrected_amount_by_sphere_radius(num : f64) -> f64 {
    num * RHO.powi(2)
}



/// Calculate the local elevation z_j for the Cremer-Pople coordinate prompted 
///
///
fn cremerpople_evelation(sphere : &SAxes, amount: usize) -> Array2<f64> {
    // spherical coordinates are by default in radians

    // 6 atomic elevations (Z_j) for any set of (r, theta, phi)
    let mut z: Array2<f64> = Array2::zeros((amount, Z_SIZE));

    // Set two constant values
    let constant1 = [0.,1.,2.,3.,4.,5.].map(|j| ((TWOPI * j) / 3.));
    let constant2 = [0, 1, 2, 3, 4, 5].map(|j| (-1_f64).powi(j));

    // Set some more constant values
    let one_over_sqrt_three: f64 = 3_f64.sqrt() ;
    let one_over_sqrt_six: f64 = 6_f64.sqrt() ;

    let mut idx_theta: usize = 0;
    for i in 0..amount { 
        // the way we generate the sphere is in layered circles.
        // every new circle, we start off again at phi == 0.0
        // if we move to a new layer; we have to the next theta value
        // NOTE :the theta and phi arrays are not of the same length
        if (sphere.phi[i] == 0.0) && i != 0 {
            idx_theta += 1 
        };

        for j in 0..Z_SIZE {
            z[[i, j]] = calculate_local_elevation(sphere.rho, sphere.theta[idx_theta], sphere.phi[i], constant1[j], constant2[j],
                                                  one_over_sqrt_three, one_over_sqrt_six);
        }
    }

    z // return local elevation of all the conformations

}

/// Theta and Phi are prompted as the spherical coordinates
/// In actually : 
/// Theta (Spherical Coordinate) => Phi_2 (Cremer-Pople)
/// Phi   (Spherical Coordinate) => Theta (Cremer-Pople)
fn calculate_local_elevation(rho : f64, theta: f64, phi: f64, c1 : f64,  c2 : f64,
                             onethree: f64, onesix: f64) -> f64 {
    
    let term1 = (theta.sin() * (phi + c1).cos()) / onethree; // first term of the equation
    let term2 = (theta.cos() * c2) / onesix ; // second term of the equation

    (term1 + term2) * rho // multiply by rho, which was pushed out by both terms to the outside of the equation
}



/// Calculate projected bondlengths and bond angles 
pub struct ProjectionPartition {
    pub rpij : Array2::<f64>,
    pub cosbpijk : Array2::<f64>,
    pub sinbpijk : Array2::<f64>,
    pub op : Array1::<f64>,
    pub qp : Array1::<f64>,
    pub oq : Array1::<f64>,

}
// Make a trait where we can implement our own function on the ArrayBase<S,D> struct.
pub trait RingPartition {
    fn projection_and_partition(&self, sphere_size : usize) -> ProjectionPartition;
}


impl<S> RingPartition for ArrayBase<S, Ix2>
where 
    S : DataOwned<Elem = f64>, // Instead of having A as a generic type
                               // we just need A to be f64 types
                               // So we just prompt it in, since we won't use the function for
                               // other type floats or integers
{
    /// The `self` parameter is actually the local_elevation matrix (z_j)
    fn projection_and_partition(&self, sphere_size : usize) -> ProjectionPartition {

        let mut rpij_arr = Array2::<f64>::zeros((sphere_size as usize, Z_SIZE));
        let mut cospb_arr = Array2::<f64>::zeros((sphere_size as usize, Z_SIZE));
        let mut sinpb_arr = Array2::<f64>::zeros((sphere_size as usize, Z_SIZE));
        let mut op_arr = Array1::<f64>::zeros(sphere_size as usize);
        let mut qp_arr = Array1::<f64>::zeros(sphere_size as usize);
        let mut oq_arr = Array1::<f64>::zeros(sphere_size as usize);

        for i in 0..sphere_size as usize {
            
            for j in 0..Z_SIZE {
                rpij_arr[[i,j]] = ( RIJSQ - 
                                    ( self[[i,j]] - self[[i, (j+1) % Z_SIZE]] ).powi(2)
                                  ).sqrt();
            }

            for j in 0..Z_SIZE {

                // sphere points are in radians
                // the values of the cosine values are abnormal
                // they all appear in values above 2PI and are often negative. This shouldnt be the
                // case, where cosine values can only be between [-1 , 1]
                cospb_arr[[i,j]] = ( (self[[i, (j+2) % Z_SIZE]] - self[[i,j]]).powi(2) // zk - zi 
                                   - (self[[i, (j+1) % Z_SIZE]] - self[[i,j]]).powi(2) // zj - zi
                                   - (self[[i, (j+2) % Z_SIZE]] - self[[i,(j+1) % Z_SIZE]]).powi(2) // zk - zj
                                   + (2. * RIJ * RIJ * COSBIJK) // 2 * rij * rjk * cos Bijk
                                   ) / (2. * rpij_arr[[i,j]] * rpij_arr[[i, (j+1) % Z_SIZE]] ); // 2 * rpij * rpjk 

                sinpb_arr[[i,j]] = (1. - &cospb_arr[[i,j]].powi(2) ).sqrt();
                
            };

            op_arr[i] = (( rpij_arr[[i,0]].powi(2) + rpij_arr[[i,1]].powi(2) ) - (2. * rpij_arr[[i,0]] * rpij_arr[[i,1]] * cospb_arr[[i, 0]])).sqrt();
            qp_arr[i] = (( rpij_arr[[i,2]].powi(2) + rpij_arr[[i,3]].powi(2) ) - (2. * rpij_arr[[i,2]] * rpij_arr[[i,3]] * cospb_arr[[i, 2]])).sqrt();
            oq_arr[i] = (( rpij_arr[[i,4]].powi(2) + rpij_arr[[i,5]].powi(2) ) - (2. * rpij_arr[[i,4]] * rpij_arr[[i,5]] * cospb_arr[[i, 4]])).sqrt();

        }

        ProjectionPartition { 
            rpij: rpij_arr,
            cosbpijk: cospb_arr,
            sinbpijk: sinpb_arr,
            op: op_arr,
            qp: qp_arr,
            oq: oq_arr,
        }

    }
}



/// Position the segments and the points on the (x',y') plane
///
///
#[allow(dead_code)] // -> fields s11, s25 and s31 are never read. Included for declarative purposes
#[derive(Debug)]
struct PointPositions {
    s11 : Coordinate,
    s12 : Coordinate,
    s13 : Coordinate,
    s23 : Coordinate,
    s24 : Coordinate,
    s25 : Coordinate,
    s35 : Coordinate,
    s36 : Coordinate,
    s31 : Coordinate,
}

#[derive(Debug)]
pub struct SixRingAtoms {
    pub p1 : Coordinate,
    pub p2 : Coordinate,
    pub p3 : Coordinate,
    pub p4 : Coordinate,
    pub p5 : Coordinate,
    pub p6 : Coordinate,

}

impl SixRingAtoms {
    fn calculate_geometric_center(&self) -> Coordinate {

    let mut p_g = [0.;3];

    for i in 0..3 { // 0 -> 1 -> 2
        p_g[i] = (self.p1[i] + self.p2[i] + self.p3[i] + self.p4[i] + self.p5[i] + self.p6[i]) / 6.
    };
    p_g
    }
    
}


/// Return a Vec of all the conformers' atoms' position in cartesian coordinates
/// Then we will derive all the alpha angles (the improper dihedrals) in a next function
pub fn reconstruct_coordinates(proj : ProjectionPartition, sphere_size : usize, z_j : Array2<f64>) -> Vec<SixRingAtoms> {
    // proj : projections and partitioning. 

    let mut pyranosecoordinates = Vec::with_capacity(sphere_size);

    for i in 0..sphere_size {


        // Add the local evelation already as the z-coordinate to the final molecule's array
        let mut sixring = SixRingAtoms {
            p1 : [0., 0., z_j[[i, 0]]],
            p2 : [0., 0., z_j[[i, 1]]],
            p3 : [0., 0., z_j[[i, 2]]],
            p4 : [0., 0., z_j[[i, 3]]],
            p5 : [0., 0., z_j[[i, 4]]],
            p6 : [0., 0., z_j[[i, 5]]],
        };

        let pyranose = PointPositions {
                s11 : 
                    [0.,
                     0.,
                     0.],
                s12 : 
                    [-proj.rpij[[i,0]],
                     0.,
                     0.],
                s13 : 
                    [(-proj.rpij[[i,0]]) + (proj.rpij[[i,1]] * proj.cosbpijk[[i,0]]),
                     proj.rpij[[i,1]] * proj.sinbpijk[[i,0]],
                     0.],
                s23 : 
                    [(proj.oq[i] + proj.rpij[[i,3]]) - (proj.rpij[[i,2]] * proj.cosbpijk[[i,2]]),
                     proj.rpij[[i,2]] * proj.sinbpijk[[i,2]],
                     0.],
                s24 : 
                    [proj.oq[i] + proj.rpij[[i,3]],
                     0.,
                     0.],
                s25 : 
                    [proj.oq[i] ,
                     0.,
                     0.],
                s35 : 
                    [proj.rpij[[i,5]] - (proj.rpij[[i,4]] * proj.cosbpijk[[i,4]]),
                     proj.rpij[[i,4]] * proj.sinbpijk[[i,4]],
                     0.],
                s36 :
                    [proj.rpij[[i,5]],
                     0.,
                     0.],
                s31 : 
                    [0.,
                     0.,
                     0.],
        };

        let rho1 = pyranose.s13[1].atan2(pyranose.s13[0]);
        let rho2 = pyranose.s23[1].atan2(pyranose.s23[0] - proj.oq[i]);
        let rho3 = pyranose.s35[1].atan2(pyranose.s35[0]);

        let p_o : Coordinate = [0.,
                                0.,
                                0.]; //pO
        let p_p : Coordinate = [(proj.op[i].powi(2) + proj.oq[i].powi(2) - proj.qp[i].powi(2))/(2. * proj.oq[i]),
                                (proj.op[i].powi(2) - ( ( (proj.op[i].powi(2) + proj.oq[i].powi(2) - proj.qp[i].powi(2)).powi(2) ) / (4. * proj.op[i].powi(2)) ) ).sqrt(),
                                0.]; //pP
        let p_q : Coordinate = [proj.oq[i],
                                0.,
                                0.]; //pQ

        let rho_ps1 = p_p[1].atan2(p_p[0]);
        let rho_ps2 = p_p[1].atan2(p_p[0] - proj.oq[i]);

        let sigma1 = rho1 - rho_ps1;
        let sigma2 = rho_ps2 - rho2;
        let sigma3 = rho3;

        // p1, p3, p5 already exist on the xy'-plane, so need only to rotate p2,p4,p6
        let tmp_sixring = SixRingAtoms {
            p1 : p_o,
            p2 : RotationMatrix::new(-sigma1).apply_rotation(pyranose.s12),
            p3 : p_p,
            p4 : RotationMatrix::new(sigma2).apply_rotation(subtract_arr(pyranose.s24, p_q)).add_arr(&p_q),
            p5 : p_q,
            p6 : RotationMatrix::new(-sigma3).apply_rotation(pyranose.s36),
        };

        // Calculate geometric center
        let p_g : Coordinate = tmp_sixring.calculate_geometric_center();
        // Derive final rotation matrix
        let rho_g = (PI / 2.) + p_g[1].atan2(p_g[0]);
        let rot4 = RotationMatrix::new(-rho_g);

        // final rotation
        sixring.p1[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p1, p_g), 0);
        sixring.p2[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p2, p_g), 0);
        sixring.p3[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p3, p_g), 0);
        sixring.p4[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p4, p_g), 0);
        sixring.p5[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p5, p_g), 0);
        sixring.p6[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p6, p_g), 0);

        sixring.p1[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p1, p_g), 1);
        sixring.p2[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p2, p_g), 1);
        sixring.p3[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p3, p_g), 1);
        sixring.p4[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p4, p_g), 1);
        sixring.p5[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p5, p_g), 1);
        sixring.p6[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p6, p_g), 1);


        pyranosecoordinates.push(sixring);
    }

    pyranosecoordinates // return Vec<SixRingAtoms>
}
