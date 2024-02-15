use crate::geometry::fundamental_ops::{LinAlg, Coordinate, RotMatrix, RotationMatrix, subtract_arr};
use std::f64::consts::PI;

use crate::conf_sampling::sixring::TWOPI;
use crate::formalism::inversion::{RIJ, RIJSQ, COSBIJK};

// Returns array of coordinates
pub fn invert_sixring(rho: f64, phi2: f64, theta: f64) -> [[f64; 3]; 6] {

    let zj = local_elevation(rho, phi2, theta);
    let projection = projection_and_partition(&zj);
    let sixring = reconstruct_coordinates(projection, zj);

    [
        sixring.p1,
        sixring.p2,
        sixring.p3,
        sixring.p4,
        sixring.p5,
        sixring.p6,
    ]
}



fn local_elevation(rho : f64, phi2: f64, theta: f64) -> [f64;6] {
    
    let phi2 = phi2.to_radians();
    let theta = theta.to_radians();
    // Calculate local elevation
    let constant1 = [0.,1.,2.,3.,4.,5.].map(|j| ((TWOPI * j) / 3.));
    let constant2 = [0,1,2,3,4,5].map(|j| (-1_f64).powi(j));

    let one_over_sqrt_three: f64 = 3_f64.sqrt() ;
    let one_over_sqrt_six: f64 = 6_f64.sqrt() ;

    // iterate over the two arrays of constants
    constant1.iter().zip(constant2.iter()).map(|(one, two)| {

    let term1 = (theta.sin() * (phi2 + one).cos()) / one_over_sqrt_three;
    let term2 = (theta.cos() * two) / one_over_sqrt_six ; // second term of the equation
    (term1 + term2) * rho // multiply by rho, which was pushed out by both terms to the outside of the equation

    }).collect::<Vec<f64>>()
      .try_into().unwrap() // we are certain that it will collect into a [f64;6] as both arrays
                           // are of the same size.
}

// Store ring partitioning in the struct
struct ProjectionPartition {
    pub rpij : [f64;6],
    pub cosbpijk : [f64;6],
    pub sinbpijk : [f64;6],
    pub op : f64,
    pub qp : f64,
    pub oq : f64,

}

fn projection_and_partition(zj: &[f64;6]) -> ProjectionPartition {

    let mut rpij_arr: [f64;6] = [0.;6];
    let mut cospb_arr: [f64;6] =  [0.;6];
    let mut sinpb_arr: [f64;6] = [0.;6];

    for j in 0..6 {
        rpij_arr[j] = ( RIJSQ - 
                            ( zj[j] - zj[(j+1) % 6] ).powi(2)
                          ).sqrt();
    }

    for j in 0..6 {

        // sphere points are in radians
        // the values of the cosine values are abnormal
        // they all appear in values above 2PI and are often negative. This shouldnt be the
        // case, where cosine values can only be between [-1 , 1]
        cospb_arr[j] = ( (zj[(j+2) % 6] - zj[j]).powi(2) // zk - zi 
                           - (zj[(j+1) % 6] - zj[j]).powi(2) // zj - zi
                           - (zj[(j+2) % 6] - zj[(j+1) % 6]).powi(2) // zk - zj
                           + (2. * RIJ * RIJ * COSBIJK) // 2 * rij * rjk * cos Bijk
                           ) / (2. * rpij_arr[j] * rpij_arr[(j+1) % 6] ); // 2 * rpij * rpjk 

        sinpb_arr[j] = (1. - &cospb_arr[j].powi(2) ).sqrt();
        
    };

    let op: f64 = (( rpij_arr[0].powi(2) + rpij_arr[1].powi(2) ) - (2. * rpij_arr[0] * rpij_arr[1] * cospb_arr[0])).sqrt();
    let qp: f64 = (( rpij_arr[2].powi(2) + rpij_arr[3].powi(2) ) - (2. * rpij_arr[2] * rpij_arr[3] * cospb_arr[2])).sqrt();
    let oq: f64 = (( rpij_arr[4].powi(2) + rpij_arr[5].powi(2) ) - (2. * rpij_arr[4] * rpij_arr[5] * cospb_arr[4])).sqrt();


    ProjectionPartition { 
        rpij: rpij_arr,
        cosbpijk: cospb_arr,
        sinbpijk: sinpb_arr,
        op,
        qp,
        oq,
    }
}

struct SixRingAtoms {
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
    }
    p_g
    }
    
}

#[allow(dead_code)] // -> fields s11, s25 and s31 are never read. Included for declarative purposes
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

fn reconstruct_coordinates(proj : ProjectionPartition, z_j : [f64;6]) -> SixRingAtoms {
    // proj : projections and partitioning. 

    // Add the local evelation already as the z-coordinate to the final molecule's array
    let mut sixring = SixRingAtoms {
        p1 : [0., 0., z_j[0]],
        p2 : [0., 0., z_j[1]],
        p3 : [0., 0., z_j[2]],
        p4 : [0., 0., z_j[3]],
        p5 : [0., 0., z_j[4]],
        p6 : [0., 0., z_j[5]],
    };

    let pyranose = PointPositions {
            s11 : 
                [0.,
                 0.,
                 0.],
            s12 : 
                [-proj.rpij[0],
                 0.,
                 0.],
            s13 : 
                [(-proj.rpij[0]) + (proj.rpij[1] * proj.cosbpijk[0]),
                 proj.rpij[1] * proj.sinbpijk[0],
                 0.],
            s23 : 
                [(proj.oq + proj.rpij[3]) - (proj.rpij[2] * proj.cosbpijk[2]),
                 proj.rpij[2] * proj.sinbpijk[2],
                 0.],
            s24 : 
                [proj.oq + proj.rpij[3],
                 0.,
                 0.],
            s25 : 
                [proj.oq ,
                 0.,
                 0.],
            s35 : 
                [proj.rpij[5] - (proj.rpij[4] * proj.cosbpijk[4]),
                 proj.rpij[4] * proj.sinbpijk[4],
                 0.],
            s36 :
                [proj.rpij[5],
                 0.,
                 0.],
            s31 : 
                [0.,
                 0.,
                 0.],
    };

    let rho1 = pyranose.s13[1].atan2(pyranose.s13[0]);
    let rho2 = pyranose.s23[1].atan2(pyranose.s23[0] - proj.oq);
    let rho3 = pyranose.s35[1].atan2(pyranose.s35[0]);

    let p_o : Coordinate = [0.,
                            0.,
                            0.]; //pO
    let p_p : Coordinate = [(proj.op.powi(2) + proj.oq.powi(2) - proj.qp.powi(2))/(2. * proj.oq),
                            (proj.op.powi(2) - ( ( (proj.op.powi(2) + proj.oq.powi(2) - proj.qp.powi(2)).powi(2) ) / (4. * proj.op.powi(2)) ) ).sqrt(),
                            0.]; //pP
    let p_q : Coordinate = [proj.oq,
                            0.,
                            0.]; //pQ

    let rho_ps1 = p_p[1].atan2(p_p[0]);
    let rho_ps2 = p_p[1].atan2(p_p[0] - proj.oq);

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
    // Rotate the xi-coordinate correctly
    sixring.p1[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p1, p_g), 0);
    sixring.p2[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p2, p_g), 0);
    sixring.p3[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p3, p_g), 0);
    sixring.p4[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p4, p_g), 0);
    sixring.p5[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p5, p_g), 0);
    sixring.p6[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p6, p_g), 0);

    // Rotate the yi-coordinate correctly
    sixring.p1[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p1, p_g), 1);
    sixring.p2[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p2, p_g), 1);
    sixring.p3[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p3, p_g), 1);
    sixring.p4[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p4, p_g), 1);
    sixring.p5[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p5, p_g), 1);
    sixring.p6[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_sixring.p6, p_g), 1);


    sixring
}
