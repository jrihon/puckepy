use crate::geometry::fundamental_ops::{Coordinate, RotMatrix, RotationMatrix, subtract_arr};
use std::f64::consts::PI;

use crate::conf_sampling::sixring::TWOPI;
use crate::formalism::inversion::{RIJ, RIJSQ, COSBIJK};

// Return the array of coordinates
pub fn invert_fivering(rho: f64, phi2: f64) -> [[f64; 3]; 5] {

    let zj = local_elevation(rho, phi2);
    let projection = projection_and_partition(&zj);
    let fivering = reconstruct_coordinates(projection, zj);

    [
        fivering.p1,
        fivering.p2,
        fivering.p3,
        fivering.p4,
        fivering.p5,
    ]
}



fn local_elevation(rho : f64, phi2: f64) -> [f64;5] {
    
    let phi2 = phi2.to_radians();
    // Calculate local elevation
    let constant1 = [0.,1.,2.,3.,4.].map(|j| ((2. * TWOPI * j) / 5.));

    let two_fifth_sqrt: f64 = (2_f64/5_f64).sqrt() ;
//    let one_over_sqrt_six: f64 = 6_f64.sqrt() ;

    // Rework of Z_j over fiverings works
    constant1.iter().map(|one| {

    ((phi2 + one).cos() * two_fifth_sqrt) * rho

    }).collect::<Vec<f64>>()
      .try_into().unwrap() // we are certain that it will collect into a [f64;6] as both arrays
                            // are of the same size

}

// Store ring partitioning in the struct
struct ProjectionPartition {
    pub rpij : [f64;5],
    pub cosbpijk : [f64;5],
    pub sinbpijk : [f64;5],
    pub op : f64,
    pub qp : f64,
    pub oq : f64,

}

fn projection_and_partition(zj: &[f64;5]) -> ProjectionPartition {

    let mut rpij_arr: [f64;5] = [0.;5];
    let mut cospb_arr: [f64;5] =  [0.;5];
    let mut sinpb_arr: [f64;5] = [0.;5];

    for j in 0..5 {
        rpij_arr[j] = ( RIJSQ - 
                            ( zj[j] - zj[(j+1) % 5] ).powi(2)
                          ).sqrt();
    }

    for j in 0..5 {

        // sphere points are in radians
        // the values of the cosine values are abnormal
        // they all appear in values above 2PI and are often negative. This shouldnt be the
        // case, where cosine values can only be between [-1 , 1]
        cospb_arr[j] = ( (zj[(j+2) % 5] - zj[j]).powi(2) // zk - zi 
                           - (zj[(j+1) % 5] - zj[j]).powi(2) // zj - zi
                           - (zj[(j+2) % 5] - zj[(j+1) % 5]).powi(2) // zk - zj
                           + (2. * RIJ * RIJ * COSBIJK) // 2 * rij * rjk * cos Bijk
                           ) / (2. * rpij_arr[j] * rpij_arr[(j+1) % 5] ); // 2 * rpij * rpjk 

        sinpb_arr[j] = (1. - cospb_arr[j].powi(2) ).sqrt();
        
    };

    // S1 : op = (r12, r23, B123)
    // S2 : qp = r34
    // S3 : oq = (r45, r51, B451)
    let op: f64 = (( rpij_arr[0].powi(2) + rpij_arr[1].powi(2) ) - (2. * rpij_arr[0] * rpij_arr[1] * cospb_arr[0])).sqrt();
    let qp: f64 = rpij_arr[2]; // S3 is just a segment of the line, since sqrt(r_51 * r_51) = r_51,
                               // and N - 3 = 2 bond angles total, so cosB234 does not exist here
    let oq: f64 = (( rpij_arr[3].powi(2) + rpij_arr[4].powi(2) ) - (2. * rpij_arr[3] * rpij_arr[4] * cospb_arr[3])).sqrt();

    ProjectionPartition { 
        rpij: rpij_arr,
        cosbpijk: cospb_arr,
        sinbpijk: sinpb_arr,
        op,
        qp,
        oq,
    }
}

struct FiveRingAtoms {
    pub p1 : Coordinate,
    pub p2 : Coordinate,
    pub p3 : Coordinate,
    pub p4 : Coordinate,
    pub p5 : Coordinate,

}

impl FiveRingAtoms {
    fn calculate_geometric_center(&self) -> Coordinate {

    let mut p_g = [0.;3];

    for i in 0..3 { // 0 -> 1 -> 2
        p_g[i] = (self.p1[i] + self.p2[i] + self.p3[i] + self.p4[i] + self.p5[i]) / 6.
    }
    p_g
    }
    
}

#[allow(dead_code)] // -> fields s11, s25 and s31 are never read. Included for completive purposes
struct PointPositions {
    s11 : Coordinate,
    s12 : Coordinate,
    s13 : Coordinate,
    s23 : Coordinate,
    s24 : Coordinate,
    s34 : Coordinate,
    s35 : Coordinate,
    s31 : Coordinate,
}

fn reconstruct_coordinates(proj : ProjectionPartition, z_j : [f64;5]) -> FiveRingAtoms {

    // Create the partitions
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
                [proj.oq + proj.qp,
                 0.,
                 0.],
            s24 : 
                [proj.oq ,
                 0.,
                 0.],
            s34 : 
                [proj.rpij[4] - (proj.rpij[3] * proj.cosbpijk[4]),
                 proj.rpij[3] * proj.sinbpijk[4],
                 0.],
            s35 :
                [proj.rpij[4],
                 0.,
                 0.],
            s31 : 
                [0.,
                 0.,
                 0.],
    };


    let rho1 = pyranose.s13[1].atan2(pyranose.s13[0]);
    let rho3 = pyranose.s34[1].atan2(pyranose.s34[0]);

    let p_o : Coordinate = [0.,
                            0.,
                            0.]; //point O in triangle, origin 
    let p_p : Coordinate = [(proj.op.powi(2) + proj.oq.powi(2) - proj.qp.powi(2))/(2. * proj.oq),
                            (proj.op.powi(2) - ( ( (proj.op.powi(2) + proj.oq.powi(2) - proj.qp.powi(2)).powi(2) ) / (4. * proj.op.powi(2)) ) ).sqrt(),
                            0.]; //point P in triangle, top
    let p_q : Coordinate = [proj.oq,
                            0.,
                            0.]; //point Q, on x axis

    let rho_ps1 = p_p[1].atan2(p_p[0]);

    let sigma1 = rho1 - rho_ps1;
    let sigma3 = rho3;

    // p1, p3, p5 already exist on the xy'-plane, so need only to rotate p2,p4,p6
    let tmp_fivering = FiveRingAtoms {
        p1 : p_o,
        p2 : RotationMatrix::new(-sigma1).apply_rotation(pyranose.s12),
        p3 : p_p,
        p4 : p_q,
        p5 : RotationMatrix::new(-sigma3).apply_rotation(pyranose.s35),
    };

    // Calculate geometric center
    let p_g : Coordinate = tmp_fivering.calculate_geometric_center();
    // Derive final rotation matrix
    let rho_g = (PI / 2.) + p_g[1].atan2(p_g[0]);
    let rot4 = RotationMatrix::new(-rho_g);



    // Add the local evelation already as the z-coordinate to the final molecule's array
    let mut fivering = FiveRingAtoms {
        p1 : [0., 0., z_j[0]],
        p2 : [0., 0., z_j[1]],
        p3 : [0., 0., z_j[2]],
        p4 : [0., 0., z_j[3]],
        p5 : [0., 0., z_j[4]],
    };

    // final rotation
    // Rotate the xi-coordinate correctly
    fivering.p1[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p1, p_g), 0);
    fivering.p2[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p2, p_g), 0);
    fivering.p3[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p3, p_g), 0);
    fivering.p4[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p4, p_g), 0);
    fivering.p5[0] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p5, p_g), 0);

    // Rotate the yi-coordinate correctly
    fivering.p1[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p1, p_g), 1);
    fivering.p2[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p2, p_g), 1);
    fivering.p3[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p3, p_g), 1);
    fivering.p4[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p4, p_g), 1);
    fivering.p5[1] = rot4.apply_rotation_around_g(subtract_arr(tmp_fivering.p5, p_g), 1);


    fivering
}
