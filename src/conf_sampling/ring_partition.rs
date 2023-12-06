use ndarray::{ArrayBase, Ix2, DataOwned, Array2, Array1};
use crate::conf_sampling::sixring::Z_SIZE;

///
/// All six (6) bond length projections are needed
/// Only three (3) angles are required to progress to reconstruction
/// We need to cosine and the sine value of this angle
///
///
///
///
///
// CONSTANTS : default bond length R_{ij}  1.54 Angstrom
// CONSTANTS : default angle size  B_{ijk} cos(a) = -1/3
const RIJ : f32 = 1.54;
const RIJSQ : f32 = 1.54*1.54;
const COSBIJK : f32 = -1./3. ;
//const BIJK : f32 = 1.9106332362490186 ; // around 109.4712206.. degrees, angle of perfect sp^3 angle


// Store ring partitioning in the struct
pub struct ProjectionPartition {
    pub rpij : Array2::<f32>,
    pub cosbpijk : Array2::<f32>,
    pub sinbpijk : Array2::<f32>,
    pub op : Array1::<f32>,
    pub qp : Array1::<f32>,
    pub oq : Array1::<f32>,

}
// Make a trait where we can implement our own function on the ArrayBase<S,D> struct.
pub trait RingPartition {
    fn projection_and_partition(&self, sphere_size : usize) -> ProjectionPartition;
}


impl<S> RingPartition for ArrayBase<S, Ix2>
where 
    S : DataOwned<Elem = f32>, // Instead of having A as a generic type
                               // we just need A to be f32 types
                               // So we just prompt it in, since we won't use the function for
                               // other type floats or integers
{
    /// The `self` parameter is actually the local_elevation matrix (z_j)
    fn projection_and_partition(&self, sphere_size : usize) -> ProjectionPartition {

        let mut rpij_arr = Array2::<f32>::zeros((sphere_size as usize, Z_SIZE));
        let mut cospb_arr = Array2::<f32>::zeros((sphere_size as usize, Z_SIZE));
        let mut sinpb_arr = Array2::<f32>::zeros((sphere_size as usize, Z_SIZE));
        let mut op_arr = Array1::<f32>::zeros(sphere_size as usize);
        let mut qp_arr = Array1::<f32>::zeros(sphere_size as usize);
        let mut oq_arr = Array1::<f32>::zeros(sphere_size as usize);

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

/* C-Code
 rpij[i]=sqrt(    rij[i] * rij[i] - 
                  (z[i]-z[(i+1)%N]) * (z[i]-z[(i+1)%N]) 
             );

 COSbetapijk[i]= (  pow((z[(i+2)%N]-z[i]),2)-
                    pow((z[(i+1)%N]-z[i]),2)-
                    pow((z[(i+2)%N]-z[(i+1)%N]),2) + 
                    2.*rij[i]*rij[(i+1)%N]*cos(betaijk[i]) )
                    /( 2.*rpij[i]*rpij[(i+1)%N] ) ;

 SINbetapijk[i]= sqrt( 1- COSbetapijk[i]*COSbetapijk[i] );
*/

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
