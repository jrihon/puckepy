pub mod cremerpople;
pub mod altonasund;
pub mod strausspickett;
pub mod moleculefile;


mod search_atomname;// match a pattern in a Vec<String>. If not found, Err(()) => panic!()
                    // Used in self.from_atomnames() methods

pub mod inversion;  // include the inversion module
                    // inversion methods on the CP*{} Structs


// Used to convert -> to radians
pub const PIS_IN_180: f64 = 57.2957795130823208767981548141051703_f64;
