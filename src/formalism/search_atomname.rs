
pub trait FindString {
    fn at_position(&self, pattern: &str) -> Result<usize, ()> ;
}

impl FindString for Vec<String> {

    fn at_position(&self, pattern: &str) -> Result<usize, ()> {

        let mut c = 0;

        for name in self {
            if name != pattern {
                c += 1
            } else { 
                return Ok(c as usize) 
            }
        };

        Err(())

//        let a = self.iter().try_fold( 0_u32, |mut prev: u32, atomname| 
//            {
//                if atomname == pattern {
//                    ControlFlow::Continue(prev += 1_u32)
//                } else {
//                    ControlFlow::Break(prev)
//                }
//            });
//
//        Ok(4)
    }
    
}
