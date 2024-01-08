
// Implement a map to find a pattern (&str) in a Vec<String> and return at which position this is,
// in linear search time .
pub trait FindString {
    fn at_position(&self, pattern: &str) -> Result<usize, ()> ;
}

// Implement this on Vec<String> . Allows the match while looping
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
    }
}
