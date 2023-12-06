use pyo3::prelude::*;

mod geometry;
mod conf_sampling;

use conf_sampling::{
    peptide::Peptide,
    fivering::Fivering,
    sixring::Sixring,
};

/// Create a function to call from the Python package
#[pyfunction]
fn helloworld() {
    println!("Hello, from Rust!");
}


#[pyfunction]
/// This is the name of the function
fn add_two(a: i32) -> i32 {
    a + 2
}

#[pyclass]
struct Number {
    nummie: i32
}

#[pymethods()]
impl Number {

    #[new]
    fn new(nummie: i32) -> Self {
        Self {nummie}
    }

    fn get_number(&self) -> i32 {
        self.nummie
    }
}
/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
/// This is the name of the module
#[pymodule]
fn puckepy(py: Python, m: &PyModule) -> PyResult<()> {

    m.add_function(wrap_pyfunction!(helloworld, m)?)?;
    m.add_function(wrap_pyfunction!(add_two, m)?)?;

    let geom_sub_module = PyModule::new(py, "geometry")?;
    geom_sub_module.add_function(wrap_pyfunction!(geometry::testing_submodule, geom_sub_module)?)?;
    m.add_submodule(geom_sub_module)?;

    // Testing
    m.add_class::<Number>()?;

    // Add module
    let cs_module = PyModule::new(py, "conf_sampling")?;
    // Add class
    cs_module.add_class::<Peptide>()?;
    cs_module.add_class::<Fivering>()?;
    cs_module.add_class::<Sixring>()?;
    // Append submodule to root module
    m.add_submodule(cs_module)?;

    Ok(())

}
