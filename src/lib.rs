use pyo3::prelude::*;

mod geometry;



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


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
/// This is the name of the module
#[pymodule]
fn puckepy(py: Python, m: &PyModule) -> PyResult<()> {

    m.add_function(wrap_pyfunction!(helloworld, m)?)?;
    m.add_function(wrap_pyfunction!(add_two, m)?)?;

    let geom_module = PyModule::new(py, "geometry")?;
    geom_module.add_function(wrap_pyfunction!(geometry::testing_submodule, geom_module)?)?;
    m.add_submodule(geom_module)?;

    Ok(())

}
