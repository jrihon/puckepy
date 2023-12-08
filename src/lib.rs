use pyo3::prelude::*;

mod geometry;
mod conf_sampling;
mod molfile;
mod formalism;

use conf_sampling::{
    peptide::Peptide,
    fivering::Fivering,
    sixring::Sixring,
};

use molfile::{
    Pdb, Xyz
};

use formalism::{
    CP, AS
};

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
/// This is the name of the module
#[pymodule]
fn puckepy(py: Python, m: &PyModule) -> PyResult<()> {

//    m.add_function(wrap_pyfunction!(helloworld, m)?)?;
//    m.add_function(wrap_pyfunction!(add_two, m)?)?;
//    // Testing
//    m.add_class::<Number>()?;

    let geom_sub_module = PyModule::new(py, "geometry")?;
    geom_sub_module.add_function(wrap_pyfunction!(geometry::testing_submodule, geom_sub_module)?)?;
    m.add_submodule(geom_sub_module)?;


    // Add conformational sampling module
    let cs_module = PyModule::new(py, "conf_sampling")?;
    // Add conf_sampling class
    cs_module.add_class::<Peptide>()?;
    cs_module.add_class::<Fivering>()?;
    cs_module.add_class::<Sixring>()?;
    //
    // Append submodule to root module
    m.add_submodule(cs_module)?;

    // Add molfile module
    let molfile_module = PyModule::new(py, "molfile")?;
    // Add function
    molfile_module.add_function(wrap_pyfunction!(molfile::from_pdb, molfile_module)?)?;
    molfile_module.add_function(wrap_pyfunction!(molfile::from_xyz, molfile_module)?)?;
    // Add molfile class
    cs_module.add_class::<Pdb>()?;
    cs_module.add_class::<Xyz>()?;
    //
    // Append submodule to root module
    m.add_submodule(molfile_module)?;

    // Add formalism module
    let form_module = PyModule::new(py, "formalism")?;
    form_module.add_class::<CP>()?;
    form_module.add_class::<AS>()?;
    //
    // Append submodule to root module
    m.add_submodule(form_module)?;
    Ok(())

}

