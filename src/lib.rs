use pyo3::prelude::*;

mod geometry;
mod conf_sampling;
mod formalism;
mod inversion;

use conf_sampling::{
    peptide::Peptide,
    fivering::Fivering,
    sixring::Sixring,
};

use formalism::{
    CP5, CP6, AS, Pdb, Xyz
};

use inversion::{
    sixring,
    fivering,
    write_pdb,
    write_xyz
};


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
/// This is the name of the module
#[pymodule]
fn puckepy(py: Python, m: &PyModule) -> PyResult<()> {

    let geom_sub_module = PyModule::new(py, "geometry")?;
    geom_sub_module.add_function(wrap_pyfunction!(geometry::testing_submodule, geom_sub_module)?)?;
    m.add_submodule(geom_sub_module)?;


    // Add conformational sampling module
    let cs_module = PyModule::new(py, "confsampling")?;
    // Add conf_sampling class
    cs_module.add_class::<Peptide>()?;
    cs_module.add_class::<Fivering>()?;
    cs_module.add_class::<Sixring>()?;
    //
    // Append submodule to root module
    m.add_submodule(cs_module)?;

    // Add formalism module
    let form_module = PyModule::new(py, "formalism")?;
    form_module.add_class::<CP5>()?;
    form_module.add_class::<CP6>()?;
    form_module.add_class::<AS>()?;
    form_module.add_class::<Pdb>()?;
    form_module.add_class::<Xyz>()?;
    //
    // Append submodule to root module
    m.add_submodule(form_module)?;

    // Add inversion module
    let inv_module = PyModule::new(py, "inversion")?;
    inv_module.add_function(wrap_pyfunction!(sixring::invert_sixring, inv_module)?)?;
    inv_module.add_function(wrap_pyfunction!(fivering::invert_fivering, inv_module)?)?;
    inv_module.add_function(wrap_pyfunction!(write_pdb, inv_module)?)?;
    inv_module.add_function(wrap_pyfunction!(write_xyz, inv_module)?)?;
    //
    // Append submodule to root module
    m.add_submodule(inv_module)?;
    Ok(())

}

