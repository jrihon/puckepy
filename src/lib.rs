use pyo3::prelude::*;

mod geometry;
use geometry::molecule_ops::{
    dihedral,
    bondangle,
    bondlength
};

mod conf_sampling;
use conf_sampling::{
    peptide::Peptide,
    fivering::Fivering,
    sixring::Sixring,
};

mod formalism;
use formalism::{
    CP5, CP6, AS, Pdb, Xyz, SP
};

mod inversion;
use inversion::{
    sixring,
    fivering,
    write_file::{write_to_pdb, write_to_xyz},
};


/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
/// This is the name of the module
#[pymodule]
fn puckepy(py: Python, m: &PyModule) -> PyResult<()> {

    // Add geometry functions to the public API
    let geom_sub_module = PyModule::new(py, "geometry")?;
    geom_sub_module.add_function(wrap_pyfunction!(dihedral, geom_sub_module)?)?;
    geom_sub_module.add_function(wrap_pyfunction!(bondangle, geom_sub_module)?)?;
    geom_sub_module.add_function(wrap_pyfunction!(bondlength, geom_sub_module)?)?;

    // Add conformational sampling methods to the public API
    let cs_module = PyModule::new(py, "confsampling")?;
    cs_module.add_class::<Peptide>()?;
    cs_module.add_class::<Fivering>()?;
    cs_module.add_class::<Sixring>()?;

    // Add formalisms to the public API
    let form_module = PyModule::new(py, "formalism")?;
    form_module.add_class::<CP5>()?;
    form_module.add_class::<CP6>()?;
    form_module.add_class::<AS>()?;
    form_module.add_class::<SP>()?;
    form_module.add_class::<Pdb>()?;
    form_module.add_class::<Xyz>()?;

    // Add inversion from formalism to molecule to the public API
    let inv_module = PyModule::new(py, "inversion")?;
    inv_module.add_function(wrap_pyfunction!(sixring::invert_sixring, inv_module)?)?;
    inv_module.add_function(wrap_pyfunction!(fivering::invert_fivering, inv_module)?)?;
    inv_module.add_function(wrap_pyfunction!(write_to_pdb, inv_module)?)?;
    inv_module.add_function(wrap_pyfunction!(write_to_xyz, inv_module)?)?;

    // Append submodule to root module
    m.add_submodule(geom_sub_module)?;
    m.add_submodule(cs_module)?;
    m.add_submodule(form_module)?;
    m.add_submodule(inv_module)?;
    Ok(())

}

