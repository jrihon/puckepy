use pyo3::prelude::*;

mod geometry;
use geometry::molecule_ops::{
    dihedral,
    bondangle,
    bondlength
};

mod conf_sampling;
use conf_sampling::{
    peptide::{Peptide, PeptideAxes},
    fivering::{Fivering, FiveringAxes},
    sixring::{Sixring, SixringAxes},
};

mod formalism;
use formalism::{
    cremerpople::{CP5, CP6},
    altonasund::AS,
    strausspickett::SP,
    moleculefile::{Pdb, 
                   Xyz,
                   write_to_pdb,
                   write_to_xyz
                    },
};

/// A Python module implemented in Rust. The name of this function must match
/// the `lib.name` setting in the `Cargo.toml`, else Python will not be able to
/// import the module.
/// This is the name of the module
#[pymodule]
fn puckepy(m: &Bound<'_, PyModule>) -> PyResult<()> {
//fn puckepy(py: Python, m: &PyModule) -> PyResult<()> {
    register_child_modules(m)?;
    Ok(())
}

fn register_child_modules(parent_module: &Bound<'_, PyModule>) -> PyResult<()> {
    //
    // Add geometry functions to the public API
    let geom_sub_module = PyModule::new_bound(parent_module.py(), "geometry")?;
    geom_sub_module.add_function(wrap_pyfunction!(dihedral, &geom_sub_module)?)?;
    geom_sub_module.add_function(wrap_pyfunction!(bondangle, &geom_sub_module)?)?;
    geom_sub_module.add_function(wrap_pyfunction!(bondlength, &geom_sub_module)?)?;

    // Add conformational sampling methods to the public API
    let cs_module = PyModule::new_bound(parent_module.py(), "confsampling")?;
    cs_module.add_class::<Peptide>()?;
    cs_module.add_class::<PeptideAxes>()?;
    cs_module.add_class::<Fivering>()?;
    cs_module.add_class::<FiveringAxes>()?;
    cs_module.add_class::<Sixring>()?;
    cs_module.add_class::<SixringAxes>()?;

    // Add formalisms to the public API
    let form_module = PyModule::new_bound(parent_module.py(), "formalism")?;
    form_module.add_class::<CP5>()?;
    form_module.add_class::<CP6>()?;
    form_module.add_class::<AS>()?;
    form_module.add_class::<SP>()?;
    form_module.add_class::<Pdb>()?;
    form_module.add_class::<Xyz>()?;
    form_module.add_function(wrap_pyfunction!(write_to_pdb, &form_module)?)?;
    form_module.add_function(wrap_pyfunction!(write_to_xyz, &form_module)?)?;

    // Append submodule to root module
    parent_module.add_submodule(&geom_sub_module)?;
    parent_module.add_submodule(&cs_module)?;
    parent_module.add_submodule(&form_module)?;
    Ok(())

}

