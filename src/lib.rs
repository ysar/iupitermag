pub mod convert;
pub mod currentsheet;
pub mod field;
pub mod internal;
pub mod legendre;
pub mod trace;

use pyo3::{pymodule, types::PyModule, PyResult, Python};

#[pymodule]
fn iupitermag<'py>(_py: Python<'py>, m: &'py PyModule) -> PyResult<()> {
    m.add_class::<internal::PyInternalField>()?;
    m.add_class::<currentsheet::PyCurrentSheetField>()?;
    m.add_class::<trace::PyFieldTrace>()?;
    Ok(())
}
