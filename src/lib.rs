pub mod convert;
pub mod currentsheet;
pub mod field;
pub mod internal;
pub mod legendre;

use pyo3::{pymodule, types::PyModule, PyResult, Python};

#[pymodule]
fn iupitermag<'py>(_py: Python<'py>, m: &'py PyModule) -> PyResult<()> {
    m.add_class::<internal::PyInternalField>()?;
    Ok(())
}
