pub mod convert;
pub mod currentsheet;
pub mod field;
pub mod internal;
pub mod legendre;
pub mod trace;

use pyo3::pymodule;

#[pymodule]
mod _core {
    #[pymodule_export]
    pub use crate::internal::PyInternalField;

    #[pymodule_export]
    pub use crate::currentsheet::PyCurrentSheetField;

    #[pymodule_export]
    pub use crate::trace::trace_field_to_planet;
}
