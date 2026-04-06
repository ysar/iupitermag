use crate::currentsheet::PyCurrentSheetField;
use crate::internal::PyInternalField;
use iupitermag::trace;
use numpy::{IntoPyArray, PyReadonlyArray2};
use pyo3::{pyfunction, types::PyList, Bound, PyResult, Python};
use std::f64;

#[pyfunction]
pub fn trace_field_to_planet<'py>(
    py: Python<'py>,
    positions: PyReadonlyArray2<f64>,
    internal_field: Bound<'py, PyInternalField>,
    currentsheet_field: Bound<'py, PyCurrentSheetField>,
) -> PyResult<Bound<'py, PyList>> {
    let pos_arr = positions.as_array();

    let mut traces = vec![];

    let internal = internal_field.borrow();
    let currentsheet = currentsheet_field.borrow();

    for pos in pos_arr.rows() {
        traces.push(
            trace::trace_field_to_planet(pos.to_owned(), &internal.field, &currentsheet.field)
                .into_pyarray(py),
        )
    }
    PyList::new(py, traces)
}
