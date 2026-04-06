use iupitermag::currentsheet::{CurrentSheetField, IntegrationType};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::{
    prelude::PyAnyMethods, pyclass, pymethods, types::IntoPyDict, types::PyDict, Bound, PyAny,
    PyResult, Python,
};
use std::collections::HashMap;

use crate::impl_field_methods;
use iupitermag::field::Field;

#[pyclass]
pub struct PyCurrentSheetField {
    pub field: CurrentSheetField,
}

#[pymethods]
impl PyCurrentSheetField {
    #[new]
    pub fn __init__(field_type: String, pyparams: Bound<'_, PyAny>, integration: String) -> Self {
        let params: HashMap<String, f64> = pyparams.extract().unwrap();

        let integration_type = match integration.to_lowercase().as_str() {
            "analytic" => IntegrationType::Analytic,
            "integral" => IntegrationType::Integral,
            _ => panic!("Unrecognized integration type. Allowed - analytic, integral ."),
        };

        PyCurrentSheetField {
            field: CurrentSheetField::new(field_type, Some(params), integration_type),
        }
    }
}

impl_field_methods!(PyCurrentSheetField);

#[pymethods]
impl PyCurrentSheetField {
    // Get parameters for a field model
    pub fn get_params<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyDict>> {
        self.field.get_params().into_py_dict(py)
    }
}
