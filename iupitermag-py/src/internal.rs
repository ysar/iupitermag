use crate::impl_field_methods;
use iupitermag::field::Field;
use iupitermag::internal::InternalField;
use numpy::ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyArrayMethods, PyReadonlyArray2};
use pyo3::{pyclass, pymethods, Bound, Python};

// #[derive(Clone)]
#[pyclass]
pub struct PyInternalField {
    pub field: InternalField,
}

#[pymethods]
impl PyInternalField {
    #[new]
    pub fn __init__(
        field_type: &str,
        g_in: Option<PyReadonlyArray2<f64>>,
        h_in: Option<PyReadonlyArray2<f64>>,
        degree_in: Option<usize>,
    ) -> Self {
        let g: Option<Array2<f64>> = g_in.map(|x| x.to_owned_array());
        let h: Option<Array2<f64>> = h_in.map(|x| x.to_owned_array());

        PyInternalField {
            field: InternalField::new(field_type, g, h, degree_in),
        }
    }

    pub fn get_coefficients<'py>(
        &self,
        py: Python<'py>,
    ) -> (Bound<'py, PyArray2<f64>>, Bound<'py, PyArray2<f64>>) {
        let (g, h) = self.field.get_coefficients();
        (g.into_pyarray(py), h.into_pyarray(py))
    }
}

impl_field_methods!(PyInternalField);
