use std::f64;

use crate::currentsheet::{CurrentSheetField, PyCurrentSheetField};
use crate::field::Field;
use crate::internal::{InternalField, PyInternalField};
use lazyivy::{RungeKutta, RungeKuttaMethod};
use ndarray::{Array1, Array2, ArrayView1};
use numpy::{IntoPyArray, PyReadonlyArray2};
use pyo3::{pyfunction, types::PyList, Bound, PyResult, Python};

const R_TRACE_MAXIMUM: f64 = 200.;

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
            _trace_field_to_planet(pos.to_owned(), &internal.field, &currentsheet.field)
                .into_pyarray(py),
        )
    }
    PyList::new(py, traces)
}

fn _trace_field_to_planet(
    start_position: Array1<f64>,
    internal_field: &InternalField,
    currentsheet_field: &CurrentSheetField,
) -> Array2<f64> {
    let planet_field = PlanetField {
        internal_field: internal_field.clone(),
        currentsheet_field: currentsheet_field.clone(),
    };

    // Also keep the predicate referenced (no capture).
    let _inside = is_inside_jupiter(start_position.view());

    let absolute_tol = Array1::from_vec(vec![1.0e-4, 1.0e-4, 1.0e-4]);
    let relative_tol = Array1::from_vec(vec![1.0e-4, 1.0e-4, 1.0e-4]);

    let integrator = RungeKutta::builder(
        |_, p, mut val| val.assign(&(calc_b_unit_vector(&planet_field, p))),
        |_, p| is_inside_jupiter(p),
    )
    .initial_condition(0., start_position.clone())
    .initial_step_size(0.025)
    .method(RungeKuttaMethod::DormandPrince, true) // `true` for adaptive step-size
    .tolerances(absolute_tol.clone(), relative_tol.clone())
    .set_max_step_size(0.25)
    .build()
    .unwrap();

    let integrator_inverse = RungeKutta::builder(
        |_, p, mut val| val.assign(&(calc_b_unit_vector_inverse(&planet_field, p))),
        |_, p| is_inside_jupiter(p),
    )
    .initial_condition(0., start_position.clone())
    .initial_step_size(0.025)
    .method(RungeKuttaMethod::DormandPrince, true) // `true` for adaptive step-size
    .tolerances(absolute_tol, relative_tol)
    .set_max_step_size(0.25)
    .build()
    .unwrap();

    let trace_pos = integrator.map(|(_, p)| p).collect::<Vec<Array1<f64>>>();
    let trace_neg = integrator_inverse
        .map(|(_, p)| p)
        .collect::<Vec<Array1<f64>>>();

    let num_points_pos = trace_pos.len();
    let num_points_neg = trace_neg.len();

    let mut result: Array2<f64> =
        Array2::from_elem((num_points_pos + num_points_neg + 1, 3), f64::NAN);

    for (i, pos) in trace_neg.iter().rev().enumerate() {
        result.row_mut(i).assign(pos)
    }

    result.row_mut(num_points_neg).assign(&start_position);

    for (i, pos) in trace_pos.iter().enumerate() {
        result.row_mut(i + num_points_neg + 1).assign(pos);
    }

    result
}

/// Calculate the magnetic field unit vector at a point.
fn calc_b_unit_vector(field: &PlanetField, pos: ArrayView1<f64>) -> Array1<f64> {
    let b = field.calc_field_xyz(pos[0], pos[1], pos[2]);
    let b_mag = b.dot(&b).sqrt();
    b / b_mag
}

/// Calculate the inverse of the magnetic field unit vector.
fn calc_b_unit_vector_inverse(field: &PlanetField, pos: ArrayView1<f64>) -> Array1<f64> {
    calc_b_unit_vector(field, pos) * (-1.0)
}

/// Check if point is inside Jupiter's ellipsoid or outside the bounds of tracing.
fn is_inside_jupiter(pos: ArrayView1<f64>) -> bool {
    let f: f64 = 1.0 / 15.4; // Polar flattening of Jupiter ellipsoid.
    let a: f64 = 1.0;
    let c: f64 = (1.0 - f) * a;
    let r_ellipsoid_norm =
        pos[0].powi(2) / a.powi(2) + pos[1].powi(2) / a.powi(2) + pos[2].powi(2) / c.powi(2);

    r_ellipsoid_norm < 1. || pos.dot(&pos).sqrt() > R_TRACE_MAXIMUM
}

struct PlanetField {
    internal_field: InternalField,
    currentsheet_field: CurrentSheetField,
}

impl Field for PlanetField {
    fn calc_field(&self, r: f64, theta: f64, phi: f64) -> Array1<f64> {
        let b_internal = self.internal_field.calc_field(r, theta, phi);
        let b_currentsheet = self.currentsheet_field.calc_field(r, theta, phi);
        b_internal + b_currentsheet
    }
}
