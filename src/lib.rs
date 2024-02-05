pub mod aux;
pub mod field;
pub mod legendre;

// use pyo3::prelude::*;

// /// Formats the sum of two numbers as string.
// #[pyfunction]
// fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
//     Ok((a + b).to_string())
// }

// /// A Python module implemented in Rust.
// #[pymodule]
// fn iupitermag(_py: Python, m: &PyModule) -> PyResult<()> {
//     m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
//     Ok(())
// }

use numpy::ndarray::{s, Array2};
use numpy::{IntoPyArray, PyArray2, PyReadonlyArray2};
use pyo3::{pyfunction, pymodule, types::PyModule, wrap_pyfunction, PyResult, Python};

#[pyfunction]
#[pyo3(name = "get_planet_field")]
fn get_planet_field_py<'py>(
    py: Python<'py>,
    pos: PyReadonlyArray2<'py, f64>,
    internal_field_type: String,
    is_in_cartesian: bool,
    is_out_cartesian: bool,
    degree: usize,
) -> &'py PyArray2<f64> {
    let positions = pos.as_array();
    let num_points = positions.nrows();
    assert!(positions.ncols() == 3);

    let internal_field = match internal_field_type.as_str() {
        "JRM09" => field::InternalField::JRM09,
        "JRM33" => field::InternalField::JRM33,
        _ => panic!("Unknown field type."),
    };
    let (g, h) = internal_field.get_coefficients();

    let mut positions_rtp: Array2<f64>;

    let positions = if is_in_cartesian {
        positions_rtp = Array2::<f64>::zeros((num_points, 3));

        let x = positions.slice(s![.., 0]);
        let y = positions.slice(s![.., 1]);
        let z = positions.slice(s![.., 2]);

        for i in 0..num_points {
            positions_rtp[[i, 0]] = (x[i].powi(2) + y[i].powi(2) + z[i].powi(2)).sqrt();
            positions_rtp[[i, 1]] = (z[i] / positions_rtp[[i, 0]]).acos();
            positions_rtp[[i, 2]] = y[i].atan2(x[i]);
        }
        positions_rtp.view()
    } else {
        positions
    };

    let mut result_arr = Array2::<f64>::zeros((num_points, 3));

    for i in 0..num_points {
        let r = positions[[i, 0]];
        let theta = positions[[i, 1]];
        let phi = positions[[i, 2]];

        let val = internal_field.calc_internal_field(r, theta, phi, g.view(), h.view(), degree);

        if is_out_cartesian {
            let val_xyz = aux::convert_vec_rtp_to_xyz(val.view(), &theta, &phi);

            (result_arr[[i, 0]], result_arr[[i, 1]], result_arr[[i, 2]]) =
                (val_xyz[0], val_xyz[1], val_xyz[2]);
        } else {
            (result_arr[[i, 0]], result_arr[[i, 1]], result_arr[[i, 2]]) = (val[0], val[1], val[2]);
        }
    }

    result_arr.into_pyarray(py)
}

#[pymodule]
fn iupitermag<'py>(_py: Python<'py>, m: &'py PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_planet_field_py, m)?)?;
    Ok(())
}
// // example using immutable borrows producing a new array
// fn axpy(a: f64, x: ArrayViewD<'_, f64>, y: ArrayViewD<'_, f64>) -> ArrayD<f64> {
//     a * &x + &y
// }

// // example using a mutable borrow to modify an array in-place
// fn mult(a: f64, mut x: ArrayViewMutD<'_, f64>) {
//     x *= a;
// }

// // wrapper of `axpy`
// #[pyfn(m)]
// #[pyo3(name = "axpy")]
// fn axpy_py<'py>(
//     py: Python<'py>,
//     a: f64,
//     x: PyReadonlyArrayDyn<'py, f64>,
//     y: PyReadonlyArrayDyn<'py, f64>,
// ) -> &'py PyArrayDyn<f64> {
//     let x = x.as_array();
//     let y = y.as_array();
//     let z = axpy(a, x, y);
//     z.into_pyarray(py)
// }

// // wrapper of `mult`
// #[pyfn(m)]
// #[pyo3(name = "mult")]
// fn mult_py<'py>(a: f64, x: &'py PyArrayDyn<f64>) {
//     let x = unsafe { x.as_array_mut() };
//     mult(a, x);
// }
