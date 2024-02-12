pub mod aux;
pub mod internal;
pub mod legendre;

use ndarray::{s, Array2};
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
        "JRM09" => internal::InternalField::JRM09,
        "JRM33" => internal::InternalField::JRM33,
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

    let mut result_arr: Array2<f64>;

    if num_points > 1000 {
        result_arr = internal::calc_arr_internal_field_parallel(
            internal_field,
            positions.view(),
            g.view(),
            h.view(),
            &degree,
        );
        if is_out_cartesian {
            result_arr = aux::convert_arr_rtp_to_xyz_parallel(result_arr.view(), positions.view());
        }
    } else {
        result_arr = internal::calc_arr_internal_field_serial(
            internal_field,
            positions.view(),
            g.view(),
            h.view(),
            &degree,
        );
        if is_out_cartesian {
            result_arr = aux::convert_arr_rtp_to_xyz_serial(result_arr.view(), positions.view());
        }
    }
    result_arr.into_pyarray(py)
}

#[pymodule]
fn iupitermag<'py>(_py: Python<'py>, m: &'py PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(get_planet_field_py, m)?)?;
    Ok(())
}
