use std::collections::HashMap;
use std::f64::consts::PI;
use ndarray::{Array1};
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::{pyclass, pymethods, Python};

use crate::field::{Field, PyField};
use crate::convert;


#[pyclass]
pub struct PyCurrentSheetField {
    _f: CurrentSheetField,
}

#[pymethods]
impl PyCurrentSheetField {
    #[new]
    pub fn __init__(field_type: &str, params: Option<HashMap<String, f64>>) -> Self {
        PyCurrentSheetField {
            _f: CurrentSheetField::new(field_type, params),
        }
    }
}

impl PyField for PyCurrentSheetField {

}


pub struct CurrentSheetField {
    r_0: f64,
    r_1: f64,
    d: f64,
    mu0_i_2: f64,
    theta_d: f64,
    phi_d: f64,
}

impl CurrentSheetField {
    pub fn new(field_type: &str, params: Option<HashMap<String, f64>>) -> Self {
        match field_type {
            "CON2020" => CurrentSheetField {
                r_0: 7.8,
                r_1: 51.4,
                d: 3.6,
                mu0_i_2: 139.6,
                theta_d: 9.3 * PI / 180.,
                phi_d: 204.2 * PI / 180.,
            },
            "Custom" => {
                let _params = params.expect(
                    "params HashMap expected for CurrentSheetField with field_type 'Custom' ",
                );

                CurrentSheetField {
                    r_0: _params["r_0"],
                    r_1: _params["r_1"],
                    d: _params["d"],
                    mu0_i_2: _params["mu0_i_2"],
                    theta_d: _params["theta_d"],
                    phi_d: _params["phi_d"],
                }
            }
            _ => panic!("Unknown field_type: Supported (CON2020, Custom)"),
        }
    }

    fn _calc_field(&self, rho: f64, z: f64) -> (f64, f64) {

        let z_star: f64;
        if z.abs() <= self.d.abs() {
            z_star = z;
        } else {
            z_star = z.signum() * self.d;
        }

        let x_neg = (rho.powi(2) + (z - self.d).powi(2)).sqrt();
        let x_pos = (rho.powi(2) + (z + self.d).powi(2)).sqrt();

        let y_neg = z - self.d;
        let y_pos = z + self.d;

        let b_rho = self.mu0_i_2
            * (1. / rho * (x_neg - x_pos)
                + rho * self.r_0.powi(2) / 4. * (1. / x_pos.powi(3) - 1. / x_neg.powi(3))
                + 2. / rho * z_star);
        let b_z = self.mu0_i_2
            * (((y_pos + x_pos) / (y_neg + x_neg)).ln()
                + self.r_0.powi(2) / 4. * (y_pos / x_pos.powi(3) - y_neg / x_neg.powi(3)));

        (b_rho, b_z)
    }
}

impl Field for CurrentSheetField {
    fn calc_field(&self, r: f64, theta: f64, phi: f64) -> Array1<f64> {

        let cosphi = phi.cos();
        let sinphi = phi.sin();

        // First we need to convert the input coordinates to cartesian
        let pos_in = Array1::from_vec(vec![r, theta, phi]);
        let pos_xyz = convert::pos_rtp_to_xyz(pos_in.view());

        // Then we convert the input coordinates from IAU to MAG frame
        let pos_xyz_mag = convert::vec_iau_to_mag(pos_xyz.view(), self.theta_d, self.phi_d);

        // Perform calculation in IAU frame and get (Brho, Bz)_MAG
        let (b_mag_rho, b_mag_z) = self._calc_field(
            (pos_xyz_mag[0].powi(2) + pos_xyz_mag[1].powi(2)).sqrt(),
            pos_xyz_mag[2]
        );

        // Convert (Brho, Bz)_MAG to (Bx, By, Bz)_MAG
        let b_mag_x = b_mag_rho * cosphi;
        let b_mag_y = b_mag_rho * sinphi;
        let b_mag = Array1::from_vec(vec![b_mag_x, b_mag_y, b_mag_z]);

        // Convert (Bx, By, Bz)_MAG to (Bx, By, Bz)_IAU
        let b_iau = convert::vec_mag_to_iau(b_mag.view(), self.theta_d, self.phi_d);

        // Convert (Bx, By, Bz)_IAU to (Br, Btheta, Bphi)_IAU and return
        convert::vec_xyz_to_rtp(b_iau.view(), &theta, &phi)
    }
}