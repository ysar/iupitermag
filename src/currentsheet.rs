use ndarray::Array1;
use numpy::{IntoPyArray, PyArray1, PyArray2, PyReadonlyArray2};
use pyo3::{
    prelude::PyAnyMethods, pyclass, pymethods, types::IntoPyDict, types::PyDict, Bound, PyAny,
    PyResult, Python,
};
use std::collections::HashMap;
use std::f64::consts::PI;

use crate::convert;
use crate::field::Field;
use crate::impl_field_methods;

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

#[derive(Clone)]
pub struct CurrentSheetField {
    r_0: f64,
    r_1: f64,
    d: f64,
    mu0_i_2: f64,
    theta_d: f64,
    phi_d: f64,
    integration_type: IntegrationType,
}

impl CurrentSheetField {
    pub fn new(
        field_type: String,
        params: Option<HashMap<String, f64>>,
        integration_type: IntegrationType,
    ) -> Self {
        match field_type.as_str() {
            "CON2020" => CurrentSheetField {
                r_0: 7.8,
                r_1: 51.4,
                d: 3.6,
                mu0_i_2: 139.6,
                theta_d: 9.3 * PI / 180.,
                phi_d: 204.2 * PI / 180.,
                integration_type,
            },
            "Custom" => {
                let _params = params.expect("params required for Custom field type.");

                let param_names = ["r_0", "r_1", "d", "mu0_i_2", "theta_d", "phi_d"];

                for key in param_names {
                    assert!(_params.contains_key(key), "Missing param - {:}", key);
                }

                CurrentSheetField {
                    r_0: _params["r_0"],
                    r_1: _params["r_1"],
                    d: _params["d"],
                    mu0_i_2: _params["mu0_i_2"],
                    theta_d: _params["theta_d"],
                    phi_d: _params["phi_d"],
                    integration_type,
                }
            }
            _ => panic!("Unknown field_type: Supported (CON2020, Custom)"),
        }
    }

    fn _calc_field(&self, rho: f64, z: f64, a: f64) -> (f64, f64) {
        match self.integration_type {
            IntegrationType::Analytic => self._calc_field_analytic(rho, z, a),
            IntegrationType::Integral => self._calc_field_integral(rho, z, a),
        }
    }

    fn _calc_field_analytic(&self, rho: f64, z: f64, a: f64) -> (f64, f64) {
        let b_rho: f64;
        let b_z: f64;

        let m_neg = z - self.d;
        let m_pos = z + self.d;

        if rho > a {
            let z_star = if z.abs() <= self.d.abs() {
                z
            } else {
                z.signum() * self.d
            };

            let n_neg = (rho.powi(2) + m_neg.powi(2)).sqrt();
            let n_pos = (rho.powi(2) + m_pos.powi(2)).sqrt();

            b_rho = self.mu0_i_2
                * (1. / rho * (n_neg - n_pos)
                    + rho * a.powi(2) / 4. * (1. / n_pos.powi(3) - 1. / n_neg.powi(3))
                    + 2. / rho * z_star);

            b_z = self.mu0_i_2
                * (((m_pos + n_pos) / (m_neg + n_neg)).ln()
                    + a.powi(2) / 4. * (m_pos / n_pos.powi(3) - m_neg / n_neg.powi(3)));
        } else {
            let n_neg = (a.powi(2) + m_neg.powi(2)).sqrt();
            let n_pos = (a.powi(2) + m_pos.powi(2)).sqrt();

            let p_neg = a.powi(2) - 2. * m_neg.powi(2);
            let p_pos = a.powi(2) - 2. * m_pos.powi(2);

            b_rho = self.mu0_i_2
                * (rho * 0.5 * (1. / n_neg - 1. / n_pos)
                    + rho.powi(3) / 16. * (p_neg / n_neg.powi(5) - p_pos / n_pos.powi(5)));

            b_z = self.mu0_i_2
                * (((m_pos + n_pos) / (m_neg + n_neg)).ln()
                    + rho.powi(2) / 4. * (m_pos / n_pos.powi(3) - m_neg / n_neg.powi(3)));
        }
        (b_rho, b_z)
    }

    fn _calc_field_integral(&self, rho: f64, z: f64, a: f64) -> (f64, f64) {
        let b_rho: f64;
        let b_z: f64;

        b_rho = 3.14;
        b_z = 3.14;

        (b_rho, b_z)
    }

    pub fn get_params(&self) -> HashMap<&str, f64> {
        HashMap::from([
            ("r_0", self.r_0),
            ("r_1", self.r_1),
            ("d", self.d),
            ("mu0_i_2", self.mu0_i_2),
            ("theta_d", self.theta_d),
            ("phi_d", self.phi_d),
        ])
    }
}

impl Field for CurrentSheetField {
    fn calc_field(&self, r: f64, theta: f64, phi: f64) -> Array1<f64> {
        // First we need to convert the input coordinates to cartesian
        let pos_in = Array1::from_vec(vec![r, theta, phi]);
        let pos_xyz = convert::pos_rtp_to_xyz(pos_in.view());

        // Then we convert the input coordinates from IAU to MAG frame
        let pos_xyz_mag = convert::vec_iau_to_mag(pos_xyz.view(), self.theta_d, self.phi_d);
        let r_mag = (pos_xyz_mag[0].powi(2) + pos_xyz_mag[1].powi(2)).sqrt();
        let z_mag = pos_xyz_mag[2];
        let phi_mag = pos_xyz_mag[1].atan2(pos_xyz_mag[0]);

        // Perform calculation in IAU frame and get (Brho, Bz)_MAG
        let (mut b_mag_rho, mut b_mag_z) = self._calc_field(r_mag, z_mag, self.r_0);

        // Calculate outer field to subtract if self.r_1 is present
        if !self.r_1.is_nan() {
            let (b_mag_rho_outer, b_mag_z_outer) = self._calc_field(r_mag, z_mag, self.r_1);

            b_mag_rho -= b_mag_rho_outer;
            b_mag_z -= b_mag_z_outer;
        }

        // Convert (Brho, Bz)_MAG to (Bx, By, Bz)_MAG
        let b_mag_x = b_mag_rho * phi_mag.cos();
        let b_mag_y = b_mag_rho * phi_mag.sin();
        let b_mag = Array1::from_vec(vec![b_mag_x, b_mag_y, b_mag_z]);

        // Convert (Bx, By, Bz)_MAG to (Bx, By, Bz)_IAU
        let b_iau = convert::vec_mag_to_iau(b_mag.view(), self.theta_d, self.phi_d);

        // Convert (Bx, By, Bz)_IAU to (Br, Btheta, Bphi)_IAU and return
        convert::vec_xyz_to_rtp(b_iau.view(), &theta, &phi)
    }
}

/// How to perform the integration to calculate the magnetic field contribution of the current sheet.
#[derive(Clone)]
pub enum IntegrationType {
    Analytic,
    Integral,
}
