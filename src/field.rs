use crate::legendre;
use nalgebra::{DMatrix, Rotation3, Vector3};
use std::f64::consts::PI;
// use ndarray::{array, Array, Ix1, Ix2};

pub fn calc_internal_field(
    r: f64,
    theta: f64,
    phi: f64,
    g: DMatrix<f64>,
    h: DMatrix<f64>,
) -> Vector3<f64> {
    //

    let mut b_r: f64 = 0.;
    let mut b_theta: f64 = 0.;
    let mut b_phi: f64 = 0.;

    let a: f64 = 1. / r;

    assert!((g.nrows() == g.ncols()) && (h.nrows() == h.ncols()));

    let n_degree = g.nrows() - 1;
    let (p, dp) = legendre::assoc_legendre_poly(theta, n_degree);

    let sintheta: f64 = theta.sin();
    let costheta: f64 = theta.cos();
    let sinphi: f64 = phi.sin();
    let cosphi: f64 = phi.cos();

    let inv_sintheta: f64 = if sintheta.abs() > 1e-9 {
        1. / sintheta
    } else {
        0.
    };

    let mut sinphi_prev: f64;
    let mut cosphi_prev: f64;
    let mut sin_mphi: f64;
    let mut cos_mphi: f64;

    for i in 0..n_degree + 1 {
        sinphi_prev = 0.;
        cosphi_prev = 1.;

        let i_i32 = i32::try_from(i).unwrap();

        b_r += a.powi(i_i32 + 2) * (i + 1) as f64 * p[(i, 0)] * g[(i, 0)];
        b_theta -= a.powi(i_i32 + 2) * dp[(i, 0)] * g[(i, 0)];

        for j in 1..i + 1 {
            sin_mphi = sinphi_prev * cosphi + cosphi_prev * sinphi;
            cos_mphi = cosphi_prev * cosphi - sinphi_prev * sinphi;

            b_r += a.powi(i_i32 + 2)
                * (i + 1) as f64
                * p[(i, j)]
                * (g[(i, j)] * cos_mphi + h[(i, j)] * sin_mphi);

            b_theta -=
                a.powi(i_i32 + 2) * dp[(i, j)] * (g[(i, j)] * cos_mphi + h[(i, j)] * sin_mphi);

            b_phi += inv_sintheta
                * a.powi(i_i32 + 2)
                * p[(i, j)]
                * j as f64
                * (g[(i, j)] * sin_mphi - h[(i, j)] * cos_mphi);

            sinphi_prev = sin_mphi;
            cosphi_prev = cos_mphi;
        }
    }

    // return array of bx, by, bz
    Vector3::new(
        (b_r * sintheta * cosphi) + (b_theta * costheta * cosphi) - b_phi * sinphi,
        (b_r * sintheta * sinphi) + (b_theta * costheta * sinphi) + b_phi * cosphi,
        (b_r * costheta) - (b_theta * sintheta),
    )
}

pub struct CurrentSheetParameters {
    r_0: f64,
    r_1: f64,
    d: f64,
    mu0_i_2: f64,
    theta_d: f64,
    phi_d: f64,
}

pub enum CurrentSheetType {
    Con2020,
}

impl CurrentSheetType {
    pub fn get_current_sheet_params(&self) -> CurrentSheetParameters {
        match &self {
            CurrentSheetType::Con2020 => CurrentSheetParameters {
                r_0: 7.8,
                r_1: 51.4,
                d: 3.6,
                mu0_i_2: 139.6,
                theta_d: 9.3 * PI / 180.,
                phi_d: 204.2 * PI / 180.,
            },
        }
    }

    pub fn calc_current_sheet_field_cylind(
        &self,
        rho: f64,
        z: f64,
        params: &CurrentSheetParameters,
    ) -> (f64, f64) {
        if ((rho - params.r_0) <= 2.) && (z.abs() <= params.d * 1.5) {
            let denom_neg: f64 = (params.r_0.powi(2) + (z - params.d).powi(2)).sqrt();
            let denom_pos: f64 = (params.r_0.powi(2) + (z + params.d).powi(2)).sqrt();
            let numer_neg: f64 = params.r_0.powi(2) - 2. * (z - params.d).powi(2);
            let numer_pos: f64 = params.r_0.powi(2) - 2. * (z + params.d).powi(2);

            // b_r
            (
                params.mu0_i_2
                    * ((rho / 2.) * (1. / denom_neg - 1. / denom_pos)
                        + rho.powi(3) / 16.
                            * (numer_neg / denom_neg.powf(2.5) - numer_pos / denom_pos.powf(2.5))),
                // b_z
                params.mu0_i_2
                    * (((z + params.d + denom_pos) / (z - params.d + denom_neg)).ln()
                        + rho.powi(2) / 4.
                            * ((z + params.d) / denom_pos.powf(1.5)
                                - (z - params.d) / denom_neg.powf(1.5))),
            )
        } else {
            let _cond_statement = if z.abs() <= params.d {
                z * (2. / rho)
            } else {
                z.signum() * params.d * (2. / rho)
            };

            let denom_neg: f64 = (rho.powi(2) + (z - params.d).powi(2)).sqrt();
            let denom_pos: f64 = (rho.powi(2) + (z + params.d).powi(2)).sqrt();

            // b_r
            (
                params.mu0_i_2 * (1. / rho) * (denom_neg - denom_pos)
                    + rho * params.r_0.powi(2) / 4. * (1. / denom_pos - 1. / denom_neg)
                    + _cond_statement,
                // b_z
                params.mu0_i_2
                    * (((z + params.d + denom_pos) / (z - params.d + denom_neg)).ln()
                        + params.r_0.powi(2) / 4.
                            * ((z + params.d) / denom_pos.powf(1.5)
                                - (z - params.d) / denom_neg.powf(1.5))),
            )
        }
    }

    pub fn calc_current_sheet_field(&self, pos: Vector3<f64>) -> Vector3<f64> {
        let params = self.get_current_sheet_params();
        let pos_mag = rotate_vector_iau2mag(pos, &params);

        let rho: f64 = (pos_mag[0].powi(2) + pos_mag[1].powi(2)).sqrt();
        let z: f64 = pos_mag[2];

        let (b_rho, b_z) = self.calc_current_sheet_field_cylind(rho, z, &params);

        let phi = pos_mag[1].atan2(pos_mag[0]);
        let b_mag = Vector3::new(b_rho * phi.cos(), b_rho * phi.sin(), b_z);
        rotate_vector_mag2iau(b_mag, &params)
    }
}

fn rotate_vector_iau2mag(v_in: Vector3<f64>, params: &CurrentSheetParameters) -> Vector3<f64> {
    let axis_y = Vector3::y_axis();
    let axis_z = Vector3::z_axis();

    Rotation3::from_axis_angle(&axis_y, -params.theta_d)
        * Rotation3::from_axis_angle(&axis_z, -params.phi_d)
        * v_in
}

fn rotate_vector_mag2iau(v_in: Vector3<f64>, params: &CurrentSheetParameters) -> Vector3<f64> {
    let axis_y = Vector3::y_axis();
    let axis_z = Vector3::z_axis();

    Rotation3::from_axis_angle(&axis_z, params.phi_d)
        * Rotation3::from_axis_angle(&axis_y, params.theta_d)
        * v_in
}
