use crate::legendre;
use nalgebra::{DMatrix, Vector3};

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
