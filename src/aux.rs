use numpy::ndarray::{Array1, ArrayView1};

pub fn convert_vec_xyz_to_rtp(arr: ArrayView1<f64>) -> Array1<f64> {
    let theta = (arr[2] / (arr[1].powi(2) + arr[0].powi(2)).sqrt()).acos();
    let phi = arr[1].atan2(arr[0]);

    let sintheta = theta.sin();
    let costheta = theta.cos();
    let sinphi = phi.sin();
    let cosphi = phi.cos();

    Array1::from_vec(vec![
        (arr[0] * sintheta * cosphi) + (arr[1] * sintheta * sinphi) + arr[2] * costheta,
        (arr[0] * costheta * cosphi) + (arr[1] * costheta * sinphi) - arr[2] * sintheta,
        -(arr[0] * sinphi) + (arr[1] * cosphi),
    ])
}

pub fn convert_vec_rtp_to_xyz(arr: ArrayView1<f64>, theta: &f64, phi: &f64) -> Array1<f64> {
    let sintheta = theta.sin();
    let costheta = theta.cos();
    let sinphi = phi.sin();
    let cosphi = phi.cos();

    Array1::from_vec(vec![
        (arr[0] * sintheta * cosphi) + (arr[1] * costheta * cosphi) - arr[2] * sinphi,
        (arr[0] * sintheta * sinphi) + (arr[1] * costheta * sinphi) + arr[2] * cosphi,
        (arr[0] * costheta) - (arr[1] * sintheta),
    ])
}
