use ndarray::{s, Array1, Array2, ArrayView1, ArrayView2, Zip};

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

pub fn convert_arr_rtp_to_xyz_serial(
    val: ArrayView2<f64>,
    positions: ArrayView2<f64>,
) -> Array2<f64> {
    let mut result = Array2::<f64>::zeros(val.dim());
    let mut _tmp = Array1::<f64>::zeros(3);
    for i in 0..val.nrows() {
        _tmp = convert_vec_rtp_to_xyz(val.slice(s![i, ..]), &positions[[i, 1]], &positions[[i, 2]]);
        (result[[i, 0]], result[[i, 1]], result[[i, 2]]) = (_tmp[0], _tmp[1], _tmp[2]);
    }
    result
}

pub fn convert_arr_rtp_to_xyz_parallel(
    val: ArrayView2<f64>,
    positions: ArrayView2<f64>,
) -> Array2<f64> {
    let mut result = Array2::<f64>::zeros(val.dim());
    Zip::from(result.rows_mut())
        .and(val.rows())
        .and(positions.rows())
        .par_for_each(|mut x, y, z| {
            let _tmp = convert_vec_rtp_to_xyz(y.view(), &z[1], &z[2]);
            (x[0], x[1], x[2]) = (_tmp[0], _tmp[1], _tmp[2]);
        });
    result
}
