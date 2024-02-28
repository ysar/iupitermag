use ndarray::{Array1, Array2, ArrayView1};

// Converts a vector of cartesian coordinates to spherical coordinates
pub fn pos_xyz_to_rtp(arr: ArrayView1<f64>) -> Array1<f64> {
    let r = (arr[0].powi(2) + arr[1].powi(2) + arr[2].powi(2)).sqrt();
    Array1::from_vec(vec![
        r,
        (arr[2] / r).acos(),
        (arr[1]).atan2(arr[0],)
    ])
}

// Converts a vector of spherical coordinates to cartesian coordinates
pub fn pos_rtp_to_xyz(arr: ArrayView1<f64>) -> Array1<f64> {
    Array1::from_vec(vec![
        arr[0] * arr[1].sin() * arr[2].cos(),
        arr[0] * arr[1].sin() * arr[2].sin(),
        arr[0] * arr[1].cos(),
    ])
}

// Converts a vector from cartesian basis to spherical basis
pub fn vec_xyz_to_rtp(arr: ArrayView1<f64>, theta: &f64, phi: &f64) -> Array1<f64> {

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

// Converts a vector from spherical basis to cartesian basis
pub fn vec_rtp_to_xyz(arr: ArrayView1<f64>, theta: &f64, phi: &f64) -> Array1<f64> {

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

// Converts a cartesian vector in IAU frame to a cartesian vector in MAG frame
pub fn vec_iau_to_mag(arr: ArrayView1<f64>, theta_d: f64, phi_d: f64) -> Array1<f64> {

    // To get from IAU to MAG we need to rotate by -Phi_d along z, then 
    // rotate by -Theta_d along Y.

    let rot_matrix_1 = rot_matrix_z(-phi_d);
    let rot_matrix_2 = rot_matrix_y(-theta_d);
    
    rot_matrix_2.dot(&rot_matrix_1).dot(&arr)
}

// Converts a cartesian vector in MAG frame to a cartesian vector in IAU frame
pub fn vec_mag_to_iau(arr: ArrayView1<f64>, theta_d: f64, phi_d: f64) -> Array1<f64> {
    // To get from IAU to MAG we need to rotate by
    // Theta_d along Y
    // Phi_d along Z

    let rot_matrix_1 = rot_matrix_y(theta_d);
    let rot_matrix_2 = rot_matrix_z(phi_d);
    
    rot_matrix_2.dot(&rot_matrix_1).dot(&arr)
}

// Creates a matrix of rotation about the X axis
pub fn rot_matrix_x(angle: f64) -> Array2<f64> {
    Array2::<f64>::from_shape_vec((3, 3), vec![
        1., 0., 0.,
        0., angle.cos(), -angle.sin(),
        0., angle.sin(), angle.cos()
    ]).unwrap()
}

// Creates a matrix of rotation about the Y axis
pub fn rot_matrix_y(angle: f64) -> Array2<f64> {
    Array2::<f64>::from_shape_vec((3, 3), vec![
        angle.cos(), 0., angle.sin(),
        0., 1., 0.,
        -angle.sin(), 0., angle.cos(),
    ]).unwrap()
}

// Creates a matrix of rotation about the Z axis
pub fn rot_matrix_z(angle: f64) -> Array2<f64> {
    Array2::<f64>::from_shape_vec((3, 3), vec![
        angle.cos(), -angle.sin(), 0.,
        angle.sin(), angle.cos(), 0.,
        0., 0., 1.,
    ]).unwrap()
}