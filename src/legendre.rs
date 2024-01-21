use nalgebra::DMatrix;

// Calculate the Gauss-normalized Legendre polynomials.
// Returns two ndarrays of size (n_degree + 1, n_degree + 1) containing the
// coefficient of the polynomial and its derivative (P(n, m) and dP(n, m))
// at order n and degree m.
pub fn assoc_legendre_poly(theta: f64, n_degree: usize) -> (DMatrix<f64>, DMatrix<f64>) {
    let n_size = n_degree + 1;

    let mut p = DMatrix::from_element(n_size, n_size, 0.);
    let mut dp = DMatrix::from_element(n_size, n_size, 0.);
    let mut k = DMatrix::from_element(n_size, n_size, 0.);

    let sintheta = theta.sin();
    let costheta = theta.cos();

    // We'd like to do the special cases manually since (i - 2) isn't defined.
    p[(0, 0)] = 1.;
    dp[(0, 0)] = 0.;

    p[(1, 1)] = sintheta * p[(0, 0)];
    p[(1, 0)] = costheta * p[(0, 0)];
    dp[(1, 1)] = sintheta * dp[(0, 0)] + costheta * p[(0, 0)];
    dp[(1, 0)] = costheta * dp[(0, 0)] - sintheta * p[(0, 0)];

    // General cases
    for i in 2..n_size {
        for j in 0..i {
            k[(i, j)] = ((i - 1).pow(2) - j.pow(2)) as f64 / ((2 * i - 1) * (2 * i - 3)) as f64;

            p[(i, j)] = costheta * p[(i - 1, j)] - k[(i, j)] * p[(i - 2, j)];

            dp[(i, j)] =
                costheta * dp[(i - 1, j)] - sintheta * p[(i - 1, j)] - k[(i, j)] * dp[(i - 2, j)];
        }

        p[(i, i)] = sintheta * p[(i - 1, i - 1)];

        dp[(i, i)] = sintheta * dp[(i - 1, i - 1)] + costheta * p[(i - 1, i - 1)];
    }
    (p, dp)
}
