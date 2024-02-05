use numpy::ndarray::Array2;

// Calculate the Gauss-normalized Legendre polynomials.
// Returns two numpy::ndarrays of size (n_degree + 1, n_degree + 1) containing the
// coefficient of the polynomial and its derivative (P(n, m) and dP(n, m))
// at order n and degree m.
pub fn assoc_legendre_poly(theta: &f64, n_degree: &usize) -> (Array2<f64>, Array2<f64>) {
    let n_size = n_degree + 1;

    let mut p = Array2::<f64>::zeros((n_size, n_size));
    let mut dp = Array2::<f64>::zeros((n_size, n_size));
    let mut k = Array2::<f64>::zeros((n_size, n_size));

    let sintheta = theta.sin();
    let costheta = theta.cos();

    // We'd like to do the special cases manually since (i - 2) isn't defined.
    p[[0, 0]] = 1.;
    dp[[0, 0]] = 0.;

    p[[1, 1]] = sintheta * p[[0, 0]];
    p[[1, 0]] = costheta * p[[0, 0]];
    dp[[1, 1]] = sintheta * dp[[0, 0]] + costheta * p[[0, 0]];
    dp[[1, 0]] = costheta * dp[[0, 0]] - sintheta * p[[0, 0]];

    // General cases
    for i in 2..n_size {
        for j in 0..i {
            k[[i, j]] = ((i - 1).pow(2) - j.pow(2)) as f64 / ((2 * i - 1) * (2 * i - 3)) as f64;

            p[[i, j]] = costheta * p[[i - 1, j]] - k[[i, j]] * p[[i - 2, j]];

            dp[[i, j]] =
                costheta * dp[[i - 1, j]] - sintheta * p[[i - 1, j]] - k[[i, j]] * dp[[i - 2, j]];
        }

        p[[i, i]] = sintheta * p[[i - 1, i - 1]];

        dp[[i, i]] = sintheta * dp[[i - 1, i - 1]] + costheta * p[[i - 1, i - 1]];
    }
    (p, dp)
}

pub fn schmidt_semi_normalization_constants(degree: &usize) -> Array2<f64> {
    let mut s = Array2::<f64>::from_elem((degree + 1, degree + 1), 1.);

    for i in 1..degree + 1 {
        s[[i, 0]] = s[[i - 1, 0]] * (2 * i - 1) as f64 / i as f64;
        s[[i, 1]] = s[[i, 0]] * (i as f64 * 2. / (i + 1) as f64).sqrt();

        for j in 2..i + 1 {
            s[[i, j]] = s[[i, j - 1]] * ((i - j + 1) as f64 / (i + j) as f64).sqrt();
        }
    }
    s // Returns here
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_assoc_legendre_poly() {
        use crate::legendre;
        use numpy::ndarray::Array2;

        let xarr = vec![0., 0.2, 0.4, 0.6, 0.8, 1.0];
        let val_check = Array2::from_shape_vec(
            (6, 2),
            vec![0., 1., 0.2, 0.9798, 0.4, 0.9165, 0.6, 0.8, 0.8, 0.6, 1., 0.],
        )
        .unwrap();
        println!("{:?}", val_check);
        let n = 1;
        for iarr in 0..6 {
            println!("x = {:}", xarr[iarr]);
            println!("val = {:?}", val_check[[iarr, 1]]);
            let (p, _dp) = legendre::assoc_legendre_poly(&((xarr[iarr] as f64).acos()), &n);
            let s = legendre::schmidt_semi_normalization_constants(&n);

            for j in 0..2 {
                assert!(
                    val_check[[iarr, j]] - p[[1, j]] * s[[1, j]] < 1e-5,
                    "Assoc. Legendre Poly. Test Failed: \n Calculated {:?}, Expected {:?}",
                    p[[1, j]] * s[[1, j]],
                    val_check[[iarr, j]]
                );
            }
        }
    }
}
