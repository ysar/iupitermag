use crate::convert;
use ndarray::{Array1, Array2, ArrayView2, Zip};

pub trait Field {
    fn calc_field(&self, r: f64, theta: f64, phi: f64) -> Array1<f64>;

    /// Return the field calculated at a point (x, y, z) in cartesian coordinates
    /// (Bx, By, Bz).
    fn calc_field_xyz(&self, x: f64, y: f64, z: f64) -> Array1<f64> {
        let pos_xyz = Array1::from(vec![x, y, z]);
        let pos_rtp = convert::pos_xyz_to_rtp(pos_xyz.view());
        let b_rtp = self.calc_field(pos_rtp[0], pos_rtp[1], pos_rtp[2]);
        convert::vec_rtp_to_xyz(b_rtp.view(), &pos_rtp[1], &pos_rtp[2])
    }

    fn map_calc_field(&self, positions: ArrayView2<f64>) -> Array2<f64> {
        let mut result = Array2::<f64>::zeros((positions.nrows(), 3));

        Zip::from(result.rows_mut())
            .and(positions.rows())
            .for_each(|mut x, y| {
                x.assign(&self.calc_field(y[0], y[1], y[2]));
            });

        result
    }

    fn parmap_calc_field(&self, positions: ArrayView2<f64>) -> Array2<f64>
    where
        Self: Sync,
    {
        let mut result = Array2::<f64>::zeros((positions.nrows(), 3));

        Zip::from(result.rows_mut())
            .and(positions.rows())
            .par_for_each(|mut x, y| {
                x.assign(&self.calc_field(y[0], y[1], y[2]));
            });

        result
    }
}

// The macro below implements pymethods for various pyclasses throughout the
// module. This is the only way that I know of to implement shared
// functionality across pyclasses. Each class contains a field of a different
// type, so it is not possible to subclass from a base class. Likewise, it was
// not possible to implement a trait since the trait definition
// has no knowledge of the struct fields beforehand.

#[macro_export]
macro_rules! impl_field_methods {
    ($cls:ident) => {
        #[pymethods]
        impl $cls {
            // Calculate the field at a point
            pub fn calc_field<'py>(
                &self,
                py: Python<'py>,
                r: f64,
                theta: f64,
                phi: f64,
            ) -> Bound<'py, PyArray1<f64>> {
                self.field.calc_field(r, theta, phi).into_pyarray(py)
            }

            // Calculate the field at a point in cartesian coordinates.
            pub fn calc_field_xyz<'py>(
                &self,
                py: Python<'py>,
                x: f64,
                y: f64,
                z: f64,
            ) -> Bound<'py, PyArray1<f64>> {
                self.field.calc_field_xyz(x, y, z).into_pyarray(py)
            }

            // Serial iterator into an array of positions to calculate field
            pub fn map_calc_field<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> Bound<'py, PyArray2<f64>> {
                self.field
                    .map_calc_field(positions.as_array())
                    .into_pyarray(py)
            }

            // Rayon iterator into an array of positions to calculate field
            pub fn parmap_calc_field<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> Bound<'py, PyArray2<f64>> {
                self.field
                    .parmap_calc_field(positions.as_array())
                    .into_pyarray(py)
            }
        }
    };
}
