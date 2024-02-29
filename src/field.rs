use ndarray::{Array1, Array2, ArrayView2, Zip};

pub trait Field {
    fn calc_field(&self, r: f64, theta: f64, phi: f64) -> Array1<f64>;

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
            ) -> &'py PyArray1<f64> {
                let result = self._f.calc_field(r, theta, phi);
                result.into_pyarray(py)
            }
            
            // Get parameters for a field model
            pub fn get_params<'py>(&self, py: Python<'py>) {
                todo!();
            }
            
            // Serial iterator into an array of positions to calculate field
            pub fn map_calc_field<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> &'py PyArray2<f64> {
                self._f
                    .map_calc_field(positions.as_array())
                    .into_pyarray(py)
            }
            
            // Rayon iterator into an array of positions to calculate field
            pub fn parmap_calc_field<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> &'py PyArray2<f64> {
                self._f
                    .parmap_calc_field(positions.as_array())
                    .into_pyarray(py)
            }
        }
    };
}
