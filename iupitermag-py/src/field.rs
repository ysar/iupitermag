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
            /// Calculate the field at a point and return to Python.
            pub fn calc_field<'py>(
                &self,
                py: Python<'py>,
                r: f64,
                theta: f64,
                phi: f64,
            ) -> Bound<'py, PyArray1<f64>> {
                self.field.calc_field(r, theta, phi).into_pyarray(py)
            }

            /// Calculate the field at a point in cartesian coordinates and return result to Python.
            pub fn calc_field_xyz<'py>(
                &self,
                py: Python<'py>,
                x: f64,
                y: f64,
                z: f64,
            ) -> Bound<'py, PyArray1<f64>> {
                self.field.calc_field_xyz(x, y, z).into_pyarray(py)
            }

            /// Serial iterator into an array of positions to calculate field.
            pub fn map_calc_field<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> Bound<'py, PyArray2<f64>> {
                self.field
                    .map_calc_field(positions.as_array())
                    .into_pyarray(py)
            }

            /// Rayon iterator into an array of positions to calculate the field.
            pub fn parmap_calc_field<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> Bound<'py, PyArray2<f64>> {
                self.field
                    .parmap_calc_field(positions.as_array())
                    .into_pyarray(py)
            }

            /// Serial iterator into an array of positions to calculate field (cartesian).
            pub fn map_calc_field_xyz<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> Bound<'py, PyArray2<f64>> {
                self.field
                    .map_calc_field_xyz(positions.as_array())
                    .into_pyarray(py)
            }

            /// Rayon iterator into an array of positions to calculate the field (cartesian).
            pub fn parmap_calc_field_xyz<'py>(
                &self,
                py: Python<'py>,
                positions: PyReadonlyArray2<f64>,
            ) -> Bound<'py, PyArray2<f64>> {
                self.field
                    .parmap_calc_field_xyz(positions.as_array())
                    .into_pyarray(py)
            }
        }
    };
}
