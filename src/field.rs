use ndarray::{Array1, Array2, ArrayView2, Zip};
use pyo3::{pyclass, pymethods, Python};

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

pub trait PyField {

    fn calc_field<'py>(
        &self,
        py: Python<'py>,
        r: f64,
        theta: f64,
        phi: f64,
    ) -> &'py PyArray1<f64> {
        let result = self._f.calc_field(r, theta, phi);
        result.into_pyarray(py)
    }

    fn get_params<'py>(&self, py: Python<'py>) {
        todo!();
    }

    fn map_calc_field<'py>(
        &self,
        py: Python<'py>,
        positions: PyReadonlyArray2<f64>,
    ) -> &'py PyArray2<f64> {
        self._f
            .map_calc_field(positions.as_array())
            .into_pyarray(py)
    }

    fn parmap_calc_field<'py>(
        &self,
        py: Python<'py>,
        positions: PyReadonlyArray2<f64>,
    ) -> &'py PyArray2<f64> {
        self._f
            .parmap_calc_field(positions.as_array())
            .into_pyarray(py)
    }

}