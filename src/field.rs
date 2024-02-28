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

// https://github.com/stm32-rs/stm32f0xx-hal/blob/master/src/gpio.rs
// Need to add trait behaviour using macros