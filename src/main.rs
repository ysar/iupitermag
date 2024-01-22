extern crate lazy_static;
extern crate nalgebra;

//use crate::legendre;
// use std::f64::consts::PI;
use nalgebra::Vector3;

pub mod coefficients;
pub mod field;
pub mod legendre;

// #[derive(Debug)]
// struct Node {
//     pos: [f64; 3],
//     b_val: [f64; 3],
// }

fn main() {
    // let some_node = Node {
    //     pos: [1., 2., 3.],
    //     b_val: [0., 2.2, -1.],
    // };

    // // Initialize a vector of nodes
    // let mut vec: Vec<Node> = Vec::new();

    // vec.push(some_node);

    // println!("{:?}", vec);

    // let theta: f64 = 16. * PI / 180.;
    // let (p, p_diff) = legendre::assoc_legendre_poly(theta, 3);

    // println!("{}", p);
    // println!("{}", p_diff);

    let type_field = coefficients::InternalField::JRM09;
    let (g, h) = type_field.get_coefficients();

    println!("{:?}", g[(1, 0)]);
    println!("{:?}", g);
    let val = field::calc_internal_field(10., 90. * 3.1415 / 180., 0., g, h);

    println!("{:?}, {:?}, {:?}", 10., 90. * 3.1415 / 180., 0.);
    println!("{:?}", val);
    // [[-48487.45378523, -12482.15605242605, -352312.72859993443]]}

    let csmodel = field::CurrentSheetType::Con2020;
    //let params = csmodel.get_current_sheet_params();
    let pos = Vector3::new(10., 0., 0.);
    let b_siii = csmodel.calc_current_sheet_field(pos);

    println!("{:?}", b_siii);
}
