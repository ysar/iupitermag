extern crate lazy_static;
extern crate ndarray;

//use crate::legendre;
use std::f64::consts::PI;

pub mod coefficients;
pub mod legendre;
pub mod mag;

#[derive(Debug)]
struct Node {
    pos: [f64; 3],
    b_val: [f64; 3],
}

fn main() {
    let some_node = Node {
        pos: [1., 2., 3.],
        b_val: [0., 2.2, -1.],
    };

    // Initialize a vector of nodes
    let mut vec: Vec<Node> = Vec::new();

    vec.push(some_node);

    println!("{:?}", vec);

    let theta: f64 = 16. * PI / 180.;
    let (p, p_diff) = legendre::assoc_legendre_poly(theta, 3);

    println!("{}", p);
    println!("{}", p_diff);

    let type_field = coefficients::InternalField::JRM33;
    let (g, h) = type_field.get_coefficients();

    println!("{:?}", g);
    // p[[0, 0]] = 1.
}
