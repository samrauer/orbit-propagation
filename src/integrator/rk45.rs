#![allow(warnings)]
#![allow(dead_code)]

use na::Vector3;

// TODO: finish RK45

const RK45_SAFETY_FACTOR: f64 = 0.88;


/// Runge-Kutta 4
pub fn integrate_rk4<F>(f: F, x0: &Vector3<f64>, dt: f64) -> Vector3<f64>
where 
    F: Fn(&Vector3<f64>) -> Vector3<f64>,

{

    todo!();
}