#![allow(dead_code)]

use na::Vector3;


/// Euler integration for 1 step
pub fn integrate_euler<F>(f: F, x0: &Vector3<f64>, dt: f64) -> Vector3<f64>
where 
    F: Fn(&Vector3<f64>) -> Vector3<f64>,

{
    // check we have valid conditions
    assert!(dt >= 0.0, "dt must be positive");
    assert!(
        !x0.iter().any(|e| e.is_nan() || e.is_infinite()),
        "x0 contains NaN or infinity"
    );

    // compute derivative
    let xdot = f(&x0);

    // euler integration
    x0 + xdot*dt
}
