use na::{Vector3};


/// Euler integration for 1 step
#[allow(dead_code)]
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


/// Euler predictor corrector (midpoint) integration for 1 step
pub fn integrate_euler_pc<F>(f: F, x0: &Vector3<f64>, dt: f64) -> Vector3<f64>
where 
    F: Fn(&Vector3<f64>) -> Vector3<f64>,

{
    // check we have valid conditions
    assert!(dt >= 0.0, "dt must be positive");
    assert!(
        !x0.iter().any(|e| e.is_nan() || e.is_infinite()),
        "x0 contains NaN or infinity"
    );

    // compute derivative at x0
    let xdot_x0 = f(&x0);
    let x1 = x0 + xdot_x0*dt;

    // compute derivative at x1
    let xdot_x1 = f(&x1);

    // compute average derivative between x0 and x1
    let xdot_avg = 0.5 * (xdot_x0 + xdot_x1);

    // use averaged derivative to propagate x0 forward
    x0 + xdot_avg * dt
}

