use std::iter::zip;

/// Euler integration for 1 step
///
/// # Arguments
///
/// * `f` - A function that computes the time derivative.
/// * `x` - The current state.
/// * `dt` - The time step (must be positive).
///
/// # Returns
///
/// A new `Vec<f64>` representing the next state.
pub fn integrate_euler<F>(f: F, x0: &Vec<f64>, dt: f64) -> Vec<f64>
where F: Fn(&Vec<f64>) -> Vec<f64>,
{
    assert!(dt >= 0.0, "dt must be positive");

    let xdot = f(&x0);

    zip(x0, xdot)
        .map(|(xi, xdoti)| xi + xdoti * dt)
        .collect()
}




/// Predictor-corrector integration for 1 step
///
/// # Arguments
///
/// * `f` - A function that computes the time derivative.
/// * `x` - The current state.
/// * `dt` - The time step (must be positive).
///
/// # Returns
///
/// A new `Vec<f64>` representing the next state.
pub fn integrate_euler_with_corrector<F>(f: F, x0: &Vec<f64>, dt: f64) -> Vec<f64>
where F: Fn(&Vec<f64>) -> Vec<f64>,
{
    assert!(dt >= 0.0, "dt must be positive");

    // x1 guess using forward euler
    let xdot_x0 = f(&x0);
    let x1: Vec<f64> = zip(x0, &xdot_x0)
        .map(|(xi, xdoti)| xi + xdoti * dt)
        .collect();

    // predictor by taking derivative at prediction spot and averaging
    let xdot_x1 = f(&x1);
    let xdot_avg: Vec<f64> = zip(&xdot_x0, &xdot_x1)
        .map(|(xdot0, xdot1)| 0.5 * (xdot0 + xdot1))
        .collect();

    // use average derivative to propagate state forward
    zip(x0, &xdot_avg)
        .map(|(xi, xdoti)| xi + xdoti * dt)
        .collect()
}
