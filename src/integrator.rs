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
pub fn integrate_euler<F>(f: F, x: &Vec<f64>, dt: f64) -> Vec<f64>
where F: Fn(&Vec<f64>) -> Vec<f64>,
{
    assert!(dt >= 0.0, "dt must be positive");

    let xdot = f(&x);

    x.iter()
        .zip(xdot.iter())
        .map(|(xi, xdoti)| xi + xdoti * dt)
        .collect()

    // collect will build vector return type because vector implements FromIterator<T>
}
