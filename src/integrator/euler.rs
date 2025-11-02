#![allow(dead_code)]

use super::state::State;


/// Euler integration
pub fn integrate_euler<F, T>(f: F, x0: &T, tspan: (f64, f64), dt: f64) -> Result<T, String>
where 
    F: Fn(&T) -> T,
    T: State,
{
    // check we have valid conditions
    if dt < 0.0 {
        return Err(format!(
            "dt ({}) must be positive", 
            dt
        ));
    }
    if tspan.1 < tspan.0 {
        return Err(format!(
            "tspan must have time moving forward, i.e. tspan.1 ({}) > tspan.0 ({})", 
            tspan.1, tspan.0
        ))
    }

    // TODO: figure out how to iterate through the state
    // assert!(
    //     !x0.iter().any(|e| e.is_nan() || e.is_infinite()),
    //     "x0 contains NaN or infinity"
    // );

    let mut t = tspan.0;
    let mut x = x0.clone();
    let mut dt = dt;

    while t < tspan.1 {
        if t + dt > tspan.1 {
            dt = tspan.1 - t;
        }

        let xdot = f(&x);
        x = x + xdot * dt;
        t = t + dt;
    }

    Ok(x)
}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn euler_const_deriv() {
        // integrating 1 should product x(t) = t
        fn f(_: &f64) -> f64 {
            1.0
        }

        let x0: f64 = 0.0;
        let tspan: (f64, f64) = (0.0, 2.0);
        let dt: f64 = (tspan.1 - tspan.0) / 10.0 + 0.01; // add 0.01 to force end logic to fire

        // expected value
        let x1_expected: f64 = 2.0;

        // true value
        let result = integrate_euler(f, &x0, tspan, dt);
        let x1_truth: f64 = result.unwrap();

        assert!((x1_expected - x1_truth).abs() < 1e-12);
    }

    #[test]
    fn euler_dependent_deriv() {
        // integrating x should product x(t) = x0 * e^t
        fn f(x: &f64) -> f64 {
            *x
        }

        let x0: f64 = 1.0;
        let tspan: (f64, f64) = (0.0, 2.0);
        let dt: f64 = (tspan.1 - tspan.0) / 1000.0; // add 0.01 to force end logic to fire

        // expected value
        let x1_expected: f64 = x0 * 2f64.exp();

        // true value
        let result = integrate_euler(f, &x0, tspan, dt);
        assert!(result.is_ok());
        let x1_truth: f64 = result.unwrap();

        assert!((x1_expected - x1_truth).abs() < 0.05);
    }
}
