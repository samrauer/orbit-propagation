#![allow(dead_code)]

use std::{ops::{Mul}};
use super::state::State;



const RK12_SAFETY_FACTOR: f64 = 0.8;
const RK12_C2: f64 = 0.5; // TODO: decide whether to change F(x) to F(t,x)
const RK12_C3: f64 = 1.0; // TODO: decide whether to change F(x) to F(t,x)
const RK12_A21: f64 = 0.5;
const RK12_A31: f64 = 1.0 / 256.0;
const RK12_A32: f64 = 255.0 / 256.0;
const RK12_B1: f64 = 1.0 / 512.0;
const RK12_B2: f64 = 255.0 / 256.0;
const RK12_B3: f64 = 1.0 / 512.0;
const RK12_B1S: f64 = 1.0 / 256.0;
const RK12_B2S: f64 = 255.0 / 256.0;
const RK12_B3S: f64 = 0.0;



pub struct RK12Options {
    // relative tolerance
    rtol: f64,
    // absolute tolerance
    atol: f64,
    // min integrator step size relative to total integration time
    rhmin: f64,
    // max integrator step size relative to total integration time
    rhmax: f64,
    scale_min: f64,
    scale_max: f64,
} 

impl Default for RK12Options {
    fn default() -> Self {
        RK12Options { 
            rtol: 1e-3, 
            atol: 1e-6,
            rhmin: 1e-7,
            rhmax: 0.34,
            scale_min: 0.2,
            scale_max: 5.0,
        }
    }
}


/// Fehlberg RK1(2) integration
///     Integrates from t in tspan
pub fn integrate_rk12<F, T>(f: F, x0: &T, tspan: (f64,f64), opts: Option<RK12Options>) -> Result<T, String>
where 
    F: Fn(&T) -> T,
    T: State,
{
    // unwrap the options
    let opts = opts.unwrap_or_default();
    
    // TODO: handle inputs that are invalid

    let n: f64 = x0.length() as f64;

    let mut t: f64 = tspan.0;
    let mut x: T = x0.clone();
    let mut result: Vec<(f64, T)> = Vec::new();
    #[allow(unused_variables)]
    let mut steps: u64 = 0;

    // starting h = min { tspan/5 , 1 }
    let mut h: f64 = ((tspan.1 - tspan.0) / 5.0).min(1.0);

    while t < tspan.1 {
        

        // make sure we don't step past dt
        if t + h > tspan.1 {
            h = tspan.1 - t;

            // TODO: add logic so if 
        } else if h < opts.rhmin * tspan.1 {
            return Err(
                format!(
                    "Step size dropped below h_min ({}) at t = {}",
                    opts.rhmin * tspan.1, t
                )
            )
        }

        // compute k values
        let k1 = f(&x);

        let xk2 = x.clone() + (k1.clone() * RK12_A21) * h;
        let k2 = f(&xk2);

        let xk3 = x.clone() + (k1.clone() * RK12_A31 + k2.clone() * RK12_A32) * h;
        let k3 = f(&xk3);


        // compute 1st and 2nd order
        let xh_1 = x.clone() + (k1.clone() * RK12_B1S + k2.clone() * RK12_B2S + k3.clone() * RK12_B3S) * h;
        let xh_2 = x.clone() + (k1.clone() * RK12_B1 + k2.clone() * RK12_B2 + k3.clone() * RK12_B3) * h;

        // compute error
        // TODO: did I compute error allowed correctly? I don't think I did
        let error_allowed = (xh_1.clone() - xh_2.clone())
            .map(|x|  opts.atol + opts.rtol * x.abs());

        // normalize the error
        // sqrt( 1/dim * (error / error allowed)^2 )
        let error_normalized = (xh_1.clone() - xh_2.clone())
            .zip_map(&error_allowed, |err, err_allow| err / err_allow)
            .map(|x| x*x)
            .sum()
            .mul(1.0 / n)
            .sqrt();
        
        
        if error_normalized <= 1.0 {
            // accept step
            t = t + h;
            x = xh_2;
            result.push((t, x.clone()));
            steps += 1;
        }

        // compute next h
        let scale = (RK12_SAFETY_FACTOR * error_normalized.powf(-1.0 / (2.0 + 1.0))).clamp(opts.scale_min, opts.scale_max);

        // keep h smaller than h upper bound
        h = (h * scale).min(opts.rhmax * tspan.1);
    }

    Ok(x)

}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rk12_const_deriv() {
        // integrating 1 should product x(t) = t
        fn f(_: &f64) -> f64 {
            1.0
        }

        let x0: f64 = 0.0;
        let tspan: (f64, f64) = (0.0, 2.0);

        // expected value
        let x1_expected: f64 = 2.0;

        // true value
        let result = integrate_rk12(f, &x0, tspan, None);
        let x1_truth: f64 = result.unwrap();

        assert!((x1_expected - x1_truth).abs() < 1e-12);
    }

    #[test]
    fn rk12_dependent_deriv() {
        // integrating x should product x(t) = x0 * e^t
        fn f(x: &f64) -> f64 {
            *x
        }

        let x0: f64 = 1.0;
        let tspan: (f64, f64) = (0.0, 2.0);

        // expected value
        let x1_expected: f64 = x0 * 2f64.exp();

        // true value
        let result = integrate_rk12(f, &x0, tspan, None);
        assert!(result.is_ok());
        let x1_truth: f64 = result.unwrap();

        assert!((x1_expected - x1_truth).abs() < 1e-3);
    }
}
