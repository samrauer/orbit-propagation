#![allow(dead_code)]

use std::{ops::{Mul}};

use na::{DefaultAllocator, Dim, OVector, Vector3, allocator::Allocator};


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


// RK12 constants
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
///     Integrates from t = 0 --> dt
pub fn integrate_rk12<F, D>(f: F, x0: &OVector<f64, D>, dt: f64, opts: Option<RK12Options>) -> Result<OVector<f64, D>, String>
where 
    F: Fn(&OVector<f64, D>) -> OVector<f64, D>,
    D: Dim,
    DefaultAllocator: Allocator<D>,
{
    // unwrap the options
    let opts = opts.unwrap_or_default();
    
    // TODO: handle inputs that are invalid

    let n = x0.nrows() as f64;

    let mut t: f64 = 0.0;
    let mut x: OVector<f64, D> = x0.clone();
    let mut result: Vec<(f64, OVector<f64, D>)> = Vec::new();

    // starting h = min { dt/5 , 1 }
    let mut h: f64 = (dt / 5.0).min(1.0);

    while t + f64::EPSILON <= dt {
        // make sure we don't step past dt
        if t + h > dt {
            h = dt - t;
        } else if h < opts.rhmin * dt {
            return Err(
                format!(
                    "Step size dropped below h_min ({}) at t = {}",
                    opts.rhmin * dt, t
                )
            )
        }

        // compute k values
        let k1 = f(&x);

        let xk2 = &x + h*(RK12_A21 * &k1);
        let k2 = f(&xk2);

        let xk3 = &x + h*(RK12_A31 * &k1 + RK12_A32 * &k2);
        let k3 = f(&xk3);


        // compute 1st and 2nd order
        let xh_1 = RK12_B1S * &k1 + RK12_B2S * &k2 + RK12_B3S * &k3;
        let xh_2 = RK12_B1 * &k1 + RK12_B2 * &k2 + RK12_B3 * &k3;

        // compute error
        let error_allowed = (&xh_1 - &xh_2).map(|x|  opts.atol + opts.rtol * x.abs());

        // normalize the error
        // sqrt( 1/dim * (error / error allowed)^2 )
        let error_normalized = (&xh_1 - &xh_2)
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
        }

        // compute next h
        let scale = (RK12_SAFETY_FACTOR * error_normalized.powf(-1.0 / (2.0 + 1.0))).clamp(opts.scale_min, opts.scale_max);

        // keep h smaller than h upper bound
        h = (h * scale).min(opts.rhmax * dt);
    }

    Ok(x)

}








// TODO: finish RK45

const RK45_SAFETY_FACTOR: f64 = 0.88;


/// Runge-Kutta 4
pub fn integrate_rk4<F>(f: F, x0: &Vector3<f64>, dt: f64) -> Vector3<f64>
where 
    F: Fn(&Vector3<f64>) -> Vector3<f64>,

{

    todo!();
}