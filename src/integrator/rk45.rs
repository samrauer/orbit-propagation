#![allow(dead_code)]

use std::ops::Mul;
use crate::integrator::state::State;


const SAFETY_FACTOR: f64 = 0.8;

// RK45 Dormand-Prince coefficients
// https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods

const C2: f64 = 1.0 / 5.0;
const C3: f64 = 3.0 / 10.0;
const C4: f64 = 4.0 / 5.0;
const C5: f64 = 8.0 / 9.0;
const C6: f64 = 1.0;
const C7: f64 = 1.0;

const A21: f64 = 1.0 / 5.0;

const A31: f64 = 3.0 / 40.0;
const A32: f64 = 9.0 / 40.0;

const A41: f64 = 44.0 / 45.0;
const A42: f64 = -56.0 / 15.0;
const A43: f64 = 32.0 / 9.0;

const A51: f64 = 19372.0 / 6561.0;
const A52: f64 = -25360.0 / 2187.0;
const A53: f64 = 64448.0 / 6561.0;
const A54: f64 = -212.0 / 729.0;

const A61: f64 = 9017.0 / 3168.0;
const A62: f64 = -355.0 / 33.0;
const A63: f64 = 46732.0 / 5247.0;
const A64: f64 = 49.0 / 176.0;
const A65: f64 = -5103.0 / 18656.0;

const A71: f64 = 35.0 / 384.0;
const A72: f64 = 0.0;
const A73: f64 = 500.0 / 1113.0;
const A74: f64 = 125.0 / 192.0;
const A75: f64 = -2187.0 / 6784.0;
const A76: f64 = 11.0 / 84.0;

// 5th order accurate solution
// NOTE: Don't even need these since it's the same as the previous function evaluation
const B1: f64 = 35.0 / 384.0;
const B2: f64 = 0.0;
const B3: f64 = 500.0 / 1113.0;
const B4: f64 = 125.0 / 192.0;
const B5: f64 = -2187.0 / 6784.0;
const B6: f64 = 11.0 / 84.0;
const B7: f64 = 0.0;

// 4th order accurate solution
const B1S: f64 = 5179.0 / 57600.0;
const B2S: f64 = 0.0;
const B3S: f64 = 7571.0 / 16695.0;
const B4S: f64 = 393.0 / 640.0;
const B5S: f64 = -92097.0 / 339200.0;
const B6S: f64 = 187.0 / 2100.0;
const B7S: f64 = 1.0 / 40.0;


pub struct RK45Options {
    // relative tolerance
    pub rtol: f64,
    // absolute tolerance
    pub atol: f64,
    // min integrator step size relative to total integration time
    pub rhmin: f64,
    // max integrator step size relative to total integration time
    pub rhmax: f64,
    pub scale_min: f64,
    pub scale_max: f64,
} 

impl Default for RK45Options {
    fn default() -> Self {
        RK45Options { 
            rtol: 1e-3, 
            atol: 1e-6,
            rhmin: 1e-7,
            rhmax: 0.34,
            scale_min: 0.2,
            scale_max: 5.0,
        }
    }
}


/// Runge-Kutta 45 Dormand-Prince
/// TODO: decide whether to change F(x) to F(t,x)
/// TODO: allow tspan to be a vector of points to evaluate at
pub fn integrate_rk45<F, T>(f: F, x0: &T, tspan: (f64,f64), opts: Option<RK45Options>) -> Result<T, String>
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

            // TODO: add logic so if h becomes super tiny, we step past tspan.1 and interpolate at tspan.1
            //       to avoid floating point edge cases
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

        let xk2 = x.clone() + (
            k1.clone() * A21
        ) * h;
        let k2 = f(&xk2);

        let xk3 = x.clone() + (
            k1.clone() * A31 + 
            k2.clone() * A32
        ) * h;
        let k3 = f(&xk3);

        let xk4 = x.clone() + (
            k1.clone() * A41 + 
            k2.clone() * A42 + 
            k3.clone() * A43
        ) * h;
        let k4 = f(&xk4);

        let xk5 = x.clone() + (
            k1.clone() * A51 + 
            k2.clone() * A52 + 
            k3.clone() * A53 + 
            k4.clone() * A54
        ) * h;
        let k5 = f(&xk5);

        let xk6 = x.clone() + (
            k1.clone() * A61 + 
            k2.clone() * A62 + 
            k3.clone() * A63 + 
            k4.clone() * A64 + 
            k5.clone() * A65
        ) * h;
        let k6 = f(&xk6);

        let xk7 = x.clone() + (
            k1.clone() * A71 + 
            k2.clone() * A72 + 
            k3.clone() * A73 + 
            k4.clone() * A74 + 
            k5.clone() * A75 +
            k6.clone() * A76
        ) * h;
        let k7 = f(&xk7);

        // compute 4th and 5th order
        let xh_5 = xk7;
        let xh_4 = x.clone() + (
            k1.clone() * B1S + 
            k2.clone() * B2S + 
            k3.clone() * B3S + 
            k4.clone() * B4S + 
            k5.clone() * B5S +
            k6.clone() * B6S +
            k7.clone() * B7S
        ) * h;

        // compute error
        // TODO: did I compute error allowed correctly? I don't think I did. Come back to this!
        let error_allowed = (xh_4.clone() - xh_5.clone())
            .map(|x|  opts.atol + opts.rtol * x.abs());

        // normalize the error
        // sqrt( 1/dim * (error / error allowed)^2 )
        let error_normalized = (xh_4.clone() - xh_5.clone())
            .zip_map(&error_allowed, |err, err_allow| err / err_allow)
            .map(|x| x*x)
            .sum()
            .mul(1.0 / n)
            .sqrt();
        
        
        if error_normalized <= 1.0 {
            // accept step
            t = t + h;
            x = xh_5;
            result.push((t, x.clone()));
            steps += 1;
        }

        // compute next h
        let scale = (SAFETY_FACTOR * error_normalized.powf(-1.0 / (2.0 + 1.0)))
            .clamp(opts.scale_min, opts.scale_max);

        // keep h smaller than h upper bound
        h = (h * scale).min(opts.rhmax * tspan.1);
    }

    Ok(x)
}


// TODO: tests!!!
