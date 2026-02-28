mod integrator;
mod numerics;
mod old;
mod orbital_mechanics;
mod state;

extern crate nalgebra as na;
use std::f64::consts::PI;

use na::{Vector3, Vector6};

use orbital_mechanics::consts::{MU_EARTH, R_EARTH};

use integrator::rk45::{RK45Options, integrate_rk45};

fn main() {
    println!("Hello, world!");

    let altitude = 250e3; // m
    let orbit_radius = R_EARTH + altitude; // m
    let r0: Vector3<f64> = Vector3::new(orbit_radius, 0.0, 0.0); // m

    let vmag0 = (orbital_mechanics::consts::MU_EARTH / r0.norm()).sqrt();
    let v0: Vector3<f64> = vmag0 * Vector3::new(0.0, 0.5 * 2.0f64.sqrt(), 0.5 * 2.0f64.sqrt()); // m/s
    // let v0: Vector3<f64> = vmag0 * Vector3::new(0.0, 1.0, 0.0); // m/s

    let x0: Vector6<f64> = stack_vectors(&r0, &v0);

    fn dynamics(x: &Vector6<f64>) -> Vector6<f64> {
        // TODO: there's gotta be a better way to do this than manually index
        let r: Vector3<f64> = Vector3::new(x[0], x[1], x[2]);
        let v: Vector3<f64> = Vector3::new(x[3], x[4], x[5]);

        let xdot = v;

        // compute gravity (F = ma => vdot = -mu/r^3 * r)
        let r_norm = r.norm();
        let vdot = (-orbital_mechanics::consts::MU_EARTH / r_norm.powi(3)) * r;

        stack_vectors(&xdot, &vdot)
    }

    let expected_period = 2.0 * PI * (r0.norm().powi(3) / MU_EARTH).sqrt();

    let tspan = (0.0, expected_period);
    let options = RK45Options {
        rtol: 1e-12,
        atol: 1e-10,
        ..Default::default()
    };
    let result = integrate_rk45(dynamics, &x0, tspan, Some(options));

    let x1 = result.unwrap();
    let r1: Vector3<f64> = Vector3::new(x1[0], x1[1], x1[2]);
    let v1: Vector3<f64> = Vector3::new(x1[3], x1[4], x1[5]);

    println!("t: {} -> {}", tspan.0, tspan.1);
    println!("r: {} -> {}", r0, r1);
    println!("v: {} -> {}", v0, v1);
}

fn stack_vectors(r: &Vector3<f64>, v: &Vector3<f64>) -> Vector6<f64> {
    Vector6::from_iterator(r.iter().chain(v.iter()).cloned())
}
