#![allow(dead_code)]


use na::Vector3;
use crate::orbital_mechanics::{orbit::Trajectory, stumpff::{stumpff_c, stumpff_s}};


/// tolerance for picking the initial guess of x.
/// chosen arbitrarily, change in future if needed.
const TOL_PARABOLA: f64 = 5e-2;

/// newton-raphson iteration parameters
const NR_TOL: f64 = 1e-9;
const NR_MAX_ITER: i64 = 100;



fn guess_x_ellipse(alpha: f64, mu: f64, dt: f64) -> f64 {
    // ellipse (BMW eq. 4.75)
    mu * dt * alpha
}

fn guess_x_parabola(mu: f64, r0: f64, dt: f64) -> f64 {
    // parabola
    // TODO: make a smarter guess
    // assume small step (linear term dominates)
    mu.sqrt() * dt / r0
}

fn guess_x_hyperbola(alpha: f64, mu: f64, r0: f64, rv0: f64, dt: f64) -> f64 {
    // hyperbola (BMW eq. 4.76)
    (-1.0/alpha).sqrt() * (
        -2.0 * mu * dt * alpha /
        (rv0 + (-mu * alpha).sqrt() * (1.0 - r0*alpha))
    ).ln()
}



/// Propagates an orbit using the universal variable formulation for Keplers problem.
/// Fundamentals of Astrodynamics (BMWS), Sec 4.4.4.
pub fn propagate_orbit(trajectory: &Trajectory, dt: f64) -> Result<Trajectory, String> {
    // technically we can propagate the orbit backwards in time, but don't want to
    //  support that right now
    assert!(dt >= 0.0, "must propagate orbit forward in time");

    // norms of the position and velocity
    let r0: f64 = trajectory.r.norm();
    let v0: f64 = trajectory.v.norm();
    let rv0: f64 = trajectory.r.dot(&trajectory.v);

    if dt < 0.1 {
        // use euler integration bc uv method has problem with very small time steps
        let rdot = trajectory.v.clone();
        let vdot = -1.0/r0.powi(3) * trajectory.r.clone();
        
        let r_new = trajectory.r.clone() + dt*rdot;
        let v_new = trajectory.v.clone() + dt*vdot;

        return Ok(Trajectory { r: r_new, v: v_new, mu: trajectory.mu })
    }


    // alpha = 1/a (a is the semi-major axis)
    let alpha: f64 = -v0.powi(2)/trajectory.mu + 2.0/r0;

    // initial guess for x (universal anomaly)
    let mut x: f64 = if alpha > TOL_PARABOLA {
        guess_x_ellipse(alpha, trajectory.mu, dt)
    } else if alpha < -TOL_PARABOLA {
        guess_x_hyperbola(alpha, trajectory.mu, r0, rv0, dt)
    } else {
        guess_x_parabola(trajectory.mu, r0, dt)
    };


    // newton-raphson iteration
    let mut counter = 0;
    while counter < NR_MAX_ITER {
        // compute z from x
        let z = alpha * x.powi(2);
        let s = stumpff_s(z);
        let c = stumpff_c(z);

        // universal kepler equation (BMW eq. 4.39 / 4.41)
        //   sqrt(mu)*dt = rv0/sqrt(mu) * chi^2 * C + (1 - r0*alpha) * chi^3 * S + r0 * chi
        //   f = sqrt(mu) * (tn - t)
        let f: f64 =
            (rv0 / trajectory.mu.sqrt()) * x.powi(2) * c
            + (1.0 - r0 * alpha) * x.powi(3) * s
            + r0 * x
            - trajectory.mu.sqrt() * dt;

        // universal variable formulation (BMW eq. 4.40 / 4.44)
        //   sqrt(mu) * dt/dx = r
        //                    = chi^2 * C + rv0/sqrt(mu) * chi * (1 - zS) + r0 * (1 - zC)
        let fp: f64 =
            c * x.powi(2)
            + (rv0 / trajectory.mu.sqrt()) * x * (1.0 - z*s)
            + r0 * (1.0 - z*c);

        // update chi (BMW eq. 4.42)
        let dx: f64 = f/fp;
        x = x - dx;

        // check tolerance
        if dx.abs() < NR_TOL {
            break
        } else {
            counter += 1;
        }
    }

    if counter >= NR_MAX_ITER {
        return Err(
            String::from("Max iterations reached when searching for x (universal anomaly)")
        )
    }

    // re-evaluate z at final x
    let z = alpha * x.powi(2);
    let s = stumpff_s(z);
    let c = stumpff_c(z);

    // evaluate f and g functions, then compute r
    // (BMW eq. 4.58)
    let f: f64 = 1.0 - (x.powi(2) / r0) * c;
    // (BMW eq. 4.61)
    let g: f64 = dt - (x.powi(3) / trajectory.mu.sqrt()) * s;
    // (BMW eq. 4.45)
    let r_new: Vector3<f64> = f * trajectory.r.clone() + g * trajectory.v.clone();
    let r1: f64 = r_new.norm();
    
    // evalute fdot and gdot functions, then compute v
    // (BMW eq. 4.63)
    let fdot: f64 = (trajectory.mu.sqrt() / (r1 * r0)) * (alpha * x.powi(3) * s - x);
    // (BMW eq. 4.62)
    let gdot: f64 = 1.0 - (x.powi(2) / r1) * c;
    // (BMW eq. 4.46)
    let v_new: Vector3<f64> = fdot * trajectory.r.clone() + gdot * trajectory.v.clone();

    Ok(Trajectory { r: r_new, v: v_new, mu: trajectory.mu })
}
