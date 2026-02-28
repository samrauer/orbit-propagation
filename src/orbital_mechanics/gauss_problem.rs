use std::f64::consts::PI;

use crate::orbital_mechanics::keplerian_propagation::propagate_orbit;
use crate::orbital_mechanics::{
    orbit::{Trajectory, TransferTrajectory},
    stumpff::{stumpff_c, stumpff_dcdz, stumpff_dsdz, stumpff_s},
};
use na::Vector3;

/// Solve the Gauss problem to get an interecpt trajectory
pub fn transfer_trajectory_solver(
    traj_target: &Trajectory,
    traj_chaser: &Trajectory,
    dt: f64,
) -> Result<(TransferTrajectory, TransferTrajectory), String> {
    assert!(dt >= 0.0, "must propagate orbit forward in time");

    // target at intercept time
    let traj_targ1 = propagate_orbit(&traj_target, dt)?;

    let r1_vec = traj_chaser.r;
    let r2_vec = traj_targ1.r;

    let r1 = r1_vec.norm();
    let r2 = r2_vec.norm();
    let mu = traj_chaser.mu;

    let nu_short = (r1_vec.cross(&r2_vec)).norm().atan2(r1_vec.dot(&r2_vec));
    let nu_long = 2.0 * PI - nu_short;

    // short way
    let (f_short, g_short, gdot_short) = gauss_uv_fg_solver(r1, r2, dt, mu, nu_short)?;
    let v1_short: Vector3<f64> = (r2_vec - f_short * r1_vec) / g_short;
    let v2_short: Vector3<f64> = (gdot_short * r2_vec - r1_vec) / g_short;
    let traj_short1 = Trajectory {
        r: r1_vec,
        v: v1_short,
        mu: mu,
    };
    let traj_short2 = Trajectory {
        r: r2_vec,
        v: v2_short,
        mu: mu,
    };

    //Trajectory { r: r1, v: v1, mu: trajectory.mu }

    // long way
    let (f_long, g_long, gdot_long) = gauss_uv_fg_solver(r1, r2, dt, mu, nu_long)?;
    let v1_long: Vector3<f64> = (r2_vec - f_long * r1_vec) / g_long;
    let v2_long: Vector3<f64> = (gdot_long * r2_vec - r1_vec) / g_long;
    let traj_long1 = Trajectory {
        r: r1_vec,
        v: v1_long,
        mu: mu,
    };
    let traj_long2 = Trajectory {
        r: r2_vec,
        v: v2_long,
        mu: mu,
    };

    Ok((
        TransferTrajectory {
            traj1: traj_short1,
            traj2: traj_short2,
        },
        TransferTrajectory {
            traj1: traj_long1,
            traj2: traj_long2,
        },
    ))
}

pub fn gauss_uv_fg_solver(
    r1: f64,
    r2: f64,
    dt: f64,
    mu: f64,
    nu: f64,
) -> Result<(f64, f64, f64), String> {
    // direction of motion
    let dm: f64 = (PI - nu).signum();

    // step 1: from r1 and r2 and "direction of motion", evaluate the constant (eq. 5-15, eq. 5-37 BMW)
    // #![ignore(non_snake_case)]
    let A: f64 = dm * (r1 * r2 * (1.0 + nu.cos())).sqrt();

    // step 2: pick a trial value for z
    //         z = deltaE^2 for elliptical
    //         z = deltaF^2 for hyperolic
    //         good initial guess is z = 0
    let mut z = 0.0;
    let mut dz = 0.0;

    // newton raphson parameters
    let dt_tol = 1e-6;
    let nr_max_iter = 500;
    let mut nr_counter = 0;

    // bisection search parameters
    let bs_max_iter = 100;
    let mut bs_counter = 0;

    loop {
        if nr_counter >= nr_max_iter {
            return Err(String::from(
                "Warning: max newton-raphson iterations reached when solving for z",
            ));
        }
        if bs_counter >= bs_max_iter {
            return Err(String::from(
                "Warning: max bisection-search iterations reached when solving for z",
            ));
        }

        bs_counter += 1;

        let z_trial = z + dz;

        // step 3: evaluate S and C for selected z (eq. 4-37 and eq. 4-38 BMW)
        let s = stumpff_s(z_trial);
        let c = stumpff_c(z_trial);

        // step 4: determine aux variable y (eq. 5-17 BMW)
        let y = r1 + r2 - A * (1.0 - z_trial * s) / c.sqrt();

        // step 4.5: check if y is negative (only can happen with short way trajectories)
        //           keep halving step size until y is no longer negative
        if y < 0.0 {
            dz = 0.5 * dz;
            continue;
        }

        // z_trial is accepted
        z = z_trial;

        // increment counters
        nr_counter += 1;
        bs_counter = 0;

        // step 5: determine x (eq. 5-18 BMW)
        let x = (y / c).sqrt();

        // step 6: check trial value of z by computing dt_trial (eq. 5-20 BMW)
        //         then compare to true dt
        //         then newton-raphson iterate

        let dt_trail = (s * x.powi(3) + A * y.sqrt()) / mu.sqrt();
        let dt_error = dt_trail - dt;

        // Debugging:
        // println!("{}, {}", z, dt_error);

        if dt_error.abs() < dt_tol {
            // Debugging:
            // println!("breaking");
            break;
        }
        // Debugging:
        // println!("did not break...");

        let sprime = stumpff_dsdz(z);
        let cprime = stumpff_dcdz(z);

        let dtdz = 1.0 / mu.sqrt()
            * (x.powi(3) * (sprime - 3.0 * s * cprime / (2.0 * c))
                + A / 8.0 * (3.0 * s * y.sqrt() / c + A / x));

        // clip to prevent huge jumps
        let mut dt_max_jump = 5.0;
        if nr_counter > 200 {
            // need bigger steps because very hyperbolic or very elliptical
            dt_max_jump = 20.0;
        }
        dz = (-dt_error / dtdz).clamp(-dt_max_jump, dt_max_jump);
    }

    // step 7: evaluate f, g, gdot (eq. 5-21, 5-22, 5-23 BMW)
    //         then compute v1 and v2 (eq. 5-24, 5-25 BMW)

    let s = stumpff_s(z);
    let c = stumpff_c(z);
    let y = r1 + r2 - A * (1.0 - z * s) / c.sqrt();

    let f = 1.0 - y / r1;
    let g = A * (y / mu).sqrt();
    let gdot = 1.0 - y / r2;

    Ok((f, g, gdot))
}

#[cfg(test)]
mod tests {
    use std::f64::consts::PI;

    use super::*;
    use crate::orbital_mechanics::consts::{MU_EARTH, R_EARTH};
    use approx::assert_relative_eq;

    #[test]
    fn test_transfer_trajectory_solver() {
        let traj_target: Trajectory = Trajectory {
            r: Vector3::new(R_EARTH + 500e3, 0.0, 0.0),
            v: Vector3::new(0.0, 0.0, (MU_EARTH / (R_EARTH + 500e3)).sqrt()),
            mu: MU_EARTH,
        };
        let period_target: f64 = 5677.0;

        let traj0_chaser: Trajectory = Trajectory {
            r: Vector3::new(R_EARTH + 600e3, 0.0, 0.0),
            v: Vector3::new(0.0, 0.0, (MU_EARTH / (R_EARTH + 600e3)).sqrt()),
            mu: MU_EARTH,
        };

        let dt: f64 = 0.5 * period_target;

        dbg!(period_target);

        let result = transfer_trajectory_solver(&traj_target, &traj0_chaser, dt);
        assert!(result.is_ok());
        let result = result.unwrap();

        // testing
        dbg!(&result);

        let traj_target_t = propagate_orbit(&traj_target, dt).unwrap();
        let traj_post_transfer_short = propagate_orbit(&result.0.traj1, dt).unwrap();
        let traj_post_transfer_long = propagate_orbit(&result.1.traj1, dt).unwrap();

        // check that position is correct
        assert_relative_eq!(
            traj_post_transfer_short.r.x,
            traj_target_t.r.x,
            epsilon = 5.0
        );
        assert_relative_eq!(
            traj_post_transfer_short.r.y,
            traj_target_t.r.y,
            epsilon = 5.0
        );
        assert_relative_eq!(
            traj_post_transfer_short.r.z,
            traj_target_t.r.z,
            epsilon = 5.0
        );

        assert_relative_eq!(
            traj_post_transfer_long.r.x,
            traj_target_t.r.x,
            epsilon = 5.0
        );
        assert_relative_eq!(
            traj_post_transfer_long.r.y,
            traj_target_t.r.y,
            epsilon = 5.0
        );
        assert_relative_eq!(
            traj_post_transfer_long.r.z,
            traj_target_t.r.z,
            epsilon = 5.0
        );
    }

    #[test]
    fn test_gauss_uv_fg_solver() {
        let r1 = 6370e3 + 800e3;
        let r2 = 6370e3 + 700e3;
        let dt = 45.0 * 60.0;
        let mu = 3.98e14;
        let nu_short = PI / 2.0;
        let nu_long = 3.0 * PI / 2.0;

        let result_short = gauss_uv_fg_solver(r1, r2, dt, mu, nu_short);
        assert!(result_short.is_ok());
        let result_short = result_short.unwrap();

        assert_relative_eq!(result_short.0, -0.5651232638486305, max_relative = 1e-5);
        assert_relative_eq!(result_short.1, 1195.5330869650286, max_relative = 1e-5);
        assert_relative_eq!(result_short.2, -0.5872607923330526, max_relative = 1e-5);

        let result_long = gauss_uv_fg_solver(r1, r2, dt, mu, nu_long);
        assert!(result_long.is_ok());
        let result_long = result_long.unwrap();

        assert_relative_eq!(result_long.0, -0.3056781728523754, max_relative = 1e-5);
        assert_relative_eq!(result_long.1, -1091.9566526520407, max_relative = 1e-5);
        assert_relative_eq!(result_long.2, -0.3241460395122393, max_relative = 1e-5);
    }

    #[test]
    fn test_gauss_uv_fg_solver_2() {
        let r1 = 6370e3 + 400e3;
        let r2 = 6370e3 + 400e3;
        let dt = 10.0 * 60.0;
        let mu = 3.98e14;
        let nu_short = PI / 3.0;
        let nu_long = 2.0 * PI - PI / 3.0;

        let result_short = gauss_uv_fg_solver(r1, r2, dt, mu, nu_short);
        assert!(result_short.is_ok());
        let result_short = result_short.unwrap();

        assert_relative_eq!(result_short.0, 0.7446003442538833, max_relative = 1e-5);
        assert_relative_eq!(result_short.1, 546.5089143188977, max_relative = 1e-5);
        assert_relative_eq!(result_short.2, 0.7446003442538833, max_relative = 1e-5);

        let result_long = gauss_uv_fg_solver(r1, r2, dt, mu, nu_long);
        assert!(result_long.is_ok());
        let result_long = result_long.unwrap();

        assert_relative_eq!(result_long.0, -8.58292, max_relative = 1e-5);
        assert_relative_eq!(result_long.1, -3347.6196, max_relative = 1e-5);
        assert_relative_eq!(result_long.2, -8.5829235, max_relative = 1e-5);
    }
}
