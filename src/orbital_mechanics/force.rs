#![allow(dead_code)]

use na::Vector3;

use super::consts::MU_EARTH;




/// Gravity force from position (m) in ECI and mass of the orbiting object (kg)
pub fn compute_force_gravity(r: &Vector3<f64>, m: f64) -> Vector3<f64> {
    let r_norm = r.norm();

    (-MU_EARTH * m / r_norm.powi(3)) * r
}



#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_abs_diff_eq};

    #[test]
    fn compute_force_gravity_along_axis() {
        // LEO is 250 - 800 km, choose 250km
        let mass = 10.0;
        let expected_fnorm = MU_EARTH * mass / 250e3f64.powi(2);

        let r: Vector3<f64> = Vector3::new(250e3, 0.0, 0.0);
        let force: Vector3<f64> = compute_force_gravity(&r, 10.0);
        let force_direction: Vector3<f64> = force / force.norm();

        assert_relative_eq!(force.norm(), expected_fnorm, max_relative = 1e-12);
        assert_abs_diff_eq!(force_direction[0], -1.0);
        assert_abs_diff_eq!(force_direction[1], 0.0);
        assert_abs_diff_eq!(force_direction[2], 0.0);

        let r: Vector3<f64> = Vector3::new(0.0, 250e3, 0.0);
        let force: Vector3<f64> = compute_force_gravity(&r, 10.0);
        let force_direction: Vector3<f64> = force / force.norm();

        assert_relative_eq!(force.norm(), expected_fnorm, max_relative = 1e-12);
        assert_abs_diff_eq!(force_direction[0], 0.0);
        assert_abs_diff_eq!(force_direction[1], -1.0);
        assert_abs_diff_eq!(force_direction[2], 0.0);

        let r: Vector3<f64> = Vector3::new(0.0, 0.0, 250e3);
        let force: Vector3<f64> = compute_force_gravity(&r, 10.0);
        let force_direction: Vector3<f64> = force / force.norm();

        assert_relative_eq!(force.norm(), expected_fnorm, max_relative = 1e-12);
        assert_abs_diff_eq!(force_direction[0], 0.0);
        assert_abs_diff_eq!(force_direction[1], 0.0);
        assert_abs_diff_eq!(force_direction[2], -1.0);
    }

    #[test]
    fn compute_force_gravity_complex_a() {
        let mass = 25.0;
        let r: Vector3<f64> = Vector3::new(150e3, 150e3, 150e3);
        let rnorm: f64 = r.norm();

        let expected_fnorm = MU_EARTH * mass / rnorm.powi(2);

        let force = compute_force_gravity(&r, mass);
        let force_direction = force / force.norm();

        assert_relative_eq!(force.norm(), expected_fnorm, max_relative = 1e-12);
        assert_abs_diff_eq!(force_direction[0], -3.0f64.sqrt() / 3.0);
        assert_abs_diff_eq!(force_direction[1], -3.0f64.sqrt() / 3.0);
        assert_abs_diff_eq!(force_direction[2], -3.0f64.sqrt() / 3.0);
    }

    #[test]
    fn compute_force_gravity_complex_b() {
        let mass = 25.0;
        let r: Vector3<f64> = Vector3::new(150e3, -120e3, -80e3);
        let rnorm: f64 = r.norm();

        let expected_fnorm = MU_EARTH * mass / rnorm.powi(2);

        let force = compute_force_gravity(&r, mass);
        let force_direction = force / force.norm();

        assert_relative_eq!(force.norm(), expected_fnorm, max_relative = 1e-12);
        assert_abs_diff_eq!(force_direction[0], -0.720853996998319, epsilon = 1e-12);
        assert_abs_diff_eq!(force_direction[1], 0.576683197598655, epsilon = 1e-12);
        assert_abs_diff_eq!(force_direction[2], 0.384455465065770, epsilon = 1e-12);
    }
}

