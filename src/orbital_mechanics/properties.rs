#![allow(dead_code)]

use na::Vector3;
use super::consts::MU_EARTH;

// specific kinetic energy
pub fn compute_energy_kinetic(v: &Vector3<f64>) -> f64 {
    0.5 * v.norm().powi(2)
}

// specific potential energy
pub fn compute_energy_potential(r: &Vector3<f64>) -> f64 {
    -MU_EARTH / r.norm()
}

// specific angular momentum
pub fn compute_angular_momentum(r: &Vector3<f64>, v: &Vector3<f64>) -> Vector3<f64> {
    r.cross(v)
}



#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_relative_eq, assert_abs_diff_eq};

    #[test]
    fn compute_energy_kinetic_tests() {
        let v: Vector3<f64> = Vector3::new(100.0, 0.0, 0.0);
        let expected_energy = 0.5 * 100.0f64.powi(2);
        let actual_energy = compute_energy_kinetic(&v);

        assert_abs_diff_eq!(expected_energy, actual_energy);


        let v: Vector3<f64> = Vector3::new(10.0, -500.0, 67.0);
        let expected_energy = 1.272945e5;
        let actual_energy = compute_energy_kinetic(&v);

        assert_relative_eq!(expected_energy, actual_energy, max_relative = 1e-12);
    }

    #[test]
    fn compute_energy_potential_tests() {
        let r: Vector3<f64> = Vector3::new(100.0e3, 0.0, 0.0);
        let expected_energy = -3.986004418e9;
        let actual_energy = compute_energy_potential(&r);

        assert_relative_eq!(expected_energy, actual_energy, max_relative = 1e-12);


        let r: Vector3<f64> = Vector3::new(10.0e3, -500.0e3, 67.0e3);
        let expected_energy = -7.899833865585763e8;
        let actual_energy = compute_energy_potential(&r);

        assert_relative_eq!(expected_energy, actual_energy, max_relative = 1e-12);
    }

    #[test]
    fn compute_angular_momentum_tests() {
        // TODO
        assert!(true);
    }
}
