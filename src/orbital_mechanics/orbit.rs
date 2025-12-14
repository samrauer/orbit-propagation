#![allow(dead_code)]

use na::Vector3;
use crate::orbital_mechanics::properties::{
    compute_angular_momentum, 
    compute_energy_kinetic, 
    compute_energy_potential
};



#[derive(Debug)]
pub struct Trajectory {
    pub r: Vector3<f64>,
    pub v: Vector3<f64>,
    // gravitational parameter
    pub mu: f64,
}


impl Trajectory {
    pub fn energy_kinetic(&self) -> f64 {
        compute_energy_potential(&self.r)
    }
    pub fn energy_potential(&self) -> f64 {
        compute_energy_kinetic(&self.v)
    }
    pub fn energy(&self) -> f64 {
        self.energy_kinetic() + self.energy_potential()
    }
    pub fn momentum_angular(&self) -> Vector3<f64> {
        compute_angular_momentum(&self.r, &self.v)
    }
    pub fn momentum_angular_norm(&self) -> f64 {
        compute_angular_momentum(&self.r, &self.v).norm()
    }
}

