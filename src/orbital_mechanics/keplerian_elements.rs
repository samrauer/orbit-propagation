#![allow(dead_code)]
#![allow(unused_variables)]

use std::f64::consts::PI;

use na::Vector3;



pub struct KeplerianElements {
    // semi-major axis
    pub sma: f64,
    // eccentricity
    pub e: f64,
    // inclination
    pub i: f64,
    // right ascension of the ascending node
    pub raan: f64,
    // argument of periapsis
    pub aop: f64,
    // true anomaly
    pub f: f64,
}


impl KeplerianElements {
    pub fn new(sma: f64, e: f64, i: f64, raan: f64, aop: f64, f: f64) -> KeplerianElements {
        // for now don't worry about error handling
        //  we can just panic and fix it later if needed
        assert_sma_valid(sma);
        assert_eccentricity_valid(e);
        assert_inclination_valid(i);
        assert_raan_valid(raan);
        assert_argument_of_periapsis_valid(aop);
        assert_true_anomaly_valid(f);

        KeplerianElements { 
            sma: sma, 
            e: e, 
            i: i, 
            raan: raan, 
            aop: aop, 
            f: f 
        }
    }

    pub fn from_circular_orbit(radius: f64, i: f64, raan: f64, aop: f64, f: f64) -> KeplerianElements {
        KeplerianElements::new(
            radius, 
            0.0, 
            i, 
            raan, 
            aop, 
            f
        )
    }

    pub fn from_rv(r: &Vector3<f64>, v: &Vector3<f64>, mu: f64) -> KeplerianElements {
        // from position and velocity
        let rmag = r.norm();
        let vmag = v.norm();

        let a = rmag / (2.0 - rmag*vmag.powi(2)/mu);

        let e = (1.0/mu) * ((vmag.powi(2) - mu/rmag)*r - (r.dot(&v))*v);

        // let x = r;

        // TODO: finish the rest of this function, then add some tests
        todo!();
    }
}



fn assert_sma_valid(sma: f64) {
    assert!(sma > 0.0, "semi-major axis > 0");
}

fn assert_eccentricity_valid(e: f64) {
    assert!(e >= 0.0, "eccentricity >= 0");
    assert!(e <= 1.0, "eccentricity <= 0");
}

fn assert_inclination_valid(i: f64) {
    assert!(i <= PI, "inclination <= PI");
    assert!(i >= 0.0, "inclination >= 0");
}

fn assert_raan_valid(raan: f64) {
    assert!(raan >= 0.0, "RAAN >= 0");
    assert!(raan < 2.0*PI, "RAAN < 2 PI");
}

fn assert_argument_of_periapsis_valid(omega: f64) {
    assert!(omega >= 0.0, "argument of periapsis >= 0");
    assert!(omega < 2.0*PI, "argument of periapsis < 2 PI");
}

fn assert_true_anomaly_valid(f: f64) {
    assert!(f >= 0.0, "true anomaly >= 0");
    assert!(f < 2.0*PI, "true anomaly < 2 PI");
}
