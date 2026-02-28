use crate::numerics::newton_raphson::newton_raphson_scalar;
use crate::orbital_mechanics::keplerian_elements::KeplerianElements;

const MA_2_EA_ATOL: f64 = 1e-12;

pub struct KeplerianOrbit {
    pub elements: KeplerianElements,
    // gravitational parameter
    pub mu: f64,
}

impl KeplerianOrbit {
    pub fn n(&self) -> f64 {
        // mean motion
        // n = sqrt(mu / a^3)
        (self.mu / self.elements.sma.powi(3)).sqrt()
    }
}

fn ea2ma(e: f64, ea: f64) -> f64 {
    // eccentric anomaly to mean anomaly
    // TODO: check for bounds on e and ea

    ea - e * ea.sin()
}

fn ma2ea(e: f64, ma: f64) -> f64 {
    // mean anomaly to eccentric anomaly
    // TODO: check for bounds on e and ma

    let f = |ea: f64| ea - e * ea.sin() - ma;
    let fp = |ea: f64| 1.0 - e * ea.cos();

    let ea0: f64 = if ma / (1.0 - e) < (6.0 * (1.0 - e) / e).sqrt() {
        ma / (1.0 - e)
    } else {
        (6.0 * ma / e).powf(1.0 / 3.0)
    };

    newton_raphson_scalar(f, fp, ea0, MA_2_EA_ATOL)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_ea2ma_circular() {
        let e = 0.0;
        let ea = 5.0 / 6.0 * PI;
        // should the the same, its a circular orbit!
        let expected_ma = ea;

        let actual_ma = ea2ma(e, ea);

        assert_abs_diff_eq!(expected_ma, actual_ma, epsilon = 1e-12);
    }

    #[test]
    fn test_ma2ea_circular() {
        let e = 0.0;
        let ma = 5.0 / 6.0 * PI;
        // should the the same, its a circular orbit!
        let expected_ea = ma;

        let actual_ea = ma2ea(e, ma);

        assert_abs_diff_eq!(expected_ea, actual_ea, epsilon = 2.0 * MA_2_EA_ATOL);
    }

    #[test]
    fn test_ea2ma_elliptical() {
        let e = 0.7;
        let ea = 9.0 / 5.0 * PI;
        let expected_ma = 6.066316453066359;

        let actual_ma = ea2ma(e, ea);

        assert_abs_diff_eq!(expected_ma, actual_ma, epsilon = 1e-12);
    }

    #[test]
    fn test_ma2ea_elliptical() {
        let e = 0.7;
        let ma = 6.066316453066359;
        // should the the same, its a circular orbit!
        let expected_ea = 9.0 / 5.0 * PI;

        let actual_ea = ma2ea(e, ma);

        assert_abs_diff_eq!(expected_ea, actual_ea, epsilon = 2.0 * MA_2_EA_ATOL);
    }

    #[test]
    fn test_ma2ea_highly_elliptical_1() {
        let e = 0.9;
        let ma = 0.021197602284445;
        // should the the same, its a circular orbit!
        let expected_ea = 0.2;

        let actual_ea = ma2ea(e, ma);

        assert_abs_diff_eq!(expected_ea, actual_ea, epsilon = 2.0 * MA_2_EA_ATOL);
    }

    #[test]
    fn test_ma2ea_highly_elliptical_2() {
        let e = 0.9;
        let ma = 0.021197602284445;
        // should the the same, its a circular orbit!
        let expected_ea = 0.2;

        let actual_ea = ma2ea(e, ma);

        assert_abs_diff_eq!(expected_ea, actual_ea, epsilon = 2.0 * MA_2_EA_ATOL);
    }
}
