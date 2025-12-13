#![allow(dead_code)]


const STUMPFF_TOL: f64 = 1e-2;

// factorials for readability
// pretty dumb but why not
const F2: f64 = 2.0;
const F3: f64 = 6.0;
const F4: f64 = 24.0;
const F5: f64 = 120.0;
const F6: f64 = 720.0;
const F7: f64 = 5040.0;
const F8: f64 = 40320.0;
const F9: f64 = 362880.0;
const F10: f64 = 3628800.0;
const F11: f64 = 39916800.0;


/// Stumpff function S(z)
pub fn stumpff_s(z: f64) -> f64 {    
    // Eq. 4-38 BMW
    if z > STUMPFF_TOL {
        let sqrt_z = z.sqrt();
        return (sqrt_z - sqrt_z.sin()) / sqrt_z.powi(3);

    }
    if z < -STUMPFF_TOL {
        let sqrt_neg_z = (-z).sqrt();
        return (sqrt_neg_z.sinh() - sqrt_neg_z) / sqrt_neg_z.powi(3)
    }
    // factorial apprixmation
    // z^0/3! - z^1/5! + z^2/7! - z^3/9!
    return 1.0/F3 - z/F5 + z.powi(2)/F7 - z.powi(3)/F9;
}

/// Stumpf function C(z)
pub fn stumpff_c(z: f64) -> f64 {
    // Eq. 4-37 BMW
    if z > STUMPFF_TOL {
        let sqrt_z = z.sqrt();
        return (1.0 - sqrt_z.cos()) / z;
    }
    if z < -STUMPFF_TOL {
        let sqrt_neg_z = (-z).sqrt();
        return (1.0 - sqrt_neg_z.cosh()) / z
    }
    // factorial approximation
    // z^0/2! - z^1/4! + z^2/6! - z^3/8!
    return 1.0/F2 - z/F4 + z.powi(2)/F6 - z.powi(3)/F8;
}

/// Stumpf function derivative S'(z)
pub fn stumpff_dsdz(z: f64) -> f64 {
    //
    if z.abs() < STUMPFF_TOL {
        // Eq. 5-35 BMW
        // use power series expansion for small z
        return -1.0/F5 + 2.0*z/F7 - 3.0*z.powi(2)/F9 + 4.0*z.powi(3)/F11;
    }
    
    let s = stumpff_s(z);
    let c = stumpff_c(z);

    // Eq. 5-30 BMW
    1.0/(2.0*z) * (c - 3.0*s)
}

/// Stumpf function derivative C'(z)
pub fn stumpff_dcdz(z: f64) -> f64 {
    if z.abs() < STUMPFF_TOL {
        // Eq. 5-34 BMW
        // use power series expansion for small z
        return -1.0/F4 + 2.0*z/F6 - 3.0*z.powi(2)/F8 + 4.0*z.powi(3)/F10;
    }
    
    let s = stumpff_s(z);
    let c = stumpff_c(z);

    // Eq. 5-31 BMW
    1.0/(2.0*z) * (1.0 - z*s - 2.0*c)
}




#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_abs_diff_eq};

    #[test]
    fn test_stumpff_s_0() {
        // z = 0 case
        let z: f64 = 0.0;
        let s_expected = 1.0/6.0;
        let s = stumpff_s(z);

        assert_abs_diff_eq!(s, s_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_s_gt0() {
        // z > 0 case
        let z: f64 = 0.0;
        let s_expected = 1.0/6.0;
        let s = stumpff_s(z);

        assert_abs_diff_eq!(s, s_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_s_lt0() {
        // z < 0 case
        let z: f64 = 0.0;
        let s_expected = 1.0/6.0;
        let s = stumpff_s(z);

        assert_abs_diff_eq!(s, s_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_c_0() {
        // z > 0 case
        let z: f64 = 0.0;
        let c_expected = 1.0/2.0;
        let c = stumpff_c(z);

        assert_abs_diff_eq!(c, c_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_c_gt0() {
        // z > 0 case
        let z: f64 = 1.0;
        let c_expected = 0.45969769413186023;
        let c = stumpff_c(z);

        assert_abs_diff_eq!(c, c_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_c_lt0() {
        // z < 0 case
        let z: f64 = -1.0;
        let c_expected = 0.5430806348152437;
        let c = stumpff_c(z);

        assert_abs_diff_eq!(c, c_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_sp_0() {
        // z = 0 case
        let z: f64 = 0.0;
        let sp_expected = -1.0/120.0;
        let sp = stumpff_dsdz(z);

        assert_abs_diff_eq!(sp, sp_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_sp_g0() {
        // z > 0 case
        let z: f64 = 1.0;
        let sp_expected = -0.007944675722225125;
        let sp = stumpff_dsdz(z);

        assert_abs_diff_eq!(sp, sp_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_cp_0() {
        // z = 0 case
        let z: f64 = 0.0;
        let cp_expected = -1.0/24.0;
        let cp = stumpff_dcdz(z);

        assert_abs_diff_eq!(cp, cp_expected, epsilon = 1e-12);
    }

    #[test]
    fn test_stumpff_cp_g0() {
        // z > 0 case
        let z: f64 = 1.0;
        let cp_expected = -0.03896220172791198;
        let cp = stumpff_dcdz(z);

        assert_abs_diff_eq!(cp, cp_expected, epsilon = 1e-12);
    }
}
