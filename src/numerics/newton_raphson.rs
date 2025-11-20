#![allow(dead_code)]


const MAX_ITER: i64 = 1000;



pub fn newton_raphson_scalar<F1,F2>(f: F1, fp: F2, x0: f64, atol: f64) -> f64 
where 
    F1: Fn(f64) -> f64,
    F2: Fn(f64) -> f64
{
    let mut x1: f64 = x0 - 2.0*atol; // force while loop condition to be false
    let mut x2: f64 = x0;

    let mut counter: i64 = MAX_ITER;

    while (x2 - x1).abs() > atol {
        if counter <= 0 {
            // panic for now, add error handling in the future when needed
            panic!("MAX_ITER reached")
        }
        counter -= 1;

        x1 = x2;

        let fx = f(x1);
        let fpx = fp(x2);

        x2 = x1 - fx / fpx;

        // TODO: check that |fp| > 1e-10 (far away from 0) before doing this divide
    }

    x2
}


#[cfg(test)]
mod tests {
    use super::*;
    use approx::{assert_abs_diff_eq};

    /// Solve x^2 - 2 = 0 -> root = sqrt(2)
    fn f_sqrt2(x: f64) -> f64 {
        x * x - 2.0
    }
    fn fp_sqrt2(x: f64) -> f64 {
        2.0 * x
    }

    #[test]
    fn finds_sqrt2() {
        let x0 = 1.0; // initial guess
        let atol = 1e-12;
        let root = newton_raphson_scalar(f_sqrt2 as fn(f64) -> f64, fp_sqrt2 as fn(f64) -> f64, x0, atol);
        let expected = 2f64.sqrt();
        
        assert_abs_diff_eq!(root, expected, epsilon = 1e-12);
    }

    /// Solve x^3 - 2 = 0 -> real root = 2^(1/3)
    fn f_cuberoot(x: f64) -> f64 {
        x * x * x - 2.0
    }
    fn fp_cuberoot(x: f64) -> f64 {
        3.0 * x * x
    }

    #[test]
    fn finds_cuberoot_of_2() {
        let x0 = 1.5;
        let atol = 1e-12;
        let root = newton_raphson_scalar(
            f_cuberoot as fn(f64) -> f64,
            fp_cuberoot as fn(f64) -> f64,
            x0,
            atol,
        );
        let expected = 2f64.powf(1.0 / 3.0);

        assert_abs_diff_eq!(root, expected, epsilon = 1e-12)
    }
}
