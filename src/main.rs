mod state;
mod old;
mod integrator;

extern crate nalgebra as na;
use na::{Vector3, Rotation3};


fn main() {
    println!("Hello, world!");

    // testing1();
    // testing2();
    testing3();
}


#[allow(dead_code)]
fn testing1() {
    let f = |x: &Vec<f64>| -> Vec<f64> {
        x.iter().map(|xi| -xi).collect()
    };

    let x0 = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    // forward euler
    let x1_euler = old::integrator::integrate_euler(f, &x0, 0.1);
    let x2_euler = old::integrator::integrate_euler(f, &x1_euler, 0.1);

    // predictor corrector
    let x1_pc = old::integrator::integrate_euler_with_corrector(f, &x0, 0.1);
    let x2_pc = old::integrator::integrate_euler_with_corrector(f, &x1_pc, 0.1);

    println!("{:?}", x0);
    println!("{:?}", x1_euler);
    println!("{:?}", x2_euler);

    println!("{:?}", x0);
    println!("{:?}", x1_pc);
    println!("{:?}", x2_pc);
}


#[allow(dead_code)]
fn testing2() {
    let axis = Vector3::<f64>::x_axis();
    let angle: f64 = 1.0;
    let b = Rotation3::from_axis_angle(&axis, angle);

    let c = Vector3::new(1.0, 2.0, 3.0);

    let d= b * c;

    println!("{}", b);
    println!("{}", d);
}


#[allow(dead_code)]
fn testing3() {
    let f = |x: &Vector3<f64>| -x;

    let x0 = Vector3::new(10.0, 1.0, -5.0);

    // forward euler
    let x1_euler = integrator::euler::integrate_euler(f, &x0, 0.1);
    let x2_euler = integrator::euler::integrate_euler(f, &x1_euler, 0.1);

    // predictor corrector
    let x1_pc = integrator::euler_pc::integrate_euler_pc(f, &x0, 0.1);
    let x2_pc = integrator::euler_pc::integrate_euler_pc(f, &x1_pc, 0.1);

    println!("{:?}", x0);
    println!("{:?}", x1_euler);
    println!("{:?}", x2_euler);

    println!("{:?}", x0);
    println!("{:?}", x1_pc);
    println!("{:?}", x2_pc);
}
