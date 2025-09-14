mod integrator;
mod state;

fn main() {
    println!("Hello, world!");

    let f = |x: &Vec<f64>| -> Vec<f64> {
        x.iter().map(|xi| -xi).collect()
    };

    let x0 = vec![1.0, 2.0, 3.0, 4.0, 5.0];

    // forward euler
    let x1_euler = integrator::integrate_euler(f, &x0, 0.1);
    let x2_euler = integrator::integrate_euler(f, &x1_euler, 0.1);

    // predictor corrector
    let x1_pc = integrator::integrate_euler_with_corrector(f, &x0, 0.1);
    let x2_pc = integrator::integrate_euler_with_corrector(f, &x1_pc, 0.1);

    println!("{:?}", x0);
    println!("{:?}", x1_euler);
    println!("{:?}", x2_euler);

    println!("{:?}", x0);
    println!("{:?}", x1_pc);
    println!("{:?}", x2_pc);

}


