mod integrator;



fn main() {
    println!("Hello, world!");

    let f = |x: &Vec<f64>| -> Vec<f64> {
        x.iter().map(|xi| -xi).collect()
    };

    let x0 = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let x1 = integrator::integrate_euler(f, &x0, 0.1);
    let x2 = integrator::integrate_euler(f, &x1, 0.1);

    println!("{:?}", x0);
    println!("{:?}", x1);
    println!("{:?}", x2);

}


