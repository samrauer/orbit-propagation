use std::ops::{Add, Mul, Sub, Div};

pub trait State:
    Clone
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    fn norm(&self) -> f64;
    fn sum(&self) -> f64;

    fn map<F>(&self, f: F) -> Self
    where
        F: FnMut(f64) -> f64;

    fn zip_map<F>(&self, rhs: &Self, f: F) -> Self
    where 
        F: FnMut(f64, f64) -> f64;
}


impl State for f64 {
    fn norm(&self) -> f64 {
        self.abs()
    }

    fn sum(&self) -> f64 {
        *self
    }
    
    fn map<F>(&self, mut f: F) -> Self 
    where
        F: FnMut(f64) -> f64
    {
        f(*self)
    }

    fn zip_map<F>(&self, rhs: &Self, mut f: F) -> Self
        where 
            F: FnMut(f64, f64) -> f64 
    {
        f(*self, *rhs)
    }
}







