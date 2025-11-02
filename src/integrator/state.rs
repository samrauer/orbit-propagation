use std::ops::{Add, Mul, Sub, Div};
use na::{DefaultAllocator, Dim, OVector, allocator::Allocator};

pub trait State:
    Clone
    + Add<Self, Output = Self>
    + Sub<Self, Output = Self>
    + Mul<f64, Output = Self>
    + Div<f64, Output = Self>
{
    fn norm(&self) -> f64;
    fn sum(&self) -> f64;
    fn length(&self) -> u64;

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

    fn length(&self) -> u64 {
        1
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



impl<D> State for OVector<f64, D> 
where
    D: Dim,
    DefaultAllocator: Allocator<D>,
{
    fn norm(&self) -> f64 {
        self.norm()
    }

    fn sum(&self) -> f64 {
        self.sum()
    }

    fn length(&self) -> u64 {
        self.nrows() as u64
    }
    
    fn map<F>(&self, mut f: F) -> Self
    where 
        F: FnMut(f64) -> f64,
    {
        self.map(|x| f(x))
    }

    fn zip_map<F>(&self, other: &Self, mut f: F) -> Self
    where 
        F: FnMut(f64, f64) -> f64,
    {
        self.zip_map(other, |x, y| f(x, y))
    }
}

