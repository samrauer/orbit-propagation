use std::ops::{Add, Mul};
use std::iter::zip;


// TODO: for all EVec operations, decide whether we want to consume
//       the two objects orjust borrow them to do our operations


#[derive(Debug, Clone, PartialEq)]
pub struct EVec<T>(pub Vec<T>);


// add self + another EVec
impl<T> Add for EVec<T>
where T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        assert_eq!(self.0.len(), rhs.0.len(), "EVecs must be same length");
        EVec(
            zip(self.0.iter(), rhs.0.iter())
                .map(|(&a, &b)| a + b)
                .collect()
        )
    }
}

// add self + scalar
impl<T> Add<T> for EVec<T>
where T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, rhs: T) -> Self::Output {
        EVec(
            self.0.into_iter()
                .map(|a| a + rhs)
                .collect()
        )
    }
}

// multiply self * another EVec
impl<T> Mul for EVec<T>
where T: Mul<Output = T> + Copy,
{
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        assert_eq!(self.0.len(), rhs.0.len(), "EVecs must be same length");
        EVec(
            zip(self.0.iter(), rhs.0.iter())
                .map(|(&a, &b)| a * b)
                .collect()
        )
    }
}

// multiply self + scalar
impl<T> Mul<T> for EVec<T>
where T: Mul<Output = T> + Copy,
{
    type Output = Self;

    fn mul(self, rhs: T) -> Self::Output {
        EVec(
            self.0.into_iter()
                .map(|a| a * rhs)
                .collect()
        )
    }
}
