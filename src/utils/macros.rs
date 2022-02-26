/// Generalization over [std::cmp::min]
/// ```
/// use wfa::min;
///
/// assert_eq!(min![45, 56, 89], 45);
/// ```
#[macro_export]
macro_rules! min {
    () => ({});
    ($x: expr) => ($x);
    ($x: expr, $($xs: expr),+) => ({ ::std::cmp::min($x, min!($($xs),*)) });
}

/// Generalization over [std::cmp::max]
/// # Example
/// ```
/// use wfa::max;
///
/// let x = 45;
/// let y = 23;
/// let z = 73;
/// assert_eq!(max![&x, &y, &z], &73);
/// ```
#[macro_export]
macro_rules! max {
    () => ({});
    ($x: expr) => ($x);
    ($x: expr, $($xs: expr),+) => ({ ::std::cmp::max($x, max!($($xs),*)) });
}

pub(crate) use {max, min};


#[cfg(test)]
mod tests {
    use super::*;
    use pretty_assertions::{assert_eq};

    #[test]
    fn test_max() {
        // max![];
        assert_eq!(max![], ());
        assert_eq!(max![23], 23);
        assert_eq!(max![45, 56], 56);
        assert_eq!(max![45, 56, 89], 89);

        assert_eq!(max![Some(5), None, None], Some(5));
    }

    #[test]
    fn test_min() {
        // min![];
        assert_eq!(min![], ());
        assert_eq!(min![23], 23);
        assert_eq!(min![45, 56], 45);
        assert_eq!(min![45, 56, 89], 45);
        assert_eq!(min![Some(5), None, None], None);

        // refs
        let x = 45;
        let y = 23;
        assert_eq!(min![&x, &y], &23);
    }

}
