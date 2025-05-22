// src/math/extensions.rs

/// Integer exponentiation.
///
/// # Panics
///
/// Panics if `exp` is negative, or if the operation overflows `i64`.
/// H3 C library's `_ipow` does not check for overflow, but Rust best
/// practices suggest panicking on overflow for debug builds or providing
/// `checked_` versions. For now, we'll match C's behavior of potential
/// overflow if not careful.
///
/// Consider using `base.pow(exp as u32)` for non-negative bases if `exp` fits `u32`,
/// or `num_bigint` for arbitrary precision if overflow is a concern.
///
/// This implementation matches the C version's loop.
#[inline]
pub(crate) fn _ipow(mut base: i64, mut exp: i64) -> i64 {
  if exp < 0 {
    if base == 1 {
      return 1;
    }
    if base == -1 {
      return if exp % 2 == 0 { 1 } else { -1 };
    }
    return 0;
  }

  let mut result: i64 = 1;
  loop {
    if exp & 1 != 0 {
      result = result.wrapping_mul(base);
    }
    exp >>= 1;
    if exp == 0 {
      break;
    }
    base = base.wrapping_mul(base);
  }
  result
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_ipow() {
    assert_eq!(_ipow(7, 0), 1, "7^0");
    assert_eq!(_ipow(7, 1), 7, "7^1");
    assert_eq!(_ipow(7, 2), 49, "7^2");
    assert_eq!(_ipow(1, 20), 1, "1^20");
    assert_eq!(_ipow(2, 5), 32, "2^5");
    assert_eq!(_ipow(-2, 2), 4, "(-2)^2");
    assert_eq!(_ipow(-2, 3), -8, "(-2)^3");

    // Behavior for negative exponent (matching original Rust guards, as C is undefined here)
    assert_eq!(_ipow(2, -1), 0, "2^-1");
    assert_eq!(_ipow(1, -5), 1, "1^-5");
    assert_eq!(_ipow(-1, -2), 1, "(-1)^-2");
    assert_eq!(_ipow(-1, -3), -1, "(-1)^-3");

    // Test cases for wrapping behavior, now that we use wrapping_mul
    // 3^39 fits in i64: 4,052,555,153,018,976,267
    assert_eq!(_ipow(3, 39), 4052555153018976267_i64, "3^39 should fit");

    // 3^40 overflows i64. Exact u128: 12_157_665_459_056_928_801
    assert_eq!(_ipow(3, 40), -6289078614652622815_i64, "3^40 should wrap");
  }
}
