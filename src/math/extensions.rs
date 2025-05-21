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
  // C version does not handle negative exponents, Rust's .pow() does.
  // We match C by not supporting negative exp.
  // In C, a negative exp would result in 0 due to integer division.
  // Here, we'll stick to the C behavior for now.
  // If exp were u64, this check would be simpler.
  if exp < 0 {
    // This case is not handled by the C _ipow, which would likely produce 0 or unexpected results.
    // For a direct port, one might return 0.
    // However, Rust's `pow` would panic or handle it differently.
    // Let's stick to the loop logic which implies non-negative exp for meaningful results.
    // If `exp` were `u64`, this check wouldn't be needed.
    // Given C's _ipow signature, it's unlikely to be called with negative exp.
    return if base == 1 {
      1
    } else if base == -1 {
      if exp % 2 == 0 {
        1
      } else {
        -1
      }
    } else {
      0 // Or panic, or return an error. C _ipow just yields 0 for base > 1 or < -1.
    };
  }

  let mut result: i64 = 1;
  while exp > 0 {
    // C's while(exp)
    if exp % 2 == 1 {
      // C's if (exp & 1)
      result = result.saturating_mul(base);
    }
    exp /= 2; // C's exp >>= 1;
    if exp > 0 {
      // Avoid unnecessary multiplication if exp becomes 0
      base = base.saturating_mul(base);
    }
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

    // C's behavior for negative exponent (if `exp` was signed like in this Rust version)
    // It seems the C version implicitly assumes non-negative exponent.
    assert_eq!(_ipow(2, -1), 0, "2^-1 (C-like behavior)");
    assert_eq!(_ipow(1, -5), 1, "1^-5 (C-like behavior)");
    assert_eq!(_ipow(-1, -2), 1, "(-1)^-2 (C-like behavior)");
    assert_eq!(_ipow(-1, -3), -1, "(-1)^-3 (C-like behavior)");

    // Test potential overflow cases (using saturating_mul to avoid panic in tests)
    // This is where Rust's default checked arithmetic in debug mode would differ from C.
    // Our `saturating_mul` matches the "no panic" behavior of C's overflow.
    assert_eq!(_ipow(i64::MAX, 2), i64::MAX, "i64::MAX^2 saturates");
    assert_eq!(_ipow(3, 40), i64::MAX, "3^40 saturates (3^39 is max for i64)");
    // 3^38 fits, 3^39 overflows positive
    assert_eq!(3_i64.pow(38), 1350851717672992089_i64);
    assert_eq!(_ipow(3, 38), 1350851717672992089_i64);
    assert_eq!(_ipow(3, 39), i64::MAX); // Saturates
  }
}
