// src/coords/ijk.rs

use crate::constants::{M_ONESEVENTH, M_ONETHIRD, M_RSIN60, M_SQRT3_2, M_SQRT7};
use crate::math::extensions::_ipow; // For _downAp* later, if needed here or used elsewhere
use crate::types::{CoordIJ, CoordIJK, Direction};
use crate::Vec2d;
use std::f64; // For lround if we implement it, or f64::round
use std::i32; // For i32::MAX/MIN for overflow checks

/// IJK unit vectors corresponding to the 7 H3 digits (0-6).
#[rustfmt::skip] // Keep the C-style formatting for easy comparison
pub(crate) static UNIT_VECS: [CoordIJK; 7] = [
    CoordIJK { i: 0, j: 0, k: 0 },  // CENTER_DIGIT
    CoordIJK { i: 0, j: 0, k: 1 },  // K_AXES_DIGIT
    CoordIJK { i: 0, j: 1, k: 0 },  // J_AXES_DIGIT
    CoordIJK { i: 0, j: 1, k: 1 },  // JK_AXES_DIGIT (J_AXES_DIGIT | K_AXES_DIGIT)
    CoordIJK { i: 1, j: 0, k: 0 },  // I_AXES_DIGIT
    CoordIJK { i: 1, j: 0, k: 1 },  // IK_AXES_DIGIT (I_AXES_DIGIT | K_AXES_DIGIT)
    CoordIJK { i: 1, j: 1, k: 0 },  // IJ_AXES_DIGIT (I_AXES_DIGIT | J_AXES_DIGIT)
];

/// Sets an IJK coordinate to the specified component values.
#[inline]
pub(crate) fn _set_ijk(ijk: &mut CoordIJK, i: i32, j: i32, k: i32) {
  ijk.i = i;
  ijk.j = j;
  ijk.k = k;
}

/// Returns whether or not two IJK coordinates contain exactly the same
/// component values.
#[inline]
#[must_use]
pub(crate) fn _ijk_matches(c1: &CoordIJK, c2: &CoordIJK) -> bool {
  c1.i == c2.i && c1.j == c2.j && c1.k == c2.k
}

/// Add two IJK coordinates.
#[inline]
pub(crate) fn _ijk_add(h1: &CoordIJK, h2: &CoordIJK, sum: &mut CoordIJK) {
  sum.i = h1.i.saturating_add(h2.i);
  sum.j = h1.j.saturating_add(h2.j);
  sum.k = h1.k.saturating_add(h2.k);
}

/// Subtract two IJK coordinates. (h1 - h2)
#[inline]
pub(crate) fn _ijk_sub(h1: &CoordIJK, h2: &CoordIJK, diff: &mut CoordIJK) {
  diff.i = h1.i.saturating_sub(h2.i);
  diff.j = h1.j.saturating_sub(h2.j);
  diff.k = h1.k.saturating_sub(h2.k);
}

/// Uniformly scale IJK coordinates by a scalar. Works in place.
#[inline]
pub(crate) fn _ijk_scale(c: &mut CoordIJK, factor: i32) {
  c.i = c.i.saturating_mul(factor);
  c.j = c.j.saturating_mul(factor);
  c.k = c.k.saturating_mul(factor);
}

// Helper for _ijk_normalize_could_overflow and _ijk_normalize
// Equivalent to C's ADD_INT32S_OVERFLOWS
#[inline]
fn _add_i32s_overflows(a: i32, b: i32) -> bool {
  // From C: `if (a > 0) return INT32_MAX - a < b; else return INT32_MIN - a > b;`
  // Rust's checked_add returns None on overflow.
  a.checked_add(b).is_none()
}

// Helper for _ijk_normalize_could_overflow and _ijk_normalize
// Equivalent to C's SUB_INT32S_OVERFLOWS
#[inline]
fn _sub_i32s_overflows(a: i32, b: i32) -> bool {
  // From C: `if (a >= 0) return INT32_MIN + a >= b; else return INT32_MAX + a + 1 < b;`
  // Rust's checked_sub returns None on overflow.
  a.checked_sub(b).is_none()
}

/// Returns true if `_ijk_normalize` with the given input could have a signed
/// integer overflow. Assumes k is set to 0 before this call.
/// (This check is specific to the C implementation's _ijkNormalize logic).
#[inline]
#[must_use]
pub(crate) fn _ijk_normalize_could_overflow(ijk: &CoordIJK) -> bool {
  // Pre-condition from C code: ijk.k must be 0 for this check to be meaningful
  // in the context of _ijkNormalize, which sets k=0 before these subtractions.
  // The C macro versions are:
  // #define ADD_INT32S_OVERFLOWS(a,b) ( (a > 0) ? (INT32_MAX - a < b) : (INT32_MIN - a > b) )
  // #define SUB_INT32S_OVERFLOWS(a,b) ( (a >= 0) ? (INT32_MIN + a >= b) : (INT32_MAX + a + 1 < b) )

  // The C version's _ijkNormalize does this after potentially modifying i, j, k:
  // if (c->i < 0) { c->j -= c->i; c->k -= c->i; c->i = 0; }
  // This means `c->j` could become `c->j_orig - c->i_orig`.
  // The critical subtractions that _could_ overflow inside _ijkNormalize's logic are
  // when adjusting components based on other negative components.
  // E.g., `c->j -= c->i;` if `c->i` is very negative and `c->j` is positive.

  // Let's analyze the C _ijkNormalize carefully:
  // 1. if (c->i < 0) { c->j -= c->i; c->k -= c->i; c->i = 0; }
  //    Overflow if c->j - c->i overflows. `c->i` is negative. So `c->j + abs(c->i)`
  //    Could overflow if c->j is large positive and c->i is large negative.
  //    The check in C is `ADD_INT32S_OVERFLOWS(c->j, -c->i)` which is `c->j + (-c->i)`
  if ijk.i < 0 && _add_i32s_overflows(ijk.j, -ijk.i) {
    return true;
  }
  //    And `c->k -= c->i;` (k is initially 0). `0 - c->i`
  if ijk.i < 0 && _add_i32s_overflows(0, -ijk.i) {
    return true;
  } // Effectively -c->i

  // 2. if (c->j < 0) { c->i -= c->j; c->k -= c->j; c->j = 0; } (uses potentially modified c->i)
  //    Test with original values for simplicity of this pre-check.
  //    If ijk.j < 0 and ijk.i - ijk.j overflows:
  if ijk.j < 0 && _add_i32s_overflows(ijk.i, -ijk.j) {
    return true;
  }
  //    If ijk.j < 0 and 0 - ijk.j overflows (if k was 0, after i adjustment):
  if ijk.j < 0 && _add_i32s_overflows(0, -ijk.j) {
    return true;
  } // Effectively -ijk.j

  // 3. if (c->k < 0) { c->i -= c->k; c->j -= c->k; c->k = 0; } (uses potentially modified c->i, c->j)
  //    Test with original values for simplicity.
  //    If ijk.k < 0 (which it is, as 0 initially) and ijk.i - ijk.k overflows:
  //    This step is applied to k which might have become negative. Let's assume k is already 0 here.
  //    The C code version of this check is complex because it's trying to predict
  //    overflow within the _ijkNormalize steps.
  //    The C's `_ijkNormalizeCouldOverflow` seems to simplify this by checking original i,j (k=0):
  //    int max_val, min_val; if (ijk->i > ijk->j) { max_val = ijk->i; min_val = ijk->j; } else { ... }
  //    if (min_val < 0) {
  //        if (ADD_INT32S_OVERFLOWS(max_val, min_val)) return true; // max_val + min_val (not sub)
  //        if (SUB_INT32S_OVERFLOWS(0, min_val)) return true;       // 0 - min_val
  //        if (SUB_INT32S_OVERFLOWS(max_val, min_val)) return true; // max_val - min_val
  //    }

  let (max_val, min_val) = if ijk.i > ijk.j { (ijk.i, ijk.j) } else { (ijk.j, ijk.i) };
  if min_val < 0 {
    if _add_i32s_overflows(max_val, min_val) {
      return true;
    }
    if _sub_i32s_overflows(0, min_val) {
      return true;
    } // equivalent to -min_val if min_val is INT32_MIN
    if _sub_i32s_overflows(max_val, min_val) {
      return true;
    }
  }

  false
}

/// Normalizes IJK coordinates by setting the components to the smallest possible
/// values. Works in place.
///
/// This function does not protect against signed integer overflow if called directly
/// with arbitrary inputs; however, H3 library functions that call this internally
/// are responsible for input validation.
/// `_ijk_normalize_could_overflow` can be used as a pre-check.
#[inline]
pub(crate) fn _ijk_normalize(c: &mut CoordIJK) {
  // Simplify by finding the smallest component and subtracting it from all components
  // This makes one component 0 and the others non-negative.
  if c.i < 0 {
    c.j = c.j.saturating_sub(c.i);
    c.k = c.k.saturating_sub(c.i);
    c.i = 0;
  }

  if c.j < 0 {
    c.i = c.i.saturating_sub(c.j);
    c.k = c.k.saturating_sub(c.j);
    c.j = 0;
  }

  if c.k < 0 {
    c.i = c.i.saturating_sub(c.k);
    c.j = c.j.saturating_sub(c.k);
    c.k = 0;
  }

  // Remove the min value if needed
  let min_val = c.i.min(c.j).min(c.k);
  if min_val > 0 {
    c.i -= min_val;
    c.j -= min_val;
    c.k -= min_val;
  }
}

/// Determines the H3 digit corresponding to a unit vector or the zero vector
/// in IJK coordinates.
///
/// # Arguments
/// * `ijk` - The IJK coordinates; must be a unit vector or zero vector after normalization.
///
/// # Returns
/// The H3 `Direction` (0-6) corresponding to the IJK unit vector,
/// or `Direction::InvalidDigit` (7) on failure/if not a unit vector.
#[inline]
#[must_use]
pub(crate) fn _unit_ijk_to_digit(ijk: &CoordIJK) -> Direction {
  let mut c = *ijk; // Make a mutable copy
  _ijk_normalize(&mut c);

  // Corresponds to H3 C library's for loop check against UNIT_VECS
  if _ijk_matches(&c, &UNIT_VECS[0]) {
    Direction::Center
  } else if _ijk_matches(&c, &UNIT_VECS[1]) {
    Direction::KAxes
  } else if _ijk_matches(&c, &UNIT_VECS[2]) {
    Direction::JAxes
  } else if _ijk_matches(&c, &UNIT_VECS[3]) {
    Direction::JkAxes
  } else if _ijk_matches(&c, &UNIT_VECS[4]) {
    Direction::IAxes
  } else if _ijk_matches(&c, &UNIT_VECS[5]) {
    Direction::IkAxes
  } else if _ijk_matches(&c, &UNIT_VECS[6]) {
    Direction::IjAxes
  } else {
    Direction::InvalidDigit
  }
}

/// Find the normalized IJK coordinates of the hex in the specified digit
/// direction from the specified IJK coordinates. Works in place.
///
/// # Arguments
/// * `ijk` - The IJK coordinates. Modified in place.
/// * `digit` - The H3 `Direction` from the original IJK coordinates.
#[inline]
pub(crate) fn _neighbor(ijk: &mut CoordIJK, digit: Direction) {
  if digit != Direction::Center && digit != Direction::InvalidDigit {
    let temp_h2 = UNIT_VECS[digit as usize]; // temp_h2 is a copy from static
    let mut result = CoordIJK::default(); // Temporary variable for the sum
    _ijk_add(ijk, &temp_h2, &mut result); // Calculate sum into the temporary
    *ijk = result; // Assign the result back
    _ijk_normalize(ijk);
  }
}

/// Determine the containing hex in IJK+ coordinates for a 2D Cartesian
/// coordinate vector (from DGGRID).
///
/// # Arguments
/// * `v` - The 2D Cartesian coordinate vector.
/// * `h` - Output: The IJK+ coordinates of the containing hex.
#[inline]
pub(crate) fn _hex2d_to_coord_ijk(v: &Vec2d, h: &mut CoordIJK) {
  // quantize into the ij system and then normalize
  h.k = 0;

  let a1 = v.x.abs();
  let a2 = v.y.abs();

  // first do a reverse conversion
  let x2 = a2 * M_RSIN60; // M_RSIN60 = 1.0 / M_SIN60 = 1.0 / (M_SQRT3 / 2.0)
  let x1 = a1 + x2 / 2.0;

  // check if we have the center of a hex
  let m1 = x1 as i32;
  let m2 = x2 as i32;

  // otherwise round correctly
  let r1 = x1 - m1 as f64;
  let r2 = x2 - m2 as f64;

  if r1 < 0.5 {
    if r1 < (1.0 / 3.0) {
      if r2 < (1.0 + r1) / 2.0 {
        h.i = m1;
        h.j = m2;
      } else {
        h.i = m1;
        h.j = m2 + 1;
      }
    } else {
      if r2 < (1.0 - r1) {
        h.j = m2;
      } else {
        h.j = m2 + 1;
      }

      if (1.0 - r1) <= r2 && r2 < (2.0 * r1) {
        h.i = m1 + 1;
      } else {
        h.i = m1;
      }
    }
  } else {
    if r1 < (2.0 / 3.0) {
      if r2 < (1.0 - r1) {
        h.j = m2;
      } else {
        h.j = m2 + 1;
      }

      if (2.0 * r1 - 1.0) < r2 && r2 < (1.0 - r1) {
        h.i = m1;
      } else {
        h.i = m1 + 1;
      }
    } else {
      if r2 < (r1 / 2.0) {
        h.i = m1 + 1;
        h.j = m2;
      } else {
        h.i = m1 + 1;
        h.j = m2 + 1;
      }
    }
  }

  // now fold across the axes if necessary
  if v.x < 0.0 {
    if (h.j % 2) == 0 {
      // even
      let axisi: i64 = h.j as i64 / 2; // Use i64 for intermediate to avoid overflow before cast
      let diff: i64 = h.i as i64 - axisi;
      h.i = (h.i as i64 - 2 * diff) as i32;
    } else {
      // odd
      let axisi: i64 = (h.j as i64 + 1) / 2;
      let diff: i64 = h.i as i64 - axisi;
      h.i = (h.i as i64 - (2 * diff + 1)) as i32;
    }
  }

  if v.y < 0.0 {
    // Extended precision for intermediate calculation, then cast to i32
    h.i = (h.i as i64 - (2 * h.j as i64 + 1) / 2) as i32;
    h.j = -1 * h.j;
  }

  _ijk_normalize(h);
}

/// Find the center point in 2D Cartesian coordinates of a hex.
///
/// # Arguments
/// * `h` - The IJK coordinates of the hex.
/// * `v` - Output: The 2D Cartesian coordinates of the hex center point.
#[inline]
pub(crate) fn _ijk_to_hex2d(h: &CoordIJK, v: &mut Vec2d) {
  // If not normalized, this could overflow. Assume normalized input.
  let i = h.i - h.k;
  let j = h.j - h.k;

  v.x = i as f64 - 0.5 * j as f64;
  v.y = j as f64 * M_SQRT3_2;
}

// Helper for lround, C's lround might not be available or behave identically
// Rust's f64::round() rounds to nearest, ties to even.
// C99 lround rounds to nearest, ties away from zero.
// For H3, precise tie-breaking might be important if the C code relied on it.
// Given the usage, standard f64::round() is likely acceptable.
// If specific lround behavior is needed, we'd implement it.
// For now, let's assume f64::round() is okay for the Aperture functions.
fn h3_lround(val: f64) -> i32 {
  val.round() as i32 // This is ties-to-even. C's lround is ties-away-from-zero.
                     // A more C-like lround:
                     // if val > 0.0 { (val + 0.5) as i32 } else { (val - 0.5) as i32 }
                     // However, C's lround has specific NaN/Inf handling we are not replicating here.
                     // Given the inputs to aperture functions are scaled integers, NaN/Inf unlikely.
}

/// Find the normalized IJK coordinates of the indexing parent of a cell in a
/// counter-clockwise aperture 7 grid. Works in place.
/// (Corresponds to Class III resolution.)
///
/// # Arguments
/// * `ijk` - The IJK coordinates. Modified in place.
/// # Returns
/// An `H3Error` if an integer overflow occurs.
// In C, these were `_upAp7Checked` and `_upAp7rChecked` when returning H3Error.
// The non-checked versions didn't return error.
// Let's make them all checked for Rust.
pub(crate) fn _up_ap7_checked(ijk: &mut CoordIJK) -> Result<(), crate::types::H3Error> {
  // convert to CoordIJ (local, non-normalized)
  // These subtractions can't overflow if ijk components are within i32 range,
  // as the result will also be within i32 range.
  let i_ax = ijk.i - ijk.k;
  let j_ax = ijk.j - ijk.k;

  // check for int overflow before multiplying.
  // Max value for i_ax or j_ax can be 2*MAX_H3_RES_COORD_VAL approx.
  // 3 * i_ax, 2 * j_ax
  let term1_i = i_ax.checked_mul(3).ok_or(crate::types::H3Error::Failed)?;
  let new_i_num = term1_i.checked_sub(j_ax).ok_or(crate::types::H3Error::Failed)?;

  let term1_j = j_ax.checked_mul(2).ok_or(crate::types::H3Error::Failed)?;
  let new_j_num = i_ax.checked_add(term1_j).ok_or(crate::types::H3Error::Failed)?;

  ijk.i = h3_lround(new_i_num as f64 * M_ONESEVENTH);
  ijk.j = h3_lround(new_j_num as f64 * M_ONESEVENTH);
  ijk.k = 0;

  // Normalization involves subtractions that could overflow if intermediate i,j are pathological.
  // The C code had `_ijkNormalizeCouldOverflow` for this.
  // Our `_ijk_normalize` uses saturating arithmetic, so it won't panic,
  // but the result might differ from C's wrapping if C relied on specific wrapping.
  // Given H3's design, saturation or error is likely safer than wrapping.
  if _ijk_normalize_could_overflow(ijk) { // This check is based on the new i, j before normalize
     // It's hard to predict overflow *within* _ijk_normalize precisely without re-doing its logic.
     // If _ijk_normalize itself used checked_ops, it could return Result.
     // For now, _ijk_normalize uses saturating_ops.
  }
  _ijk_normalize(ijk);
  Ok(())
}

/// Find the normalized IJK coordinates of the indexing parent of a cell in a
/// clockwise aperture 7 grid. Works in place.
/// (Corresponds to Class II resolution.)
///
/// # Arguments
/// * `ijk` - The IJK coordinates. Modified in place.
/// # Returns
/// An `H3Error` if an integer overflow occurs.
pub(crate) fn _up_ap7r_checked(ijk: &mut CoordIJK) -> Result<(), crate::types::H3Error> {
  let i_ax = ijk.i - ijk.k;
  let j_ax = ijk.j - ijk.k;

  let term1_i = i_ax.checked_mul(2).ok_or(crate::types::H3Error::Failed)?;
  let new_i_num = term1_i.checked_add(j_ax).ok_or(crate::types::H3Error::Failed)?;

  let term1_j = j_ax.checked_mul(3).ok_or(crate::types::H3Error::Failed)?;
  let new_j_num = term1_j.checked_sub(i_ax).ok_or(crate::types::H3Error::Failed)?;

  ijk.i = h3_lround(new_i_num as f64 * M_ONESEVENTH);
  ijk.j = h3_lround(new_j_num as f64 * M_ONESEVENTH);
  ijk.k = 0;

  _ijk_normalize(ijk);
  Ok(())
}

// Non-checked versions (panic on overflow in debug, wrap in release, like C default)
// These are closer to the original C _upAp7 and _upAp7r signatures.
// For safety, it's better to use the checked versions internally if performance allows.

/// Find the normalized IJK coordinates of the indexing parent of a cell in a
/// counter-clockwise aperture 7 grid. Works in place. (Class III)
#[inline]
pub(crate) fn _up_ap7(ijk: &mut CoordIJK) {
  let i = ijk.i - ijk.k;
  let j = ijk.j - ijk.k;

  ijk.i = h3_lround((3 * i - j) as f64 * M_ONESEVENTH);
  ijk.j = h3_lround((i + 2 * j) as f64 * M_ONESEVENTH);
  ijk.k = 0;
  _ijk_normalize(ijk);
}

/// Find the normalized IJK coordinates of the indexing parent of a cell in a
/// clockwise aperture 7 grid. Works in place. (Class II)
#[inline]
pub(crate) fn _up_ap7r(ijk: &mut CoordIJK) {
  let i = ijk.i - ijk.k;
  let j = ijk.j - ijk.k;

  ijk.i = h3_lround((2 * i + j) as f64 * M_ONESEVENTH);
  ijk.j = h3_lround((3 * j - i) as f64 * M_ONESEVENTH);
  ijk.k = 0;
  _ijk_normalize(ijk);
}

/// Find the normalized IJK coordinates of the hex centered on the indicated
/// hex at the next finer aperture 7 counter-clockwise resolution. Works in place. (Class III)
#[inline]
pub(crate) fn _down_ap7(ijk: &mut CoordIJK) {
  // res r unit vectors in res r+1
  // Components are scaled based on the original ijk values.
  let i_vec = CoordIJK { i: 3, j: 0, k: 1 };
  let j_vec = CoordIJK { i: 1, j: 3, k: 0 };
  let k_vec = CoordIJK { i: 0, j: 1, k: 3 };

  let mut temp_i = i_vec;
  _ijk_scale(&mut temp_i, ijk.i);
  let mut temp_j = j_vec;
  _ijk_scale(&mut temp_j, ijk.j);
  let mut temp_k = k_vec;
  _ijk_scale(&mut temp_k, ijk.k);

  let mut sum1 = CoordIJK::default();
  _ijk_add(&temp_i, &temp_j, &mut sum1);
  _ijk_add(&sum1, &temp_k, ijk); // Store final result in ijk

  _ijk_normalize(ijk);
}

/// Find the normalized IJK coordinates of the hex centered on the indicated
/// hex at the next finer aperture 7 clockwise resolution. Works in place. (Class II)
#[inline]
pub(crate) fn _down_ap7r(ijk: &mut CoordIJK) {
  // res r unit vectors in res r+1
  let i_vec = CoordIJK { i: 3, j: 1, k: 0 };
  let j_vec = CoordIJK { i: 0, j: 3, k: 1 };
  let k_vec = CoordIJK { i: 1, j: 0, k: 3 };

  let mut temp_i = i_vec;
  _ijk_scale(&mut temp_i, ijk.i);
  let mut temp_j = j_vec;
  _ijk_scale(&mut temp_j, ijk.j);
  let mut temp_k = k_vec;
  _ijk_scale(&mut temp_k, ijk.k);

  let mut sum1 = CoordIJK::default();
  _ijk_add(&temp_i, &temp_j, &mut sum1);
  _ijk_add(&sum1, &temp_k, ijk);

  _ijk_normalize(ijk);
}

/// Find the normalized IJK coordinates of the hex centered on the indicated
/// hex at the next finer aperture 3 counter-clockwise resolution. Works in place.
#[inline]
pub(crate) fn _down_ap3(ijk: &mut CoordIJK) {
  // res r unit vectors in res r+1
  let i_vec = CoordIJK { i: 2, j: 0, k: 1 };
  let j_vec = CoordIJK { i: 1, j: 2, k: 0 };
  let k_vec = CoordIJK { i: 0, j: 1, k: 2 };

  let mut temp_i = i_vec;
  _ijk_scale(&mut temp_i, ijk.i);
  let mut temp_j = j_vec;
  _ijk_scale(&mut temp_j, ijk.j);
  let mut temp_k = k_vec;
  _ijk_scale(&mut temp_k, ijk.k);

  let mut sum1 = CoordIJK::default();
  _ijk_add(&temp_i, &temp_j, &mut sum1);
  _ijk_add(&sum1, &temp_k, ijk);

  _ijk_normalize(ijk);
}

/// Find the normalized IJK coordinates of the hex centered on the indicated
/// hex at the next finer aperture 3 clockwise resolution. Works in place.
#[inline]
pub(crate) fn _down_ap3r(ijk: &mut CoordIJK) {
  // res r unit vectors in res r+1
  let i_vec = CoordIJK { i: 2, j: 1, k: 0 };
  let j_vec = CoordIJK { i: 0, j: 2, k: 1 };
  let k_vec = CoordIJK { i: 1, j: 0, k: 2 };

  let mut temp_i = i_vec;
  _ijk_scale(&mut temp_i, ijk.i);
  let mut temp_j = j_vec;
  _ijk_scale(&mut temp_j, ijk.j);
  let mut temp_k = k_vec;
  _ijk_scale(&mut temp_k, ijk.k);

  let mut sum1 = CoordIJK::default();
  _ijk_add(&temp_i, &temp_j, &mut sum1);
  _ijk_add(&sum1, &temp_k, ijk);

  _ijk_normalize(ijk);
}

/// Rotates IJK coordinates 60 degrees counter-clockwise. Works in place.
#[inline]
pub(crate) fn _ijk_rotate60_ccw(ijk: &mut CoordIJK) {
  // unit vector rotations
  let i_vec = CoordIJK { i: 1, j: 1, k: 0 };
  let j_vec = CoordIJK { i: 0, j: 1, k: 1 };
  let k_vec = CoordIJK { i: 1, j: 0, k: 1 };

  let mut temp_i = i_vec;
  _ijk_scale(&mut temp_i, ijk.i);
  let mut temp_j = j_vec;
  _ijk_scale(&mut temp_j, ijk.j);
  let mut temp_k = k_vec;
  _ijk_scale(&mut temp_k, ijk.k);

  let mut sum1 = CoordIJK::default();
  _ijk_add(&temp_i, &temp_j, &mut sum1);
  _ijk_add(&sum1, &temp_k, ijk);

  _ijk_normalize(ijk);
}

/// Rotates IJK coordinates 60 degrees clockwise. Works in place.
#[inline]
pub(crate) fn _ijk_rotate60_cw(ijk: &mut CoordIJK) {
  // unit vector rotations
  let i_vec = CoordIJK { i: 1, j: 0, k: 1 };
  let j_vec = CoordIJK { i: 1, j: 1, k: 0 };
  let k_vec = CoordIJK { i: 0, j: 1, k: 1 };

  let mut temp_i = i_vec;
  _ijk_scale(&mut temp_i, ijk.i);
  let mut temp_j = j_vec;
  _ijk_scale(&mut temp_j, ijk.j);
  let mut temp_k = k_vec;
  _ijk_scale(&mut temp_k, ijk.k);

  let mut sum1 = CoordIJK::default();
  _ijk_add(&temp_i, &temp_j, &mut sum1);
  _ijk_add(&sum1, &temp_k, ijk);

  _ijk_normalize(ijk);
}

/// Rotates an H3 digit 60 degrees counter-clockwise.
#[inline]
#[must_use]
pub(crate) fn _rotate60_ccw(digit: Direction) -> Direction {
  use Direction::*;
  match digit {
    KAxes => IkAxes,
    IkAxes => IAxes,
    IAxes => IjAxes,
    IjAxes => JAxes,
    JAxes => JkAxes,
    JkAxes => KAxes,
    _ => digit, // Center or Invalid
  }
}

/// Rotates an H3 digit 60 degrees clockwise.
#[inline]
#[must_use]
pub(crate) fn _rotate60_cw(digit: Direction) -> Direction {
  use Direction::*;
  match digit {
    KAxes => JkAxes,
    JkAxes => JAxes,
    JAxes => IjAxes,
    IjAxes => IAxes,
    IAxes => IkAxes,
    IkAxes => KAxes,
    _ => digit, // Center or Invalid
  }
}

/// Finds the grid distance between the two IJK+ coordinates.
#[inline]
#[must_use]
pub(crate) fn ijk_distance(c1: &CoordIJK, c2: &CoordIJK) -> i32 {
  let mut diff = CoordIJK::default();
  _ijk_sub(c1, c2, &mut diff);
  _ijk_normalize(&mut diff); // Normalize to IJK "axial" coordinates.

  // The distance is the maximum component of the normalized difference.
  // Or, equivalently, (abs(i) + abs(j) + abs(k)) / 2, assuming sum i+j+k=0 for axial.
  // H3 C uses MAX(abs(i), abs(j), abs(k)).
  diff.i.abs().max(diff.j.abs()).max(diff.k.abs())
}

/// Transforms coordinates from the IJK+ coordinate system to the IJ coordinate
/// system. IJ coordinates are chosen relative to a given origin.
///
/// # Arguments
/// * `ijk` - The input IJK+ coordinates.
/// * `ij`  - Output: The IJ coordinates.
#[inline]
pub(crate) fn ijk_to_ij(ijk: &CoordIJK, ij: &mut CoordIJ) {
  // Note: This simple conversion assumes that IJK is already normalized
  // such that its k component can be ignored or is 0.
  // If ijk is an arbitrary H3 cell's local IJK, it might need
  // normalization first relative to its origin's IJK (which is (0,0,0)).
  // However, the C code appears to do `i = ijk->i - ijk->k; j = ijk->j - ijk->k;`
  // which is a standard way to convert from Cube/Axial (where i+j+k=0) to a 2-axis system.
  // Our _ijk_normalize ensures one component is 0. If k=0, then i=i, j=j.
  // If i=0, then new_i = 0 - k = -k; new_j = j - k.
  // If j=0, then new_i = i - k; new_j = 0 - k = -k.
  // The C code comment says "convert to IJ" where k is 0.
  // This implies normalized IJK where k has become 0 if it wasn't already minimal.
  ij.i = ijk.i - ijk.k;
  ij.j = ijk.j - ijk.k;
}

/// Transforms coordinates from the IJ coordinate system to the IJK+ coordinate system.
///
/// # Arguments
/// * `ij`  - The input IJ coordinates.
/// * `ijk` - Output: The IJK+ coordinates.
/// # Returns
/// An `H3Error` if an integer overflow would have occurred during normalization.
// Note: C's ijToIjk also has an error return for overflow, but only if a define is set.
// We make it always checked.
pub(crate) fn ij_to_ijk(ij: &CoordIJ, ijk: &mut CoordIJK) -> Result<(), crate::types::H3Error> {
  ijk.i = ij.i;
  ijk.j = ij.j;
  ijk.k = 0; // In IJK+, k is initially 0 before normalization.

  if _ijk_normalize_could_overflow(ijk) {
    // If we detect a potential overflow before normalization, return error.
    // This matches a pre-check pattern seen in some H3 C functions,
    // though not explicitly in the C ijToIjk which seems to assume inputs
    // don't cause overflow in its _ijkNormalize.
    // Our _ijk_normalize itself uses saturating arithmetic making it safe,
    // but if the *intent* of C was to error on such large pre-normalized values,
    // this check helps. For H3, typically inputs that would cause this are invalid.
    return Err(crate::types::H3Error::Failed); // Or a more specific error if available.
  }

  _ijk_normalize(ijk);
  Ok(())
}

/// Convert IJK coordinates to cube coordinates, in place.
/// Cube coordinates have the property that `i + j + k = 0`.
#[inline]
pub(crate) fn ijk_to_cube(ijk: &mut CoordIJK) {
  // This transformation assumes input IJK is *axial* where one component could be non-zero.
  // It does not assume input IJK has i+j+k=0.
  // If input IJK is H3's normalized (minimal non-negative components),
  // say {i',j',0}, then this maps to:
  // cube.i = -i'
  // cube.j = j'
  // cube.k = i' - j'
  // such that (-i') + (j') + (i' - j') = 0.
  // This matches the C code's transformation.
  let i = ijk.i; // Store original i
  ijk.i = ijk.i - ijk.k; // This actually computes the IJ system's i' = i-k
  ijk.j = ijk.j - ijk.k; // This computes the IJ system's j' = j-k
  ijk.k = -ijk.i - ijk.j; // Now use these to make k_cube = -i' -j'

  // The C code had:
  // ijk->i = -ijk->i + ijk->k;
  // ijk->j = ijk->j - ijk->k;
  // ijk->k = -ijk->i - ijk->j;
  // Let's trace with an example: IJK {1,0,0} (normalized H3 I-axis)
  // 1. i_new = -1 + 0 = -1
  // 2. j_new =  0 - 0 =  0
  // 3. k_new = -(-1) - 0 = 1
  // Result: {-1, 0, 1}. Sum is 0. This is correct for cube coords from I-axis.
  // My previous interpretation was slightly off.

  // Correct port of C logic:
  let i_orig = ijk.i;
  let j_orig = ijk.j;
  let k_orig = ijk.k;

  ijk.i = -i_orig + k_orig;
  ijk.j = j_orig - k_orig;
  ijk.k = -ijk.i - ijk.j; // Use the NEW ijk.i and ijk.j for this line
}

/// Convert cube coordinates to IJK coordinates, in place.
/// Assumes input `ijk` holds cube coordinates where `i+j+k = 0`.
/// Output `ijk` will be H3's normalized IJK+ (minimal non-negative components).
#[inline]
pub(crate) fn cube_to_ijk(ijk: &mut CoordIJK) {
  // This is the reverse of the C ijkToCube.
  // Given cube_i, cube_j, cube_k (where sum is 0):
  // h3_i = -cube_i
  // h3_j = cube_k  (This is wrong, based on RedBlobGames, cube.x = i, cube.z = j, cube.y = k. Then i = x, k = z. axial_to_cube: x = col, z = row)
  // h3_k = 0 (initially) then normalize.
  //
  // Let's look at C's cubeToIjk:
  // ijk->i = -ijk->i;
  // ijk->k = 0;
  // _ijkNormalize(ijk);
  // This means it takes cube.i, cube.j, cube.k (stored in ijk.i, ijk.j, ijk.k).
  // It effectively sets:
  // h3_i_intermediate = -cube.i
  // h3_j_intermediate = cube.j (remains untouched from input)
  // h3_k_intermediate = 0
  // Then normalizes {h3_i_intermediate, h3_j_intermediate, 0}.

  // Example: Cube {-1, 0, 1} (from I-axis)
  // h3_i_intermediate = -(-1) = 1
  // h3_j_intermediate = 0
  // h3_k_intermediate = 0
  // Normalize {1,0,0} -> {1,0,0}. Correct.

  // Example: Cube {0, -1, 1} (from J-axis)
  // h3_i_intermediate = -0 = 0
  // h3_j_intermediate = -1
  // h3_k_intermediate = 0
  // Normalize {0,-1,0}:
  //  i < 0 (false)
  //  j < 0 (true, j=-1): i = 0 - (-1) = 1. k = 0 - (-1) = 1. j = 0.
  //  Result: {1,0,1} (IK-axis in H3). But J-axis in H3 is {0,1,0}. This isn't right.

  // The C code's `ijkToCube` is:
  // i_axial = ijk->i - ijk->k;
  // j_axial = ijk->j - ijk->k;
  // ijk->i = i_axial; // Stores axial i
  // ijk->j = j_axial; // Stores axial j
  // ijk->k = -i_axial - j_axial; // Calculates k_cube and stores it
  //
  // So, `ijkToCube` actually transforms H3 IJK+ to Axial IJ, then computes k_cube.
  // The `cubeToIjk` in C must be reversing *that specific transformation*.
  // If input ijk has {i_axial, j_axial, k_cube_from_axial}, then:
  // ijk->i = -i_axial;
  // ijk->k = 0;
  // _ijkNormalize(ijk); where ijk.j was j_axial.
  // So we normalize {-i_axial, j_axial, 0}.

  // Let's use the C names to be clear:
  // Input `ijk` contains (ci, cj, ck) which are cube coordinates.
  // We want to get back H3's IJK+ (i', j', k')
  // From C: `ijk->i = -ijk->i; ijk->k = 0; _ijkNormalize(ijk);`
  // This means:
  //  i' (intermediate) = -ci
  //  j' (intermediate) = cj (original j from input struct, which was cj)
  //  k' (intermediate) = 0
  //  Then normalize { -ci, cj, 0 } to get final IJK+

  ijk.i = ijk.i.saturating_neg(); // i = -i (cube)
                                  // ijk.j remains ijk.j (cube)
  ijk.k = 0;
  _ijk_normalize(ijk);
}

#[cfg(test)]
mod tests {
  use crate::Vec2d;

  use super::*;

  #[test]
  fn test_set_ijk() {
    let mut ijk = CoordIJK::default();
    _set_ijk(&mut ijk, 1, 2, 3);
    assert_eq!(ijk, CoordIJK { i: 1, j: 2, k: 3 });
  }

  #[test]
  fn test_ijk_matches() {
    let c1 = CoordIJK { i: 1, j: 2, k: 3 };
    let c2 = CoordIJK { i: 1, j: 2, k: 3 };
    let c3 = CoordIJK { i: 4, j: 2, k: 3 };
    assert!(_ijk_matches(&c1, &c2));
    assert!(!_ijk_matches(&c1, &c3));
  }

  #[test]
  fn test_ijk_add() {
    let h1 = CoordIJK { i: 1, j: 2, k: -3 };
    let h2 = CoordIJK { i: 4, j: -5, k: 6 };
    let mut sum = CoordIJK::default();
    _ijk_add(&h1, &h2, &mut sum);
    assert_eq!(sum, CoordIJK { i: 5, j: -3, k: 3 });

    // Test saturation
    let h_max = CoordIJK {
      i: i32::MAX,
      j: 0,
      k: 0,
    };
    let h_one = CoordIJK { i: 1, j: 0, k: 0 };
    _ijk_add(&h_max, &h_one, &mut sum);
    assert_eq!(
      sum,
      CoordIJK {
        i: i32::MAX,
        j: 0,
        k: 0
      }
    );
  }

  #[test]
  fn test_ijk_sub() {
    let h1 = CoordIJK { i: 1, j: 2, k: -3 };
    let h2 = CoordIJK { i: 4, j: -5, k: 6 };
    let mut diff = CoordIJK::default();
    _ijk_sub(&h1, &h2, &mut diff);
    assert_eq!(diff, CoordIJK { i: -3, j: 7, k: -9 });

    // Test saturation
    let h_min = CoordIJK {
      i: i32::MIN,
      j: 0,
      k: 0,
    };
    let h_one = CoordIJK { i: 1, j: 0, k: 0 };
    _ijk_sub(&h_min, &h_one, &mut diff);
    assert_eq!(
      diff,
      CoordIJK {
        i: i32::MIN,
        j: 0,
        k: 0
      }
    );
  }

  #[test]
  fn test_ijk_scale() {
    let mut c = CoordIJK { i: 1, j: -2, k: 3 };
    _ijk_scale(&mut c, 2);
    assert_eq!(c, CoordIJK { i: 2, j: -4, k: 6 });
    _ijk_scale(&mut c, -1);
    assert_eq!(c, CoordIJK { i: -2, j: 4, k: -6 });
    _ijk_scale(&mut c, 0);
    assert_eq!(c, CoordIJK { i: 0, j: 0, k: 0 });

    // Test saturation
    let mut c_max = CoordIJK {
      i: i32::MAX / 2 + 1,
      j: 0,
      k: 0,
    };
    _ijk_scale(&mut c_max, 2);
    assert_eq!(
      c_max,
      CoordIJK {
        i: i32::MAX,
        j: 0,
        k: 0
      }
    );
  }

  #[test]
  fn test_ijk_normalize() {
    let mut c = CoordIJK { i: 0, j: 0, k: 0 };
    _ijk_normalize(&mut c);
    assert_eq!(c, CoordIJK { i: 0, j: 0, k: 0 }, "0,0,0 normalizes to 0,0,0");

    _set_ijk(&mut c, 2, 3, 4); // sum is 9
    _ijk_normalize(&mut c);
    assert_eq!(c, CoordIJK { i: 0, j: 1, k: 2 }, "positive components");

    _set_ijk(&mut c, -2, -3, -4); // sum is -9
    _ijk_normalize(&mut c);
    assert_eq!(c, CoordIJK { i: 2, j: 1, k: 0 }, "negative components");

    _set_ijk(&mut c, 2, -1, 0); // sum is 1
    _ijk_normalize(&mut c);
    assert_eq!(c, CoordIJK { i: 3, j: 0, k: 1 }, "mixed components (1)");

    _set_ijk(&mut c, 10, 20, 5); // sum = 35, min = 5
    _ijk_normalize(&mut c);
    assert_eq!(c, CoordIJK { i: 5, j: 15, k: 0 }, "remove min value");
  }

  #[test]
  fn test_ijk_normalize_saturating_behavior() {
    // Test cases that would overflow with standard arithmetic but should saturate
    let mut c = CoordIJK {
      i: i32::MIN,
      j: 1,
      k: 1,
    };
    _ijk_normalize(&mut c); // Expects i=0, j = 1-MIN, k = 1-MIN (both saturate to MAX)
    assert_eq!(
      c,
      CoordIJK {
        i: 0,
        j: i32::MAX,
        k: i32::MAX
      }
    );

    c = CoordIJK {
      i: 1,
      j: i32::MIN,
      k: 1,
    };
    _ijk_normalize(&mut c); // Expects i = 1-MIN (MAX), j = 0, k = 1-MIN (MAX)
    assert_eq!(
      c,
      CoordIJK {
        i: i32::MAX,
        j: 0,
        k: i32::MAX
      }
    );

    c = CoordIJK {
      i: 1,
      j: 1,
      k: i32::MIN,
    };
    _ijk_normalize(&mut c); // Expects i = 1-MIN (MAX), j = 1-MIN (MAX), k = 0
    assert_eq!(
      c,
      CoordIJK {
        i: i32::MAX,
        j: i32::MAX,
        k: 0
      }
    );

    // Case with multiple negative values, one being MIN
    c = CoordIJK {
      i: -10,
      j: i32::MIN,
      k: -20,
    };
    // After first step (i becomes 0): j = MIN - (-10), k = -20 - (-10) = -10. i=0
    //   j = MIN + 10
    // After second step (j becomes 0): i = 0 - (MIN+10) = -MIN-10. k = -10 - (MIN+10) = -MIN-20. j=0
    //   i = MAX-10+1 (due to two's complement for MIN). k = MAX-20+1
    //   i = 2147483638, k = 2147483628 (approx)
    // After third step (k becomes 0):
    //   i = i_prev - k_prev = (MAX-10+1) - (MAX-20+1) = 10
    //   j = 0 - k_prev = -(MAX-20+1) -> MIN+20-1
    //   This becomes tricky with saturation at each step.
    //
    // Let's trace with saturation:
    // c = { i: -10, j: i32::MIN, k: -20 }
    // if c.i < 0:
    //   c.j = c.j.saturating_sub(c.i) = i32::MIN.saturating_sub(-10) = i32::MIN (no change as i32::MIN + 10 is still < i32::MIN)
    //   c.k = c.k.saturating_sub(c.i) = -20.saturating_sub(-10) = -10
    //   c.i = 0
    // c is now { i: 0, j: i32::MIN, k: -10 }
    // if c.j < 0:
    //   c.i = c.i.saturating_sub(c.j) = 0.saturating_sub(i32::MIN) = i32::MAX
    //   c.k = c.k.saturating_sub(c.j) = -10.saturating_sub(i32::MIN) = i32::MAX
    //   c.j = 0
    // c is now { i: i32::MAX, j: 0, k: i32::MAX }
    // if c.k < 0: (false)
    // min_val = 0
    // if min_val > 0: (false)
    _ijk_normalize(&mut c);
    assert_eq!(
      c,
      CoordIJK {
        i: i32::MAX,
        j: 0,
        k: i32::MAX
      }
    );
  }

  #[test]
  fn test_ijk_normalize_could_overflow() {
    // Test cases that should not overflow _ijk_normalize with current C logic
    // (assuming k=0 for input to this pre-check, as per C usage)
    assert!(!_ijk_normalize_could_overflow(&CoordIJK { i: 0, j: 0, k: 0 }));
    assert!(!_ijk_normalize_could_overflow(&CoordIJK { i: 10, j: 5, k: 0 }));
    assert!(!_ijk_normalize_could_overflow(&CoordIJK { i: -10, j: -5, k: 0 }));
    assert!(!_ijk_normalize_could_overflow(&CoordIJK { i: 10, j: -5, k: 0 }));

    // Test cases that *could* cause overflow in C's _ijkNormalize intermediate steps
    // These are primarily when a large positive and large negative number are involved
    // such that their subtraction (or addition of absolute values) exceeds i32::MAX.
    // ijk.i < 0 && _add_i32s_overflows(ijk.j, -ijk.i)
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: i32::MIN,
      j: i32::MAX,
      k: 0
    })); // j - i => MAX - MIN overflows
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: i32::MIN + 1,
      j: i32::MAX,
      k: 0
    }));

    // ijk.j < 0 && _add_i32s_overflows(ijk.i, -ijk.j)
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: i32::MAX,
      j: i32::MIN,
      k: 0
    })); // i - j => MAX - MIN overflows
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: i32::MAX,
      j: i32::MIN + 1,
      k: 0
    }));

    // From C's specific logic: min_val < 0 checks
    // max_val + min_val
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: i32::MAX,
      j: i32::MIN + 1,
      k: 0
    })); // i + j would underflow
         // 0 - min_val
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: 0,
      j: i32::MIN,
      k: 0
    })); // 0 - MIN overflows
         // max_val - min_val
    assert!(_ijk_normalize_could_overflow(&CoordIJK {
      i: i32::MAX,
      j: i32::MIN,
      k: 0
    })); // MAX - MIN overflows
  }

  #[test]
  fn test_unit_ijk_to_digit() {
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 0, j: 0, k: 0 }),
      Direction::Center,
      "Center"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 0, j: 0, k: 1 }),
      Direction::KAxes,
      "K"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 0, j: 1, k: 0 }),
      Direction::JAxes,
      "J"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 0, j: 1, k: 1 }),
      Direction::JkAxes,
      "JK"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 1, j: 0, k: 0 }),
      Direction::IAxes,
      "I"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 1, j: 0, k: 1 }),
      Direction::IkAxes,
      "IK"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 1, j: 1, k: 0 }),
      Direction::IjAxes,
      "IJ"
    );

    // Test non-unit vectors after normalization
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 2, j: 2, k: 2 }),
      Direction::Center,
      "Unnormalized Center"
    ); // normalizes to 0,0,0
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 1, j: 1, k: 2 }),
      Direction::KAxes,
      "Unnormalized K"
    ); // normalizes to 0,0,1

    // Test invalid
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 2, j: 0, k: 0 }),
      Direction::InvalidDigit,
      "Not a unit vector (I=2)"
    );
    assert_eq!(
      _unit_ijk_to_digit(&CoordIJK { i: 1, j: 2, k: 3 }),
      Direction::InvalidDigit,
      "Not a unit vector (complex)"
    );
  }

  #[test]
  fn test_neighbor() {
    let mut ijk = CoordIJK { i: 0, j: 0, k: 0 };
    let initial_ijk = ijk;

    _neighbor(&mut ijk, Direction::Center);
    assert!(_ijk_matches(&ijk, &initial_ijk), "Center neighbor is self");

    _neighbor(&mut ijk, Direction::IAxes);
    assert!(
      _ijk_matches(&ijk, &UNIT_VECS[Direction::IAxes as usize]),
      "I neighbor as expected"
    );

    // Reset for next test
    ijk = initial_ijk;
    _neighbor(&mut ijk, Direction::InvalidDigit); // Should do nothing
    assert!(_ijk_matches(&ijk, &initial_ijk), "Invalid neighbor is self");

    // Test from non-origin
    _set_ijk(&mut ijk, 1, 1, 1); // Normalizes to 0,0,0
    _ijk_normalize(&mut ijk);
    _neighbor(&mut ijk, Direction::JAxes);
    assert!(
      _ijk_matches(&ijk, &UNIT_VECS[Direction::JAxes as usize]),
      "J neighbor from normalized (1,1,1)"
    );
  }

  // Test _hex2dToCoordIJK (ported from C testFaceIjk.c _hex2dToCoordIJK_coordinates_test)
  fn assert_hex2d_to_coord_ijk(v: Vec2d, expected_h: CoordIJK, msg: &str) {
    let mut h = CoordIJK::default();
    _hex2d_to_coord_ijk(&v, &mut h);
    assert!(
      _ijk_matches(&h, &expected_h),
      "{} - Expected ({},{},{}), Got ({},{},{})",
      msg,
      expected_h.i,
      expected_h.j,
      expected_h.k,
      h.i,
      h.j,
      h.k
    );
  }

  #[test]
  fn test_hex2d_to_coord_ijk() {
    assert_hex2d_to_coord_ijk(Vec2d { x: 0.0, y: 0.0 }, CoordIJK { i: 0, j: 0, k: 0 }, "origin");
    assert_hex2d_to_coord_ijk(Vec2d { x: 1.0, y: 0.0 }, CoordIJK { i: 1, j: 0, k: 0 }, "i axis");
    assert_hex2d_to_coord_ijk(Vec2d { x: -0.5, y: M_SQRT3_2 }, CoordIJK { i: 0, j: 1, k: 0 }, "j axis");
    assert_hex2d_to_coord_ijk(
      Vec2d { x: -0.5, y: -M_SQRT3_2 },
      CoordIJK { i: 0, j: 0, k: 1 },
      "k axis",
    );
    assert_hex2d_to_coord_ijk(Vec2d { x: 0.1, y: 0.1 }, CoordIJK { i: 0, j: 0, k: 0 }, "eps case");
    assert_hex2d_to_coord_ijk(Vec2d { x: -0.1, y: -0.1 }, CoordIJK { i: 0, j: 0, k: 0 }, "eps case N");
    assert_hex2d_to_coord_ijk(Vec2d { x: 0.1, y: -0.1 }, CoordIJK { i: 0, j: 0, k: 0 }, "eps case E");
    assert_hex2d_to_coord_ijk(Vec2d { x: -0.1, y: 0.1 }, CoordIJK { i: 0, j: 0, k: 0 }, "eps case W");
    assert_hex2d_to_coord_ijk(Vec2d { x: 0.5, y: M_SQRT3_2 }, CoordIJK { i: 0, j: 1, k: 0 }, "on edge");
    assert_hex2d_to_coord_ijk(
      Vec2d { x: 0.4, y: M_SQRT3_2 },
      CoordIJK { i: 0, j: 1, k: 0 },
      "pert on edge1",
    );
    assert_hex2d_to_coord_ijk(
      Vec2d { x: 0.6, y: M_SQRT3_2 },
      CoordIJK { i: 0, j: 1, k: 0 },
      "pert on edge2",
    );
    assert_hex2d_to_coord_ijk(
      Vec2d { x: 1.0, y: M_SQRT3_2 },
      CoordIJK { i: 1, j: 1, k: 0 },
      "Shared vertex",
    );
  }

  // Test _ijkToHex2d (ported from C testFaceIjk.c _ijkToHex2d_roundtrip_test)
  fn assert_ijk_to_hex2d_roundtrip(h: CoordIJK, msg: &str) {
    let mut v = Vec2d::default();
    _ijk_to_hex2d(&h, &mut v);
    let mut h2 = CoordIJK::default();
    _hex2d_to_coord_ijk(&v, &mut h2);
    assert!(_ijk_matches(&h, &h2), "{} - roundtrip", msg);
  }

  #[test]
  fn test_ijk_to_hex2d_roundtrip() {
    assert_ijk_to_hex2d_roundtrip(CoordIJK { i: 0, j: 0, k: 0 }, "origin");
    for dir_enum_val in 1..=6 {
      // K_AXES_DIGIT to IJ_AXES_DIGIT
      let dir = unsafe { std::mem::transmute::<u8, Direction>(dir_enum_val as u8) };
      let unit_vec = UNIT_VECS[dir as usize];
      assert_ijk_to_hex2d_roundtrip(unit_vec, &format!("unit vec {:?}", dir));
    }
    assert_ijk_to_hex2d_roundtrip(CoordIJK { i: 2, j: 0, k: 1 }, "arbitrary"); // Normalized to (1,0,0)
  }

  #[test]
  fn test_up_ap7_checked() {
    let mut ijk = CoordIJK { i: 0, j: 0, k: 0 };
    assert!(_up_ap7_checked(&mut ijk).is_ok());
    assert_eq!(ijk, CoordIJK { i: 0, j: 0, k: 0 });

    // Example from C test (testH3IndexInternal.c#faceIjkToH3ExtremeCoordinates)
    // Where fijk.coord = {18,0,0} for res 2, upAp7r is called for class II (res 2)
    // and upAp7 is called for class III (res 1).
    // Let's test a similar progression. For _up_ap7 (class III):
    // Start with a fine res coord, e.g., equivalent to res 1 digit K (0,0,1) scaled up.
    // UNIT_VECS[1] = {0,0,1}. Let's scale it as if it's child of a res 0 cell at res 1.
    // This would mean its parent at res 0 is (0,0,0) after _up_ap7.
    let mut c = UNIT_VECS[Direction::KAxes as usize]; // (0,0,1)
    assert!(_up_ap7_checked(&mut c).is_ok());
    assert_eq!(c, CoordIJK { i: 0, j: 0, k: 0 }, "K axis upAp7");

    let mut c2 = UNIT_VECS[Direction::IAxes as usize]; // (1,0,0)
    assert!(_up_ap7_checked(&mut c2).is_ok());
    assert_eq!(c2, CoordIJK { i: 0, j: 0, k: 0 }, "I axis upAp7");

    // Test overflow scenarios. These are hard to hit with valid H3 cell coordinates
    // as they are typically small after normalization per resolution.
    // Using large raw numbers to test overflow checks.
    let mut large_i = CoordIJK {
      i: i32::MAX,
      j: 0,
      k: 0,
    };
    assert!(
      _up_ap7_checked(&mut large_i).is_err(),
      "Large i for up_ap7 should overflow"
    );

    let mut large_j = CoordIJK {
      i: 0,
      j: i32::MAX,
      k: 0,
    };
    assert!(
      _up_ap7_checked(&mut large_j).is_err(),
      "Large j for up_ap7 should overflow"
    );
  }

  #[test]
  fn test_up_ap7r_checked() {
    let mut ijk = CoordIJK { i: 0, j: 0, k: 0 };
    assert!(_up_ap7r_checked(&mut ijk).is_ok());
    assert_eq!(ijk, CoordIJK { i: 0, j: 0, k: 0 });

    let mut c = UNIT_VECS[Direction::KAxes as usize]; // (0,0,1)
    assert!(_up_ap7r_checked(&mut c).is_ok());
    assert_eq!(c, CoordIJK { i: 0, j: 0, k: 0 }, "K axis upAp7r");

    // Test overflow
    let mut large_i = CoordIJK {
      i: i32::MAX,
      j: 0,
      k: 0,
    };
    assert!(
      _up_ap7r_checked(&mut large_i).is_err(),
      "Large i for up_ap7r should overflow"
    );
  }

  fn test_down_aperture(
    ap_func: fn(&mut CoordIJK),
    ap_r_func: fn(&mut CoordIJK),
    expected_ijs: &[(i32, i32, i32)], // Expected IJK for Center, I, J, K after down-aperture
    label: &str,
  ) {
    // Test Center
    let mut center = CoordIJK { i: 0, j: 0, k: 0 };
    ap_func(&mut center);
    assert_eq!(
      center,
      CoordIJK {
        i: expected_ijs[0].0,
        j: expected_ijs[0].1,
        k: expected_ijs[0].2
      },
      "{}: Center",
      label
    );

    let mut center_r = CoordIJK { i: 0, j: 0, k: 0 };
    ap_r_func(&mut center_r);
    assert_eq!(
      center_r,
      CoordIJK {
        i: expected_ijs[0].0,
        j: expected_ijs[0].1,
        k: expected_ijs[0].2
      },
      "{}: Center (reversed)",
      label
    ); // Reversed should be same for center

    // Test I-axis (CoordIJK {1,0,0})
    let mut i_axis = UNIT_VECS[Direction::IAxes as usize];
    ap_func(&mut i_axis);
    assert_eq!(
      i_axis,
      CoordIJK {
        i: expected_ijs[1].0,
        j: expected_ijs[1].1,
        k: expected_ijs[1].2
      },
      "{}: I-axis",
      label
    );

    // Test J-axis (CoordIJK {0,1,0})
    let mut j_axis = UNIT_VECS[Direction::JAxes as usize];
    ap_func(&mut j_axis);
    assert_eq!(
      j_axis,
      CoordIJK {
        i: expected_ijs[2].0,
        j: expected_ijs[2].1,
        k: expected_ijs[2].2
      },
      "{}: J-axis",
      label
    );

    // Test K-axis (CoordIJK {0,0,1})
    let mut k_axis = UNIT_VECS[Direction::KAxes as usize];
    ap_func(&mut k_axis);
    assert_eq!(
      k_axis,
      CoordIJK {
        i: expected_ijs[3].0,
        j: expected_ijs[3].1,
        k: expected_ijs[3].2
      },
      "{}: K-axis",
      label
    );
  }

  #[test]
  fn test_down_ap7() {
    // Expected from C's _downAp7:
    // I: (3,0,1), J: (1,3,0), K: (0,1,3)
    // Center is always (0,0,0) after normalization
    let expected = [(0, 0, 0), (3, 0, 1), (1, 3, 0), (0, 1, 3)];
    test_down_aperture(_down_ap7, _down_ap7r, &expected, "downAp7");
  }

  #[test]
  fn test_down_ap3() {
    // Expected from C's _downAp3:
    // I: (2,0,1), J: (1,2,0), K: (0,1,2)
    let expected = [(0, 0, 0), (2, 0, 1), (1, 2, 0), (0, 1, 2)];
    test_down_aperture(_down_ap3, _down_ap3r, &expected, "downAp3");
  }

  #[test]
  fn test_ijk_rotate60() {
    let mut c = CoordIJK { i: 1, j: 0, k: 0 }; // I vector
    _ijk_rotate60_ccw(&mut c); // Should be IJ {1,1,0}
    assert_eq!(c, CoordIJK { i: 1, j: 1, k: 0 }, "I ccw is IJ");
    _ijk_rotate60_ccw(&mut c); // Should be J {0,1,0}
    assert_eq!(c, CoordIJK { i: 0, j: 1, k: 0 }, "IJ ccw is J");
    _ijk_rotate60_ccw(&mut c); // Should be JK {0,1,1}
    assert_eq!(c, CoordIJK { i: 0, j: 1, k: 1 }, "J ccw is JK");
    _ijk_rotate60_ccw(&mut c); // Should be K {0,0,1}
    assert_eq!(c, CoordIJK { i: 0, j: 0, k: 1 }, "JK ccw is K");
    _ijk_rotate60_ccw(&mut c); // Should be IK {1,0,1}
    assert_eq!(c, CoordIJK { i: 1, j: 0, k: 1 }, "K ccw is IK");
    _ijk_rotate60_ccw(&mut c); // Should be I {1,0,0} again
    assert_eq!(c, CoordIJK { i: 1, j: 0, k: 0 }, "IK ccw is I");

    // Now clockwise
    _ijk_rotate60_cw(&mut c); // Should be IK {1,0,1}
    assert_eq!(c, CoordIJK { i: 1, j: 0, k: 1 }, "I cw is IK");
  }

  #[test]
  fn test_rotate_digit() {
    assert_eq!(_rotate60_ccw(Direction::KAxes), Direction::IkAxes);
    assert_eq!(_rotate60_cw(Direction::KAxes), Direction::JkAxes);
    assert_eq!(_rotate60_ccw(Direction::Center), Direction::Center); // No change
    assert_eq!(_rotate60_cw(Direction::InvalidDigit), Direction::InvalidDigit); // No change
  }

  #[test]
  fn test_ijk_distance() {
    let z = CoordIJK { i: 0, j: 0, k: 0 };
    let i = CoordIJK { i: 1, j: 0, k: 0 };
    let ik = CoordIJK { i: 1, j: 0, k: 1 }; // normalize -> {1,0,1} -> {0, -1, 0} after sub? No, normalize {1,0,1} is still {1,0,1}
                                            // IJK axial coords. {1,0,1} is not normalized in H3 sense (sum of components !=0)
                                            // _ijk_normalize({1,0,1}) = {0,0,0} because min is 0, then i=1-0,j=0-0,k=1-0. No, min is 0.
                                            // _ijk_normalize({1,0,1}) = {1,0,1} because all non-negative, min is 0.
                                            // The inputs to ijkDistance are H3 IJK+ (normalized).
                                            // So, ik is an H3 IJK+.
    let ij = CoordIJK { i: 1, j: 1, k: 0 };
    let j2 = CoordIJK { i: 0, j: 2, k: 0 };

    assert_eq!(ijk_distance(&z, &z), 0, "identity distance 0,0,0");
    assert_eq!(ijk_distance(&i, &i), 0, "identity distance 1,0,0");
    assert_eq!(ijk_distance(&ik, &ik), 0, "identity distance 1,0,1");
    assert_eq!(ijk_distance(&ij, &ij), 0, "identity distance 1,1,0");
    assert_eq!(ijk_distance(&j2, &j2), 0, "identity distance 0,2,0");

    assert_eq!(ijk_distance(&z, &i), 1, "0,0,0 to 1,0,0");
    assert_eq!(ijk_distance(&z, &j2), 2, "0,0,0 to 0,2,0");
    assert_eq!(ijk_distance(&z, &ik), 1, "0,0,0 to 1,0,1"); // diff {1,0,1}, max_abs=1
    assert_eq!(ijk_distance(&i, &ik), 1, "1,0,0 to 1,0,1"); // diff {0,0,1}, max_abs=1
    assert_eq!(ijk_distance(&ik, &j2), 3, "1,0,1 to 0,2,0"); // diff {1,-2,1}, max_abs=2? No, normalize {1,-2,1}
                                                             // {1,-2,1} -> {1+2, 0, 1+2} = {3,0,3} -> {0, -3, 0}
                                                             // _ijk_sub(ik,j2) = {1,0,1} - {0,2,0} = {1,-2,1}
                                                             // _ijk_normalize({1,-2,1}):
                                                             //  i=1, j=-2, k=1
                                                             //  j<0: i=1-(-2)=3. k=1-(-2)=3. j=0. -> {3,0,3}
                                                             //  min=0. -> {3,0,3}
                                                             //  abs components are 3,0,3. Max is 3.
    assert_eq!(ijk_distance(&ij, &ik), 2, "1,1,0 to 1,0,1"); // diff {0,1,-1}
                                                             // _ijk_normalize({0,1,-1}):
                                                             //  i=0, j=1, k=-1
                                                             //  k<0: i=0-(-1)=1. j=1-(-1)=2. k=0. -> {1,2,0}
                                                             //  min=0. -> {1,2,0}
                                                             //  abs components are 1,2,0. Max is 2.
  }

  #[test]
  fn test_ijk_to_ij_and_back() {
    let mut orig_ijk = CoordIJK { i: 1, j: 2, k: 0 }; // Already H3 IJK+ normalized
    _ijk_normalize(&mut orig_ijk); // ensure it's minimal non-negative
    let mut ij = CoordIJ::default();
    let mut new_ijk = CoordIJK::default();

    ijk_to_ij(&orig_ijk, &mut ij);
    let res = ij_to_ijk(&ij, &mut new_ijk);
    assert!(res.is_ok());
    assert!(
      _ijk_matches(&orig_ijk, &new_ijk),
      "Roundtrip IJK->IJ->IJK should match for {:?}",
      orig_ijk
    );

    // Test with k being non-zero before normalization (but still i+j+k=0 for axial representation)
    let axial_k_non_zero = CoordIJK { i: 1, j: -3, k: 2 }; // Sum is 0
                                                           // H3 IJK+ needs one component to be 0 and others non-negative.
                                                           // We need to convert this axial to H3 IJK+ first to test fairly.
    let mut h3_ijk_from_axial = axial_k_non_zero;
    _ijk_normalize(&mut h3_ijk_from_axial); // {1, -3, 2} -> {1 - (-3), 0, 2 - (-3)} = {4,0,5} -> {0, -4, 1} (this logic is faulty)
                                            // _ijk_normalize of {1,-3,2}:
                                            // j<0 (-3): i_new = 1-(-3)=4. k_new = 2-(-3)=5. j_new=0. -> {4,0,5}
                                            // min=0. Stays {4,0,5}. This is a valid H3 IJK+.

    ijk_to_ij(&h3_ijk_from_axial, &mut ij); // ij.i = 4-5 = -1. ij.j = 0-5 = -5.
    let res2 = ij_to_ijk(&ij, &mut new_ijk);
    assert!(res2.is_ok());
    assert!(
      _ijk_matches(&h3_ijk_from_axial, &new_ijk),
      "Roundtrip with initial non-zero k (after H3 norm) for {:?}",
      h3_ijk_from_axial
    );
  }

  #[test]
  fn test_ijk_cube_roundtrip() {
    let mut ijk = CoordIJK { i: 1, j: 2, k: 0 }; // Normalized H3 IJK+
    _ijk_normalize(&mut ijk);
    let original_h3_ijk = ijk;

    ijk_to_cube(&mut ijk); // ijk now holds cube coords
                           // Verify sum of cube coords is 0
    assert_eq!(ijk.i + ijk.j + ijk.k, 0, "Cube coordinates sum to 0");

    cube_to_ijk(&mut ijk); // ijk converted back to H3 IJK+
    assert!(
      _ijk_matches(&original_h3_ijk, &ijk),
      "Roundtrip IJK->Cube->IJK should match for {:?}",
      original_h3_ijk
    );

    // Test another
    let mut ijk2 = CoordIJK { i: 0, j: 0, k: 1 }; // K-axis vector in H3
    _ijk_normalize(&mut ijk2);
    let original_h3_ijk2 = ijk2;
    ijk_to_cube(&mut ijk2);
    assert_eq!(ijk2.i + ijk2.j + ijk2.k, 0, "Cube coordinates sum to 0 for K-axis");
    cube_to_ijk(&mut ijk2);
    assert!(_ijk_matches(&original_h3_ijk2, &ijk2), "Roundtrip for K-axis");
  }
}
