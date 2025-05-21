// src/h3_index/string_conv.rs

use crate::types::{H3Error, H3Index, H3_NULL};

/// Converts a string representation of an H3 index into an `H3Index`.
///
/// # Arguments
/// * `s` - The string representation of an H3 index (hexadecimal).
///
/// # Returns
/// `Ok(H3Index)` on success, or an `H3Error` if parsing fails.
pub fn string_to_h3(s: &str) -> Result<H3Index, H3Error> {
  if s.is_empty() {
    return Err(H3Error::Failed); // Match C behavior of failing on empty string
  }
  match u64::from_str_radix(s, 16) {
    Ok(val) => Ok(H3Index(val)),
    Err(_) => Err(H3Error::Failed), // C sscanf returns 0 or EOF for bad parse
  }
}

/// Converts an `H3Index` into its string representation.
///
/// # Arguments
/// * `h` - The `H3Index` to convert.
/// * `buffer` - A mutable byte slice to write the hexadecimal string into.
///              The string will be null-terminated if space allows.
///
/// # Returns
/// `Ok(())` on success, or `H3Error::MemoryBounds` if the buffer is too small.
/// The H3 C library requires a buffer of at least 17 bytes (16 hex chars + null).
pub fn h3_to_string(h: H3Index, buffer: &mut [u8]) -> Result<(), H3Error> {
  const MIN_BUF_SIZE_FOR_CONTENT: usize = 16; // Max 16 hex chars for u64
  const MIN_BUF_SIZE_WITH_NULL: usize = MIN_BUF_SIZE_FOR_CONTENT + 1;

  if buffer.len() < MIN_BUF_SIZE_WITH_NULL {
    return Err(H3Error::MemoryBounds);
  }

  let s = format!("{:x}", h.0);
  let bytes = s.as_bytes();

  // We need to check if the formatted string + null terminator will fit.
  // bytes.len() is the length of the hex string.
  // We need bytes.len() + 1 for the null terminator.
  if bytes.len() + 1 > buffer.len() {
    // This case should ideally be caught by the MIN_BUF_SIZE_WITH_NULL check,
    // but it's good to be explicit if format! somehow produced a longer string
    // than expected (though for {:x} on u64, it's max 16).
    return Err(H3Error::MemoryBounds);
  }

  // Use the original `buffer` mutable reference directly.
  // No need to rebind to `writer` if we are just going to use `buffer`.
  buffer[..bytes.len()].copy_from_slice(bytes);
  buffer[bytes.len()] = 0; // Null terminate

  Ok(())
}

/// Converts an `H3Index` into a `String`.
/// This is a Rust-idiomatic alternative to `h3_to_string`.
#[must_use]
pub fn h3_to_string_alloc(h: H3Index) -> String {
  format!("{:x}", h.0)
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn test_string_to_h3() {
    assert_eq!(string_to_h3("8928308280fffff"), Ok(H3Index(0x8928308280fffff)));
    assert_eq!(string_to_h3("0"), Ok(H3Index(0)));
    assert_eq!(string_to_h3("ffffffffffffffff"), Ok(H3Index(0xffffffffffffffff)));

    assert_eq!(string_to_h3(""), Err(H3Error::Failed));
    assert_eq!(string_to_h3("invalid"), Err(H3Error::Failed));
    assert_eq!(string_to_h3("123INVALID"), Err(H3Error::Failed));
    // Test string too long to be a u64 hex (more than 16 chars)
    assert_eq!(string_to_h3("10000000000000000"), Err(H3Error::Failed));
  }

  #[test]
  fn test_h3_to_string() {
    let mut buffer = [0u8; 17];

    assert_eq!(h3_to_string(H3Index(0x8928308280fffff), &mut buffer), Ok(()));
    assert_eq!(&buffer[0..15], b"8928308280fffff");
    assert_eq!(buffer[15], 0); // Null terminated

    assert_eq!(h3_to_string(H3Index(0), &mut buffer), Ok(()));
    assert_eq!(&buffer[0..1], b"0");
    assert_eq!(buffer[1], 0);

    assert_eq!(h3_to_string(H3Index(0xffffffffffffffff), &mut buffer), Ok(()));
    assert_eq!(&buffer[0..16], b"ffffffffffffffff");
    assert_eq!(buffer[16], 0);

    // Test buffer too small
    let mut small_buffer = [0u8; 16]; // exactly 16, no room for null if string is 16 chars
    assert_eq!(
      h3_to_string(H3Index(0xffffffffffffffff), &mut small_buffer),
      Err(H3Error::MemoryBounds)
    );

    let mut too_small_buffer = [0u8; 5];
    assert_eq!(
      h3_to_string(H3Index(0x8928308280fffff), &mut too_small_buffer),
      Err(H3Error::MemoryBounds)
    );
  }

  #[test]
  fn test_h3_to_string_alloc() {
    assert_eq!(h3_to_string_alloc(H3Index(0x8928308280fffff)), "8928308280fffff");
    assert_eq!(h3_to_string_alloc(H3Index(0)), "0");
    assert_eq!(h3_to_string_alloc(H3Index(0xffffffffffffffff)), "ffffffffffffffff");
  }
}
