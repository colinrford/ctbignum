//
// This file is part of
//
// CTBignum
//
// C++ Library for Compile-Time and Run-Time Multi-Precision and Modular Arithmetic
//
//
// This file is distributed under the Apache License, Version 2.0. See the LICENSE
// file for details.

export module lam.ctbignum:decimal_literals;

import std;

import :bigint;
import :addition;
import :mult;

namespace lam::cbn
{

/**
 * Runtime (and constexpr) decimal string to big_int conversion.
 *
 * Companion to chars_to_big_int which only works at compile time
 * via the UDL integer_sequence<char,...> path.
 *
 * Accepted input:
 *   - Decimal digit characters '0'..'9', optionally separated by '_'.
 *   - Leading zeros are permitted.
 *   - Empty strings and strings with no digit characters return std::nullopt.
 *   - Any non-digit, non-underscore character returns std::nullopt.
 *
 * Overflow detection: N+1 limbs are used internally; if the accumulated value
 * exceeds N limbs the extra limb is non-zero and std::nullopt is returned.
 */
export template<std::size_t N, typename T = std::uint64_t>
constexpr auto big_int_from_string(std::string_view s) -> std::optional<big_int<N, T>>
{
  constexpr std::size_t max_digits = (N * std::numeric_limits<T>::digits * 301) / 1000 + 2;

  std::size_t digit_count = 0;
  for (char c : s)
  {
    if (c >= '0' && c <= '9')
      ++digit_count;
    else if (c != '_')
      return std::nullopt;
  }
  if (digit_count == 0)
    return std::nullopt;
  if (digit_count > max_digits)
    return std::nullopt;

  big_int<N + 1, T> num{};
  big_int<N + 1, T> power_of_ten{};
  power_of_ten[0] = T{1};

  for (std::size_t idx = s.size(); idx > 0; --idx)
  {
    if (s[idx - 1] == '_')
      continue;
    const T digit = static_cast<T>(s[idx - 1] - '0');
    num = add_ignore_carry(num, partial_mul<N + 1>(big_int<1, T>{digit}, power_of_ten));
    power_of_ten = partial_mul<N + 1>(big_int<1, T>{T{10}}, power_of_ten);
  }

  if (num[N] != T{0})
    return std::nullopt;

  big_int<N, T> result{};
  for (std::size_t i = 0; i < N; ++i)
    result[i] = num[i];
  return result;
}

} // namespace lam::cbn
