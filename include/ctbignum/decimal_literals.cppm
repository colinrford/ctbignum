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
 * Parse a big_int from a runtime string. Supports '_' digit separators.
 * Uses N+1 limbs internally to detect overflow.
 */
export template<std::size_t N, typename T = std::uint64_t>
constexpr auto big_int_from_string(std::string_view s) -> big_int<N, T>
{
  big_int<N + 1, T> num{0};
  big_int<N + 1, T> power_of_ten{1};

  for (std::size_t idx = s.size(); idx > 0; --idx)
  {
    char c = s[idx - 1];
    if (c == '_') continue;

    num = add_ignore_carry(num, partial_mul<N + 1>(big_int<1, T>{static_cast<T>(c - '0')}, power_of_ten));
    power_of_ten = partial_mul<N + 1>(big_int<1, T>{static_cast<T>(10)}, power_of_ten);
  }

  if (num[N] != 0)
    throw std::overflow_error("big_int_from_string: value overflows N limbs");

  big_int<N, T> result{};
  for (std::size_t i = 0; i < N; ++i)
    result[i] = num[i];
  return result;
}

} // namespace lam::cbn
