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

export module lam.ctbignum:literals;

import std;

import :bigint;
import :slicing;
import :utility;
import :addition;
import :mult;
import :relational;

namespace lam::cbn
{
namespace detail
{

export template <typename T = std::uint64_t, std::size_t L = 0, char... Chars> //, std::size_t... Is>
constexpr auto chars_to_big_int(std::integer_sequence<char, Chars...>)
{
  // might return a 'non-tight' representation, meaning that there could be
  // leading zero-limbs
  constexpr std::size_t len = sizeof...(Chars);
  constexpr std::size_t N = std::max(L, 1 + (10 * len) / (3 * std::numeric_limits<T>::digits));
  std::array<char, len> digits{Chars...};
  big_int<N, T> num{0};
  big_int<N, T> power_of_ten{1};

  for (int i = len - 1; i >= 0; --i)
  {
    num = add_ignore_carry(num, partial_mul<N>(big_int<1, T>{static_cast<T>(digits[i]) - 48}, power_of_ten));
    power_of_ten = partial_mul<N>(big_int<1, T>{static_cast<T>(10)}, power_of_ten);
  }
  return num;
}

export template <typename T = std::uint64_t, char... Chars, std::size_t... Is>
constexpr auto chars_to_integer_seq(std::integer_sequence<char, Chars...>, std::index_sequence<Is...>)
{
  constexpr auto num = detail::chars_to_big_int<T, sizeof...(Chars)>(std::integer_sequence<char, Chars...>{});
  return std::integer_sequence<T, num[Is]...>{};
}

} // namespace detail

namespace literals
{
export template <typename T, char... Chars> constexpr auto generic_limb_literal()
{
  constexpr std::size_t len = sizeof...(Chars);
  constexpr std::size_t N = 1 + (10 * len) / (3 * std::numeric_limits<T>::digits);
  auto num = detail::chars_to_integer_seq<T>(std::integer_sequence<char, Chars...>{}, std::make_index_sequence<N>{});
  constexpr auto L = detail::tight_length(num) + (to_big_int(num) == big_int<1, T>{});
  return detail::take_first(num, std::make_index_sequence<L>{});
}

export template <char... Chars> constexpr auto operator""_Z() // for backwards compatibility
{
  return generic_limb_literal<std::uint64_t, Chars...>();
}

export template <char... Chars> constexpr auto operator""_Z64()
{
  return generic_limb_literal<std::uint64_t, Chars...>();
}

export template <char... Chars> constexpr auto operator""_Z32()
{
  return generic_limb_literal<std::uint32_t, Chars...>();
}

} // namespace literals
} // namespace lam::cbn
