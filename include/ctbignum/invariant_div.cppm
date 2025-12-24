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

export module lam.ctbignum:invariant_div;

import std;

import :bigint;
import :slicing;
import :division;
import :mult;
import :addition;
import :utility;
import :bitshift;

namespace lam::cbn
{
namespace detail
{

export template <std::size_t N, typename T = std::uint64_t, T... Divisor, std::size_t... Is>
constexpr auto precompute_m_prime_nontight(std::integer_sequence<T, Divisor...>, std::index_sequence<Is...>)
{
  constexpr auto D = sizeof...(Divisor);
  constexpr big_int<D, T> d{Divisor...};
  constexpr auto ell = bit_length(d - big_int<1, T>{1}); // TODO: should this indeed be d-1 ???
  constexpr auto w = std::numeric_limits<T>::digits;
  constexpr auto limb_shifts = ell / w;
  constexpr auto bit_shifts = ell % w;
  constexpr auto pow2ell = place_at<std::max(D, limb_shifts + 1), T>(static_cast<T>(1) << bit_shifts, limb_shifts);
  constexpr auto pow2N = unary_encoding<N, N + 1, T>();
  constexpr auto divrem = div(mul(pow2N, subtract(pow2ell, d)), d);
  constexpr auto mp = to_length<N>(add(divrem.quotient, big_int<1, T>{static_cast<T>(1)}));
  return std::integer_sequence<T, mp[Is]...>{};
}

export template <std::size_t N, typename T = std::uint64_t, T... Divisor>
constexpr auto precompute_m_prime(std::integer_sequence<T, Divisor...>)
{
  auto m = precompute_m_prime_nontight<N, T>(std::integer_sequence<T, Divisor...>{}, std::make_index_sequence<N>{});
  return take_first(m, std::make_index_sequence<tight_length(m)>{});
}

} // namespace detail

export template <typename T, std::size_t N, T... Divisor>
constexpr auto quotient(big_int<N, T> n, std::integer_sequence<T, Divisor...>)
{
  // Integer division with compile-time divisor
  // as described in "Division by Invariant Integers using Multiplication",
  // by Granlund and Montgomery, 1994
  // https://gmplib.org/~tege/divcnst-pldi94.pdf
  //
  // inputs:
  //  n           divident
  //  Divisor...  compile-time divisor
  //

  using detail::skip;
  using detail::to_length;

  constexpr big_int<sizeof...(Divisor), T> d{Divisor...};
  if constexpr (sizeof...(Divisor) > N)
    return big_int<1, T>{static_cast<T>(0)};
  else if constexpr (d == big_int<1, T>{static_cast<T>(1)})
    return n;
  else
  {
    // Compile-time precomputation of m_prime
    constexpr auto ell = detail::bit_length(d - big_int<1, T>{1});
    constexpr auto w = std::numeric_limits<T>::digits;
    constexpr auto m_prime = to_big_int(detail::precompute_m_prime<N>(std::integer_sequence<T, Divisor...>{}));
    // end of pre-computation

    // Perform the division
    auto t1 = skip<N>(mul(m_prime, n));
    auto q = shift_right(skip<(ell - 1) / w>(add(t1, shift_right(subtract_ignore_carry(n, to_length<N>(t1)), 1))),
                         (ell - 1) % w); // n >= t1
    return to_length<N>(q);
  }
}

export template <typename T, std::size_t N, T... Modulus>
constexpr auto mod(big_int<N, T> n, std::integer_sequence<T, Modulus...>)
// Constant-time modulo operation with a fixed modulus
{
  auto d = quotient(n, std::integer_sequence<T, Modulus...>{});
  constexpr auto M = sizeof...(Modulus);
  return detail::to_length<M>(subtract_ignore_carry(n, partial_mul<N>(big_int<M, T>{Modulus...}, d)));
}

export template <typename T, std::size_t N, T... Modulus>
constexpr DivisionResult<big_int<N, T>, big_int<sizeof...(Modulus), T>> div(big_int<N, T> n,
                                                                            std::integer_sequence<T, Modulus...>)
{
  auto quot = quotient(n, std::integer_sequence<T, Modulus...>{});
  constexpr auto M = sizeof...(Modulus);
  auto rem = detail::to_length<M>(subtract_ignore_carry(n, partial_mul<N>(big_int<M, T>{Modulus...}, quot)));

  return {quot, rem};
}

} // namespace lam::cbn
