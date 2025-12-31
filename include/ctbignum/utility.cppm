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

export module lam.ctbignum:utility;

import std;

import :bigint;

namespace lam::cbn
{
namespace detail
{

export 
template <std::size_t N, typename T> 
constexpr auto tight_length(big_int<N, T> num)
{
  // count the effective number of limbs
  // (ignoring zero-limbs at the most-significant-limb side)
  std::size_t L = N;
  while (L > 0 && num[L - 1] == 0)
    --L;

  return L;
}

export 
template <typename T, T... Is> 
constexpr auto tight_length(std::integer_sequence<T, Is...>)

{
  // count the effective number of limbs
  // (ignoring zero-limbs at the most-significant-limb side)

  std::size_t L = sizeof...(Is);
  std::array<T, sizeof...(Is)> num{Is...};
  while (L > 0 && num[L - 1] == 0)
    --L;

  return L;
}

export 
template <std::size_t N, typename T> 
constexpr auto bit_length(big_int<N, T> num)
{
  // we define bit_length(0) := 1

  auto L = tight_length(num);
  L += (L == 0U); // ensure L > 0
  std::size_t bitlen = L * std::numeric_limits<T>::digits;
  T msb = num[L - 1];
  while (bitlen > 1 && (msb & (static_cast<T>(1) << (std::numeric_limits<T>::digits - 1))) == 0)
  {
    msb <<= 1;
    --bitlen;
  }
  return bitlen;
}

export 
template <std::size_t N1, std::size_t N2, typename T>
constexpr void assign(big_int<N1, T> &dst, big_int<N2, T> src)
{
  // assignment for the scenario where N1 >= N2

  static_assert(N1 >= N2, "cannot assign: destination has smaller size than source");

  for (std::size_t i = 0; i < N1; ++i)
    dst[i] = src[i];
  for (auto i = N1; i < N2; ++i)
    dst[i] = 0;
}

} // namespace detail

export 
template <std::size_t ExplicitLength = 0, typename T, T... Limbs>
constexpr auto to_big_int(std::integer_sequence<T, Limbs...>)
{ return big_int<ExplicitLength ? ExplicitLength : sizeof...(Limbs), T>{Limbs...}; }

} // namespace lam::cbn
