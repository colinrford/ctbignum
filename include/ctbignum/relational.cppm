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

export module lam.ctbignum:relational;

import std;

import :bigint;

namespace lam::cbn
{

export template <typename T, std::size_t N1, std::size_t N2>
constexpr bool operator<(big_int<N1, T> a, big_int<N2, T> b)
{
  constexpr auto L = std::max(N1, N2);
  for (auto i = L; i > 0; --i)
  {
    T ai = (i <= N1) ? a[i - 1] : 0;
    T bi = (i <= N2) ? b[i - 1] : 0;
    if (ai < bi)
      return true;
    if (ai > bi)
      return false;
  }
  return false;
}

export template <typename T, std::size_t N1, std::size_t N2>
constexpr bool operator>(big_int<N1, T> a, big_int<N2, T> b)
{
  return b < a;
}

export template <typename T, std::size_t N1, std::size_t N2>
constexpr bool operator<=(big_int<N1, T> a, big_int<N2, T> b)
{
  return !(b < a);
}

export template <typename T, std::size_t N1, std::size_t N2>
constexpr bool operator>=(big_int<N1, T> a, big_int<N2, T> b)
{
  return !(a < b);
}

export template <typename T, std::size_t N1, std::size_t N2>
constexpr bool operator==(big_int<N1, T> a, big_int<N2, T> b)
{
  constexpr auto L = std::max(N1, N2);
  for (std::size_t i = 0; i < L; ++i)
  {
    T ai = (i < N1) ? a[i] : 0;
    T bi = (i < N2) ? b[i] : 0;
    if (ai != bi)
      return false;
  }
  return true;
}

export template <typename T, std::size_t N1, std::size_t N2>
constexpr bool operator!=(big_int<N1, T> a, big_int<N2, T> b)
{
  return !(a == b);
}

} // end namespace lam::cbn
