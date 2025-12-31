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

export module lam.ctbignum:addition;

import std;

import :bigint;
import :slicing;
import :utility;

namespace lam::cbn
{

export 
template <typename T, std::size_t M, std::size_t N> 
constexpr auto add(big_int<M, T> a, big_int<N, T> b)
{
  constexpr auto L = std::max(M, N);
  return add_same(detail::pad<L - M>(a), detail::pad<L - N>(b));
}

export 
template <typename T, std::size_t M, std::size_t N> 
constexpr auto subtract(big_int<M, T> a, big_int<N, T> b)
{
  constexpr auto L = std::max(M, N);
  return subtract_same(detail::pad<L - M>(a), detail::pad<L - N>(b));
}

export 
template <typename T, std::size_t N> 
constexpr auto add_same(big_int<N, T> a, big_int<N, T> b)
{
  T carry{};
  big_int<N + 1, T> r{};

  for (auto i = 0U; i < N; ++i)
  {
    auto aa = a[i];
    auto sum = aa + b[i];
    auto res = sum + carry;
    carry = (sum < aa) | (res < sum);
    r[i] = res;
  }

  r[N] = carry;
  return r;
}

export 
template <typename T, std::size_t N> 
constexpr auto subtract_same(big_int<N, T> a, big_int<N, T> b)
{
  T carry{};
  big_int<N + 1, T> r{};

  for (auto i = 0U; i < N; ++i)
  {
    auto aa = a[i];
    auto diff = aa - b[i];
    auto res = diff - carry;
    carry = (diff > aa) | (res > diff);
    r[i] = res;
  }

  r[N] = carry * static_cast<T>(-1); // sign extension
  return r;
}

export 
template <typename T, std::size_t N> 
constexpr auto add_ignore_carry(big_int<N, T> a, big_int<N, T> b)
{
  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i)
  {
    T aa = a[i];
    T sum = aa + b[i];
    T res = sum + carry;
    carry = (sum < aa) | (res < sum);
    r[i] = res;
  }

  return r;
}

export 
template <typename T, std::size_t N> 
constexpr auto subtract_ignore_carry(big_int<N, T> a, big_int<N, T> b)
{
  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i)
  {
    auto aa = a[i];
    auto diff = aa - b[i];
    auto res = diff - carry;
    carry = (diff > aa) | (res > diff);
    r[i] = res;
  }

  return r;
}

export 
template <typename T, std::size_t N>
constexpr auto mod_add(big_int<N, T> a, big_int<N, T> b, big_int<N, T> modulus)
{
  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i)
  {
    auto aa = a[i];
    auto sum = aa + b[i];
    auto res = sum + carry;
    carry = (sum < aa) | (res < sum);
    r[i] = res;
  }

  auto reduced = subtract(r, modulus);
  T r_geq_modulus = reduced[N] ? 0 : 1;
  big_int<N, T> res = (carry + r_geq_modulus != 0) ? detail::first<N>(reduced) : r;
  return res;
}

export 
template <typename T, std::size_t N>
constexpr auto mod_sub(big_int<N, T> a, big_int<N, T> b, big_int<N, T> modulus)
{
  T carry{};
  big_int<N, T> r{};

  for (auto i = 0U; i < N; ++i)
  {
    auto aa = a[i];
    auto diff = aa - b[i];
    auto res = diff - carry;
    carry = (diff > aa) | (res > diff);
    r[i] = res;
  }

  auto adjusted_r = add_ignore_carry(r, modulus);
  big_int<N, T> res = carry ? adjusted_r : r;
  return res;
}

export 
template <typename T, std::size_t N, T... Modulus>
constexpr auto mod_add(big_int<N, T> a, big_int<N, T> b, std::integer_sequence<T, Modulus...>)
{
  big_int<sizeof...(Modulus), T> modulus{{Modulus...}};
  return mod_add(a, b, modulus);
}

export 
template <typename T, std::size_t N1, std::size_t N2>
constexpr auto operator+(big_int<N1, T> a, big_int<N2, T> b)
{
  return add(a, b);
}

export 
template <typename T, std::size_t N1, std::size_t N2>
constexpr auto operator-(big_int<N1, T> a, big_int<N2, T> b)
{
  return subtract(a, b);
}

export 
template <typename T, std::size_t N, T... Modulus>
constexpr auto mod_sub(big_int<N, T> a, big_int<N, T> b, std::integer_sequence<T, Modulus...>)
{
  big_int<sizeof...(Modulus), T> modulus{{Modulus...}};
  return mod_sub(a, b, modulus);
}

} // end namespace lam::cbn
