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

export module lam.ctbignum:roots;
import std;

import :bigint;
import :field;
import :mod_exp;
import :addition;
import :bitshift;
import :utility;

namespace lam::cbn
{

namespace detail
{

// Factor p - 1 = Q * 2^S where Q is odd
template<std::size_t N, typename T>
constexpr auto factor_out_twos(big_int<N, T> n)
{
  std::size_t S = 0;
  while ((n[0] & 1) == 0)
  {
    n = shift_right(n, 1);
    ++S;
  }
  return std::pair{n, S};
}

// Factor n = d * val^S
template<std::size_t N, typename T>
constexpr auto factor_out_val(big_int<N, T> n, T val)
{
  std::size_t S = 0;
  auto val_bi = big_int<N, T>{val};
  auto zero = big_int<N, T>{};

  while (true)
  {
    auto division = div(n, val_bi);
    if (division.remainder == zero)
    {
      n = division.quotient;
      ++S;
    }
    else
      break;
  }
  return std::pair{n, S};
}

// Miller-Rabin with the first 12 primes as bases: exact for n < 3.18e23,
// probabilistic above
template<typename T, T... Modulus>
constexpr bool is_prime(std::integer_sequence<T, Modulus...>)
{
  constexpr std::size_t N = sizeof...(Modulus);
  constexpr auto n = big_int<N, T>{Modulus...};
  constexpr auto one = big_int<N, T>{1};
  constexpr auto two = big_int<N, T>{2};

  if (n < two)
    return false;
  if (n == two || n == big_int<N, T>{3})
    return true;
  if ((n[0] & 1) == 0)
    return false;

  // n - 1 = d * 2^r
  auto n_minus_1 = subtract_ignore_carry(n, one);
  auto [d, r] = factor_out_twos(n_minus_1);

  constexpr unsigned bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};

  for (auto a_val : bases)
  {
    auto a = big_int<N, T>{static_cast<T>(a_val)};
    if (n <= a)
      break;

    auto x = mod_exp(a, d, std::integer_sequence<T, Modulus...>{});

    if (x == one || x == n_minus_1)
      continue;

    bool witness_found = false;
    for (std::size_t i = 1; i < r; ++i)
    {
      x = mod_exp(x, two, std::integer_sequence<T, Modulus...>{});
      if (x == n_minus_1)
      {
        witness_found = true;
        break;
      }
    }

    if (!witness_found)
      return false; // Composite
  }
  return true; // Likely prime
}

} // namespace detail

// Check if n is a quadratic residue mod p using Euler's criterion
// n^((p - 1) / 2) ≡ 1 (mod p) iff n is a quadratic residue
export template<typename T, T... Modulus>
constexpr bool is_quadratic_residue(ZqElement<T, Modulus...> n)
{
  constexpr auto p = big_int<sizeof...(Modulus), T>{Modulus...};
  constexpr auto one = big_int<sizeof...(Modulus), T>{1};
  constexpr auto p_minus_1 = subtract_ignore_carry(p, one);
  constexpr auto exp = shift_right(p_minus_1, 1); // (p - 1) / 2

  if (n.data == big_int<sizeof...(Modulus), T>{})
    return true; // 0 is considered a residue (sqrt(0) = 0)

  auto result = mod_exp(n.data, exp, std::integer_sequence<T, Modulus...>{});
  return result == one;
}

// Tonelli-Shanks algorithm for computing modular square roots
// Returns std::optional<ZqElement> - nullopt if n is not a quadratic residue
export template<typename T, T... Modulus>
constexpr auto sqrt(ZqElement<T, Modulus...> n) -> std::optional<ZqElement<T, Modulus...>>
{
  constexpr std::size_t N = sizeof...(Modulus);
  constexpr auto p = big_int<N, T>{Modulus...};
  constexpr auto one = big_int<N, T>{1};
  constexpr auto two = big_int<N, T>{2};
  constexpr auto p_minus_1 = subtract_ignore_carry(p, one);
  if (n.data == big_int<N, T>{})
    return ZqElement<T, Modulus...>{};

  // Check primality to prevent infinite loops (and ensure field properties)
  if (!detail::is_prime(std::integer_sequence<T, Modulus...>{}))
    return std::nullopt;
  // Check if n is a quadratic residue to prevent infinite loops
  if (!is_quadratic_residue(n))
    return std::nullopt;
  // Factor p - 1 = Q * 2^S
  constexpr auto QS = detail::factor_out_twos(p_minus_1);
  constexpr auto Q = QS.first;
  constexpr std::size_t S = QS.second;
  // Special case: p ≡ 3 (mod 4), i.e., S == 1
  // sqrt(n) = n^((p + 1) / 4)
  if constexpr (S == 1)
  {
    constexpr auto exp = shift_right(add_ignore_carry(p, one), 2);
    auto result = mod_exp(n.data, exp, std::integer_sequence<T, Modulus...>{});
    return ZqElement<T, Modulus...>{result};
  }
  else
  { // Find a quadratic non-residue z
    constexpr auto neg_one = subtract_ignore_carry(p, one);
    constexpr auto legendre_exp = shift_right(p_minus_1, 1);
    big_int<N, T> z{2};
    while (mod_exp(z, legendre_exp, std::integer_sequence<T, Modulus...>{}) != neg_one)
      z = add_ignore_carry(z, one);
    std::size_t M = S;
    auto c = mod_exp(z, Q, std::integer_sequence<T, Modulus...>{});
    auto t = mod_exp(n.data, Q, std::integer_sequence<T, Modulus...>{});
    // R = n^((Q + 1) / 2)
    auto Q_plus_1_div_2 = shift_right(add_ignore_carry(Q, one), 1);
    auto R = mod_exp(n.data, Q_plus_1_div_2, std::integer_sequence<T, Modulus...>{});
    while (t != one)
    { // Find the least i such that t^(2^i) = 1
      std::size_t i = 1;
      auto temp = mod_exp(t, two, std::integer_sequence<T, Modulus...>{});
      while (temp != one && i < M)
      {
        temp = mod_exp(temp, two, std::integer_sequence<T, Modulus...>{});
        ++i;
      }
      // b = c^(2^(M - i - 1))
      auto b = c;
      for (std::size_t j = 0; j < M - i - 1; ++j)
        b = mod_exp(b, two, std::integer_sequence<T, Modulus...>{});
      M = i;
      c = mod_exp(b, two, std::integer_sequence<T, Modulus...>{});
      // t = t * c, R = R * b
      t = mod(mul(t, c), std::integer_sequence<T, Modulus...>{});
      R = mod(mul(R, b), std::integer_sequence<T, Modulus...>{});
    }
    return ZqElement<T, Modulus...>{R};
  }
}

// Cube root in finite field
// Uses the formula: cbrt(a) = a^((2 * (p - 1)) / 3) when p ≡ 2 (mod 3)
// Returns nullopt if n is not a cubic residue (when p ≡ 1 mod 3)
export template<typename T, T... Modulus>
constexpr auto cbrt(ZqElement<T, Modulus...> n) -> std::optional<ZqElement<T, Modulus...>>
{
  constexpr std::size_t N = sizeof...(Modulus);
  constexpr auto p = big_int<N, T>{Modulus...};
  constexpr auto one = big_int<N, T>{1};
  constexpr auto two = big_int<N, T>{2};
  constexpr auto three = big_int<N, T>{3};
  constexpr auto p_minus_1 = subtract_ignore_carry(p, one);

  if (n.data == big_int<N, T>{})
    return ZqElement<T, Modulus...>{};

  // Check primality
  if (!detail::is_prime(std::integer_sequence<T, Modulus...>{}))
    return std::nullopt;

  // For prime p, gcd(3, p - 1) is either 1 (p ≡ 2 mod 3) or 3 (p ≡ 1 mod 3)
  constexpr auto p_mod_3 = div(p, three).remainder[0];

  if constexpr (p_mod_3 == 2)
  { // p ≡ 2 (mod 3): every element is a cubic residue, unique cube root
    // cbrt(n) = n^((2p - 1) / 3) = n^(2k + 1) with p = 3k + 2; avoids computing 2p
    constexpr auto k = div(subtract_ignore_carry(p, two), three).quotient;
    constexpr auto exp = add_ignore_carry(add_ignore_carry(k, k), one);
    auto result = mod_exp(n.data, exp, std::integer_sequence<T, Modulus...>{});
    return ZqElement<T, Modulus...>{result};
  }
  else
  {
    // p ≡ 1 (mod 3): Adleman-Manders-Miller
    constexpr auto residue_exp = div(p_minus_1, three).quotient;
    if (mod_exp(n.data, residue_exp, std::integer_sequence<T, Modulus...>{}) != one)
    {
      return std::nullopt;
    }

    // p - 1 = 3^S * t with 3 not dividing t
    constexpr auto Qt = detail::factor_out_val(p_minus_1, static_cast<T>(3));
    constexpr auto t = Qt.first;
    constexpr std::size_t S = Qt.second;

    // k = 1/3 mod t = (1 + qt) / 3, with q in {1, 2} chosen so 3 divides 1 + qt
    auto t_rem_3 = div(t, three).remainder;
    auto k = t;
    if (t_rem_3 == big_int<sizeof...(Modulus), T>{1})
    {
      auto two_t_full = mul(two, t);
      // 1 + 2t < p, so 2t fits in N limbs
      constexpr auto N_limbs = sizeof...(Modulus);
      big_int<N_limbs, T> two_t;
      for (std::size_t i = 0; i < N_limbs; ++i)
        two_t[i] = two_t_full[i];

      auto num = add_ignore_carry(one, two_t);
      k = div(num, three).quotient;
    }
    else
    {
      // t ≡ 2 (mod 3)
      auto num = add_ignore_carry(one, t);
      k = div(num, three).quotient;
    }
    // cubic non-residue z: z^((p - 1) / 3) != 1
    auto z = two;
    while (true)
    {
      if (mod_exp(z, residue_exp, std::integer_sequence<T, Modulus...>{}) != one)
        break;
      z = add_ignore_carry(z, one);
    }

    auto c = mod_exp(z, t, std::integer_sequence<T, Modulus...>{});
    auto r = mod_exp(n.data, k, std::integer_sequence<T, Modulus...>{});
    auto n_inv = mod_inv(n.data, p);

    auto h = r;
    auto h_cubed = mod_exp(h, three, std::integer_sequence<T, Modulus...>{});
    auto prod = mul(n_inv, h_cubed);
    auto b = mod(prod, std::integer_sequence<T, Modulus...>());

    auto x = r;
    auto g = c;
    auto r_val = S;
    auto b_term = b;
    std::size_t safety_counter = 0;
    while (b_term != one)
    {
      if (++safety_counter > S + 10) // converges within S steps for prime p
        return std::nullopt;

      // least m with b^(3^m) = 1
      auto temp = b_term;
      std::size_t m = 0;
      for (; m < r_val; ++m)
      {
        if (temp == one)
          break;
        temp = mod_exp(temp, three, std::integer_sequence<T, Modulus...>{});
      }

      if (m == r_val)
        return std::nullopt;

      // g = g^(3^(r - m - 1))
      auto g2 = g;
      for (std::size_t j = 0; j < r_val - m - 1; ++j)
        g2 = mod_exp(g2, three, std::integer_sequence<T, Modulus...>{});

      g = g2;

      // x = x * g
      // b = b * g^3
      x = mod(mul(x, g), std::integer_sequence<T, Modulus...>());
      g = mod_exp(g, three, std::integer_sequence<T, Modulus...>{});
      b_term = mod(mul(b_term, g), std::integer_sequence<T, Modulus...>());
      r_val = m;
    }

    return ZqElement<T, Modulus...>{x};
  }
}

} // namespace lam::cbn
