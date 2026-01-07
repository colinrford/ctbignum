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
#include "catch.hpp"

import std;
import lam.ctbignum;

TEST_CASE("Modular square root - Tonelli-Shanks")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("sqrt(0) = 0")
  {
    using GF = decltype(Zq(17_Z));
    constexpr GF zero{0};
    auto result = sqrt(zero);
    REQUIRE(result.has_value());
    REQUIRE(result->data == zero.data);
  }

  SECTION("sqrt(1) = 1")
  {
    using GF = decltype(Zq(17_Z));
    constexpr GF one{1};
    auto result = sqrt(one);
    REQUIRE(result.has_value());
    REQUIRE(result->data == one.data);
  }

  SECTION("Known square: sqrt(4) mod 17 = ±2")
  {
    using GF = decltype(Zq(17_Z));
    constexpr GF four{4};
    constexpr GF two{2};
    constexpr GF neg_two{-2}; // which is 15 mod 17

    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == neg_two.data));
  }

  SECTION("Known square: sqrt(9) mod 17 = ±3")
  {
    using GF = decltype(Zq(17_Z));
    constexpr GF nine{9};
    constexpr GF three{3};
    constexpr GF neg_three{-3}; // 14 mod 17

    auto result = sqrt(nine);
    REQUIRE(result.has_value());
    REQUIRE((result->data == three.data || result->data == neg_three.data));
  }

  SECTION("Prime p ≡ 3 (mod 4): p = 7")
  {
    // 7 ≡ 3 (mod 4), uses simple formula
    using GF = decltype(Zq(7_Z));
    constexpr GF four{4};
    constexpr GF two{2};
    constexpr GF neg_two{-2}; // 5 mod 7

    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == neg_two.data));
  }

  SECTION("Prime p ≡ 1 (mod 4): p = 13")
  {
    // 13 ≡ 1 (mod 4), uses full Tonelli-Shanks
    using GF = decltype(Zq(13_Z));
    constexpr GF four{4};
    constexpr GF two{2};
    constexpr GF neg_two{-2}; // 11 mod 13

    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == neg_two.data));
  }

  SECTION("Large prime: curve25519 prime")
  {
    // p = 2^255 - 19
    constexpr auto p25519 = 57896044618658097711785492504343953926634992332820282019728792003956564819949_Z;
    using GF = decltype(Zq(p25519));

    // Test sqrt(4) = ±2
    constexpr GF four{4};
    constexpr GF two{2};
    constexpr GF neg_two{-2};

    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == neg_two.data));

    // Verify by squaring
    auto squared = (*result) * (*result);
    REQUIRE(squared.data == four.data);
  }

  SECTION("Large prime: secp256k1 prime")
  {
    // p = 2^256 - 2^32 - 977
    constexpr auto secp256k1_p = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    using GF = decltype(Zq(secp256k1_p));

    // sqrt(4) = 2 (Trivial check)
    constexpr GF four{4};
    constexpr GF two{2};
    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == (GF{0} - two).data));

    // Non-trivial check: Pick a large number, square it, then take sqrt
    // Let's us a large arbitrary number < p
    // 2^200 + 12345
    constexpr auto large_val = 1606938044258990275541962092341162602522202993782792835301376_Z; // approx 2^200
    constexpr GF large_elem{large_val};

    // Calculate square: s = x^2
    constexpr auto square = large_elem * large_elem;

    // Recover root: r = sqrt(s)
    auto root = sqrt(square);
    REQUIRE(root.has_value());

    // Check: r == x OR r == -x
    // Note: -x in field is (0 - x)
    constexpr GF neg_large_elem = GF{0} - large_elem;

    bool correct_root = (root->data == large_elem.data) || (root->data == neg_large_elem.data);
    REQUIRE(correct_root);

    // Verify square property explicitly
    REQUIRE(((*root) * (*root)).data == square.data);
  }

  SECTION("Compile-Time Execution (constexpr)")
  {
    using GF = decltype(Zq(17_Z));

    // This MUST compile if sqrt is truly constexpr
    constexpr GF four{4};
    constexpr auto root = sqrt(four);

    static_assert(root.has_value());
    static_assert(root->data == to_big_int(2_Z) || root->data == to_big_int(15_Z));
  }

  SECTION("Composite modulus (safety check): mod 15")
  {
    // 15 = 3 * 5, not a prime. sqrt should return nullopt.
    // Even for valid squares like 4 (2^2 = 4), we reject the operation because
    // Tonelli-Shanks is only for primes.
    using GF = decltype(Zq(15_Z));
    constexpr GF four{4};

    auto result = sqrt(four);
    REQUIRE_FALSE(result.has_value());
  }

  SECTION("Carmichael number 1729 (safety check)")
  {
    // 1729 is a Carmichael number (pseudoprime).
    // Miller-Rabin should correctly identify it as composite and return nullopt.
    using GF = decltype(Zq(1729_Z));
    constexpr GF four{4};
    // 4 is a square mod 1729 (2^2), but we reject it because 1729 is composite.

    auto result = sqrt(four);
    REQUIRE_FALSE(result.has_value());
  }
}

// Helper for fuzzing
template <typename T>
class Randomizer
{
public:
  template <size_t N>
  void operator()(lam::cbn::big_int<N, T> &a)
  {
    for (auto &limb : a) limb = distribution(generator);
  }

private:
  std::mt19937 generator{std::random_device{}()};
  std::uniform_int_distribution<T> distribution;
};

TEST_CASE("Stress Tests and Edge Cases")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("High valuation of 2: p = 65537 (2^16 + 1)")
  {
    // 65537 is F4. p-1 = 2^16. S = 16.
    // This exercises the Tonelli-Shanks loop depth.
    using GF = decltype(Zq(65537_Z));

    // 3 is a primitive root mod 65537. 3^2 = 9.
    constexpr GF nine{9};
    auto result = sqrt(nine);
    REQUIRE(result.has_value());
    REQUIRE((result->data == to_big_int(3_Z) || result->data == (to_big_int(65537_Z) - to_big_int(3_Z))));

    // Check various powers
    constexpr GF val{123};
    auto sq = val * val;
    auto root = sqrt(sq);
    REQUIRE(root.has_value());
    REQUIRE(((*root * *root).data == sq.data));
  }

  SECTION("Randomized Fuzzing on secp256k1")
  {
    constexpr auto p_seq = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    using GF = decltype(Zq(p_seq));
    constexpr auto p_val = to_big_int(p_seq);

    Randomizer<uint64_t> randomizer;
    big_int<4, uint64_t> rand_data;

    // Run 50 random iterations
    for (int i = 0; i < 50; ++i)
    {
      // Generate random valid field element
      do { randomizer(rand_data); } while (rand_data >= p_val); // Ensure < p

      GF r{rand_data};

      // Square it
      auto sq = r * r;

      // Compute sqrt
      auto root = sqrt(sq);

      // Verify
      REQUIRE(root.has_value());
      auto root_sq = (*root) * (*root);
      REQUIRE(root_sq.data == sq.data);

      bool is_orig = (root->data == r.data);
      bool is_neg = (root->data == (GF{0} - r).data);
      REQUIRE((is_orig || is_neg));
    }
  }

  SECTION("Cube Root Fuzzing on secp256k1 (p = 1 mod 3 - Hard Case)")
  {
    // secp256k1 p is = 1 mod 3.
    // This means only 1/3 of elements have cube roots.
    // And there are 3 cube roots for any cubic residue.

    constexpr auto p_seq = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    using GF = decltype(Zq(p_seq));
    constexpr auto p_val = to_big_int(p_seq);

    Randomizer<uint64_t> randomizer;
    big_int<4, uint64_t> rand_data;

    for (int i = 0; i < 50; ++i)
    {
      // 1. Generate valid random r
      do { randomizer(rand_data); } while (rand_data >= p_val);
      GF r{rand_data};

      // 2. Compute cubic residue x = r^3
      auto x = r * r * r;
      // 3. Verify cbrt(x) finds a root
      auto root = cbrt(x);
      REQUIRE(root.has_value());
      // Verify root^3 == x
      auto val = *root;
      auto cubed = val * val * val;
      REQUIRE(cubed.data == x.data);
      // 4. Test random element (likely non-residue)
      do { randomizer(rand_data); } while (rand_data >= p_val);
      GF z{rand_data};

      auto z_root = cbrt(z);
      if (z_root.has_value())
      {
        auto z3 = (*z_root) * (*z_root) * (*z_root);
        REQUIRE(z3.data == z.data);
      }

      else
      {
        // Verify it was correctly rejected (Euler criterion for cubes)
        // z^((p - 1) / 3) != 1
        constexpr auto one = big_int<4, uint64_t>{1};
        constexpr auto three = big_int<4, uint64_t>{3};
        auto p_minus_1 = subtract_ignore_carry(p_val, one);
        auto exp = div(p_minus_1, three).quotient;
        auto res = mod_exp(z.data, exp, p_seq);
        REQUIRE(res != one);
      }
    }
  }
}

TEST_CASE("is_quadratic_residue")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  using GF17 = decltype(Zq(17_Z));

  SECTION("Quadratic residues mod 17: 1, 2, 4, 8, 9, 13, 15, 16")
  {
    REQUIRE(is_quadratic_residue(GF17{1}));
    REQUIRE(is_quadratic_residue(GF17{2}));
    REQUIRE(is_quadratic_residue(GF17{4}));
    REQUIRE(is_quadratic_residue(GF17{8}));
    REQUIRE(is_quadratic_residue(GF17{9}));
    REQUIRE(is_quadratic_residue(GF17{13}));
    REQUIRE(is_quadratic_residue(GF17{15}));
    REQUIRE(is_quadratic_residue(GF17{16}));
  }

  SECTION("Non-residues mod 17: 3, 5, 6, 7, 10, 11, 12, 14")
  {
    using GF = decltype(Zq(17_Z));

    // Check residue property
    REQUIRE_FALSE(is_quadratic_residue(GF17{3}));
    REQUIRE_FALSE(is_quadratic_residue(GF17{5}));

    // Check safe failure (returns nullopt)
    REQUIRE_FALSE(sqrt(GF17{3}).has_value());
    REQUIRE_FALSE(sqrt(GF17{5}).has_value());
    REQUIRE_FALSE(sqrt(GF17{6}).has_value());
  }
}

TEST_CASE("Cube root")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("cbrt(0) = 0")
  {
    // p = 7 ≡ 1 (mod 3), but 0 is always 0
    using GF = decltype(Zq(7_Z));
    constexpr GF zero{0};
    auto result = cbrt(zero);
    REQUIRE(result.has_value());
    REQUIRE(result->data == zero.data);
  }

  SECTION("cbrt(1) = 1 for p ≡ 2 (mod 3)")
  {
    // p = 5 ≡ 2 (mod 3): every element has unique cube root
    using GF = decltype(Zq(5_Z));
    constexpr GF one{1};
    auto result = cbrt(one);
    REQUIRE(result.has_value());
    REQUIRE(result->data == one.data);
  }

  SECTION("cbrt(8) = 2 for p ≡ 2 (mod 3)")
  {
    // p = 11 ≡ 2 (mod 3)
    using GF = decltype(Zq(11_Z));
    constexpr GF eight{8};
    constexpr GF two{2};
    auto result = cbrt(eight);
    REQUIRE(result.has_value());
    // Verify: result^3 = 8
    auto cubed = (*result) * (*result) * (*result);
    REQUIRE(cubed.data == eight.data);
  }

  SECTION("Composite modulus (safety check): mod 15")
  {
    using GF = decltype(Zq(15_Z));
    constexpr GF one{1};
    auto result = cbrt(one);
    REQUIRE_FALSE(result.has_value());
  }
}
