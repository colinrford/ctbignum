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
    
    // sqrt(4) = 2
    constexpr GF four{4};
    constexpr GF two{2};
    
    // This implicitly tests is_prime(secp256k1_p)
    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == (GF{0} - two).data));
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
