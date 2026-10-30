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

    constexpr GF four{4};
    constexpr GF two{2};
    constexpr GF neg_two{-2};

    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == neg_two.data));

    auto squared = (*result) * (*result);
    REQUIRE(squared.data == four.data);
  }

  SECTION("Large prime: secp256k1 prime")
  {
    // p = 2^256 - 2^32 - 977
    constexpr auto secp256k1_p = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    using GF = decltype(Zq(secp256k1_p));

    constexpr GF four{4};
    constexpr GF two{2};
    auto result = sqrt(four);
    REQUIRE(result.has_value());
    REQUIRE((result->data == two.data || result->data == (GF{0} - two).data));

    // 2^200
    constexpr auto large_val = 1606938044258990275541962092341162602522202993782792835301376_Z;
    constexpr GF large_elem{large_val};

    constexpr auto square = large_elem * large_elem;

    auto root = sqrt(square);
    REQUIRE(root.has_value());

    constexpr GF neg_large_elem = GF{0} - large_elem;

    bool correct_root = (root->data == large_elem.data) || (root->data == neg_large_elem.data);
    REQUIRE(correct_root);

    REQUIRE(((*root) * (*root)).data == square.data);
  }

  SECTION("Compile-Time Execution (constexpr)")
  {
    using GF = decltype(Zq(17_Z));

    constexpr GF four{4};
    constexpr auto root = sqrt(four);

    static_assert(root.has_value());
    static_assert(root->data == to_big_int(2_Z) || root->data == to_big_int(15_Z));
  }

}

template<typename E>
concept has_sqrt = requires(E e) { sqrt(e); };

template<typename E>
concept has_cbrt = requires(E e) { cbrt(e); };

TEST_CASE("Composite moduli are rejected at compile time")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  using GF15 = decltype(Zq(15_Z));
  using GF1729 = decltype(Zq(1729_Z)); // Carmichael number
  using GF17 = decltype(Zq(17_Z));

  STATIC_REQUIRE_FALSE(has_sqrt<GF15>);
  STATIC_REQUIRE_FALSE(has_cbrt<GF15>);
  STATIC_REQUIRE_FALSE(has_sqrt<GF1729>);
  STATIC_REQUIRE_FALSE(has_cbrt<GF1729>);
  STATIC_REQUIRE(has_sqrt<GF17>);
  STATIC_REQUIRE(has_cbrt<GF17>);
  STATIC_REQUIRE(std::is_same_v<decltype(Fp(17_Z)), GF17>);
}

template<typename T>
class Randomizer
{
public:
  template<size_t N>
  void operator()(lam::cbn::big_int<N, T>& a)
  {
    for (auto& limb : a)
      limb = distribution(generator);
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
    // p - 1 = 2^16, so S = 16: the deepest Tonelli-Shanks loop here
    using GF = decltype(Zq(65537_Z));

    constexpr GF nine{9};
    auto result = sqrt(nine);
    REQUIRE(result.has_value());
    REQUIRE((result->data == to_big_int(3_Z) || result->data == (to_big_int(65537_Z) - to_big_int(3_Z))));

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

    for (int i = 0; i < 50; ++i)
    {
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p_val);

      GF r{rand_data};

      auto sq = r * r;

      auto root = sqrt(sq);

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
    // secp256k1 p ≡ 1 (mod 3): a third of the nonzero elements are cubes, each with 3 cube roots
    constexpr auto p_seq = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    using GF = decltype(Zq(p_seq));
    constexpr auto p_val = to_big_int(p_seq);

    Randomizer<uint64_t> randomizer;
    big_int<4, uint64_t> rand_data;

    for (int i = 0; i < 50; ++i)
    {
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p_val);
      GF r{rand_data};

      auto x = r * r * r;
      auto root = cbrt(x);
      REQUIRE(root.has_value());
      auto val = *root;
      auto cubed = val * val * val;
      REQUIRE(cubed.data == x.data);
      // random element, a non-residue with probability 2/3
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p_val);
      GF z{rand_data};

      auto z_root = cbrt(z);
      if (z_root.has_value())
      {
        auto z3 = (*z_root) * (*z_root) * (*z_root);
        REQUIRE(z3.data == z.data);
      }

      else
      {
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
    REQUIRE_FALSE(is_quadratic_residue(GF17{3}));
    REQUIRE_FALSE(is_quadratic_residue(GF17{5}));

    REQUIRE_FALSE(sqrt(GF17{3}).has_value());
    REQUIRE_FALSE(sqrt(GF17{5}).has_value());
    REQUIRE_FALSE(sqrt(GF17{6}).has_value());
  }
}

TEST_CASE("Cube root")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  // Two-limb primes where p mod 3 differs from the low limb mod 3, and one
  // with the top bit set so that 2p overflows two limbs.
  SECTION("Two-limb p ≡ 2 (mod 3), low limb ≡ 1 (mod 3)")
  {
    using GF = decltype(Zq(1267650600228229401496703205653_Z));
    constexpr GF x{125};
    auto root = cbrt(x);
    REQUIRE(root.has_value());
    REQUIRE(*root == GF{5});
  }

  SECTION("Two-limb p ≡ 1 (mod 3), low limb ≡ 2 (mod 3)")
  {
    using GF = decltype(Zq(2535301200456458802993406410901_Z));
    constexpr GF x{125};
    auto root = cbrt(x);
    REQUIRE(root.has_value());
    REQUIRE((*root) * (*root) * (*root) == x);
  }

  SECTION("128-bit p ≡ 2 (mod 3) with top bit set")
  {
    using GF = decltype(Zq(170141183460469231750134047789593657433_Z));
    constexpr GF x{125};
    auto root = cbrt(x);
    REQUIRE(root.has_value());
    REQUIRE(*root == GF{5});
  }

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
    // TODO: decide whether two is meant to be used
    [[maybe_unused]] constexpr GF two{2};
    auto result = cbrt(eight);
    REQUIRE(result.has_value());
    // Verify: result^3 = 8
    auto cubed = (*result) * (*result) * (*result);
    REQUIRE(cubed.data == eight.data);
  }

  SECTION("Compile-Time Execution (constexpr)")
  {
    using GF = decltype(Zq(5_Z)); // 5 ≡ 2 (mod 3): unique cube roots
    constexpr GF two{2};
    // 3^3 = 27 ≡ 2 (mod 5)
    constexpr auto res = cbrt(two);
    static_assert(res.has_value());
    constexpr GF three{3};
    static_assert(res->data == three.data);

    // 7 ≡ 1 (mod 3): Adleman-Manders-Miller
    using GF7 = decltype(Zq(7_Z));
    constexpr GF7 six{6}; // -1
    constexpr auto res_amm = cbrt(six);
    static_assert(res_amm.has_value());
    static_assert(res_amm->data == six.data);
  }
}

TEST_CASE("Comprehensive Primes Matrix (Sqrt + Cbrt)", "[roots]")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  // 5 (1 mod 4, 2 mod 3)
  SECTION("Prime 5: Sqrt(TS) + Cbrt(Unique)")
  {
    using GF = decltype(Zq(5_Z));
    constexpr GF four{4};

    // Runtime
    auto r_sqrt = sqrt(four);
    REQUIRE(r_sqrt.has_value());
    REQUIRE((r_sqrt->data == to_big_int(2_Z) || r_sqrt->data == to_big_int(3_Z))); // 2^2=4, 3^2=9=4

    auto r_cbrt = cbrt(four); // 4 = -1. (-1)^3 = -1 = 4.
    REQUIRE(r_cbrt.has_value());
    REQUIRE(r_cbrt->data == to_big_int(4_Z));

    // Compile-Time
    constexpr auto c_sqrt = sqrt(four);
    static_assert(c_sqrt.has_value());
    static_assert(c_sqrt->data == to_big_int(2_Z) || c_sqrt->data == to_big_int(3_Z));

    constexpr auto c_cbrt = cbrt(four);
    static_assert(c_cbrt.has_value());
    static_assert(c_cbrt->data == to_big_int(4_Z));
  }

  // 7 (3 mod 4, 1 mod 3)
  SECTION("Prime 7: Sqrt(Simple) + Cbrt(AMM)")
  {
    using GF = decltype(Zq(7_Z));
    constexpr GF two{2};
    constexpr GF four{4};

    // Runtime
    auto r_sqrt = sqrt(two); // 3^2=2, 4^2=2
    REQUIRE(r_sqrt.has_value());
    REQUIRE((r_sqrt->data == to_big_int(3_Z) || r_sqrt->data == to_big_int(4_Z)));

    auto r_sqrt4 = sqrt(four); // 2^2=4, 5^2=4
    REQUIRE(r_sqrt4.has_value());
    REQUIRE((r_sqrt4->data == to_big_int(2_Z) || r_sqrt4->data == to_big_int(5_Z)));

    // cbrt(1) -> 1, 2, 4 (roots of unity)
    // cbrt(6) -> 6
    constexpr GF six{6};
    auto r_cbrt = cbrt(six);
    REQUIRE(r_cbrt.has_value());
    REQUIRE(r_cbrt->data == to_big_int(6_Z));

    // Compile-Time
    constexpr auto c_sqrt = sqrt(two);
    static_assert(c_sqrt.has_value());
    static_assert(c_sqrt->data == to_big_int(3_Z) || c_sqrt->data == to_big_int(4_Z));

    constexpr auto c_cbrt = cbrt(six);
    static_assert(c_cbrt.has_value());
    static_assert(c_cbrt->data == to_big_int(6_Z));
  }

  // 11 (3 mod 4, 2 mod 3)
  SECTION("Prime 11: Sqrt(Simple) + Cbrt(Unique)")
  {
    using GF = decltype(Zq(11_Z));
    constexpr GF five{5}; // 4^2=16=5. 7^2=49=5.

    // Runtime
    auto r_sqrt = sqrt(five);
    REQUIRE(r_sqrt.has_value());
    REQUIRE((r_sqrt->data == to_big_int(4_Z) || r_sqrt->data == to_big_int(7_Z)));

    // cbrt(2) -> ? 7^3 = 343 = 31*11 + 2. So 7.
    constexpr GF two{2};
    auto r_cbrt = cbrt(two);
    REQUIRE(r_cbrt.has_value());
    REQUIRE(r_cbrt->data == to_big_int(7_Z));

    // Compile-Time
    constexpr auto c_sqrt = sqrt(five);
    static_assert(c_sqrt.has_value());
    static_assert(c_sqrt->data == to_big_int(4_Z) || c_sqrt->data == to_big_int(7_Z));

    constexpr auto c_cbrt = cbrt(two);
    static_assert(c_cbrt.has_value());
    static_assert(c_cbrt->data == to_big_int(7_Z));
  }

  // 13 (1 mod 4, 1 mod 3)
  SECTION("Prime 13: Sqrt(TS) + Cbrt(AMM)")
  {
    using GF = decltype(Zq(13_Z));

    // Sqrt: 3 (9 is sq). sqrt(9) = 3, 10.
    constexpr GF nine{9};
    auto r_sqrt = sqrt(nine);
    REQUIRE(r_sqrt.has_value());
    REQUIRE((r_sqrt->data == to_big_int(3_Z) || r_sqrt->data == to_big_int(10_Z)));

    // Cbrt: 8 = 2^3. Roots of 8 are 2, 2*3=6, 2*9=18=5.
    // 2^3=8. 6^3=216=16*13+8. 5^3=125=9*13+8.
    // So cbrt(8) must be one of {2, 5, 6}.
    constexpr GF eight{8};
    auto r_cbrt = cbrt(eight);
    REQUIRE(r_cbrt.has_value());
    bool is_valid_cbrt =
      (r_cbrt->data == to_big_int(2_Z) || r_cbrt->data == to_big_int(5_Z) || r_cbrt->data == to_big_int(6_Z));
    REQUIRE(is_valid_cbrt);

    // Compile-Time
    constexpr auto c_sqrt = sqrt(nine);
    static_assert(c_sqrt.has_value());
    static_assert(c_sqrt->data == to_big_int(3_Z) || c_sqrt->data == to_big_int(10_Z));

    constexpr auto c_cbrt = cbrt(eight);
    static_assert(c_cbrt.has_value());
    static_assert(c_cbrt->data == to_big_int(2_Z) || c_cbrt->data == to_big_int(5_Z) ||
                  c_cbrt->data == to_big_int(6_Z));
  }

  // 17 (1 mod 4, 2 mod 3)
  SECTION("Prime 17: Sqrt(TS) + Cbrt(Unique)")
  {
    using GF = decltype(Zq(17_Z));

    // Sqrt: 2 (6^2=36=2, 11^2=121=7*17+2).
    constexpr GF two{2};
    auto r_sqrt = sqrt(two);
    REQUIRE(r_sqrt.has_value());
    REQUIRE((r_sqrt->data == to_big_int(6_Z) || r_sqrt->data == to_big_int(11_Z)));

    // Cbrt: 1. cbrt(1)=1. (Unique because 2 mod 3).
    constexpr GF one{1};
    auto r_cbrt = cbrt(one);
    REQUIRE(r_cbrt.has_value());
    REQUIRE(r_cbrt->data == to_big_int(1_Z));

    // Compile-Time
    constexpr auto c_sqrt = sqrt(two);
    static_assert(c_sqrt.has_value());
    static_assert(c_sqrt->data == to_big_int(6_Z) || c_sqrt->data == to_big_int(11_Z));

    constexpr auto c_cbrt = cbrt(one);
    static_assert(c_cbrt.has_value());
    static_assert(c_cbrt->data == to_big_int(1_Z));
  }
}
