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
  // Catch's seed: fixed by default so failures reproduce; vary it with --rng-seed
  std::mt19937 generator{Catch::rngSeed()};
  std::uniform_int_distribution<T> distribution;
};

TEST_CASE("Cross-validate: NTTP sqrt vs runtime sqrt")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("Small primes")
  {
    // p = 7 (3 mod 4, S==1 branch)
    {
      constexpr auto p_seq = 7_Z;
      constexpr auto p = to_big_int(p_seq);
      using GF = decltype(Zq(p_seq));

      for (long v = 1; v < 7; ++v)
      {
        GF elem{v};
        auto nttp_result = sqrt(elem);
        auto rt_result = sqrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

        REQUIRE(nttp_result.has_value() == rt_result.has_value());
        if (nttp_result.has_value())
          REQUIRE(nttp_result->data == *rt_result);
      }
    }

    // p = 13 (1 mod 4, full Tonelli-Shanks)
    {
      constexpr auto p_seq = 13_Z;
      constexpr auto p = to_big_int(p_seq);
      using GF = decltype(Zq(p_seq));

      for (long v = 1; v < 13; ++v)
      {
        GF elem{v};
        auto nttp_result = sqrt(elem);
        auto rt_result = sqrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

        REQUIRE(nttp_result.has_value() == rt_result.has_value());
        if (nttp_result.has_value())
          REQUIRE(nttp_result->data == *rt_result);
      }
    }

    // p = 17 (1 mod 4, 2 mod 3)
    {
      constexpr auto p_seq = 17_Z;
      constexpr auto p = to_big_int(p_seq);
      using GF = decltype(Zq(p_seq));

      for (long v = 1; v < 17; ++v)
      {
        GF elem{v};
        auto nttp_result = sqrt(elem);
        auto rt_result = sqrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

        REQUIRE(nttp_result.has_value() == rt_result.has_value());
        if (nttp_result.has_value())
          REQUIRE(nttp_result->data == *rt_result);
      }
    }
  }

  SECTION("p = 2 (x^2 = x)")
  {
    constexpr auto p_seq = 2_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    REQUIRE(is_quadratic_residue(GF{1}));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{1}, p));
    auto nttp_result = sqrt(GF{1});
    auto rt_result = sqrt(big_int<1, std::uint64_t>{1}, p);
    REQUIRE(nttp_result.has_value());
    REQUIRE(rt_result.has_value());
    REQUIRE(nttp_result->data == big_int<1, std::uint64_t>{1});
    REQUIRE(*rt_result == nttp_result->data);
  }

  SECTION("Fermat prime 65537 (high 2-valuation)")
  {
    constexpr auto p_seq = 65537_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    // Test a selection of values
    for (long v : {1L, 2L, 3L, 4L, 9L, 16L, 100L, 123L, 999L, 65536L})
    {
      GF elem{v};
      auto nttp_result = sqrt(elem);
      auto rt_result = sqrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value() == rt_result.has_value());
      if (nttp_result.has_value())
        REQUIRE(nttp_result->data == *rt_result);
    }
  }

  SECTION("Large prime: secp256k1 randomized")
  {
    constexpr auto p_seq = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    Randomizer<std::uint64_t> randomizer;
    big_int<4, std::uint64_t> rand_data;

    for (int i = 0; i < 20; ++i)
    {
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p);

      GF elem{rand_data};
      auto sq = elem * elem; // guaranteed quadratic residue

      auto nttp_result = sqrt(sq);
      auto rt_result = sqrt(sq.data, p);

      REQUIRE(nttp_result.has_value());
      REQUIRE(rt_result.has_value());
      REQUIRE(nttp_result->data == *rt_result);
    }
  }
}

TEST_CASE("Cross-validate: NTTP cbrt vs runtime cbrt")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("p = 2 (x^3 = x)")
  {
    constexpr auto p_seq = 2_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    auto nttp_result = cbrt(GF{1});
    auto rt_result = cbrt(big_int<1, std::uint64_t>{1}, p);
    REQUIRE(nttp_result.has_value());
    REQUIRE(rt_result.has_value());
    REQUIRE(nttp_result->data == big_int<1, std::uint64_t>{1});
    REQUIRE(*rt_result == nttp_result->data);
  }

  SECTION("p = 3 (x^3 = x)")
  {
    constexpr auto p_seq = 3_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    for (long v = 1; v < 3; ++v)
    {
      GF elem{v};
      auto nttp_result = cbrt(elem);
      auto rt_result = cbrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value());
      REQUIRE(rt_result.has_value());
      REQUIRE(nttp_result->data == big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)});
      REQUIRE(*rt_result == nttp_result->data);
    }
  }

  SECTION("p = 5 (2 mod 3, unique roots)")
  {
    constexpr auto p_seq = 5_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    for (long v = 1; v < 5; ++v)
    {
      GF elem{v};
      auto nttp_result = cbrt(elem);
      auto rt_result = cbrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value() == rt_result.has_value());
      if (nttp_result.has_value())
        REQUIRE(nttp_result->data == *rt_result);
    }
  }

  SECTION("p = 7 (1 mod 3, AMM)")
  {
    constexpr auto p_seq = 7_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    for (long v = 1; v < 7; ++v)
    {
      GF elem{v};
      auto nttp_result = cbrt(elem);
      auto rt_result = cbrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value() == rt_result.has_value());
      if (nttp_result.has_value())
        REQUIRE(nttp_result->data == *rt_result);
    }
  }

  SECTION("p = 13 (1 mod 3, 1 mod 4, AMM + TS)")
  {
    constexpr auto p_seq = 13_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    for (long v = 1; v < 13; ++v)
    {
      GF elem{v};
      auto nttp_result = cbrt(elem);
      auto rt_result = cbrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value() == rt_result.has_value());
      if (nttp_result.has_value())
        REQUIRE(nttp_result->data == *rt_result);
    }
  }

  SECTION("Large prime: secp256k1 randomized (p = 1 mod 3)")
  {
    constexpr auto p_seq = 115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    Randomizer<std::uint64_t> randomizer;
    big_int<4, std::uint64_t> rand_data;

    for (int i = 0; i < 20; ++i)
    {
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p);

      GF elem{rand_data};
      auto cubed = elem * elem * elem; // guaranteed cubic residue

      auto nttp_result = cbrt(cubed);
      auto rt_result = cbrt(cubed.data, p);

      REQUIRE(nttp_result.has_value());
      REQUIRE(rt_result.has_value());
      REQUIRE(nttp_result->data == *rt_result);
    }
  }
}

TEST_CASE("Runtime sqrt: modulus as big_int value")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("Basic correctness")
  {
    auto p = to_big_int(17_Z);
    auto four = big_int<1, std::uint64_t>{4};
    auto result = sqrt(four, p);

    REQUIRE(result.has_value());
    // Verify: result^2 = 4 mod 17
    auto sq = mod(mul(*result, *result), p);
    REQUIRE(sq == four);
  }

  SECTION("Zero")
  {
    auto p = to_big_int(17_Z);
    auto zero = big_int<1, std::uint64_t>{};
    auto result = sqrt(zero, p);
    REQUIRE(result.has_value());
    REQUIRE(*result == zero);
  }

  SECTION("Non-residue returns nullopt")
  {
    auto p = to_big_int(17_Z);
    auto three = big_int<1, std::uint64_t>{3};
    REQUIRE_FALSE(sqrt(three, p).has_value());
  }

  SECTION("p = 3 mod 4 (simple formula)")
  {
    auto p = to_big_int(7_Z);
    auto two = big_int<1, std::uint64_t>{2};
    auto result = sqrt(two, p);
    REQUIRE(result.has_value());
    auto sq = mod(mul(*result, *result), p);
    REQUIRE(sq == two);
  }

  SECTION("p = 1 mod 4 (full Tonelli-Shanks)")
  {
    auto p = to_big_int(13_Z);
    auto four = big_int<1, std::uint64_t>{4};
    auto result = sqrt(four, p);
    REQUIRE(result.has_value());
    auto sq = mod(mul(*result, *result), p);
    REQUIRE(sq == four);
  }

  SECTION("Large prime: curve25519")
  {
    auto p = to_big_int(57896044618658097711785492504343953926634992332820282019728792003956564819949_Z);
    auto four = big_int<4, std::uint64_t>{4};
    auto result = sqrt(four, p);
    REQUIRE(result.has_value());
    auto sq = mod(mul(*result, *result), p);
    REQUIRE(sq == four);
  }

  SECTION("Fuzz: square-then-sqrt on secp256k1")
  {
    auto p = to_big_int(115792089237316195423570985008687907853269984665640564039457584007908834671663_Z);

    Randomizer<std::uint64_t> randomizer;
    big_int<4, std::uint64_t> rand_data;

    for (int i = 0; i < 30; ++i)
    {
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p);

      auto sq = mod(mul(rand_data, rand_data), p);
      auto root = sqrt(sq, p);
      REQUIRE(root.has_value());
      auto check = mod(mul(*root, *root), p);
      REQUIRE(check == sq);
    }
  }
}

TEST_CASE("Runtime cbrt: modulus as big_int value")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("p = 2 mod 3 (unique roots)")
  {
    auto p = to_big_int(11_Z);
    auto eight = big_int<1, std::uint64_t>{8};
    auto result = cbrt(eight, p);
    REQUIRE(result.has_value());
    auto cubed = mod(mul(mod(mul(*result, *result), p), *result), p);
    REQUIRE(cubed == eight);
  }

  SECTION("p = 1 mod 3 (AMM)")
  {
    auto p = to_big_int(13_Z);
    auto eight = big_int<1, std::uint64_t>{8};
    auto result = cbrt(eight, p);
    REQUIRE(result.has_value());
    auto cubed = mod(mul(mod(mul(*result, *result), p), *result), p);
    REQUIRE(cubed == eight);
  }

  SECTION("Non-cubic-residue returns nullopt (p = 1 mod 3)")
  {
    auto p = to_big_int(7_Z);
    // 3 is not a cubic residue mod 7
    // Cubic residues mod 7: {0, 1, 6} (since 1^3=1, 2^3=1, 3^3=6, 4^3=1, 5^3=6, 6^3=6)
    auto three = big_int<1, std::uint64_t>{3};
    REQUIRE_FALSE(cbrt(three, p).has_value());
  }

  SECTION("Fuzz: cube-then-cbrt on secp256k1 (p = 1 mod 3)")
  {
    auto p = to_big_int(115792089237316195423570985008687907853269984665640564039457584007908834671663_Z);

    Randomizer<std::uint64_t> randomizer;
    big_int<4, std::uint64_t> rand_data;

    for (int i = 0; i < 20; ++i)
    {
      do
      {
        randomizer(rand_data);
      } while (rand_data >= p);

      // r^3 is guaranteed cubic residue
      auto r_sq = mod(mul(rand_data, rand_data), p);
      auto r_cubed = mod(mul(r_sq, rand_data), p);

      auto root = cbrt(r_cubed, p);
      REQUIRE(root.has_value());

      // Verify root^3 == r^3
      auto root_sq = mod(mul(*root, *root), p);
      auto root_cubed = mod(mul(root_sq, *root), p);
      REQUIRE(root_cubed == r_cubed);
    }
  }
}

TEST_CASE("Runtime is_quadratic_residue")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  auto p = to_big_int(17_Z);

  SECTION("Known residues mod 17")
  {
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{1}, p));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{2}, p));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{4}, p));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{8}, p));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{9}, p));
  }

  SECTION("Known non-residues mod 17")
  {
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{3}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{5}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{6}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{7}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{10}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{11}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{12}, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{14}, p));
  }

  SECTION("Remaining residues mod 17: 13, 15, 16")
  {
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{13}, p));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{15}, p));
    REQUIRE(is_quadratic_residue(big_int<1, std::uint64_t>{16}, p));
  }
}

TEST_CASE("Runtime sqrt/cbrt edge cases")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("sqrt(1) = 1")
  {
    auto p = to_big_int(17_Z);
    auto one = big_int<1, std::uint64_t>{1};
    auto result = sqrt(one, p);
    REQUIRE(result.has_value());
    REQUIRE(*result == one);
  }

  SECTION("Values >= p are reduced first")
  {
    auto p = to_big_int(17_Z);
    auto four_plus_p = big_int<1, std::uint64_t>{4 + 17};
    REQUIRE(is_quadratic_residue(four_plus_p, p));
    REQUIRE_FALSE(is_quadratic_residue(big_int<1, std::uint64_t>{3 + 17}, p));
    auto root = sqrt(four_plus_p, p);
    REQUIRE(root.has_value());
    REQUIRE(mod(mul(*root, *root), p) == big_int<1, std::uint64_t>{4});
    REQUIRE(*sqrt(p, p) == big_int<1, std::uint64_t>{});

    auto q = to_big_int(11_Z);
    auto eight_plus_q = big_int<1, std::uint64_t>{8 + 11};
    auto croot = cbrt(eight_plus_q, q);
    REQUIRE(croot.has_value());
    REQUIRE(mod(mul(mod(mul(*croot, *croot), q), *croot), q) == big_int<1, std::uint64_t>{8});
    REQUIRE(*cbrt(q, q) == big_int<1, std::uint64_t>{});
  }

  SECTION("cbrt(0) = 0")
  {
    auto p = to_big_int(7_Z);
    auto zero = big_int<1, std::uint64_t>{};
    auto result = cbrt(zero, p);
    REQUIRE(result.has_value());
    REQUIRE(*result == zero);
  }

  SECTION("cbrt(1) = 1 for p = 2 mod 3")
  {
    auto p = to_big_int(5_Z);
    auto one = big_int<1, std::uint64_t>{1};
    auto result = cbrt(one, p);
    REQUIRE(result.has_value());
    REQUIRE(*result == one);
  }

  SECTION("cbrt(1) for p = 1 mod 3 (three cube roots of unity)")
  {
    auto p = to_big_int(7_Z);
    auto one = big_int<1, std::uint64_t>{1};
    auto result = cbrt(one, p);
    REQUIRE(result.has_value());
    // result^3 must equal 1
    auto cubed = mod(mul(mod(mul(*result, *result), p), *result), p);
    REQUIRE(cubed == one);
  }
}

TEST_CASE("Cross-validate: NTTP cbrt vs runtime cbrt (p = 2 mod 3)")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("p = 11 exhaustive")
  {
    constexpr auto p_seq = 11_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    for (long v = 1; v < 11; ++v)
    {
      GF elem{v};
      auto nttp_result = cbrt(elem);
      auto rt_result = cbrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value() == rt_result.has_value());
      if (nttp_result.has_value())
        REQUIRE(nttp_result->data == *rt_result);
    }
  }

  SECTION("p = 17 exhaustive")
  {
    constexpr auto p_seq = 17_Z;
    constexpr auto p = to_big_int(p_seq);
    using GF = decltype(Zq(p_seq));

    for (long v = 1; v < 17; ++v)
    {
      GF elem{v};
      auto nttp_result = cbrt(elem);
      auto rt_result = cbrt(big_int<1, std::uint64_t>{static_cast<std::uint64_t>(v)}, p);

      REQUIRE(nttp_result.has_value() == rt_result.has_value());
      if (nttp_result.has_value())
        REQUIRE(nttp_result->data == *rt_result);
    }
  }
}

TEST_CASE("Constexpr: runtime API sqrt (dual-use)")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("sqrt(4) mod 17 (full Tonelli-Shanks, p = 1 mod 4)")
  {
    constexpr auto p = to_big_int(17_Z);
    constexpr auto four = big_int<1, std::uint64_t>{4};
    constexpr auto result = sqrt(four, p);

    static_assert(result.has_value());
    // 2^2 = 4, 15^2 = 225 = 13*17 + 4
    static_assert(*result == to_big_int(2_Z) || *result == to_big_int(15_Z));
  }

  SECTION("sqrt(4) mod 7 (S == 1 branch, p = 3 mod 4)")
  {
    constexpr auto p = to_big_int(7_Z);
    constexpr auto four = big_int<1, std::uint64_t>{4};
    constexpr auto result = sqrt(four, p);

    static_assert(result.has_value());
    static_assert(*result == to_big_int(2_Z) || *result == to_big_int(5_Z));
  }

  SECTION("sqrt of non-residue returns nullopt at compile time")
  {
    constexpr auto p = to_big_int(17_Z);
    constexpr auto three = big_int<1, std::uint64_t>{3};
    constexpr auto result = sqrt(three, p);

    static_assert(!result.has_value());
  }

  SECTION("sqrt(0) = 0 at compile time")
  {
    constexpr auto p = to_big_int(17_Z);
    constexpr auto zero = big_int<1, std::uint64_t>{};
    constexpr auto result = sqrt(zero, p);

    static_assert(result.has_value());
    static_assert(*result == zero);
  }

  SECTION("is_quadratic_residue at compile time")
  {
    constexpr auto p = to_big_int(17_Z);
    static_assert(is_quadratic_residue(big_int<1, std::uint64_t>{4}, p));
    static_assert(!is_quadratic_residue(big_int<1, std::uint64_t>{3}, p));
  }
}

TEST_CASE("Constexpr: runtime API cbrt (dual-use)")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  SECTION("cbrt(2) mod 5 (p = 2 mod 3, unique root)")
  {
    constexpr auto p = to_big_int(5_Z);
    constexpr auto two = big_int<1, std::uint64_t>{2};
    constexpr auto result = cbrt(two, p);

    static_assert(result.has_value());
    // 3^3 = 27 = 2 mod 5
    static_assert(*result == to_big_int(3_Z));
  }

  SECTION("cbrt(6) mod 7 (p = 1 mod 3, AMM)")
  {
    constexpr auto p = to_big_int(7_Z);
    constexpr auto six = big_int<1, std::uint64_t>{6};
    constexpr auto result = cbrt(six, p);

    static_assert(result.has_value());
    // 6^3 = 216 = 6 mod 7
    static_assert(*result == to_big_int(6_Z));
  }

  SECTION("cbrt of non-residue returns nullopt at compile time")
  {
    constexpr auto p = to_big_int(7_Z);
    constexpr auto three = big_int<1, std::uint64_t>{3};
    constexpr auto result = cbrt(three, p);

    static_assert(!result.has_value());
  }

  SECTION("cbrt(0) = 0 at compile time")
  {
    constexpr auto p = to_big_int(7_Z);
    constexpr auto zero = big_int<1, std::uint64_t>{};
    constexpr auto result = cbrt(zero, p);

    static_assert(result.has_value());
    static_assert(*result == zero);
  }
}
