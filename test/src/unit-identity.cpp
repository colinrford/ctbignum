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

using namespace lam::cbn::literals;

TEST_CASE("Field Identity Elements") {
    // Curve25519 prime: 2^255 - 19
    constexpr auto p25519 = 57896044618658097711785492504343953926634992332820282019728792003956564819949_Z;
    using Field25519 = decltype(lam::cbn::Zq(p25519));

    SECTION("Identity Elements") {
        constexpr auto zero = Field25519::additive_identity();
        constexpr auto one = Field25519::multiplicative_identity();

        bool is_zero = (zero.data == lam::cbn::to_big_int(0_Z));
        bool is_one = (one.data == lam::cbn::to_big_int(1_Z));

        REQUIRE(is_zero);
        REQUIRE(is_one);

        // Verify properties
        constexpr Field25519 a{12345};
        constexpr auto a_plus_zero = a + zero;
        constexpr auto a_times_one = a * one;

        bool add_id_holds = (a_plus_zero == a);
        bool mul_id_holds = (a_times_one == a);

        REQUIRE(add_id_holds);
        REQUIRE(mul_id_holds);

        // Verify aliases
        bool is_zero_alias = (Field25519::zero() == zero);
        bool is_one_alias = (Field25519::one() == one);
        REQUIRE(is_zero_alias);
        REQUIRE(is_one_alias);
    }
}
