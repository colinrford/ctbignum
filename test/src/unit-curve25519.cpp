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

TEST_CASE("Curve25519 Field Operations") {
    // Curve25519 prime: 2^255 - 19
    // Decimal: 57896044618658097711785492504343953926634992332820282019728792003956564819949
    constexpr auto p25519 = 57896044618658097711785492504343953926634992332820282019728792003956564819949_Z;

    using Field25519 = decltype(lam::cbn::Zq(p25519));

    SECTION("Basic Multiplication") {
        constexpr Field25519 a{12345};
        constexpr Field25519 b{67890};
        constexpr auto c = a * b;
        
        // Expected result: 12345 * 67890 = 838102050
        constexpr auto expected = 838102050_Z;
        
        bool is_equal = (c.data == lam::cbn::to_big_int(expected));
        REQUIRE(is_equal);
    }

    SECTION("Modular Reduction") {
        // p = 2^255 - 19.
        // Let's try x = p + 1. It should reduce to 1.
        // 57896044618658097711785492504343953926634992332820282019728792003956564819950 = p + 1
        constexpr auto p_plus_1 = 57896044618658097711785492504343953926634992332820282019728792003956564819950_Z;
        constexpr Field25519 reduced{p_plus_1};
        
        bool is_p_plus_1 = (reduced.data == lam::cbn::to_big_int(1_Z));
        REQUIRE(is_p_plus_1);

    }
}
