#include "catch.hpp"

import std;
import lam.ctbignum;


TEST_CASE("std::format support for big_int", "[print]")
{
  using namespace lam::cbn::literals;
  // Construct big_int via field element to handle literal conversion
  // using a large enough dummy modulus
  using GF_Large = decltype(lam::cbn::Zq(1000000000_Z)); 
  auto x = GF_Large(123456789_Z).data;
  
  // Test basic formatting
  auto s = std::format("{}", x);
  REQUIRE(s == "123456789");
  
  // Test inside a sentence
  auto msg = std::format("Value is {}", x);
  REQUIRE(msg == "Value is 123456789");
}

TEST_CASE("std::format support for ZqElement", "[print]")
{
  using namespace lam::cbn::literals;
  using GF = decltype(lam::cbn::Zq(100_Z)); // dummy modulus
  GF z{42_Z};

  // Test ZqElement formatting
  auto s = std::format("{}", z);
  REQUIRE(s == "42");
  
  // Test with surrounding text
  auto msg = std::format("Element: {}", z);
  REQUIRE(msg == "Element: 42");
}

TEST_CASE("std::print compilation check", "[print]")
{
    // This is mostly compile-time check to ensure std::print accepts our types
    // We can't easily capture stdout in this environment without extra work, 
    // so we trust std::format tests for correctness and just ensure this compiles.
    using namespace lam::cbn::literals;
    using GF_Large = decltype(lam::cbn::Zq(1000000000_Z)); 
    auto x = GF_Large(10_Z).data;
    std::print("Testing printing: {}\n", x);
}
