#include "catch.hpp"

import std;
import lam.ctbignum;

TEST_CASE("std::format support for big_int", "[print]")
{
  using namespace lam::cbn::literals;
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
    using namespace lam::cbn::literals;
    using GF_Large = decltype(lam::cbn::Zq(1000000000_Z)); 
    auto x = GF_Large(10_Z).data;
    std::print("Testing printing: {}\n", x);
}
