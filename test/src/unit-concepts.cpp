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

namespace 
{
  template<typename T>
  concept has_additive_identity_c_weak = requires {
    { T(0) } -> std::same_as<T>;
  } || requires {
    { T::zero() } -> std::same_as<T>;
  };
  template<typename T>
  concept has_multiplicative_identity_c_weak = requires {
    { T(1) } -> std::same_as<T>;
  } || requires {
    { T::one() } -> std::same_as<T>;
  };
  /*
   *  operator overloaded weak syntactical 'requirements' modeling group
   *  elements (does not check values)
   */
  template<typename G>
  concept additive_group_element_c_weak = has_additive_identity_c_weak<G> and requires(G g, G h, G k)
  {
    { -g } -> std::same_as<decltype(g)>;
    { g + h } -> std::same_as<decltype(h + g)>;
    { g - h } -> std::same_as<decltype(h - g)>;
    { (g + h) + k } -> std::same_as<decltype(g + (h + k))>;
  };

  template<typename G>
  concept multiplicative_group_element_c_weak = has_multiplicative_identity_c_weak<G> and requires(G g, G h, G k)
  {
    { g * h } -> std::same_as<decltype(h * g)>;
    { (g * h) * k } -> std::same_as<decltype(g * (h * k))>;
    { g / h } -> std::same_as<decltype(g)>;
    requires requires { { 1 / g } -> std::same_as<decltype(g)>; } 
          or requires { { G::one() / g } -> std::same_as<decltype(g)>; }
          or requires { { inv(g) } -> std::same_as<decltype(g)>; }; 
  };

  template<typename G>
  concept group_element_c_weak = additive_group_element_c_weak<G> or multiplicative_group_element_c_weak<G>;

  /*
   *  operator overloaded weak syntactical 'requirements' modeling group
   *  elements (does not check values)
   */
  template<typename R>
  concept ring_element_c_weak = additive_group_element_c_weak<R> 
                             and has_multiplicative_identity_c_weak<R>
                             and requires(R r, R s, R t)
  {
    { r * s } -> std::same_as<decltype(s * r)>;
    { (r * s) * t } -> std::same_as<decltype(r * (s * t))>;
    { r * (s + t) } -> std::same_as<decltype(r * s + r * t)>;
    { (r + s) * t } -> std::same_as<decltype(r * t + s * t)>;
  };

  /*
   *  operator overloaded weak syntactical 'requirements' modeling field
   *  elements (does not check values)
   */
  template<typename K>
  concept field_element_c_weak = ring_element_c_weak<K> 
                              and multiplicative_group_element_c_weak<K>;

} // namespace

TEST_CASE("Concept checks for Finite Field elements")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  auto modulus = 1267650600228229401496703205653_Z;
  using GF = decltype(Zq(modulus));
  
  // Static assertions to verify concepts
  static_assert(has_additive_identity_c_weak<GF>, "GF must have additive identity");
  static_assert(has_multiplicative_identity_c_weak<GF>, "GF must have multiplicative identity");
  static_assert(additive_group_element_c_weak<GF>, "GF must be an additive group element");
  static_assert(multiplicative_group_element_c_weak<GF>, "GF must be a multiplicative group element");
  static_assert(ring_element_c_weak<GF>, "GF must be a ring element");
  static_assert(field_element_c_weak<GF>, "GF must be a field element");
  
  // Runtime check just to ensure the test runs and passes if static_asserts pass
  REQUIRE(field_element_c_weak<GF>);
}

TEST_CASE("Zq elements in std::array")
{
  using namespace lam::cbn;
  using namespace lam::cbn::literals;

  auto modulus = 1267650600228229401496703205653_Z;
  using GF = decltype(Zq(modulus));

  // Verify can be used in std::array
  std::array<GF, 3> arr = { GF(1_Z), GF(2_Z), GF(3_Z) };
  
  REQUIRE(arr[0] == GF(1_Z));
  REQUIRE(arr[1] == GF(2_Z));
  REQUIRE(arr[2] == GF(3_Z));

  // Verify compilation of default construction in array
  std::array<GF, 5> arr2;
  arr2[0] = GF(10_Z);
  REQUIRE(arr2[0] == GF(10_Z));
}

TEST_CASE("ZqElement data access")
{
  using namespace lam::cbn::literals;
  using GF = decltype(lam::cbn::Zq(100_Z)); // dummy modulus
  GF z{42_Z};

  // 1. Access via .data
  REQUIRE(z.data == GF(42_Z).data);

  // 2. Access via explicit cast (due to explicit operator auto())
  // method 1: static_cast
  auto b1 = static_cast<decltype(z.data)>(z);
  REQUIRE(b1 == GF(42_Z).data);

  // method 2: functional cast
  auto b2 = decltype(z.data)(z);
  REQUIRE(b2 == GF(42_Z).data);
}
