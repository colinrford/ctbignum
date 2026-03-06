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

export module lam.ctbignum:type_traits;

import std;
import :config;

namespace lam::cbn
{

export template<typename T>
struct dbl_bitlen
{
  using type = void;
};
export template<>
struct dbl_bitlen<std::uint8_t>
{
  using type = std::uint16_t;
};
export template<>
struct dbl_bitlen<std::uint16_t>
{
  using type = std::uint32_t;
};
export template<>
struct dbl_bitlen<std::uint32_t>
{
  using type = std::uint64_t;
};
export template<>
struct dbl_bitlen<std::uint64_t>
{
  using type = lam::ctbignum::config::uint128_or_void;
};

} // namespace lam::cbn
