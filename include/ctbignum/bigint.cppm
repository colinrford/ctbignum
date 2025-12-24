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

export module lam.ctbignum:bigint;

import std;

namespace lam::cbn
{

export template <std::size_t N, std::unsigned_integral T = std::uint64_t> struct big_int : std::array<T, N>
{
};

} // namespace lam::cbn
