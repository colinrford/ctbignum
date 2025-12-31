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

export module lam.ctbignum:pow;

import std;

namespace lam::cbn
{

export 
template <std::size_t N1, typename T> 
constexpr auto pow(big_int<N1, T> base, T exp)
{
  big_int<N1, T> result{1};

  if (exp == 0)
    return result;

  while (true)
  {
    auto lsb = exp & 1;
    exp >>= 1;
    if (lsb)
    {
      result = partial_mul<N1>(base, result);
      if (exp == 0)
        break;
    }
    base = partial_mul<N1>(base, base);
  }

  return result;
}
} // namespace lam::cbn
