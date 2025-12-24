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

#include <NTL/ZZ.h>
#include <NTL/ZZ_p.h>
#include <benchmark/benchmark.h>
#include <random>

import std;
import lam.ctbignum;

using namespace lam::cbn::literals;

static void modexp_ntl(benchmark::State &state)
{
  using NTL::ZZ;
  using NTL::ZZ_p;

  // Use RandomBits_ZZ instead of conv to avoid libc++/libstdc++ ABI issues
  auto modulus = NTL::RandomBits_ZZ(200);
  if (IsZero(modulus))
    modulus = 1 + NTL::RandomBits_ZZ(200);
  ZZ_p::init(modulus);

  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution(0);

  ZZ_p z;
  ZZ exp = NTL::RandomBits_ZZ(100);

  ZZ_p a = NTL::random_ZZ_p();

  for (auto _ : state)
  {
    NTL::power(z, a, exp);
    benchmark::DoNotOptimize(z);
  }
}

template <size_t Len> static void modexp_cbn(benchmark::State &state)
{

  using namespace lam::cbn;

  size_t total_sz = 2 * Len * 1000;

  std::vector<uint64_t> data(total_sz);
  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution(0);
  for (auto &limb : data)
    limb = distribution(generator);

  size_t i = 0;
  auto base_ptr = data.data();

  for (auto _ : state)
  {

    auto x = reinterpret_cast<big_int<Len> *>(base_ptr + i);
    auto y = reinterpret_cast<big_int<Len> *>(base_ptr + i + Len);
    auto j = lam::cbn::mod_exp(*x, *y, 14474011154664524427946373126085988481658748083205070504932198000989141205031_Z);
    benchmark::DoNotOptimize(j);

    i += 2 * Len;
    if (i == total_sz)
      i = 0;
  }
}

BENCHMARK(modexp_ntl);
BENCHMARK_TEMPLATE(modexp_cbn, 4);

BENCHMARK_MAIN();
