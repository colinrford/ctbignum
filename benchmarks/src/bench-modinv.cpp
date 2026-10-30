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

#include <gmp.h>
#include <benchmark/benchmark.h>

import std;
import lam.ctbignum;

using namespace lam::cbn::literals;

template<std::size_t N, typename T>
void modinv_gmp(benchmark::State& state, lam::cbn::big_int<N, T> const& prime)
{

  constexpr std::size_t n = N;
  std::size_t k = static_cast<std::size_t>(mpn_sec_invert_itch(n));
  std::vector<mp_limb_t> scratch(k, 0);

  std::array<mp_limb_t, n> modulus{};
  memcpy(modulus.data(), prime.data(), N * sizeof(T));

  std::size_t total_sz = n * 1000;
  std::vector<mp_limb_t> data(total_sz);
  std::default_random_engine generator;
  std::uniform_int_distribution<mp_limb_t> distribution(0);
  for (auto& limb : data)
    limb = distribution(generator);

  mp_limb_t* base_ptr = data.data();

  mp_limb_t dummy_quotient[n];
  for (std::size_t i = 0; i < 1000; ++i)
  {
    mpn_tdiv_qr(dummy_quotient, base_ptr + i * n, 0, base_ptr + i * n, n, modulus.data(), n);
    // modular reduction
  }

  std::size_t i = 0;
  mp_limb_t result[n];

  for (auto _ : state)
  {

    mpn_sec_invert(result, base_ptr + i, modulus.data(), n, 2 * n * GMP_NUMB_BITS, scratch.data());
    benchmark::DoNotOptimize(result);

    i += n;
    if (i == total_sz)
      i = 0;
  }
}

template<std::size_t N, typename T>
void modinv_cbn(benchmark::State& state, lam::cbn::big_int<N, T> const& prime)
{

  constexpr std::size_t n = N;

  std::array<mp_limb_t, n> modulus{};
  memcpy(modulus.data(), prime.data(), N * sizeof(T));

  std::size_t total_sz = n * 1000;
  std::vector<mp_limb_t> data(total_sz);
  std::default_random_engine generator;
  std::uniform_int_distribution<mp_limb_t> distribution(0);
  for (auto& limb : data)
    limb = distribution(generator);

  mp_limb_t* base_ptr = data.data();

  mp_limb_t dummy_quotient[n];
  for (std::size_t i = 0; i < 1000; ++i)
  {
    mpn_tdiv_qr(dummy_quotient, base_ptr + i * n, 0, base_ptr + i * n, n, modulus.data(), n);
    // modular reduction
  }

  std::size_t i = 0;

  for (auto _ : state)
  {

    auto result = lam::cbn::mod_inv(*reinterpret_cast<lam::cbn::big_int<N, T>*>(base_ptr + i), prime);
    benchmark::DoNotOptimize(result);

    i += n;
    if (i == total_sz)
      i = 0;
  }
}

// Registers a benchmark named "BM_takes_args/int_string_test" that passes
// the specified values to `extra_args`.
BENCHMARK_CAPTURE(
  modinv_gmp, gmp_modular_inverse,
  lam::cbn::to_big_int(115792089237316195423570985008687907853269984665640564039457584007908834671663_Z));
BENCHMARK_CAPTURE(
  modinv_cbn, cbn_modular_inverse,
  lam::cbn::to_big_int(115792089237316195423570985008687907853269984665640564039457584007908834671663_Z));

BENCHMARK_MAIN();
