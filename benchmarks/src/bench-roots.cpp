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

#include <benchmark/benchmark.h>

import std;
import lam.ctbignum;

using namespace lam::cbn;
using namespace lam::cbn::literals;

template<std::size_t N>
static big_int<N, std::uint64_t> random_element(big_int<N, std::uint64_t> p)
{
  static std::mt19937_64 rng{42}; // fixed so runs are comparable
  std::uniform_int_distribution<std::uint64_t> dist;
  // masking the top limb to p's bit length keeps the draw uniform and accepts
  // more than half the time; unmasked, p = 65537 would almost never accept
  const auto top_mask = ~std::uint64_t{0} >> std::countl_zero(p[N - 1]);

  big_int<N, std::uint64_t> r;
  do
  {
    for (auto& limb : r)
      limb = dist(rng);
    r[N - 1] &= top_mask;
  } while (r >= p);
  return r;
}

// Random squares (quadratic residues)
template<std::size_t Limbs, typename Modulus>
static auto make_squares_nttp(Modulus, std::size_t count)
{
  using GF = decltype(Zq(Modulus{}));
  constexpr auto p = to_big_int(Modulus{});
  std::vector<GF> out;
  out.reserve(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    GF r{random_element<Limbs>(p)};
    out.push_back(r * r);
  }
  return out;
}

template<std::size_t Limbs>
static auto make_squares_rt(big_int<Limbs, std::uint64_t> p, std::size_t count)
{
  std::vector<big_int<Limbs, std::uint64_t>> out;
  out.reserve(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    auto r = random_element<Limbs>(p);
    out.push_back(mod(mul(r, r), p));
  }
  return out;
}

// Random cubes (cubic residues)
template<std::size_t Limbs, typename Modulus>
static auto make_cubes_nttp(Modulus, std::size_t count)
{
  using GF = decltype(Zq(Modulus{}));
  constexpr auto p = to_big_int(Modulus{});
  std::vector<GF> out;
  out.reserve(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    GF r{random_element<Limbs>(p)};
    out.push_back(r * r * r);
  }
  return out;
}

template<std::size_t Limbs>
static auto make_cubes_rt(big_int<Limbs, std::uint64_t> p, std::size_t count)
{
  std::vector<big_int<Limbs, std::uint64_t>> out;
  out.reserve(count);
  for (std::size_t i = 0; i < count; ++i)
  {
    auto r = random_element<Limbs>(p);
    auto sq = mod(mul(r, r), p);
    out.push_back(mod(mul(sq, r), p));
  }
  return out;
}

// secp256k1: p ≡ 3 mod 4 (S == 1 sqrt), p ≡ 1 mod 3 (AMM cbrt)

constexpr auto secp256k1_seq =
  115792089237316195423570985008687907853269984665640564039457584007908834671663_Z;
constexpr auto secp256k1_bi = to_big_int(secp256k1_seq);

static void BM_sqrt_nttp_secp256k1(benchmark::State& state)
{
  auto inputs = make_squares_nttp<4>(secp256k1_seq, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = sqrt(inputs[i & 63]);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

static void BM_sqrt_runtime_secp256k1(benchmark::State& state)
{
  auto inputs = make_squares_rt<4>(secp256k1_bi, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = sqrt(inputs[i & 63], secp256k1_bi);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

static void BM_cbrt_nttp_secp256k1(benchmark::State& state)
{
  auto inputs = make_cubes_nttp<4>(secp256k1_seq, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = cbrt(inputs[i & 63]);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

static void BM_cbrt_runtime_secp256k1(benchmark::State& state)
{
  auto inputs = make_cubes_rt<4>(secp256k1_bi, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = cbrt(inputs[i & 63], secp256k1_bi);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

// curve25519: p ≡ 1 mod 4 (S = 2), so sqrt runs the Tonelli-Shanks loop

constexpr auto curve25519_seq =
  57896044618658097711785492504343953926634992332820282019728792003956564819949_Z;
constexpr auto curve25519_bi = to_big_int(curve25519_seq);

static void BM_sqrt_nttp_curve25519(benchmark::State& state)
{
  auto inputs = make_squares_nttp<4>(curve25519_seq, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = sqrt(inputs[i & 63]);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

static void BM_sqrt_runtime_curve25519(benchmark::State& state)
{
  auto inputs = make_squares_rt<4>(curve25519_bi, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = sqrt(inputs[i & 63], curve25519_bi);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

// 65537 = 2^16 + 1: S = 16, the deepest Tonelli-Shanks loop here

constexpr auto fermat_seq = 65537_Z;
constexpr auto fermat_bi = to_big_int(fermat_seq);

static void BM_sqrt_nttp_fermat(benchmark::State& state)
{
  auto inputs = make_squares_nttp<1>(fermat_seq, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = sqrt(inputs[i & 63]);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

static void BM_sqrt_runtime_fermat(benchmark::State& state)
{
  auto inputs = make_squares_rt<1>(fermat_bi, 64);
  std::size_t i = 0;
  for (auto _ : state)
  {
    auto result = sqrt(inputs[i & 63], fermat_bi);
    benchmark::DoNotOptimize(result);
    ++i;
  }
}

BENCHMARK(BM_sqrt_nttp_secp256k1);
BENCHMARK(BM_sqrt_runtime_secp256k1);
BENCHMARK(BM_cbrt_nttp_secp256k1);
BENCHMARK(BM_cbrt_runtime_secp256k1);
BENCHMARK(BM_sqrt_nttp_curve25519);
BENCHMARK(BM_sqrt_runtime_curve25519);
BENCHMARK(BM_sqrt_nttp_fermat);
BENCHMARK(BM_sqrt_runtime_fermat);

BENCHMARK_MAIN();
