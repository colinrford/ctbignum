//
// Scaling benchmark for ctbignum
//

#include <NTL/ZZ.h>
#include <benchmark/benchmark.h>

import std;
import lam.ctbignum;

using namespace lam::cbn;

// GMP Comparison
template<size_t Len>
static void mul_gmp(benchmark::State& state)
{
  size_t total_sz = 2 * Len * 1000;
  std::vector<uint64_t> data(total_sz);
  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution(0);
  for (auto& limb : data)
    limb = distribution(generator);

  size_t i = 0;
  // Assumption: mp_limb_t is uint64_t compatible on this platform (macOS x86_64)
  auto base_ptr = reinterpret_cast<mp_limb_t*>(data.data());
  mp_limb_t result[2 * Len];

  for (auto _ : state)
  {
    mpn_mul(result, base_ptr + i, Len, base_ptr + i + Len, Len);
    benchmark::DoNotOptimize(result);

    i += 2 * Len;
    if (i == total_sz)
      i = 0;
  }
}

// ctbignum
template<size_t Len>
static void mul_cbn(benchmark::State& state)
{
  size_t total_sz = 2 * Len * 1000;
  std::vector<uint64_t> data(total_sz);
  std::default_random_engine generator;
  std::uniform_int_distribution<uint64_t> distribution(0);
  for (auto& limb : data)
    limb = distribution(generator);

  size_t i = 0;
  auto base_ptr = data.data();

  for (auto _ : state)
  {
    auto x = reinterpret_cast<big_int<Len>*>(base_ptr + i);
    auto y = reinterpret_cast<big_int<Len>*>(base_ptr + i + Len);
    auto j = mul(*x, *y);
    benchmark::DoNotOptimize(j);

    i += 2 * Len;
    if (i == total_sz)
      i = 0;
  }
}

// 4 limbs = 256 bits
BENCHMARK_TEMPLATE(mul_cbn, 4);
BENCHMARK_TEMPLATE(mul_gmp, 4);

// 8 limbs = 512 bits
BENCHMARK_TEMPLATE(mul_cbn, 8);
BENCHMARK_TEMPLATE(mul_gmp, 8);

// 16 limbs = 1024 bits
BENCHMARK_TEMPLATE(mul_cbn, 16);
BENCHMARK_TEMPLATE(mul_gmp, 16);

// 32 limbs = 2048 bits
BENCHMARK_TEMPLATE(mul_cbn, 32);
BENCHMARK_TEMPLATE(mul_gmp, 32);

// 64 limbs = 4096 bits
BENCHMARK_TEMPLATE(mul_cbn, 64);
BENCHMARK_TEMPLATE(mul_gmp, 64);

BENCHMARK_MAIN();
