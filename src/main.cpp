#include "hash.hpp"
#include <algorithm>
#include <bit>
#include <cassert>
#include <chrono>
#include <cstdint>
#include <immintrin.h>
#include <iostream>
#include <random>

typedef __m256i Row;
typedef __m256i Poly256;
typedef __m128i Poly128;

inline void xoroshiro_next(uint64_t state[2]) {
  uint64_t seed_1 = state[1] ^ state[0];
  const uint64_t seed_0 = std::rotl(state[0], 24) ^ seed_1 ^ (seed_1 << 16ULL);
  seed_1 = std::rotl(seed_1, 37);
  state[0] = seed_0;
  state[1] = seed_1;
}

inline bool check_bit(Row row, int n) {
  uint64_t w = 0;
  if ((n >> 6) == 0) {
    w = _mm256_extract_epi64(row, 0);
  } else {
    w = _mm256_extract_epi64(row, 1);
  }

  return (w >> (n & 63)) & 1;
}

inline void solve_system(Row matrix[129], uint64_t out[2]) {
  for (int j = 0; j < 128; j++) {
    int pivot = -1;
    for (int i = j; i < 128; i++) {
      if (check_bit(matrix[i], j)) {
        pivot = i;
        break;
      }
    }
    if (pivot == -1)
      continue;
    if (pivot != j) {
      std::swap(matrix[j], matrix[pivot]);
    }

    for (int i = 0; i < 129; i++) {
      if (i != j && check_bit(matrix[i], j)) {
        matrix[i] = _mm256_xor_si256(matrix[i], matrix[j]);
      }
    }
  }

  out[0] = _mm256_extract_epi64(matrix[128], 2);
  out[1] = _mm256_extract_epi64(matrix[128], 3);
}

inline void compute_jump_polynomial(const uint64_t state_a[2],
                                    const uint64_t state_b[2],
                                    uint64_t jump_polynomial[2]) {
  Row matrix[129] = {0};

  uint64_t state[2] = {state_a[0], state_a[1]};
  for (int i = 0; i < 128; ++i) {
    uint64_t identity_0 = 0;
    uint64_t identity_1 = 0;
    if (i < 64) {
      identity_0 = 1ULL << i;
    } else {
      identity_1 = 1ULL << (i - 64);
    }
    matrix[i] = _mm256_set_epi64x(identity_1, identity_0, state[1], state[0]);
    xoroshiro_next(state);
  }
  matrix[128] = _mm256_set_epi64x(0, 0, state_b[1], state_b[0]);

  solve_system(matrix, jump_polynomial);
}

inline Poly256 gf2_128_multiply(Poly128 a, Poly128 b) {
  // karatsuba multiplication
  // a_0 * b_0
  Poly128 a_0b_0 = _mm_clmulepi64_si128(a, b, 0x00);
  // a_1 * b_1
  Poly128 a_1b_1 = _mm_clmulepi64_si128(a, b, 0x11);

  // a_0 + a_1
  Poly128 a_0_p_a_1 = _mm_xor_si128(a, _mm_srli_si128(a, 8));
  // b_0 + b_1
  Poly128 b_0_p_b_1 = _mm_xor_si128(b, _mm_srli_si128(b, 8));

  // (a_0 + a_1) * (b_0 + b_1) = a_1 * b_0 + a_0 * b_1 + a_0 * b_0 + a_1 * b_1
  Poly128 middle_product = _mm_clmulepi64_si128(a_0_p_a_1, b_0_p_b_1, 0x00);

  // middle_product - a_0 * b_0 = a_1 * b_0 + a_0 * b_1 + a_1 * b_1
  middle_product = _mm_xor_si128(middle_product, a_0b_0);
  // - a_1 * b_1 = a_1 * b_0 + a_0 * b_1
  middle_product = _mm_xor_si128(middle_product, a_1b_1);

  Poly128 middle_product_low = _mm_slli_si128(middle_product, 8);
  Poly128 middle_product_high = _mm_srli_si128(middle_product, 8);

  Poly128 result_low = _mm_xor_si128(a_0b_0, middle_product_low);
  Poly128 result_high = _mm_xor_si128(a_1b_1, middle_product_high);

  Poly256 product = _mm256_set_m128i(result_high, result_low);

  return product;
}

inline Poly256 gf2_128_square(Poly128 current) {
  // (A + B)^2 = A^2 + B^2
  // C = A + Bx^64 => C^2 = A^2 + B^2 x^128
  Poly128 a_squared = _mm_clmulepi64_si128(current, current, 0x00);
  Poly128 b_squared = _mm_clmulepi64_si128(current, current, 0x11);
  return _mm256_setr_m128i(a_squared, b_squared);
}

inline Poly128 gf2_128_mod(Poly256 a) {
  // barrett reduction
  const Poly128 mu = _mm_set_epi64x(0x000882CE131D9096, 0x72F817C07CFCAEAE);
  const Poly128 p = _mm_set_epi64x(0x0008828E513B43D5, 0x095B8F76579AA001);
  Poly128 q1 = _mm256_extracti128_si256(a, 1);
  Poly256 q2 = gf2_128_multiply(q1, mu);
  Poly128 q3 = _mm256_extracti128_si256(q2, 1);
  q3 = _mm_xor_si128(q3, q1);
  Poly256 q4 = gf2_128_multiply(q3, p);
  Poly128 r = _mm_xor_si128(_mm256_extracti128_si256(a, 0),
                            _mm256_extracti128_si256(q4, 0));
  return r;
}

inline Poly128 gf2_128_exponentiate(Poly128 base, uint64_t exp_low,
                                    uint64_t exp_high) {
  Poly128 result = _mm_set_epi64x(0, 1);
  Poly128 current = base;

  for (int i = 0; i < 128; i++) {
    uint64_t exp_i;
    if (i < 64) {
      exp_i = (exp_low >> i) & 1;
    } else {
      exp_i = (exp_high >> (i - 64)) & 1;
    }
    if (exp_i) {
      Poly256 product = gf2_128_multiply(result, current);
      result = gf2_128_mod(product);
    }
    Poly256 square = gf2_128_square(current);
    current = gf2_128_mod(square);
  }

  return result;
}

inline uint64_t discrete_logarithm(Poly128 element) {
  uint64_t distance = 0;
  constexpr uint64_t product = 3 * 5 * 17 * 257 * 641;
  // p = product / mod
  // multiplier = pow(p, -1, mod) * p
  constexpr uint64_t multipliers[] = {14002645, 25204761, 14826330, 14057130,
                                      15925005};
  // (2**128 - 1) / mod
  // upper and lower bits of the exponent are the same so only one is stored
  constexpr uint64_t cofactor_exponents[] = {
      0x5555555555555555, 0x3333333333333333, 0x0F0F0F0F0F0F0F0F,
      0x00FF00FF00FF00FF, 0x00663D80FF99C27F};

  Poly128 subgroup_generator;

  // p=3
  subgroup_generator = gf2_128_exponentiate(element, cofactor_exponents[0],
                                            cofactor_exponents[0]);
  distance += lookup_p3(_mm_extract_epi64(subgroup_generator, 0),
                        _mm_extract_epi64(subgroup_generator, 1)) *
              multipliers[0];
  // p=5
  subgroup_generator = gf2_128_exponentiate(element, cofactor_exponents[1],
                                            cofactor_exponents[1]);
  distance += lookup_p5(_mm_extract_epi64(subgroup_generator, 0),
                        _mm_extract_epi64(subgroup_generator, 1)) *
              multipliers[1];
  // p=17
  subgroup_generator = gf2_128_exponentiate(element, cofactor_exponents[2],
                                            cofactor_exponents[2]);
  distance += lookup_p17(_mm_extract_epi64(subgroup_generator, 0),
                         _mm_extract_epi64(subgroup_generator, 1)) *
              multipliers[2];
  // p=257
  subgroup_generator = gf2_128_exponentiate(element, cofactor_exponents[3],
                                            cofactor_exponents[3]);
  distance += lookup_p257(_mm_extract_epi64(subgroup_generator, 0),
                          _mm_extract_epi64(subgroup_generator, 1)) *
              multipliers[3];
  // p=641
  subgroup_generator = gf2_128_exponentiate(element, cofactor_exponents[4],
                                            cofactor_exponents[4]);
  distance += lookup_p641(_mm_extract_epi64(subgroup_generator, 0),
                          _mm_extract_epi64(subgroup_generator, 1)) *
              multipliers[4];
  return distance % product;
}

uint64_t compute_distance(const uint64_t state_a[2],
                          const uint64_t state_b[2]) {
  uint64_t jump_polynomial_int[2] = {0};
  compute_jump_polynomial(state_a, state_b, jump_polynomial_int);

  Poly128 jump_poly =
      _mm_set_epi64x(jump_polynomial_int[1], jump_polynomial_int[0]);

  return discrete_logarithm(jump_poly);
}

uint64_t bruteforce_distance(const uint64_t state_a[2],
                             const uint64_t state_b[2]) {
  uint64_t state[2] = {state_a[0], state_a[1]};
  for (uint64_t distance = 0;; distance++) {
    if (state[0] == state_b[0] && state[1] == state_b[1]) {
      return distance;
    }
    xoroshiro_next(state);
  }
}

uint64_t random_advance() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint64_t> distr(1, 3 * 5 * 17 * 257 * 641 - 1);
  return distr(gen);
}

std::tuple<uint64_t, uint64_t> random_state() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<uint64_t> distr(
      std::numeric_limits<std::uint64_t>::min(),
      std::numeric_limits<std::uint64_t>::max());
  return {distr(gen), distr(gen)};
}

template <bool bruteforce = false>
std::tuple<std::chrono::duration<double>, uint64_t>
benchmark(uint64_t test_distance, uint64_t num_trials) {
  auto [state_a_0, state_a_1] = random_state();
  uint64_t state_a[2] = {state_a_0, state_a_1};
  uint64_t state_b[2] = {state_a_0, state_a_1};

  for (uint64_t i = 0; i < test_distance; i++) {
    xoroshiro_next(state_b);
  }

  auto start = std::chrono::steady_clock::now();
  uint64_t result_distance;
  for (uint64_t i = 0; i < num_trials; i++) {
    if constexpr (bruteforce) {
      result_distance = bruteforce_distance(state_a, state_b);
    } else {
      result_distance = compute_distance(state_a, state_b);
    }
  }
  auto end = std::chrono::steady_clock::now();
  assert(result_distance != 0);
  return {end - start, test_distance};
}

int main() {
  for (int i = 0; i < 1000; i++) {
    constexpr int num_trials = 1000;
    uint64_t test_distance = random_advance();
    auto [optimized_time, optimized_distance] =
        benchmark(test_distance, num_trials);
    auto [bruteforce_time, bruteforce_distance] =
        benchmark<true>(test_distance, num_trials);
    std::cout << test_distance << "\t" << optimized_time.count() << "\t"
              << bruteforce_time.count() << "\n";
  }
  return 0;
}