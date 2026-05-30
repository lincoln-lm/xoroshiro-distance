"""Generate constants for c++ code + lookup json for hash generation"""

import json
import galois

# characteristic polynomial of the xoroshiro state transition matrix
char_poly_int = 0x10008828E513B43D5095B8F76579AA001
char_poly = galois.Poly(
    [(char_poly_int >> exp) & 1 for exp in range(128, -1, -1)], field=galois.GF(2)
)
x = galois.Poly([1, 0], field=galois.GF(2))

# barrett reduction constant floor(x^(2 * 128) / char_poly)
mu, _ = divmod(x**256, char_poly)
mu_coeffs = list(map(int, reversed(mu.coeffs)))
mu_int = 0
for bit_i, bit in enumerate(mu_coeffs):
    if bit == 1:
        mu_int |= 1 << bit_i
print("// gf2_128_mod")
# cut off 128th bit
print(
    f"const Poly128 mu = _mm_set_epi64x(0x{(mu_int >> 64) & 0xFFFFFFFFFFFFFFFF:016X}, 0x{mu_int & 0xFFFFFFFFFFFFFFFF:016X});"
)
print(
    f"const Poly128 p = _mm_set_epi64x(0x{(char_poly_int >> 64) & 0xFFFFFFFFFFFFFFFF:016X}, 0x{char_poly_int & 0xFFFFFFFFFFFFFFFF:016X});"
)

period = 2**128 - 1
# 3, 5, 17, 257, 641, 65537, 274177, 6700417, 67280421310721
primes = (3, 5, 17, 257, 641)
product = 3 * 5 * 17 * 257 * 641

print("// discrete_logarithm")
multipliers = []
cofactor_exponents = []
for mod in primes:
    p = product // mod
    multiplier = pow(p, -1, mod) * p
    multipliers.append(multiplier)
    cofactor_exponent = period // mod
    assert (cofactor_exponent >> 64) == (cofactor_exponent & 0xFFFFFFFFFFFFFFFF)
    cofactor_exponents.append(cofactor_exponent >> 64)

print(f"constexpr uint64_t multipliers[] = {{{', '.join(map(str, multipliers))}}};")
print(
    f"constexpr uint64_t cofactor_exponents[] = {{{', '.join((f"0x{exp:016X}" for exp in cofactor_exponents))}}};"
)

lookup = {}
for p in primes:
    stride = period // p
    lookup_p = {}

    for i in range(p):
        generator = list(map(int, reversed(pow(x, i * stride, char_poly).coeffs)))
        generator_int = 0
        for bit_i, bit in enumerate(generator):
            if bit == 1:
                generator_int |= 1 << bit_i
        lookup_p[generator_int] = i
    lookup[p] = lookup_p

with open("lookup.json", "w", encoding="utf-8") as f:
    json.dump(lookup, f)
