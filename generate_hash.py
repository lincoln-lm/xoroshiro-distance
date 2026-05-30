"""Build a minimal perfect hash function/table for each prime table"""

import json
import random

random.seed(0)


def mix64(x, seed):
    x = (x ^ seed) ^ 0x9E3779B97F4A7C15
    x ^= x >> 30
    x *= 0xBF58476D1CE4E5B9
    x &= 0xFFFFFFFFFFFFFFFF
    x ^= x >> 27
    return x


def build_mphf(data):
    n = len(data)
    keys = [hi ^ lo for hi, lo, _ in data]

    num_buckets = max(1, n // 4)
    base_seed = 0
    seed_limit = 1 << 8
    while True:
        seed_limit_reached = False
        buckets = [[] for _ in range(num_buckets)]

        for i, k in enumerate(keys):
            bucket_index = mix64(k, base_seed) % num_buckets
            buckets[bucket_index].append(i)

        buckets = sorted(enumerate(buckets), key=lambda x: len(x[1]), reverse=True)
        bucket_seeds = [0] * num_buckets
        table = [-1] * n
        for bucket_index, bucket in buckets:
            seed = 0
            while True:
                seen = []
                for i in bucket:
                    key = keys[i]
                    idx = mix64(key, seed) % n
                    if idx in seen or table[idx] != -1:
                        break
                    seen.append(idx)
                if len(seen) == len(bucket):
                    bucket_seeds[bucket_index] = seed
                    for j, idx in enumerate(seen):
                        table[idx] = bucket[j]
                    break
                seed += 1
                if seed >= seed_limit:
                    seed_limit_reached = True
                    break
            if seed_limit_reached:
                break

        if not seed_limit_reached:
            return num_buckets, base_seed, bucket_seeds, table, data
        base_seed += 1
        if base_seed >= 4096 * 3:
            base_seed = 0
            seed_limit <<= 8


def build_mphf_c(data, name):
    num_buckets, base_seed, bucket_seeds, table, data = build_mphf(data)
    seed_type = "uint8_t" if max(bucket_seeds) < 256 else "uint16_t"
    data_type = "uint8_t" if max(data, key=lambda x: x[2])[2] < 256 else "uint16_t"
    hash_table = [str(data[i][2]) for i in table]
    for hi, lo, value in data:
        seed = bucket_seeds[mix64(hi ^ lo, base_seed) % num_buckets]
        hash_ = mix64(hi ^ lo, seed) % len(table)
        assert data[table[hash_]][2] == value
    print(f"""
inline {seed_type} bucket_seeds_{name}[] = {{{', '.join(map(str, bucket_seeds))}}};
inline {data_type} hash_table_{name}[] = {{{', '.join(hash_table)}}};
inline {data_type} lookup_{name}(uint64_t key_a, uint64_t key_b) {{
    uint64_t key = key_a ^ key_b;
    uint64_t x = (key ^ {base_seed}) ^ 0x9E3779B97F4A7C15;
    x ^= x >> 30;
    x *= 0xBF58476D1CE4E5B9;
    x ^= x >> 27;
    int bucket_index = x % {num_buckets};
    x = (key ^ bucket_seeds_{name}[bucket_index]) ^ 0x9E3779B97F4A7C15;
    x ^= x >> 30;
    x *= 0xBF58476D1CE4E5B9;
    x ^= x >> 27;
    return hash_table_{name}[x % {len(table)}];
}}""")


with open("lookup.json", "r", encoding="utf-8") as f:
    lookup = json.load(f)

print("#include <cstdint>")
for p, l in lookup.items():
    print(f"// {p=}")
    data = []
    for k, v in l.items():
        k = int(k)
        data.append((k >> 64, k & 0xFFFFFFFFFFFFFFFF, v))
    build_mphf_c(data, f"p{p}")
