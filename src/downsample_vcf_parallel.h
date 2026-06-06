#ifndef _CELLBOUNCER_DOWNSAMPLE_VCF_PARALLEL_H
#define _CELLBOUNCER_DOWNSAMPLE_VCF_PARALLEL_H

#include <bitset>
#include <cstdint>
#include <functional>
#include <utility>

// std::hash specializations for bitset<NBITS> and
// pair<bitset<NBITS>, bitset<NBITS>>. Required for unordered_map keying
// in downsample_vcf_parallel.cpp. Encoding-agnostic: hashes raw bits,
// independent of what those bits represent.

namespace std {
    template <>
    struct hash<std::bitset<NBITS>> {
        size_t operator()(const std::bitset<NBITS>& bs) const {
            const size_t numBlocks = NBITS / 64 + 1;
            size_t h = 0;
            for (size_t i = 0; i < numBlocks; i++) {
                uint64_t block = 0;
                size_t base = i * 64;
                for (size_t bit = 0; bit < 64; ++bit) {
                    size_t idx = base + bit;
                    if (idx >= NBITS) break;
                    if (bs.test(idx))
                        block |= (uint64_t(1) << bit);
                }
                h ^= std::hash<uint64_t>{}(block)
                     + 0x9e3779b97f4a7c15ULL
                     + (h << 6)
                     + (h >> 2);
            }
            return h;
        }
    };

    template <>
    struct hash<std::pair<std::bitset<NBITS>, std::bitset<NBITS>>> {
        size_t operator()(const std::pair<std::bitset<NBITS>, std::bitset<NBITS>>& p) const {
            size_t h1 = std::hash<std::bitset<NBITS>>{}(p.first);
            size_t h2 = std::hash<std::bitset<NBITS>>{}(p.second);
            return h1 ^ (h2 + 0x9e3779b97f4a7c15ULL + (h1 << 6) + (h1 >> 2));
        }
    };
}

#endif
