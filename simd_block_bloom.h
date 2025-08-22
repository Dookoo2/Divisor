
#pragma once
#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <new>
#include <stdexcept>
#include <vector>

#include <x86intrin.h>        
#include "hashutil.h"       

namespace simd_bloom {

using uint32 = std::uint32_t;
using uint64 = std::uint64_t;

[[gnu::always_inline]] inline uint64 rotl64(uint64 x, int r) noexcept {
  return (x << r) | (x >> (64 - r));
}
[[gnu::always_inline]] inline uint32 fastRange(uint32 x, uint32 n) noexcept {
  return static_cast<uint32>((uint64_t(x) * n) >> 32);
}

#ifdef __AVX2__

template <typename HashFamily = ::hashing::SimpleMixSplit>
class SimdBlockFilterFixed final {
  using Bucket = uint32[8];        

  static constexpr int kBlockShift = 14;   
  static constexpr int kBlockLen   = 1 << kBlockShift;

  const uint32 bucket_count_;
  std::unique_ptr<Bucket, void (*)(void*)> directory_;
  HashFamily hasher_;
  std::atomic<uint64> keys_{0};

  SimdBlockFilterFixed(const SimdBlockFilterFixed&)            = delete;
  SimdBlockFilterFixed& operator=(const SimdBlockFilterFixed&) = delete;
  SimdBlockFilterFixed(SimdBlockFilterFixed&&)                 = delete;
  SimdBlockFilterFixed& operator=(SimdBlockFilterFixed&&)      = delete;

  [[gnu::always_inline]] static __m256i MakeMask(uint32 h) noexcept {
    const __m256i ones = _mm256_set1_epi32(1);
    const __m256i mul  = _mm256_setr_epi32(
        0x47b6137bU, 0x44974d91U, 0x8824ad5bU, 0xa2b7289dU,
        0x705495c7U, 0x2df1424bU, 0x9efc4947U, 0x5c6bfb31U);

    __m256i x = _mm256_set1_epi32(h);
    x = _mm256_mullo_epi32(x, mul);
    x = _mm256_srli_epi32(x, 27);        
    return _mm256_sllv_epi32(ones, x);
  }

public:
  explicit SimdBlockFilterFixed(uint64 bits)
      : bucket_count_(std::max<uint32>(1, static_cast<uint32>(bits / 256))),
        directory_(nullptr, &std::free)
  {
    if (!__builtin_cpu_supports("avx2"))
      throw std::runtime_error("AVX2 не поддерживается CPU");

    const size_t bytes = static_cast<size_t>(bucket_count_) * sizeof(Bucket);
    Bucket* raw = nullptr;
    if (posix_memalign(reinterpret_cast<void**>(&raw), 64, bytes) != 0)
      throw std::bad_alloc();
    std::memset(raw, 0, bytes);
    directory_.reset(raw);
  }

  [[nodiscard]] uint64 SizeInBytes() const noexcept {
    return static_cast<uint64>(bucket_count_) * sizeof(Bucket);
  }
  [[nodiscard]] uint64 Keys() const noexcept { return keys_.load(); }

  [[gnu::always_inline]] void Add(uint64 key) noexcept {
    const uint64 h   = hasher_(key);
    const uint32 idx = fastRange(static_cast<uint32>(rotl64(h, 32)), bucket_count_);
    const __m256i mask = MakeMask(static_cast<uint32>(h));

    auto* bucket = &reinterpret_cast<__m256i*>(directory_.get())[idx];
    const __m256i cur = _mm256_load_si256(bucket);
    _mm256_store_si256(bucket, _mm256_or_si256(cur, mask));
    keys_.fetch_add(1, std::memory_order_relaxed);
  }

  void AddAll(const uint64* keys, size_t start, size_t end) {
    const size_t blocks = 1 + bucket_count_ / kBlockLen;
    std::vector<uint64> buf(blocks * kBlockLen);
    std::vector<size_t> used(blocks, 0);

    uint64 added = 0;
    for (size_t i = start; i < end; ++i) {
      const uint64 h   = hasher_(keys[i]);
      const uint32 idx = fastRange(static_cast<uint32>(rotl64(h, 32)), bucket_count_);
      const size_t blk = idx >> kBlockShift;

      size_t& len = used[blk];
      const size_t base = (blk << kBlockShift) + len;
      buf[base]     = h;
      buf[base + 1] = idx;
      len += 2;
      ++added;

      if (len == kBlockLen) { ApplyBlock(buf.data(), blk, len); len = 0; }
    }
    for (size_t blk = 0; blk < blocks; ++blk)
      if (used[blk]) ApplyBlock(buf.data(), blk, used[blk]);

    keys_.fetch_add(added, std::memory_order_relaxed);
  }
  void AddAll(const std::vector<uint64>& v) { AddAll(v.data(), 0, v.size()); }

  [[gnu::always_inline]] bool Find(uint64 key) const noexcept {
    const uint64 h   = hasher_(key);
    const uint32 idx = fastRange(static_cast<uint32>(rotl64(h, 32)), bucket_count_);
    const __m256i mask   = MakeMask(static_cast<uint32>(h));
    const __m256i bucket = reinterpret_cast<const __m256i*>(directory_.get())[idx];
    return _mm256_testc_si256(bucket, mask) != 0;
  }

private:
  [[gnu::noinline]] void ApplyBlock(const uint64* buf, size_t blk, size_t len) noexcept {
    const size_t base = blk << kBlockShift;
    for (size_t i = 0; i < len; i += 2) {
      const uint64 h   = buf[base + i];
      const uint32 idx = static_cast<uint32>(buf[base + i + 1]);
      const __m256i mask = MakeMask(static_cast<uint32>(h));

      auto* bucket = &reinterpret_cast<__m256i*>(directory_.get())[idx];
      const __m256i cur = _mm256_load_si256(bucket);
      _mm256_store_si256(bucket, _mm256_or_si256(cur, mask));
    }
  }
};

#endif
} 
