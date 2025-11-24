#ifndef BITMASK_H
#define BITMASK_H

/*
 * Bitmask is intended and optimized for 1) random reads and writes, and 2) bulk reads
 *   and writes, e.g., reading and writing 64 bits at a tiem.
 *   For heavy use of streaming reads and writes, please use the Bitstream class instead.
 *
 * Methods rlong(idx) and wlong(idx, value) are used for bulk reads and writes of
 *   64 bits. The idx parameter is the position of bits, specifying where the long integer is
 *   read or written. As a result, rlong(0) and rlong(63) will return the same integer,
 *   and wlong(64, val) and wlong(127, val) will write val to the same location.
 *
 * Bitmask does not automatically adjust its size. The size of a Bitmask is initialized
 *   at construction time, and is only changed by users calling the resize() method.
 *   The current size of a Bitmask can be queried by the size() method.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <algorithm>
#include <cassert>
#include <limits>

#if __cplusplus >= 202002L
#include <bit>
#endif
namespace sperr {

class Bitmask {
 public:
  // Constructor
  //
  Bitmask(size_t nbits = 0);  // How many bits does it hold initially?

  // Size related functions
  //
  auto size() const -> size_t;  // Num. of useful bits in this mask.
  void resize(size_t nbits);    // Resize to hold n bits.

  // Functions for read
  //
  auto rlong(size_t idx) const -> uint64_t;  // `idx` of the bit, not the long.
  auto rbit(size_t idx) const -> bool;

  // Functions to perform bulk tests.
  //
  // Two versions of the `has_true()` function. Both versions return -1 in case of no true found.
  //   1) Position == false: it returns 1 indicating finding a true.
  //   2) Position == true:  it returns the offset relative to `start` of the first true.
  template <bool Position>
  auto has_true(size_t start, size_t len) const -> int64_t;
  auto count_true() const -> size_t;  // How many 1's in this mask?

  // Functions for write
  //
  void wbit(size_t idx, bool bit);
  void wlong(size_t idx, uint64_t value);  // `idx` of the bit, not the long.
  void wtrue(size_t idx);                  // This is faster than `wbit(idx, true)`.
  void wfalse(size_t idx);                 // This is faster than `wbit(idx, false)`.
  void reset();                            // Set the current bitmask to be all 0's.
  void reset_true();                       // Set the current bitmask to be all 1's.

  // Functions for direct access of the underlying data buffer
  // Note: `use_bitstream()` reads the number of values (uint64_t type) that provide
  //       enough bits for the specified size of this mask.
  //
  auto view_buffer() const -> const std::vector<uint64_t>&;
  void use_bitstream(const void* p);

#if __cplusplus >= 202002L && defined __cpp_lib_three_way_comparison
  auto operator<=>(const Bitmask& rhs) const noexcept;
  auto operator==(const Bitmask& rhs) const noexcept -> bool;
#endif

 private:
  size_t m_num_bits = 0;
  std::vector<uint64_t> m_buf;
};

};  // namespace sperr


sperr::Bitmask::Bitmask(size_t nbits)
{
  auto num_longs = (nbits + 63) / 64;
  m_buf.assign(num_longs, 0);
  m_num_bits = nbits;
}

auto sperr::Bitmask::size() const -> size_t
{
  return m_num_bits;
}

void sperr::Bitmask::resize(size_t nbits)
{
  auto num_longs = (nbits + 63) / 64;
  m_buf.resize(num_longs, 0);
  m_num_bits = nbits;
}

auto sperr::Bitmask::rlong(size_t idx) const -> uint64_t
{
  return m_buf[idx >> 6];
}

auto sperr::Bitmask::rbit(size_t idx) const -> bool
{
  auto div = idx >> 6;  // idx / 64
  auto rem = idx & 63;  // idx % 64
  auto word = m_buf[div];
  word &= uint64_t{1} << rem;
  return word;
}

template <bool Position>
auto sperr::Bitmask::has_true(size_t start, size_t len) const -> int64_t
{
  auto long_idx = start >> 6;
  auto processed_bits = int64_t{0};
  auto word = m_buf[long_idx];
  auto answer = uint64_t{0};

  // Collect the remaining bits from the start long.
  auto begin_idx = start & 63;
  auto nbits = std::min(size_t{64}, begin_idx + len);
  for (auto i = begin_idx; i < nbits; i++) {
    answer |= word & (uint64_t{1} << i);
    if constexpr (Position) {
      if (answer != 0)
        return processed_bits;
    }
    processed_bits++;
  }
  if constexpr (!Position) {
    if (answer != 0)
      return 1;
  }

  // Examine the subsequent full longs.
  while (processed_bits + 64 <= len) {
    word = m_buf[++long_idx];
    if (word) {
      if constexpr (Position) {
#if __cplusplus >= 202002L
        int64_t i = std::countr_zero(word);
        return processed_bits + i;
#else
        for (int64_t i = 0; i < 64; i++)
          if (word & (uint64_t{1} << i))
            return processed_bits + i;
#endif
      }
      else
        return 1;
    }
    processed_bits += 64;
  }

  // Examine the remaining bits
  if (processed_bits < len) {
    nbits = len - processed_bits;
    assert(nbits < 64);
    word = m_buf[++long_idx];
    answer = 0;
    for (int64_t i = 0; i < nbits; i++) {
      answer |= word & (uint64_t{1} << i);
      if constexpr (Position) {
        if (answer != 0)
          return processed_bits + i;
      }
    }
    if constexpr (!Position) {
      if (answer != 0)
        return 1;
    }
  }

  return -1;
}
template auto sperr::Bitmask::has_true<true>(size_t, size_t) const -> int64_t;
template auto sperr::Bitmask::has_true<false>(size_t, size_t) const -> int64_t;

auto sperr::Bitmask::count_true() const -> size_t
{
  size_t counter = 0;
  if (m_buf.empty())
    return counter;

  // Note that unused bits in the last long are not guaranteed to be all 0's.
  for (size_t i = 0; i < m_buf.size() - 1; i++) {
    auto val = m_buf[i];
#if __cplusplus >= 202002L
    counter += std::popcount(val);
#else
    if (val != 0) {
      for (size_t j = 0; j < 64; j++)
        counter += ((val >> j) & uint64_t{1});
    }
#endif
  }
  const auto val = m_buf.back();
  if (val != 0) {
    for (size_t j = 0; j < m_num_bits - (m_buf.size() - 1) * 64; j++)
      counter += ((val >> j) & uint64_t{1});
  }

  return counter;
}

void sperr::Bitmask::wlong(size_t idx, uint64_t value)
{
  m_buf[idx >> 6] = value;
}

void sperr::Bitmask::wbit(size_t idx, bool bit)
{
  const auto wstart = idx >> 6;
  auto word = m_buf[wstart];

  auto mask1 = uint64_t{1} << (idx & 63);
  word &= ~mask1;

  auto mask2 = uint64_t{bit} << (idx & 63);
  word |= mask2;

  m_buf[wstart] = word;
}

void sperr::Bitmask::wtrue(size_t idx)
{
  const auto wstart = idx >> 6;
  const auto mask = uint64_t{1} << (idx & 63);

  auto word = m_buf[wstart];
  word |= mask;
  m_buf[wstart] = word;
}

void sperr::Bitmask::wfalse(size_t idx)
{
  const auto wstart = idx >> 6;
  const auto mask = uint64_t{1} << (idx & 63);

  auto word = m_buf[wstart];
  word &= ~mask;
  m_buf[wstart] = word;
}

void sperr::Bitmask::reset()
{
  std::fill(m_buf.begin(), m_buf.end(), 0);
}

void sperr::Bitmask::reset_true()
{
  std::fill(m_buf.begin(), m_buf.end(), std::numeric_limits<uint64_t>::max());
}

auto sperr::Bitmask::view_buffer() const -> const std::vector<uint64_t>&
{
  return m_buf;
}

void sperr::Bitmask::use_bitstream(const void* p)
{
  const auto* pu64 = static_cast<const uint64_t*>(p);
  std::copy(pu64, pu64 + m_buf.size(), m_buf.begin());
}

#if __cplusplus >= 202002L && defined __cpp_lib_three_way_comparison
auto sperr::Bitmask::operator<=>(const Bitmask& rhs) const noexcept
{
  auto cmp = m_num_bits <=> rhs.m_num_bits;
  if (cmp != 0)
    return cmp;

  if (m_num_bits % 64 == 0)
    return m_buf <=> rhs.m_buf;
  else {
    // Compare each fully used long.
    for (size_t i = 0; i < m_buf.size() - 1; i++) {
      cmp = m_buf[i] <=> rhs.m_buf[i];
      if (cmp != 0)
        return cmp;
    }
    // Compare the last partially used long.
    auto mylast = m_buf.back();
    auto rhslast = rhs.m_buf.back();
    for (size_t i = m_num_bits % 64; i < 64; i++) {
      auto mask = uint64_t{1} << i;
      mylast &= ~mask;
      rhslast &= ~mask;
    }
    return mylast <=> rhslast;
  }
}
auto sperr::Bitmask::operator==(const Bitmask& rhs) const noexcept -> bool
{
  return (operator<=>(rhs) == 0);
}
#endif




#endif
