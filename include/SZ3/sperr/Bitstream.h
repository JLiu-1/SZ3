#ifndef BITSTREAM_H
#define BITSTREAM_H

/*
 * Bitstream is intended and optimized for streaming reads and writes.
 *   It is heavily modeled after the bitstream structure in ZFP:
 *   (https://github.com/LLNL/zfp/blob/develop/include/zfp/bitstream.inl)
 *   For heavy use of random reads and writes, please use the Bitmask class instead.
 *
 * Bitstream uses words of 64 bits in its storage. A buffered word of 64 bits is also
 *   used to speed up stream reads and writes.
 *
 * A few caveats:
 *   1. Random reads CAN be achieved via repeated rseek() and rbit() calls.
 *      However, it will be much less efficient than the random reads in Bitmask.
 *   2. A function call of wseek() will erase the remaining bits in a buffered word, i.e.,
 *      from the wseek() position to the next word boundary, though the bits up to the wseek()
 *      position will be preserved. This design is for better efficiency of wbit().
 *   3. Because of 2, true random writes is not possible; it's only possible at the end of
 *      each word, e.g., positions of 63, 127, 191.
 *   4. A function call of flush() will align the writing position to the beginning of the
 *      next word, i.e., the number of truly useful bits is lost!
 *      One wants to call wtell() to retrieve and keep that info.
 *   5. Functions write_bitstream() and parse_bitstream() take in a raw pointer and the
 *      number of bits to write/read. The memory pointed to by the raw pointer needs to
 *      be big enough to hold the number of bits specified.
 *   6. get_bitstream() and write_bitstream() need to be supplied a number of bits because
 *      a Bitstream itself will lose track of how many useful bits are there after flush().
 *   7. Unlike std::vector, a bitstream does NOT have an equivalent concept of "size."
 *      Thus, capacity change brought by `reserve()` can be immediately used to read/write.
 */

#include <cstddef>
#include <cstdint>
#include <vector>
#include <cassert>
#include <cstring>
#include <iterator>  // std::distance()
#ifdef __SSE2__
#include <immintrin.h>
#endif

namespace sperr {

class Bitstream {
 public:
  // Constructor
  //
  Bitstream(size_t nbits = 0);  // How many bits does it hold initially?

  // Functions for both read and write
  //
  void rewind();
  auto capacity() const -> size_t;
  void reserve(size_t nbits);
  void reset();  // Reset the bitstream to be all 0's.

  // Functions for read
  //
  auto rtell() const -> size_t;
  void rseek(size_t offset);
  auto rbit() -> bool;

  // Functions for write
  //
  auto wtell() const -> size_t;
  void wseek(size_t offset);
  void wbit(bool bit);
  void flush();

  // Functions that provide or parse a compact bitstream
  //
  void write_bitstream(void* p, size_t num_bits) const;
  void parse_bitstream(const void* p, size_t num_bits);
  auto get_bitstream(size_t num_bits) const -> std::vector<std::byte>;

 private:
  uint64_t m_buffer = 0;  // incoming/outgoing bits
  size_t m_bits = 0;      // number of buffered bits

  std::vector<uint64_t>::iterator m_itr;  // Iterator to the next word to be read/written.
  std::vector<uint64_t> m_buf;
};

};  // namespace sperr

// Constructor
sperr::Bitstream::Bitstream(size_t nbits)
{
  m_itr = m_buf.begin();
  this->reserve(nbits);
}

// Functions for both read and write
void sperr::Bitstream::rewind()
{
  m_itr = m_buf.begin();
  m_buffer = 0;
  m_bits = 0;
}

auto sperr::Bitstream::capacity() const -> size_t
{
  return m_buf.size() * 64;
}

void sperr::Bitstream::reserve(size_t nbits)
{
  if (nbits > m_buf.size() * 64) {
#ifdef __SSE2__
    _mm_sfence();
#endif
    auto num_longs = (nbits + 63) / 64;
    const auto dist = std::distance(m_buf.begin(), m_itr);
    m_buf.resize(num_longs);         // trigger a memroy allocation.
    m_buf.resize(m_buf.capacity());  // be able to make use of all available capacity.
    m_itr = m_buf.begin() + dist;
  }
}

void sperr::Bitstream::reset()
{
  std::fill(m_buf.begin(), m_buf.end(), 0);
}

// Functions for read
auto sperr::Bitstream::rtell() const -> size_t
{
  // Stupid C++ insists that `m_buf.begin()` gives me a const iterator...
  std::vector<uint64_t>::const_iterator itr2 = m_itr;  // NOLINT
  return std::distance(m_buf.begin(), itr2) * 64 - m_bits;
}

void sperr::Bitstream::rseek(size_t offset)
{
  size_t div = offset >> 6;
  size_t rem = offset & 63;
  m_itr = m_buf.begin() + div;
  if (rem) {
    m_buffer = *m_itr >> rem;
    ++m_itr;
    m_bits = 64 - rem;
  }
  else {
    m_buffer = 0;
    m_bits = 0;
  }
}

auto sperr::Bitstream::rbit() -> bool
{
  if (m_bits == 0) {
    m_buffer = *m_itr;
    ++m_itr;
    m_bits = 64;
  }
  --m_bits;
  bool bit = m_buffer & uint64_t{1};
  m_buffer >>= 1;
  return bit;
}

// Functions for write
auto sperr::Bitstream::wtell() const -> size_t
{
  // Stupid C++ insists that `m_buf.begin()` gives me a const iterator...
  std::vector<uint64_t>::const_iterator itr2 = m_itr;  // NOLINT
  return std::distance(m_buf.begin(), itr2) * 64 + m_bits;
}

void sperr::Bitstream::wseek(size_t offset)
{
  size_t div = offset >> 6;
  size_t rem = offset & 63;
  m_itr = m_buf.begin() + div;
  if (rem) {
    m_buffer = *m_itr;
    m_buffer &= (uint64_t{1} << rem) - 1;
    m_bits = rem;
  }
  else {
    m_buffer = 0;
    m_bits = 0;
  }
}

void sperr::Bitstream::wbit(bool bit)
{
  m_buffer |= uint64_t{bit} << m_bits;

  if (++m_bits == 64) {
    if (m_itr == m_buf.end()) {  // allocate memory if necessary.
#ifdef __SSE2__
      _mm_sfence();
#endif
      auto dist = m_buf.size();
      m_buf.resize(std::max(size_t{1}, dist) * 2);
      m_itr = m_buf.begin() + dist;
    }
#ifdef __SSE2__
    auto dist = m_itr - m_buf.begin();
    long long int* ptr = reinterpret_cast<long long int*>(&m_buf[dist]);
    _mm_stream_si64(ptr, m_buffer);
#else
    *m_itr = m_buffer;
#endif
    ++m_itr;
    m_buffer = 0;
    m_bits = 0;
  }
}

void sperr::Bitstream::flush()
{
#ifdef __SSE2__
  _mm_sfence();
#endif
  if (m_bits) {  // only really flush when there are remaining bits.
    if (m_itr == m_buf.end()) {
      auto dist = m_buf.size();
      m_buf.resize(m_buf.size() + 1);
      m_itr = m_buf.begin() + dist;
    }
    *m_itr = m_buffer;
    ++m_itr;
    m_buffer = 0;
    m_bits = 0;
  }
}

// Functions to provide or parse a compact bitstream
void sperr::Bitstream::write_bitstream(void* p, size_t num_bits) const
{
  assert(num_bits <= m_buf.size() * 64);

  const auto num_longs = num_bits / 64;
  auto rem_bytes = num_bits / 8 - num_longs * sizeof(uint64_t);
  if (num_bits % 8 != 0)
    rem_bytes++;

  if (num_longs > 0)
    std::memcpy(p, m_buf.data(), num_longs * sizeof(uint64_t));

  if (rem_bytes > 0) {
    uint64_t value = m_buf[num_longs];
    auto* const p_byte = static_cast<std::byte*>(p);
    std::memcpy(p_byte + num_longs * sizeof(uint64_t), &value, rem_bytes);
  }
}

auto sperr::Bitstream::get_bitstream(size_t num_bits) const -> std::vector<std::byte>
{
  assert(num_bits <= m_buf.size() * 64);

  auto num_bytes = (num_bits + 7) / 8;
  auto tmp = std::vector<std::byte>(num_bytes);
  this->write_bitstream(tmp.data(), num_bits);

  return tmp;
}

void sperr::Bitstream::parse_bitstream(const void* p, size_t num_bits)
{
  this->reserve(num_bits);

  const auto num_longs = num_bits / 64;
  auto rem_bytes = num_bits / 8 - num_longs * sizeof(uint64_t);
  if (num_bits % 8 != 0)
    rem_bytes++;

  const auto* const p_byte = static_cast<const std::byte*>(p);

  if (num_longs > 0)
    std::memcpy(m_buf.data(), p_byte, num_longs * sizeof(uint64_t));

  if (rem_bytes > 0)
    std::memcpy(m_buf.data() + num_longs, p_byte + num_longs * sizeof(uint64_t), rem_bytes);

  this->rewind();
}




#endif
