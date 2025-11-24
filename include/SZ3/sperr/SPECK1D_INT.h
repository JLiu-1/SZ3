#ifndef SPECK1D_INT_H
#define SPECK1D_INT_H

#include "SPECK_INT.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>

#include <cstring>  // std::memcpy()

namespace sperr {

class Set1D {
  // In an effort to reduce the object size of Set1D, I choose to store only the first 7 bytes of
  //    both the `start` and `length` information of a Set1D. Using another 2 bytes to store the
  //    `part_level` info, this object fits in 16 bytes nicely.
  //
 private:
  std::array<uint8_t, 16> m_16 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

 public:
  void set_start(uint64_t val) { std::memcpy(m_16.data(), &val, 7); };
  void set_length(uint64_t val) { std::memcpy(m_16.data() + 7, &val, 7); };
  void set_level(uint16_t val) { std::memcpy(m_16.data() + 14, &val, 2); };
  auto get_start() const -> uint64_t
  {
    auto val = uint64_t{0};
    std::memcpy(&val, m_16.data(), 7);
    return val;
  }
  auto get_length() const -> uint64_t
  {
    auto val = uint64_t{0};
    std::memcpy(&val, m_16.data() + 7, 7);
    return val;
  }
  auto get_level() const -> uint16_t
  {
    auto val = uint16_t{0};
    std::memcpy(&val, m_16.data() + 14, 2);
    return val;
  }
};

//
// Main SPECK1D_INT class; intended to be the base class of both encoder and decoder.
//
template <typename T>
class SPECK1D_INT : public SPECK_INT<T> {
 protected:
  //
  // Bring members from the base class to this derived class.
  //
  using SPECK_INT<T>::m_LIP_mask;
  using SPECK_INT<T>::m_LSP_new;
  using SPECK_INT<T>::m_dims;
  using SPECK_INT<T>::m_coeff_buf;

  // The 1D case is different from 3D and 2D cases in that it implements additional logic that
  //    infers the significance of subsets by where the significant point is. With this
  //    consideration, functions such as m_process_S() and m_process_P() have different signatures
  //    during decoding/encoding, so they're implemented in their respective subclasses.
  //
  void m_clean_LIS() final;
  void m_initialize_lists() final;

  auto m_partition_set(Set1D) const -> std::array<Set1D, 2>;

  //
  // SPECK1D_INT specific data members
  //
  std::vector<std::vector<Set1D>> m_LIS;
};

};  // namespace sperr

template <typename T>
void sperr::SPECK1D_INT<T>::m_clean_LIS()
{
  for (auto& list : m_LIS) {
    auto it =
        std::remove_if(list.begin(), list.end(), [](const auto& s) { return s.get_length() == 0; });
    list.erase(it, list.end());
  }
}

template <typename T>
void sperr::SPECK1D_INT<T>::m_initialize_lists()
{
  const auto total_len = m_dims[0];
  auto num_of_parts = sperr::num_of_partitions(total_len);
  auto num_of_lists = num_of_parts + 1;
  if (m_LIS.size() < num_of_lists)
    m_LIS.resize(num_of_lists);
  std::for_each(m_LIS.begin(), m_LIS.end(), [](auto& list) { list.clear(); });

  // Put in two sets, each representing a half of the long array.
  Set1D set;
  set.set_length(total_len);  // Set represents the whole 1D array.
  auto sets = m_partition_set(set);
  m_LIS[sets[0].get_level()].emplace_back(sets[0]);
  m_LIS[sets[1].get_level()].emplace_back(sets[1]);
}

template <typename T>
auto sperr::SPECK1D_INT<T>::m_partition_set(Set1D set) const -> std::array<Set1D, 2>
{
  const auto start = set.get_start();
  const auto length = set.get_length();
  const auto level = set.get_level();
  std::array<Set1D, 2> subsets;

  // Prepare the 1st set
  auto& set1 = subsets[0];
  set1.set_start(start);
  set1.set_length(length - length / 2);
  set1.set_level(level + 1);
  // Prepare the 2nd set
  auto& set2 = subsets[1];
  set2.set_start(start + length - length / 2);
  set2.set_length(length / 2);
  set2.set_level(level + 1);

  return subsets;
}

template class sperr::SPECK1D_INT<uint64_t>;
template class sperr::SPECK1D_INT<uint32_t>;
template class sperr::SPECK1D_INT<uint16_t>;
template class sperr::SPECK1D_INT<uint8_t>;



#endif
