#ifndef SPECK2D_INT_DEC_H
#define SPECK2D_INT_DEC_H

#include "SPECK2D_INT.h"

#include <cassert>


namespace sperr {

//
// Main SPECK2D_INT_DEC class
//
template <typename T>
class SPECK2D_INT_DEC final : public SPECK2D_INT<T> {
 private:
  //
  // Bring members from parent classes to this derived class.
  //
  using SPECK_INT<T>::m_LIP_mask;
  using SPECK_INT<T>::m_dims;
  using SPECK_INT<T>::m_LSP_new;
  using SPECK_INT<T>::m_bit_buffer;
  using SPECK_INT<T>::m_sign_array;
  using SPECK2D_INT<T>::m_LIS;
  using SPECK2D_INT<T>::m_I;
  using SPECK2D_INT<T>::m_code_S;
  using SPECK2D_INT<T>::m_code_I;

  void m_process_S(size_t idx1, size_t idx2, size_t& counter, bool need_decide) final;
  void m_process_P(size_t idx, size_t& counter, bool need_decide) final;
  void m_process_I(bool need_decide) final;
};

};  // namespace sperr

template <typename T>
void sperr::SPECK2D_INT_DEC<T>::m_process_S(size_t idx1,
                                            size_t idx2,
                                            size_t& counter,
                                            bool need_decide)
{
  auto& set = m_LIS[idx1][idx2];
  assert(!set.is_pixel());
  bool is_sig = true;

  if (need_decide)
    is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    counter++;
    m_code_S(idx1, idx2);
    set.make_empty();
  }
}

template <typename T>
void sperr::SPECK2D_INT_DEC<T>::m_process_P(size_t idx, size_t& counter, bool need_decide)
{
  bool is_sig = true;

  if (need_decide)
    is_sig = m_bit_buffer.rbit();

  if (is_sig) {
    counter++;
    m_sign_array.wbit(idx, m_bit_buffer.rbit());

    m_LSP_new.push_back(idx);
    m_LIP_mask.wfalse(idx);
  }
}

template <typename T>
void sperr::SPECK2D_INT_DEC<T>::m_process_I(bool need_decide)
{
  if (m_I.part_level > 0) {  // Only process `m_I` when it's not empty.
    bool is_sig = true;
    if (need_decide)
      is_sig = m_bit_buffer.rbit();
    if (is_sig)
      m_code_I();
  }
}

template class sperr::SPECK2D_INT_DEC<uint8_t>;
template class sperr::SPECK2D_INT_DEC<uint16_t>;
template class sperr::SPECK2D_INT_DEC<uint32_t>;
template class sperr::SPECK2D_INT_DEC<uint64_t>;


#endif
