// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <exception>
#include <sstream>

#define THROW_EXCEPTION(MESSAGE)                          \
  {                                                       \
    std::stringstream except_str{};                       \
    except_str << __FILE__ << ":" << __LINE__ << " --- "; \
    except_str << MESSAGE;                                \
    throw std::runtime_error(except_str.str());           \
  }
namespace Acts {
/// Call the changeStripes method in the constructor
template <class T, size_t D>
template <typename... args>
DynamicArray<T, D>::DynamicArray(args... jks) {
  changeStripes(jks...);
}

template <class T, size_t D>
template <typename... args>
void DynamicArray<T, D>::changeStripes(args... jks) {
  /// Make sure that the correct number of arguments is given to the constructor
  constexpr size_t n_args = sizeof...(args);
  static_assert(n_args == D,
                "changeStripes()  -- Wrong number of arguments was given ");
  /// Use variadic templates to map the i-th argument to the i-th index in the
  /// m_strip_lengths array
  assign_stripes<0>(jks...);
  allocate();
}
template <class T, size_t D>
void DynamicArray<T, D>::allocate() {
  /// Recalculate the new allocated size
  size_t new_size = 1;
  for (size_t i : m_stripe_lengths)
    new_size *= i;
  /// Cache the offsets for each dimension. The last element does not have any
  m_stripe_offsets[D - 1] = 1;
  for (int off = static_cast<int>(D - 2); off >= 0; --off) {
    m_stripe_offsets[off] =
        m_stripe_offsets[off + 1] * m_stripe_lengths[off + 1];
  }
  /// Adapt for the new size of the vector
  if (m_size != new_size) {
    m_size = new_size;
    if (new_size == 0)
      m_data_ptr.clear();
    else
      m_data_ptr.resize(new_size);
  }
}
template <class T, size_t D>
void DynamicArray<T, D>::changeStripes(const std::array<size_t, D>& stripes) {
  m_stripe_lengths = stripes;
  allocate();
}
template <class T, size_t D>
template <size_t pos, typename... args>
void DynamicArray<T, D>::assign_stripes(size_t val, args... jks) {
  assign_stripes<pos>(val);
  assign_stripes<pos + 1>(jks...);
}
template <class T, size_t D>
template <size_t pos>
void DynamicArray<T, D>::assign_stripes(size_t val) {
  m_stripe_lengths[pos] = val;
}

/// Iterators
template <class T, size_t D>
inline typename std::vector<T>::iterator DynamicArray<T, D>::begin() {
  return m_data_ptr.begin();
}
template <class T, size_t D>
inline typename std::vector<T>::const_iterator DynamicArray<T, D>::begin()
    const {
  return m_data_ptr.begin();
}
template <class T, size_t D>
inline typename std::vector<T>::iterator DynamicArray<T, D>::end() {
  return m_data_ptr.end();
}
template <class T, size_t D>
inline typename std::vector<T>::const_iterator DynamicArray<T, D>::end() const {
  return m_data_ptr.end();
}
template <class T, size_t D>
inline typename std::vector<T>::reverse_iterator DynamicArray<T, D>::rbegin() {
  return m_data_ptr.rbegin();
}
template <class T, size_t D>
inline typename std::vector<T>::const_reverse_iterator
DynamicArray<T, D>::rbegin() const {
  return m_data_ptr.rbegin();
}
template <class T, size_t D>
inline typename std::vector<T>::reverse_iterator DynamicArray<T, D>::rend() {
  return m_data_ptr.rend();
}
template <class T, size_t D>
inline typename std::vector<T>::const_reverse_iterator
DynamicArray<T, D>::rend() const {
  return m_data_ptr.rend();
}
template <class T, size_t D>
inline size_t DynamicArray<T, D>::size() const {
  return m_size;
}
///
template <class T, size_t D>
inline DynamicArray<T, D>::operator bool() const {
  return !empty();
}
template <class T, size_t D>
inline bool DynamicArray<T, D>::operator!() const {
  return m_data_ptr.empty();
}
template <class T, size_t D>
inline bool DynamicArray<T, D>::empty() const {
  return m_data_ptr.empty();
}

/// Information of the underlying layout
template <class T, size_t D>
inline const std::array<size_t, D>& DynamicArray<T, D>::getStripes() const {
  return m_stripe_lengths;
}
template <class T, size_t D>
inline const std::array<size_t, D>& DynamicArray<T, D>::getMemoryOffsets()
    const {
  return m_stripe_offsets;
}

/// Value getters using the global index
template <class T, size_t D>
inline T& DynamicArray<T, D>::get_val(size_t idx) {
  if (idx >= size())
    THROW_EXCEPTION("get_val() -- Index "
                    << idx << " is out of range. Allowed is " << size() << ".");
  return m_data_ptr[idx];
}
template <class T, size_t D>
inline const T& DynamicArray<T, D>::get_val(size_t idx) const {
  if (idx >= size())
    THROW_EXCEPTION("get_val() -- Index "
                    << idx << " is out of range. Allowed is " << size() << ".");
  return m_data_ptr[idx];
}

//// Value getter the dimensional indices
template <class T, size_t D>
template <typename... args>
inline T& DynamicArray<T, D>::get(args... jks) {
  /// Ensure that the right number of arguments is given
  constexpr size_t n_args = sizeof...(args);
  static_assert(n_args == D, "get()  -- Wrong number of arguments was given ");
  /// Calculate the global index
  size_t idx = index(jks...);
  return get_val(idx);
}
template <class T, size_t D>
template <typename... args>
inline const T& DynamicArray<T, D>::get(args... jks) const {
  /// Ensure that the right number of arguments is given
  constexpr size_t n_args = sizeof...(args);
  static_assert(n_args == D, "get()  -- Wrong number of arguments was given ");
  size_t idx = index(jks...);
  return get_val(idx);
}
/// Calculate the global index from the array
//// Methods to calculate  the global index
template <class T, size_t D>
inline size_t DynamicArray<T, D>::index(
    const std::array<size_t, D>& indices) const {
  size_t idx{0};
  for (size_t i = 0; i < D; ++i) {
    idx += indices[i] * m_stripe_offsets[i];
    if (indices[i] >= m_stripe_lengths[i]) {
      THROW_EXCEPTION("index() --  The i-th "
                      << i << " index is out of range " << indices[i]
                      << ". Allowed is maximum " << m_stripe_lengths[i] << ".");
    }
  }
  return idx;
}
//// Calculate the global index from the local dimensional indices
template <class T, size_t D>
template <typename... args>
inline size_t DynamicArray<T, D>::index(args... jks) const {
  /// Ensure that the correct number of arguments is given
  constexpr size_t n_args = sizeof...(args);
  static_assert(n_args == D,
                "index()  -- Wrong number of arguments was given ");
  return index<0>(jks...);
}
template <class T, size_t D>
template <size_t pos, typename... args>
inline size_t DynamicArray<T, D>::index(size_t arg_val, args... jks) const {
  static_assert(pos < D, "index() -- Dimension exceeded");
  return index<pos>(arg_val) + index<pos + 1>(jks...);
}
template <class T, size_t D>
template <size_t pos>
inline size_t DynamicArray<T, D>::index(size_t arg_val) const {
  if (arg_val >= m_stripe_lengths[pos]) {
    THROW_EXCEPTION("index() --  The i-th "
                    << pos << " index is out of range " << arg_val
                    << ". Allowed is maximum " << m_stripe_lengths[pos] << ".");
  }
  return arg_val * m_stripe_offsets[pos];
}
template <class T, size_t D>
std::array<size_t, D> DynamicArray<T, D>::axisIndices(size_t globIdx) const {
  std::array<size_t, D> indices{};
  /// global idx = i*S1 + j*S2 + k*S3 + l*S4
  for (size_t d = 0; d < D; ++d) {
    const int axisIdx = globIdx / m_stripe_offsets[d];
    indices[d] = axisIdx;
    globIdx -= axisIdx * m_stripe_offsets[d];
  }
  return indices;
}

//// Assignemt
template <class T, size_t D>
void DynamicArray<T, D>::assign(const T& value) {
  for (size_t i = 0; i < size(); ++i)
    m_data_ptr[i] = value;
}

/// Access via [][][] operators

template <class T, size_t D>
inline typename DynamicArray<T, D>::array_type DynamicArray<T, D>::operator[](
    size_t i) {
  if constexpr (D != 1)
    return array_type(this, i);
  else
    return get(i);
}
template <class T, size_t D>
inline typename DynamicArray<T, D>::const_array_type
DynamicArray<T, D>::operator[](size_t i) const {
  if constexpr (D != 1)
    return const_array_type(this, i);
  else
    return get(i);
}

///
///                 ArrayNavigator
///
template <class T, size_t D>
template <size_t N>
inline DynamicArray<T, D>::ArrayNavigator<N>::ArrayNavigator(
    parent_type* parent, size_t idx)
    : m_parent{parent}, m_idx{idx} {}

template <class T, size_t D>
template <size_t N>
inline typename DynamicArray<T, D>::ArrayNavigator<N>::value_type
DynamicArray<T, D>::ArrayNavigator<N>::operator[](size_t i) {
  if constexpr (N == 1)
    return m_parent->get_val(index(i));
  else
    return value_type(this, i);
}

template <class T, size_t D>
template <size_t N>
template <class... args>
inline size_t DynamicArray<T, D>::ArrayNavigator<N>::index(
    args... indices) const {
  return m_parent->index(m_idx, indices...);
}
template <class T, size_t D>
template <size_t N>
inline T& DynamicArray<T, D>::ArrayNavigator<N>::get_val(size_t idx) {
  return m_parent->get_val(idx);
}

///
///                 ConstArrayNavigator
///
template <class T, size_t D>
template <size_t N>
inline DynamicArray<T, D>::ConstArrayNavigator<N>::ConstArrayNavigator(
    const parent_type* parent, size_t idx)
    : m_parent{parent}, m_idx{idx} {}

template <class T, size_t D>
template <size_t N>
inline typename DynamicArray<T, D>::ConstArrayNavigator<N>::value_type
DynamicArray<T, D>::ConstArrayNavigator<N>::operator[](size_t i) const {
  if constexpr (N == 1)
    return m_parent->get_val(index(i));
  else
    return value_type(this, i);
}

template <class T, size_t D>
template <size_t N>
template <class... args>
inline size_t DynamicArray<T, D>::ConstArrayNavigator<N>::index(
    args... indices) const {
  return m_parent->index(m_idx, indices...);
}
template <class T, size_t D>
template <size_t N>
inline const T& DynamicArray<T, D>::ConstArrayNavigator<N>::get_val(
    size_t idx) const {
  return m_parent->get_val(idx);
}
}  // namespace Acts
