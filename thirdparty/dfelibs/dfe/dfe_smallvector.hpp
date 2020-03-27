// SPDX-License-Identifier: MIT
// Copyright 2019 Moritz Kiehn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file
/// \brief   Minimal small vector implementation
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2019-04-00, Initial version

#pragma once

#include <iterator>
#include <memory>

namespace dfe {

/// An continous container that stores some elements in-place.
///
/// \tparam T Stored element type
/// \tparam N Maximum number of elements stored in-place.
/// \tparam Allocator Allocator for elements of type T
///
/// If the vector contains less or equal than N elements, they are stored in
/// the vector itself without the need to allocate additional memory.
///
/// Supports access by index, iteration over elements, deleting all existing
/// elements from the vector, and adding elements at a specified location or at
/// the back.
template<typename T, std::size_t N, typename Allocator = std::allocator<T>>
class SmallVector {
public:
  using value_type = T;
  using size_type = std::size_t;
  using iterator = T*;
  using const_iterator = const T*;

  SmallVector() = default;
  ~SmallVector() { clear(); }

  value_type& operator[](size_type idx);
  const value_type& operator[](size_type idx) const;

  iterator begin();
  iterator end() { return begin() + m_size; }
  const_iterator begin() const;
  const_iterator end() const { return begin() + m_size; }

  /// Return true if there are no elements in the vector.
  bool empty() const { return (m_size == 0); }
  /// Return the number of elements in the vector.
  size_type size() const { return m_size; }
  /// Return the number of elements that can be stored in the available memory.
  size_type capacity() const { return (m_size <= N) ? N : m_onheap.capacity; }

  /// Remove all elements.
  ///
  /// This will release allocated memory if the vector contains more elements
  /// than can be stored in-place.
  void clear();
  /// Construct an element directly before the given position in the vector.
  template<typename... Args>
  iterator emplace(const_iterator pos, Args&&... args);
  /// Construct an element at the back of the vector and return its reference.
  template<typename... Args>
  T& emplace_back(Args&&... args);

private:
  struct AllocatedStorage {
    size_type capacity;
    T* data;
  };

  AllocatedStorage allocate_storage(size_type capacity);
  void destruct_inplace();
  void destruct_deallocate_onheap();

  size_type m_size = 0;
  union {
    AllocatedStorage m_onheap;
    // use 'raw' memory to have full control over constructor/destructor calls
    alignas(T) char m_inplace[N * sizeof(T)];
  };
  Allocator m_alloc;
};

// implementation

template<typename T, std::size_t N, typename Allocator>
inline typename SmallVector<T, N, Allocator>::AllocatedStorage
SmallVector<T, N, Allocator>::allocate_storage(size_type capacity) {
  AllocatedStorage s;
  s.capacity = capacity;
  s.data = std::allocator_traits<Allocator>::allocate(m_alloc, capacity);
  return s;
}

// Destruct elements in in-place storage assuming they are valid.
template<typename T, std::size_t N, typename Allocator>
inline void
SmallVector<T, N, Allocator>::destruct_inplace() {
  T* ptr = reinterpret_cast<T*>(m_inplace);
  for (T* end = ptr + m_size; ptr != end; ++ptr) {
    ptr->~T();
  }
}

// Destruct and deallocate elements in heap-allocated storage.
template<typename T, std::size_t N, typename Allocator>
inline void
SmallVector<T, N, Allocator>::destruct_deallocate_onheap() {
  T* ptr = m_onheap.data;
  for (T* end = ptr + m_size; ptr != end; ++ptr) {
    std::allocator_traits<Allocator>::destroy(m_alloc, ptr);
  }
  std::allocator_traits<Allocator>::deallocate(
    m_alloc, m_onheap.data, m_onheap.capacity);
  m_onheap.capacity = 0;
  m_onheap.data = nullptr;
}

template<typename T, std::size_t N, typename Allocator>
inline T& SmallVector<T, N, Allocator>::operator[](size_type idx) {
  if (m_size <= N) {
    return *(reinterpret_cast<T*>(m_inplace) + idx);
  } else {
    return m_onheap.data[idx];
  }
}

template<typename T, std::size_t N, typename Allocator>
inline const T& SmallVector<T, N, Allocator>::operator[](size_type idx) const {
  if (m_size <= N) {
    return *(reinterpret_cast<const T*>(m_inplace) + idx);
  } else {
    return m_onheap.data[idx];
  }
}

template<typename T, std::size_t N, typename Allocator>
inline typename SmallVector<T, N, Allocator>::iterator
SmallVector<T, N, Allocator>::begin() {
  return (m_size <= N) ? reinterpret_cast<T*>(m_inplace) : m_onheap.data;
}

template<typename T, std::size_t N, typename Allocator>
inline typename SmallVector<T, N, Allocator>::const_iterator
SmallVector<T, N, Allocator>::begin() const {
  return (m_size <= N) ? reinterpret_cast<const T*>(m_inplace) : m_onheap.data;
}

template<typename T, std::size_t N, typename Allocator>
inline void
SmallVector<T, N, Allocator>::clear() {
  if (m_size <= N) {
    destruct_inplace();
  } else {
    destruct_deallocate_onheap();
  }
  m_size = 0;
}

template<typename T, std::size_t N, typename Allocator>
template<typename... Args>
typename SmallVector<T, N, Allocator>::iterator
SmallVector<T, N, Allocator>::emplace(const_iterator pos, Args&&... args) {
  using AllocatorTraits = std::allocator_traits<Allocator>;

  // TODO how, when to check iterator validity?

  // available storage is sufficient to hold one extra element.
  // existing data after the insertion point needs to be shifted by 1 to the
  // right so the new element can be constructed at the given position.
  if (size() < capacity()) {
    T* i = const_cast<T*>(pos);
    T* e = end();

    // the underlying storage is raw memory. in order for the move assignment
    // to be well-defined, the memory for the additiona element needs to be
    // (default-)initialized using placement-new first.
    (void)new (e) T();
    std::move_backward(i, e, std::next(e));
    // insertion points contains moved-from object. move-assignment should be
    // valid. placement-new construction would double-initialize.
    *i = T(std::forward<Args>(args)...);

    m_size += 1;
    return i;
  }

  // available storage is in-sufficient. move to larger heap-allocated storage.
  auto storage = allocate_storage(1.3f * (m_size + 1));
  T* source = begin();
  T* target = storage.data;

  // move data before insertion point to the new storage
  for (T* e = const_cast<T*>(pos); source != e; ++source, ++target) {
    AllocatorTraits::construct(m_alloc, target, std::move(*source));
  }
  // construct element directly in the new storage
  T* insert = target++;
  AllocatorTraits::construct(m_alloc, insert, std::forward<Args>(args)...);
  // move data after the insertion point to the new storage
  for (T* e = end(); source != e; ++source, ++target) {
    AllocatorTraits::construct(m_alloc, target, std::move(*source));
  }

  // clear previous data before replacing it with the next storage
  if (m_size == N) {
    destruct_inplace();
  } else {
    destruct_deallocate_onheap();
  }
  m_onheap = storage;

  m_size += 1;
  return insert;
}

template<typename T, std::size_t N, typename Allocator>
template<typename... Args>
inline typename SmallVector<T, N, Allocator>::value_type&
SmallVector<T, N, Allocator>::emplace_back(Args&&... args) {
  return *emplace(end(), std::forward<Args>(args)...);
}

} // namespace dfe
