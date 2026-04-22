// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/detail/qualifiers.hpp"

// Vecmem include(s)
#include <vecmem/containers/vector.hpp>

namespace detray {

/// Device vector class with a minimized memory footprint
template <typename value_t>
struct compact_device_vector {
  using value_type = value_t;
  using size_type = unsigned int;
  using difference_type = std::ptrdiff_t;
  using reference = std::add_lvalue_reference_t<value_type>;
  using const_reference =
      std::add_lvalue_reference_t<std::add_const_t<value_type>>;
  using pointer = std::add_pointer_t<value_type>;
  using const_pointer = std::add_pointer_t<std::add_const_t<value_type>>;

  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;

  DETRAY_HOST_DEVICE explicit compact_device_vector(
      const vecmem::data::vector_view<value_type>& data)
      : m_ptr(data.ptr()), m_size(data.size()) {}

  DETRAY_HOST_DEVICE
  reference at(size_type pos) {
    assert(pos < size());
    return m_ptr[pos];
  }
  DETRAY_HOST_DEVICE
  const_reference at(size_type pos) const {
    assert(pos < size());
    return m_ptr[pos];
  }

  DETRAY_HOST_DEVICE
  reference operator[](size_type pos) { return m_ptr[pos]; }
  DETRAY_HOST_DEVICE
  const_reference operator[](size_type pos) const { return m_ptr[pos]; }

  DETRAY_HOST_DEVICE
  reference front() {
    assert(!empty());
    return m_ptr[0];
  }
  DETRAY_HOST_DEVICE
  const_reference front() const {
    assert(!empty());
    return m_ptr[0];
  }

  DETRAY_HOST_DEVICE
  reference back() {
    assert(!empty());
    return m_ptr[size() - 1];
  }
  DETRAY_HOST_DEVICE
  const_reference back() const {
    assert(!empty());
    return m_ptr[size() - 1];
  }

  DETRAY_HOST_DEVICE
  pointer data() { return m_ptr; }
  DETRAY_HOST_DEVICE
  const_pointer data() const { return m_ptr; }

  DETRAY_HOST_DEVICE
  iterator begin() { return iterator(m_ptr); }
  DETRAY_HOST_DEVICE
  const_iterator begin() const { return const_iterator(m_ptr); }
  DETRAY_HOST_DEVICE
  const_iterator cbegin() const { return begin(); }

  DETRAY_HOST_DEVICE
  iterator end() { return iterator(m_ptr + size()); }
  DETRAY_HOST_DEVICE
  const_iterator end() const { return const_iterator(m_ptr + size()); }
  DETRAY_HOST_DEVICE
  const_iterator cend() const { return end(); }

  DETRAY_HOST_DEVICE
  reverse_iterator rbegin() { return reverse_iterator(end()); }
  DETRAY_HOST_DEVICE
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  DETRAY_HOST_DEVICE
  const_reverse_iterator crbegin() const { return rbegin(); }

  DETRAY_HOST_DEVICE
  reverse_iterator rend() { return reverse_iterator(begin()); }
  DETRAY_HOST_DEVICE
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  DETRAY_HOST_DEVICE
  const_reverse_iterator crend() const { return rend(); }

  DETRAY_HOST_DEVICE
  bool empty() const { return (size() == 0); }
  DETRAY_HOST_DEVICE
  size_type size() const { return m_size; }
  DETRAY_HOST_DEVICE
  size_type capacity() const { return m_size; }

 private:
  pointer m_ptr;
  size_type m_size;
};

}  // namespace detray
