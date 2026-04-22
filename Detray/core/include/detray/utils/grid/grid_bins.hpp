// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/invalid_values.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <algorithm>

namespace detray::bins {

/// @brief Bin with a single entry
template <typename entry_t>
class single : public detray::ranges::single_view<entry_t> {
  using base_type = detray::ranges::single_view<entry_t>;

 public:
  using entry_type = entry_t;
  using base_type::base_type;
  using base_type::operator=;

  /// Default constructor initializer the bin with an invalid value
  DETRAY_HOST_DEVICE constexpr single()
      : base_type{detail::invalid_value<entry_t>()} {}

  /// @returns the storage capacity of this bin
  DETRAY_HOST_DEVICE
  constexpr dindex capacity() noexcept { return 1u; }

  /// Add a new entry to the bin
  template <typename E = entry_t>
  DETRAY_HOST_DEVICE constexpr void push_back(E&& entry) noexcept {
    (*this).ref() = std::forward<E>(entry);
  }

  /// @returns an initialized bin in the backend storage
  DETRAY_HOST_DEVICE
  constexpr auto init(entry_t entry = detail::invalid_value<entry_t>())
      -> single& {
    (*this).ref() = entry;
    return *this;
  }

  /// Equality operator
  ///
  /// @param rhs the single view to compare with
  ///
  /// @returns true if the single value is equal
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const single& rhs) const {
    return (*this).value() == rhs.value();
  }
};

/// @brief Bin that holds a collection of entries.
///
/// Keeps an additional counter to track the number of entries per bin.
/// The bin capacity is static.
template <typename entry_t, std::size_t N>
class static_array
    : public detray::ranges::view_interface<static_array<entry_t, N>> {
  using bin_view_t = detray::ranges::subrange<darray<entry_t, N>>;
  using const_bin_view_t = detray::ranges::subrange<const darray<entry_t, N>>;
  using bin_iterator_t = typename detray::ranges::iterator_t<bin_view_t>;
  using const_bin_iterator_t =
      typename detray::ranges::const_iterator_t<const_bin_view_t>;

 public:
  using entry_type = entry_t;

  /// Default constructor initializer the bin with an invalid value
  DETRAY_HOST_DEVICE constexpr static_array() { init(); };
  constexpr static_array(const static_array& other) = default;
  constexpr static_array(static_array&& other) noexcept = default;
  static_array& operator=(const static_array& other) noexcept = default;

  /// @returns view iterator over bin content in start or end position
  /// @{
  DETRAY_HOST_DEVICE bin_iterator_t begin() {
    bin_view_t bv{view()};
    return detray::ranges::begin(bv);
  }
  DETRAY_HOST_DEVICE bin_iterator_t end() {
    bin_view_t bv{view()};
    return detray::ranges::end(bv);
  }
  DETRAY_HOST_DEVICE
  const_bin_iterator_t begin() const { return detray::ranges::cbegin(view()); }
  DETRAY_HOST_DEVICE
  const_bin_iterator_t end() const { return detray::ranges::cend(view()); }
  /// @}

  /// @returns the number of entries in this bin - const
  DETRAY_HOST_DEVICE
  constexpr dindex size() const { return m_size; }

  /// The storage capacity of this bin
  DETRAY_HOST_DEVICE
  constexpr dindex capacity() const noexcept {
    return static_cast<dindex>(m_content.size());
  }

  /// Add a new entry to the bin
  template <typename E = entry_t>
  DETRAY_HOST_DEVICE constexpr void push_back(E&& entry) noexcept {
    m_content[m_size] = std::forward<E>(entry);
    ++m_size;
  }

  /// Initialize with a single entry filled
  ///
  /// @returns Access to the initialized bin
  DETRAY_HOST_DEVICE
  constexpr auto init(entry_t entry = detail::invalid_value<entry_t>())
      -> static_array& {
    // Initialize the storage element
    for (auto& e : m_content) {
      e = detail::invalid_value<entry_t>();
    }

    if (entry == detail::invalid_value<entry_t>()) {
      m_size = 0u;
      return *this;
    }

    // The bin has at least a capacity of 1
    m_content[0] = entry;
    m_size = 1u;

    return *this;
  }

  /// Initialize from an entire bin content @param content.
  ///
  /// @returns Access to the initialized bin
  DETRAY_HOST_DEVICE
  constexpr auto init(darray<entry_t, N> content) -> static_array& {
    m_content = content;
    m_size = 0u;
    for (const auto& entry : m_content) {
      if (entry != detail::invalid_value<entry_t>()) {
        ++m_size;
      }
    }

    return *this;
  }

  /// Equality operator
  ///
  /// @param rhs the bin entry to compare with
  ///
  /// @returns true if the content is equal
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const static_array& rhs) const {
    return m_content == rhs.m_content;
  }

 private:
  /// @returns the subrange on the valid bin content - const
  DETRAY_HOST_DEVICE constexpr auto view() const {
    return const_bin_view_t{m_content, dindex_range{0u, m_size}};
  }

  /// @returns the subrange on the valid bin content
  DETRAY_HOST_DEVICE constexpr auto view() {
    return bin_view_t{m_content, dindex_range{0u, m_size}};
  }

  /// Number of valid elements in the bin
  dindex m_size{0u};
  /// Bin entry container
  darray<entry_t, N> m_content{};
};

/// @brief Bin that views a collection of entries it does not own.
///
/// Used if the bin capacity is not static.
template <typename entry_t>
class dynamic_array
    : public detray::ranges::view_interface<dynamic_array<entry_t>> {
  using container_t = device_container_types::template vector_type<entry_t>;
  using bin_view_t = detray::ranges::subrange<container_t>;
  using const_bin_view_t = detray::ranges::subrange<const container_t>;

 public:
  struct data {
    dindex offset{0u};
    dindex size{0u};
    dindex capacity{0u};

    constexpr bool operator==(const data& rhs) const = default;

    DETRAY_HOST_DEVICE
    constexpr void update_offset(std::size_t shift) {
      offset += static_cast<dindex>(shift);
    }
  };

  using entry_type = entry_t;
  using entry_ptr_t = const entry_type*;
  using data_ptr_t = const data*;

  /// Default constructor initializer the bin with an invalid value
  DETRAY_HOST_DEVICE constexpr dynamic_array() { init(); };

  /// Construct from an externally owned container of bin content
  /// @param bin_storage and access to an offset, size and capacity
  /// in @param bin_data
  DETRAY_HOST_DEVICE
  dynamic_array(entry_type* bin_storage, data& bin_data)
      : m_data{&bin_data}, m_capacity{bin_data.capacity} {
    // Prevent null-dereference warning
    if (bin_storage) {
      m_global_storage = bin_storage + bin_data.offset;
    }
  }

  /// Construct from an externally owned container of bin content
  /// @param bin_storage and access to an offset, size and capacity
  /// in @param bin_data - const
  DETRAY_HOST_DEVICE
  dynamic_array(const entry_type* bin_storage, const data& bin_data)
      : m_data{&bin_data}, m_capacity{bin_data.capacity} {
    // Prevent null-dereference warning
    if (bin_storage) {
      m_global_storage = bin_storage + bin_data.offset;
    }
  }

  /// @returns view iterator over bin content in start or end position
  /// @{
  DETRAY_HOST_DEVICE auto begin() {
    bin_view_t bv{view()};
    return detray::ranges::begin(bv);
  }
  DETRAY_HOST_DEVICE auto end() {
    bin_view_t bv{view()};
    return detray::ranges::end(bv);
  }
  DETRAY_HOST_DEVICE
  auto begin() const { return detray::ranges::cbegin(view()); }
  DETRAY_HOST_DEVICE auto end() const { return detray::ranges::cend(view()); }
  /// @}

  /// @returns the number of entries in this bin - const
  DETRAY_HOST_DEVICE
  constexpr dindex size() const { return m_data->size; }

  /// The storage capacity of this bin
  DETRAY_HOST_DEVICE
  constexpr dindex capacity() const noexcept { return m_capacity; }

  /// Add a new entry to the bin
  /// @note This does not check the state of the container it points to!!!
  template <typename E = entry_type>
  DETRAY_HOST_DEVICE constexpr void push_back(E&& entry) {
    assert(m_capacity > 0);
    assert(m_data->size < m_capacity);
    assert(m_global_storage);

    if (m_data->size < m_capacity) {
      *(const_cast<entry_type*>(m_global_storage) + m_data->size) =
          std::forward<E>(entry);
      ++(const_cast<data*>(m_data)->size);
    }
  }

  /// @note The bin capacity has to be set correctly before calling this
  /// method
  /// @returns Access to an initialized bin in the backend storage
  DETRAY_HOST_DEVICE
  constexpr auto init(entry_type entry = detail::invalid_value<entry_type>())
      -> dynamic_array& {
    if (capacity() == 0u) {
      return *this;
    }

    // Initialize the storage element
    std::ranges::fill(this, detail::invalid_value<entry_type>());

    if (entry == detail::invalid_value<entry_type>()) {
      const_cast<data*>(m_data)->size = 0u;
      return *this;
    }

    // The bin has at least a capacity of 1
    this->front() = entry;
    const_cast<data*>(m_data)->size = 1u;

    return *this;
  }

  /// Initialize from an entire bin content @param content.
  ///
  /// @returns Access to the initialized bin
  template <typename storage_t>
  DETRAY_HOST_DEVICE constexpr auto init(const storage_t& content)
      -> dynamic_array& {
    const_cast<data*>(m_data)->size = 0u;
    for (dindex i{0u}; i < m_capacity; ++i) {
      if (content[i] != detail::invalid_value<entry_type>()) {
        const_cast<entry_type*>(m_global_storage)[i] = content[i];
        ++(const_cast<data*>(m_data)->size);
      }
    }
    return *this;
  }

  /// Equality operator
  ///
  /// @param rhs the bin to be compared with
  ///
  /// @returns true if the view is identical
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const dynamic_array& rhs) const {
    // Check if the bin points to the same data
    if (m_data == rhs.m_data || *m_data == *rhs.m_data) {
      return true;
    }
    if (m_data->size != rhs.m_data->size) {
      return false;
    }
    // It could still point to different data, but the
    // content is the same
    auto this_view = view();
    auto rhs_view = rhs.view();
    // Loop over the size of the bin and compare
    for (dindex i{0u}; i < m_data->size; ++i) {
      if (this_view[i] != rhs_view[i]) {
        return false;
      }
    }
    return true;
  }

 private:
  /// @returns the subrange on the valid bin content - const
  DETRAY_HOST_DEVICE auto view() const {
    return const_bin_view_t{m_global_storage, m_global_storage + m_data->size};
  }

  /// @returns the subrange on the valid bin content
  DETRAY_HOST_DEVICE auto view() {
    return bin_view_t{const_cast<entry_type*>(m_global_storage),
                      const_cast<entry_type*>(m_global_storage) + m_data->size};
  }

  /// Pointer to the global bin storage that is not owned by this class
  /// Includes the offset when part of a larger collection
  entry_ptr_t m_global_storage{nullptr};
  /// Access to bin data in the global storage
  data_ptr_t m_data{nullptr};
  /// Current bin capacity
  dindex m_capacity{0u};
};

template <typename entry_t, typename data_t>
dynamic_array(entry_t* bin_storage, data_t& bin_data) -> dynamic_array<entry_t>;

template <typename entry_t, typename data_t>
dynamic_array(const entry_t* bin_storage, const data_t& bin_data)
    -> dynamic_array<entry_t>;
/// @}

}  // namespace detray::bins
