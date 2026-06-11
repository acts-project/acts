// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/utils/grid/grid_bins.hpp"
#include "detray/utils/ranges.hpp"

namespace detray::detail {

/// @brief bin data state of a grid
///
/// Can be data-owning or not. Does not contain the data of the axes,
/// as that is managed by the multi-axis type directly.
template <bool is_owning, typename bin_t, typename containers>
class bin_storage : public detray::ranges::view_interface<
                        bin_storage<is_owning, bin_t, containers>> {
  template <typename T>
  using vector_t = typename containers::template vector_type<T>;
  using bin_range_t =
      std::conditional_t<is_owning, vector_t<bin_t>,
                         detray::ranges::subrange<vector_t<bin_t>>>;

 public:
  /// Bin type: single or static_array
  using bin_type = bin_t;
  /// Backend storage type for the grid
  using bin_container_type = vector_t<bin_t>;

  // Vecmem based view type
  using view_type = dvector_view<bin_type>;
  using const_view_type = dvector_view<const bin_type>;

  // Vecmem based buffer type
  using buffer_type = dvector_buffer<bin_type>;

  /// Default constructor
  bin_storage() = default;
  /// Copy constructor
  bin_storage(const bin_storage&) noexcept = default;
  /// Move constructor
  bin_storage(bin_storage&&) noexcept = default;

  /// Copy assignment
  bin_storage& operator=(const bin_storage&) noexcept = default;
  /// Move assignment
  bin_storage& operator=(bin_storage&&) noexcept = default;

  /// Construct containers using a memory resources
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST explicit bin_storage(vecmem::memory_resource& resource)
      : m_bin_data(&resource) {}

  /// Construct grid data from containers - move
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST_DEVICE explicit bin_storage(bin_container_type&& bin_data)
      : m_bin_data(std::move(bin_data)) {}

  /// Construct the non-owning type from the @param offset into the global
  /// container @param bin_data and the number of bins @param size
  template <bool owner = is_owning>
    requires(!owner)
  DETRAY_HOST_DEVICE bin_storage(bin_container_type& bin_data, dindex offset,
                                 dindex size)
      : m_bin_data(bin_data, dindex_range{offset, offset + size}) {}

  /// Construct the non-owning type from the @param offset into the global
  /// container @param bin_data and the number of bins @param size
  template <bool owner = is_owning>
    requires(!owner)
  DETRAY_HOST_DEVICE bin_storage(const bin_container_type& bin_data,
                                 dindex offset, dindex size)
      : m_bin_data(bin_data, dindex_range{offset, offset + size}) {}

  /// Construct bin storage from its vecmem view
  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE explicit bin_storage(const view_t& view)
      : m_bin_data(view) {}

  /// begin and end of the bin range
  /// @{
  DETRAY_HOST_DEVICE
  auto begin() { return detray::ranges::begin(m_bin_data); }
  DETRAY_HOST_DEVICE
  auto begin() const { return detray::ranges::cbegin(m_bin_data); }
  DETRAY_HOST_DEVICE
  auto end() { return detray::ranges::end(m_bin_data); }
  DETRAY_HOST_DEVICE
  auto end() const { return detray::ranges::cend(m_bin_data); }
  /// @}

  /// @returns the vecmem view of the bin storage
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST auto get_data() -> view_type {
    return detray::get_data(m_bin_data);
  }

  /// @returns the vecmem view of the bin storage - const
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST auto get_data() const -> const_view_type {
    return detray::get_data(m_bin_data);
  }

  /// Equality operator
  ///
  /// @param rhs bin storage to compare with
  ///
  /// @returns true if the bin data is equal
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const bin_storage& rhs) const {
    return m_bin_data == rhs.m_bin_data;
  }

 private:
  /// Container that holds all bin data when owning or a view into an
  /// externally owned container
  bin_range_t m_bin_data{};
};

/// Facade/wrapper for the data containers of the dynamic bin storage to fit in
/// the grid collection
template <typename bin_t, typename containers>
struct dynamic_bin_container {
  template <typename T>
  using vector_t = typename containers::template vector_type<T>;
  using bin_data_t = typename bin_t::data;

  vector_t<bin_data_t> bins{};
  vector_t<typename bin_t::entry_type> entries{};

  // Vecmem based view type
  using view_type = dmulti_view<dvector_view<bin_data_t>,
                                dvector_view<typename bin_t::entry_type>>;
  using const_view_type =
      dmulti_view<dvector_view<const bin_data_t>,
                  dvector_view<const typename bin_t::entry_type>>;

  // Vecmem based buffer type
  using buffer_type = dmulti_buffer<dvector_buffer<bin_data_t>,
                                    dvector_buffer<typename bin_t::entry_type>>;

  constexpr dynamic_bin_container() = default;
  DETRAY_HOST
  explicit dynamic_bin_container(vecmem::memory_resource* resource)
      : bins{resource}, entries{resource} {}
  dynamic_bin_container(const dynamic_bin_container& other) = default;
  dynamic_bin_container(dynamic_bin_container&& other) noexcept = default;

  dynamic_bin_container& operator=(const dynamic_bin_container&) noexcept =
      default;
  dynamic_bin_container& operator=(dynamic_bin_container&&) noexcept = default;

  /// Device-side construction from a vecmem based view type
  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE explicit dynamic_bin_container(view_t& view)
      : bins(detail::get<0>(view.m_view)),
        entries(detail::get<1>(view.m_view)) {}

  /// Insert bin data at the end
  template <typename grid_bin_range_t>
  DETRAY_HOST void append(const grid_bin_range_t& grid_bins) {
    const auto& g_bins = grid_bins.bin_data();
    const auto& g_entries = grid_bins.entry_data();

    bins.insert(bins.end(), g_bins.begin(), g_bins.end());

    // Update the bin offsets
    dindex_range new_bins{static_cast<dindex>(bins.size() - grid_bins.size()),
                          static_cast<dindex>(bins.size())};
    for (auto& bin : detray::ranges::subrange{bins, new_bins}) {
      bin.update_offset(entries.size());
    }

    entries.insert(entries.end(), g_entries.begin(), g_entries.end());
  }

  /// @returns a vecmem view on the bin data - non-const
  DETRAY_HOST auto get_data() -> view_type {
    return view_type{detray::get_data(bins), detray::get_data(entries)};
  }

  /// @returns a vecmem view on the bin data - const
  DETRAY_HOST
  auto get_data() const -> const_view_type {
    return const_view_type{detray::get_data(bins), detray::get_data(entries)};
  }

  /// @returns the number of bins
  DETRAY_HOST_DEVICE
  std::size_t size() const { return bins.size(); }

  /// Clear out all data
  DETRAY_HOST
  void clear() {
    bins.clear();
    entries.clear();
  }
};

/// @brief bin data state of a grid with dynamic bin capacities
///
/// Can be data-owning or not. Does not contain the data of the axes,
/// as that is managed by the multi-axis type directly.
template <bool is_owning, typename entry_t, typename containers>
class bin_storage<is_owning, detray::bins::dynamic_array<entry_t>, containers>
    : public detray::ranges::view_interface<bin_storage<
          is_owning, detray::bins::dynamic_array<entry_t>, containers>> {
  template <typename T>
  using vector_t = typename containers::template vector_type<T>;
  using bin_t = detray::bins::dynamic_array<entry_t>;
  using bin_data_t = typename bin_t::data;
  using bin_range_t =
      std::conditional_t<is_owning, vector_t<bin_data_t>,
                         detray::ranges::subrange<vector_t<bin_data_t>>>;
  using bin_iterator_t = typename detray::ranges::iterator_t<bin_range_t>;
  using const_bin_iterator_t =
      typename detray::ranges::const_iterator_t<bin_range_t>;

  using entry_range_t =
      std::conditional_t<is_owning, vector_t<entry_t>,
                         detray::ranges::subrange<vector_t<entry_t>>>;

  /// Iterator adapter that makes sure the bin storage returns a correctly
  /// initialized bin instance
  template <typename bin_itr_t>
  struct iterator_adapter {
    using difference_type = std::iter_difference_t<bin_itr_t>;
    using value_type = bin_t;
    using pointer = bin_t*;
    using reference = bin_t;
    using iterator_category =
        typename std::iterator_traits<bin_itr_t>::iterator_category;

    using data_ptr_t =
        std::conditional_t<std::is_same_v<bin_itr_t, const_bin_iterator_t>,
                           const entry_t*, entry_t*>;

    /// Default constructor required by LegacyIterator trait
    constexpr iterator_adapter() = default;

    DETRAY_HOST_DEVICE
    iterator_adapter(bin_itr_t&& itr, data_ptr_t entry_data)
        : m_entry_data{entry_data}, m_itr{std::move(itr)} {}

    /// Wrap iterator functionality
    /// @{
    DETRAY_HOST_DEVICE iterator_adapter& operator++() {
      ++m_itr;
      return *this;
    }
    DETRAY_HOST_DEVICE constexpr iterator_adapter operator++(int) {
      auto tmp(*this);
      ++(*this);
      return tmp;
    }
    DETRAY_HOST_DEVICE iterator_adapter& operator--()
      requires std::bidirectional_iterator<bin_itr_t>
    {
      --m_itr;
      return *this;
    }
    DETRAY_HOST_DEVICE constexpr iterator_adapter operator--(int)
      requires std::bidirectional_iterator<bin_itr_t>
    {
      auto tmp(*this);
      --(*this);
      return tmp;
    }
    DETRAY_HOST_DEVICE constexpr iterator_adapter& operator+=(
        const difference_type j)
      requires std::random_access_iterator<bin_itr_t>
    {
      m_itr += j;
      return *this;
    }
    DETRAY_HOST_DEVICE constexpr iterator_adapter& operator-=(
        const difference_type j)
      requires std::random_access_iterator<bin_itr_t>
    {
      m_itr -= j;
      return *this;
    }
    DETRAY_HOST_DEVICE
    constexpr decltype(auto) operator[](const difference_type i) const
      requires std::random_access_iterator<bin_itr_t>
    {
      return *(*this + i);
    }
    /// @}

    /// Construct and @returns a bin on the fly
    /// @{
    DETRAY_HOST_DEVICE
    constexpr auto operator*() const {
      return detray::bins::dynamic_array{m_entry_data, *m_itr};
    }
    /// @}

   private:
    DETRAY_HOST_DEVICE friend constexpr bool operator==(
        const iterator_adapter& lhs, const iterator_adapter& rhs) {
      return lhs.m_itr == rhs.m_itr;
    }
    DETRAY_HOST_DEVICE friend constexpr auto operator<=>(
        const iterator_adapter& lhs, const iterator_adapter& rhs)
      requires detray::ranges::random_access_iterator<bin_itr_t>
    {
#if defined(__apple_build_version__)
      const auto l{lhs.m_itr};
      const auto r{rhs.m_itr};
      if (l < r || (l == r && l < r)) {
        return std::strong_ordering::less;
      }
      if (l > r || (l == r && l > r)) {
        return std::strong_ordering::greater;
      }
      return std::strong_ordering::equivalent;
#else
      return lhs.m_itr <=> rhs.m_itr;
#endif
    }
    DETRAY_HOST_DEVICE
    friend difference_type operator-(const iterator_adapter& lhs,
                                     const iterator_adapter& rhs)
      requires detray::ranges::random_access_iterator<bin_itr_t>
    {
      return lhs.m_itr - rhs.m_itr;
    }
    DETRAY_HOST_DEVICE
    friend iterator_adapter operator-(const iterator_adapter& itr,
                                      difference_type i)
      requires std::random_access_iterator<bin_itr_t>
    {
      return {itr.m_itr - i, itr.m_entry_data};
    }
    DETRAY_HOST_DEVICE
    friend iterator_adapter operator+(const iterator_adapter& itr,
                                      difference_type i)
      requires std::random_access_iterator<bin_itr_t>
    {
      return {itr.m_itr + i, itr.m_entry_data};
    }
    DETRAY_HOST_DEVICE
    friend iterator_adapter operator+(difference_type i,
                                      const iterator_adapter& itr)
      requires detray::ranges::random_access_iterator<bin_itr_t>
    {
      return itr + i;
    }

    /// Access to the bin content
    data_ptr_t m_entry_data{};
    /// Iterator over bin data from which to construct a bin instance
    bin_itr_t m_itr{};
  };

 public:
  /// Bin type: dynamic_array
  using bin_type = bin_t;
  /// Backend storage type for the grid
  using bin_container_type = dynamic_bin_container<bin_t, containers>;

  // Vecmem based view type
  using view_type =
      dmulti_view<dvector_view<bin_data_t>, dvector_view<entry_t>>;
  using const_view_type =
      dmulti_view<dvector_view<const bin_data_t>, dvector_view<const entry_t>>;

  // Vecmem based buffer type
  using buffer_type =
      dmulti_buffer<dvector_buffer<bin_data_t>, dvector_buffer<entry_t>>;

  /// Default constructor
  bin_storage() = default;
  /// Copy constructor
  bin_storage(const bin_storage&) noexcept = default;
  /// Move constructor
  bin_storage(bin_storage&&) noexcept = default;

  /// Construct containers using a memory resources
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST explicit bin_storage(vecmem::memory_resource& resource)
      : m_bin_data(&resource), m_entry_data(&resource) {}

  /// Construct grid data from containers - move
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST_DEVICE explicit bin_storage(bin_container_type&& bin_data)
      : m_bin_data(std::move(bin_data.bins)),
        m_entry_data(std::move(bin_data.entries)) {}

  /// Construct the non-owning type from the @param offset into the global
  /// containers @param bin_data and the number of bins @param size
  template <bool owner = is_owning>
    requires(!owner)
  DETRAY_HOST_DEVICE bin_storage(bin_container_type& bin_data, dindex offset,
                                 dindex size)
      : m_bin_data(bin_data.bins, dindex_range{offset, offset + size}),
        m_entry_data(
            bin_data.entries,
            dindex_range{0u, static_cast<dindex>(bin_data.entries.size())}) {}

  /// Construct bin storage from its vecmem view
  template <concepts::device_view view_t>
  DETRAY_HOST_DEVICE explicit bin_storage(const view_t& view)
      : m_bin_data(detray::detail::get<0>(view.m_view)),
        m_entry_data(detray::detail::get<1>(view.m_view)) {}

  /// Copy assignment
  bin_storage& operator=(const bin_storage&) noexcept = default;
  /// Move assignment
  bin_storage& operator=(bin_storage&&) noexcept = default;

  const bin_range_t& bin_data() const { return m_bin_data; }
  const entry_range_t& entry_data() const { return m_entry_data; }

  /// begin and end of the bin range
  /// @{
  DETRAY_HOST_DEVICE
  auto begin() {
    return iterator_adapter<bin_iterator_t>{detray::ranges::begin(m_bin_data),
                                            m_entry_data.data()};
  }
  DETRAY_HOST_DEVICE
  auto begin() const {
    return iterator_adapter<const_bin_iterator_t>{
        detray::ranges::cbegin(m_bin_data), m_entry_data.data()};
  }
  DETRAY_HOST_DEVICE
  auto end() {
    return iterator_adapter<bin_iterator_t>{detray::ranges::end(m_bin_data),
                                            m_entry_data.data()};
  }
  DETRAY_HOST_DEVICE
  auto end() const {
    return iterator_adapter<const_bin_iterator_t>{
        detray::ranges::cend(m_bin_data), m_entry_data.data()};
  }
  /// @}

  /// @returns the vecmem view of the bin storage
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST auto get_data() -> view_type {
    return view_type{detray::get_data(m_bin_data),
                     detray::get_data(m_entry_data)};
  }

  /// @returns the vecmem view of the bin storage - const
  template <bool owner = is_owning>
    requires owner
  DETRAY_HOST auto get_data() const -> const_view_type {
    return const_view_type{detray::get_data(m_bin_data),
                           detray::get_data(m_entry_data)};
  }

  /// Equality operator
  ///
  /// @param rhs bin storage to compare with
  ///
  /// @returns true if the bin data is equal
  DETRAY_HOST_DEVICE
  constexpr bool operator==(const bin_storage& rhs) const {
    return m_bin_data == rhs.m_bin_data && m_entry_data == rhs.m_entry_data;
  }

 private:
  /// Container that holds all bin data when owning or a view into an
  /// externally owned container
  bin_range_t m_bin_data{};
  /// Container that holds all bin entries when owning or a view into an
  /// externally owned container
  entry_range_t m_entry_data{};
};

}  // namespace detray::detail
