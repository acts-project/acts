// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <type_traits>
#include <utility>

#include <traccc/definitions/hints.hpp>
#include <traccc/definitions/qualifiers.hpp>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>

namespace traccc {

/**
 * @brief Compile-time description of a tiled struct-of-arrays layout.
 *
 * Every type in `Ts...` becomes one column, and a one-dimensional array type
 * such as `float[3]` becomes one column per element, so that a warp-wide
 * access to any single scalar is contiguous in memory. The rows are grouped
 * into tiles of `TileSize` rows. Within a tile, the columns are laid out
 * column after column; the tiles follow each other in one allocation, which
 * holds `padded_size(n)` rows, a whole number of tiles.
 *
 * Because the tile size and the row size are compile-time constants, the
 * address of a cell is the tile start, which depends only on the row, plus a
 * constant offset per column. The capacity of the buffer never enters it.
 *
 * @tparam TileSize The number of rows in one tile.
 * @tparam Ts The column types, either scalars or one-dimensional arrays.
 */
template <unsigned int TileSize, typename... Ts>
struct basic_soa_layout {
  static_assert(sizeof...(Ts) > 0, "An SoA layout needs at least one column.");
  static_assert(((std::rank_v<Ts> <= 1) && ...),
                "SoA columns must be scalars or one-dimensional arrays.");
  static_assert((std::is_trivially_copyable_v<Ts> && ...),
                "SoA columns must be trivially copyable.");
  static_assert(TileSize > 0u, "An SoA tile needs at least one row.");

  /**
   * @brief Row count and index type. All sizes are 32-bit so that the device
   * side computes column addresses with 32-bit arithmetic.
   */
  using size_type = unsigned int;

  /**
   * @brief Number of rows in one tile.
   */
  static constexpr size_type tile_size = TileSize;

  /**
   * @brief Alignment of every column start, in bytes; at least one sector.
   */
  static constexpr std::size_t max_alignment =
      std::max({alignof(Ts)..., std::size_t{32}});

  /**
   * @brief Number of bytes that one row occupies across all columns.
   */
  static constexpr std::size_t total_size = (sizeof(Ts) + ... + 0);

  using types = std::tuple<Ts...>;

  /**
   * @brief The type of column `I`, an array type for an array column.
   */
  template <std::size_t I>
  using type_at = std::tuple_element_t<I, types>;

  /**
   * @brief The scalar type stored in column `I`.
   */
  template <std::size_t I>
  using flat_type_at = std::remove_all_extents_t<type_at<I>>;

  template <std::size_t I>
  static constexpr std::size_t size_at = sizeof(type_at<I>);

  template <std::size_t I>
  static constexpr std::size_t flat_size_at = sizeof(flat_type_at<I>);

  template <std::size_t I>
  static constexpr std::size_t extent_at = std::extent_v<type_at<I>>;

  /**
   * @brief Byte offset of column `I` within one row.
   */
  template <std::size_t I>
  static constexpr std::size_t offset_at =
      []<std::size_t... Is>(std::index_sequence<Is...>) {
        return (size_at<Is> + ... + 0);
      }(std::make_index_sequence<I>());

  /**
   * @brief Round a row count up to the next multiple of `TileSize`.
   */
  static constexpr size_type padded_size(size_type n) {
    return n % TileSize == 0 ? n : n + (TileSize - n % TileSize);
  }

  // Every column within a tile, and so every tile, starts at a multiple of
  // `max_alignment` bytes.
  static_assert(((TileSize * sizeof(std::remove_all_extents_t<Ts>) %
                      max_alignment ==
                  0) &&
                 ...),
                "SoA tile size times every scalar size must be a multiple "
                "of max_alignment.");
  // The device computes byte offsets within a tile in 32 bits.
  static_assert(static_cast<std::size_t>(TileSize) * total_size <=
                    std::numeric_limits<size_type>::max(),
                "SoA tile exceeds the 32-bit address range.");
};

/**
 * @brief Struct-of-arrays layout with the default tile size of 32 rows,
 * which is one warp.
 */
template <typename... Ts>
using soa_layout = basic_soa_layout<32u, Ts...>;

template <typename T>
struct is_soa_layout_helper {
  static constexpr bool value = false;
};

template <unsigned int TileSize, typename... Ts>
struct is_soa_layout_helper<basic_soa_layout<TileSize, Ts...>> {
  static constexpr bool value = true;
};

namespace concepts {
template <typename T>
concept soa_layout = is_soa_layout_helper<T>::value;
}

template <concepts::soa_layout layout_t>
class soa_buffer;

/**
 * @brief Non-owning description of one struct-of-arrays buffer.
 *
 * This is the type that a kernel payload carries. It is trivially copyable
 * and has no accessors; a `soa_device` is built from it inside the kernel.
 * The layout of this type must never depend on the build configuration,
 * because the host code that fills a payload and the device code that reads
 * it are compiled by different compilers.
 *
 * @tparam layout_t The `soa_layout` of the buffer.
 * @tparam is_const True to hold the allocation as a `const void*`. A
 *         `soa_device` over such a view cannot form a mutable column pointer,
 *         so each write accessor fails to compile where it is used.
 */
template <concepts::soa_layout layout_t, bool is_const = false>
class soa_view {
 public:
  using layout = layout_t;
  using size_type = unsigned int;
  using pointer_type = std::conditional_t<is_const, const void*, void*>;

  /// Default constructor, leaves the members uninitialised
  soa_view() = default;

  // Constructor from a buffer
  soa_view(const soa_buffer<layout>& buff)
      : m_ptr(buff.ptr()), m_size(buff.size()) {}

  /// Constructor from a mutable view, which drops the write permission.
  template <bool other_is_const>
    requires(is_const && !other_is_const)
  TRACCC_HOST_DEVICE soa_view(const soa_view<layout_t, other_is_const>& v)
      : m_ptr(v.ptr()), m_size(v.size()) {}

  /// Number of rows that the user asked for
  TRACCC_HOST_DEVICE size_type size() const { return m_size; }

  /// Number of rows that was allocated
  TRACCC_HOST_DEVICE size_type capacity() const {
    return layout::padded_size(m_size);
  }

  /// Start of the allocation
  TRACCC_HOST_DEVICE pointer_type ptr() const { return m_ptr; }

 private:
  pointer_type m_ptr;
  size_type m_size;
};

/**
 * @brief Read-only view over a struct-of-arrays buffer.
 */
template <concepts::soa_layout layout_t>
using soa_const_view = soa_view<layout_t, true>;

/**
 * @brief Owning struct-of-arrays buffer, allocated from a memory resource.
 *
 * Mirrors `vecmem::data::vector_buffer`: the buffer owns the memory, and it
 * converts to a `soa_view` over that memory.
 *
 * @tparam layout_t The `soa_layout` of the buffer.
 */
template <concepts::soa_layout layout_t>
class soa_buffer {
 public:
  using layout = layout_t;
  using view_type = soa_view<layout_t>;
  using size_type = unsigned int;

  /// Default constructor, an empty buffer
  soa_buffer() : m_ptr(nullptr), m_size(0) {}

  /// Standard constructor
  soa_buffer(size_type size, vecmem::memory_resource& resource) {
    m_size = size;

    // Exit early for null-capacity buffers.
    if (size == 0u) {
      return;
    }

    const size_type padded = layout::padded_size(size);
    const std::size_t bytes =
        layout::total_size * static_cast<std::size_t>(padded);
    if (padded < size || bytes > std::numeric_limits<size_type>::max()) {
      throw std::length_error(
          "soa_buffer: size exceeds the 32-bit address range");
    }

    void* ptr = resource.allocate(bytes, layout::max_alignment);
    m_ptr = vecmem::unique_alloc_ptr<char[]>(
        static_cast<char*>(ptr), vecmem::details::unique_alloc_deleter<char[]>(
                                     resource, bytes, layout::max_alignment));
  }

  /// Move constructor
  soa_buffer(soa_buffer&&) noexcept = default;

  /// Move assignment
  soa_buffer& operator=(soa_buffer&&) noexcept = default;

  /// Number of rows that the user asked for
  size_type size() const { return m_size; }

  /// Number of rows that was allocated
  size_type capacity() const { return layout::padded_size(m_size); }

  /// Start of the allocation
  void* ptr() const { return m_ptr.get(); }

 private:
  /// Data object owning the allocated memory
  vecmem::unique_alloc_ptr<char[]> m_ptr;
  size_type m_size;
};

/**
 * @brief Device-side accessor over a struct-of-arrays view.
 *
 * Built from a view inside a kernel, in the way a `vecmem::device_vector` is
 * built from a `vecmem::data::vector_view`. Types that name the columns are
 * expected to derive privately from this class and to inherit its
 * constructor, so that column indices do not leak to call sites.
 *
 * @tparam layout_t The `soa_layout` of the buffer.
 * @tparam is_in_global_memory False only for a view over shared or local
 *         memory.
 * @tparam is_const True to read from a `soa_const_view`.
 */
template <concepts::soa_layout layout_t, bool is_in_global_memory = true,
          bool is_const = false>
class soa_device {
  static constexpr unsigned int tile_size = layout_t::tile_size;

 public:
  using view_type = soa_view<layout_t, is_const>;
  using layout_type = layout_t;
  using size_type = unsigned int;

  template <std::size_t I>
  using type_at = typename layout_type::template type_at<I>;

  template <std::size_t I>
  using flat_type_at = typename layout_type::template flat_type_at<I>;

  /// Constructor from a view
  TRACCC_HOST_DEVICE explicit soa_device(const view_type& v)
      : m_ptr(v.ptr())
#ifndef NDEBUG
        ,
        m_size(v.size())
#endif
  {
  }

  /**
   * @brief Access element `i` of scalar column `I`.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 0 && !is_const)
  TRACCC_HOST_DEVICE type_at<I>& at(size_type i) {
    using value_t = type_at<I>;
#ifndef NDEBUG
    assert(i < m_size);
#endif
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I>), i);
    return ptr[i % tile_size];
  }

  /**
   * @brief Read element `i` of scalar column `I`, by value.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 0)
  TRACCC_HOST_DEVICE type_at<I> at(size_type i) const {
    using value_t = const type_at<I>;
#ifndef NDEBUG
    assert(i < m_size);
#endif
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I>), i);
    return ptr[i % tile_size];
  }

  /**
   * @brief Access element `j` of array column `I` in row `i`.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 1 && !is_const)
  TRACCC_HOST_DEVICE flat_type_at<I>& at(size_type i, size_type j) {
    using value_t = flat_type_at<I>;
#ifndef NDEBUG
    assert(i < m_size);
#endif
    assert(j < layout_type::template extent_at<I>);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               j * layout_type::template flat_size_at<I>),
        i);
    return ptr[i % tile_size];
  }

  /**
   * @brief Read element `j` of array column `I` in row `i`, by value.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 1)
  TRACCC_HOST_DEVICE flat_type_at<I> at(size_type i, size_type j) const {
    using value_t = const flat_type_at<I>;
#ifndef NDEBUG
    assert(i < m_size);
#endif
    assert(j < layout_type::template extent_at<I>);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               j * layout_type::template flat_size_at<I>),
        i);
    return ptr[i % tile_size];
  }

  /**
   * @brief Access element `J` of array column `I` in row `i`, with the
   * element index fixed at compile time.
   */
  template <std::size_t I, std::size_t J>
    requires(std::rank_v<type_at<I>> == 1 &&
             J < layout_type::template extent_at<I> && !is_const)
  TRACCC_HOST_DEVICE flat_type_at<I>& at(size_type i) {
    using value_t = flat_type_at<I>;
#ifndef NDEBUG
    assert(i < m_size);
#endif
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               J * layout_type::template flat_size_at<I>),
        i);
    return ptr[i % tile_size];
  }

  /**
   * @brief Read element `J` of array column `I` in row `i`, by value, with
   * the element index fixed at compile time.
   */
  template <std::size_t I, std::size_t J>
    requires(std::rank_v<type_at<I>> == 1 &&
             J < layout_type::template extent_at<I>)
  TRACCC_HOST_DEVICE flat_type_at<I> at(size_type i) const {
    using value_t = const flat_type_at<I>;
#ifndef NDEBUG
    assert(i < m_size);
#endif
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               J * layout_type::template flat_size_at<I>),
        i);
    return ptr[i % tile_size];
  }

 private:
  /**
   * @brief Start of the column that begins `byte_offset` bytes into a row,
   * within the tile that holds row `row`.
   */
  template <typename value_t>
  TRACCC_HOST_DEVICE value_t* column(size_type byte_offset,
                                     size_type row) const {
    // A const view gives a `const char*` base, so a write accessor cannot
    // reach its column, and only that accessor fails to compile.
    std::conditional_t<is_const, const char, char>* base =
        static_cast<std::conditional_t<is_const, const char, char>*>(m_ptr);
    value_t* ptr = reinterpret_cast<value_t*>(
        base + (row / tile_size) * tile_size * layout_type::total_size +
        byte_offset * tile_size);
#if defined(__CUDA_ARCH__)
    if constexpr (is_in_global_memory) {
      TRACCC_ASSUME_GLOBAL(base);
    }
    return ptr;
#elif defined(__GNUC__)
    return static_cast<value_t*>(
        __builtin_assume_aligned(ptr, layout_type::max_alignment));
#else
    return ptr;
#endif
  }

  typename view_type::pointer_type m_ptr;
#ifndef NDEBUG
  size_type m_size;
#endif
};

/**
 * @brief Read-only `soa_device`, built from a `soa_const_view`.
 */
template <concepts::soa_layout layout_t, bool is_in_global_memory = true>
using soa_const_device = soa_device<layout_t, is_in_global_memory, true>;

}  // namespace traccc
