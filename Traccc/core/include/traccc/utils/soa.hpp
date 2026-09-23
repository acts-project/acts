/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

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
 * @brief Compile-time description of a struct-of-arrays layout.
 *
 * Every type in `Ts...` becomes one column, and a one-dimensional array type
 * such as `float[3]` becomes one column per element, so that a warp-wide
 * access to any single scalar is contiguous in memory. All columns of one
 * buffer live in a single allocation, laid out column after column, and each
 * column holds `padded_size(n)` rows. That row count is a multiple of
 * `max_alignment`, so every column starts at a multiple of `max_alignment`
 * bytes, which is at least one 32-byte sector.
 *
 * @tparam Ts The column types, either scalars or one-dimensional arrays.
 */
template <typename... Ts>
struct soa_layout {
  static_assert(sizeof...(Ts) > 0, "An SoA layout needs at least one column.");
  static_assert(((std::rank_v<Ts> <= 1) && ...),
                "SoA columns must be scalars or one-dimensional arrays.");
  static_assert((std::is_trivially_copyable_v<Ts> && ...),
                "SoA columns must be trivially copyable.");

  /**
   * @brief Row count and index type. All sizes are 32-bit so that the device
   * side computes column addresses with 32-bit arithmetic.
   */
  using size_type = unsigned int;

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
   * @brief Round a row count up to the next multiple of `max_alignment`.
   */
  static constexpr size_type padded_size(size_type n) {
    constexpr size_type a = static_cast<size_type>(max_alignment);
    return n % a == 0 ? n : n + (a - n % a);
  }
};

template <typename T>
struct is_soa_layout_helper {
  static constexpr bool value = false;
};

template <typename... Ts>
struct is_soa_layout_helper<soa_layout<Ts...>> {
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
      : m_ptr(buff.ptr()), m_size(buff.size()), m_capacity(buff.capacity()) {}

  /// Constructor from a mutable view, which drops the write permission.
  ///
  /// This is a template so that it is never a copy constructor, because a
  /// user-declared copy constructor would stop the type being trivially
  /// copyable, and a kernel payload carries this type.
  template <bool other_is_const>
    requires(is_const && !other_is_const)
  TRACCC_HOST_DEVICE soa_view(const soa_view<layout_t, other_is_const>& v)
      : m_ptr(v.ptr()), m_size(v.size()), m_capacity(v.capacity()) {}

  /// Number of rows that the user asked for
  TRACCC_HOST_DEVICE size_type size() const { return m_size; }

  /// Number of rows that was allocated
  TRACCC_HOST_DEVICE size_type capacity() const { return m_capacity; }

  /// Start of the allocation
  TRACCC_HOST_DEVICE pointer_type ptr() const { return m_ptr; }

 private:
  pointer_type m_ptr;
  size_type m_size;
  size_type m_capacity;
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
 * converts to a `soa_view` over that memory. The device side computes column
 * addresses in 32 bits, so a buffer is bounded at 4 GB; a larger request
 * throws `std::length_error` before any allocation is made.
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
  soa_buffer() : m_ptr(nullptr), m_size(0), m_capacity(0) {}

  /// Standard constructor
  soa_buffer(size_type size, vecmem::memory_resource& resource) {
    m_size = size;
    m_capacity = layout::padded_size(size);

    // Exit early for null-capacity buffers.
    if (size == 0u) {
      return;
    }

    // The first test catches a wrap-around in padded_size, the second the
    // 32-bit bound on the whole allocation.
    const std::size_t bytes =
        layout::total_size *
        static_cast<std::size_t>(layout::padded_size(size));
    if (m_capacity < size || bytes > std::numeric_limits<size_type>::max()) {
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
  size_type capacity() const { return m_capacity; }

  /// Start of the allocation
  void* ptr() const { return m_ptr.get(); }

 private:
  /// Data object owning the allocated memory
  vecmem::unique_alloc_ptr<char[]> m_ptr;
  size_type m_size;
  size_type m_capacity;
};

/**
 * @brief Device-side accessor over a struct-of-arrays view.
 *
 * Built from a view inside a kernel, in the way a `vecmem::device_vector` is
 * built from a `vecmem::data::vector_view`. Types that name the columns are
 * expected to derive privately from this class and to inherit its
 * constructor, so that column indices do not leak to call sites.
 *
 * Each accessor puts the column pointer in a `TRACCC_RESTRICT` local. The
 * device compiler uses that to move loads across stores to other buffers,
 * which is only correct if a kernel reaches every buffer through one pointer:
 * never hand two views, or a view and a raw pointer, of the same buffer to
 * one kernel. Two columns of one buffer are not separated by the qualifier.
 *
 * @tparam layout_t The `soa_layout` of the buffer.
 * @tparam is_in_global_memory False only for a view over shared or local
 *         memory, which is built on the device and never comes from a host
 *         view. A wrong value here is a silent miscompile.
 * @tparam is_const True to read from a `soa_const_view`. The write accessors
 *         then cannot form their column pointer, so each one fails to compile
 *         where it is used, and nowhere else.
 */
template <concepts::soa_layout layout_t, bool is_in_global_memory = true,
          bool is_const = false>
class soa_device {
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
      : m_ptr(v.ptr()), m_size(v.size()), m_capacity(v.capacity()) {}

  /**
   * @brief Access element `i` of scalar column `I`.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 0 && !is_const)
  TRACCC_HOST_DEVICE type_at<I>& at(size_type i) {
    using value_t = type_at<I>;
    assert(i < m_size);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I>));
    return ptr[i];
  }

  /**
   * @brief Read element `i` of scalar column `I`, by value.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 0)
  TRACCC_HOST_DEVICE type_at<I> at(size_type i) const {
    using value_t = const type_at<I>;
    assert(i < m_size);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I>));
    return ptr[i];
  }

  /**
   * @brief Access element `j` of array column `I` in row `i`.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 1 && !is_const)
  TRACCC_HOST_DEVICE flat_type_at<I>& at(size_type i, size_type j) {
    using value_t = flat_type_at<I>;
    assert(i < m_size);
    assert(j < layout_type::template extent_at<I>);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               j * layout_type::template flat_size_at<I>));
    return ptr[i];
  }

  /**
   * @brief Read element `j` of array column `I` in row `i`, by value.
   */
  template <std::size_t I>
    requires(std::rank_v<type_at<I>> == 1)
  TRACCC_HOST_DEVICE flat_type_at<I> at(size_type i, size_type j) const {
    using value_t = const flat_type_at<I>;
    assert(i < m_size);
    assert(j < layout_type::template extent_at<I>);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               j * layout_type::template flat_size_at<I>));
    return ptr[i];
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
    assert(i < m_size);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               J * layout_type::template flat_size_at<I>));
    return ptr[i];
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
    assert(i < m_size);
    value_t* TRACCC_RESTRICT ptr = column<value_t>(
        static_cast<size_type>(layout_type::template offset_at<I> +
                               J * layout_type::template flat_size_at<I>));
    return ptr[i];
  }

 private:
  /**
   * @brief Start of the column that begins `byte_offset` bytes into a row.
   *
   * The address-space hint has to run here, on the callee side of any
   * non-inlined boundary; in the constructor it would not reach a function
   * that receives this object by reference. The alignment hint is host-only:
   * on the device it buys nothing, and a pointer that passes through
   * `__builtin_assume_aligned` loses the effect of the restrict qualifier.
   */
  template <typename value_t>
  TRACCC_HOST_DEVICE value_t* column(size_type byte_offset) const {
    // A const view gives a `const char*` base, so a write accessor cannot
    // reach its column, and only that accessor fails to compile.
    std::conditional_t<is_const, const char, char>* base =
        static_cast<std::conditional_t<is_const, const char, char>*>(m_ptr);
    value_t* ptr = reinterpret_cast<value_t*>(base + byte_offset * m_capacity);
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

  /// Only the bounds asserts read `m_size`, so a release build leaves it
  /// unused. This object never crosses a translation unit and never escapes a
  /// kernel, so the compiler is free to drop it.
  typename view_type::pointer_type m_ptr;
  size_type m_size;
  size_type m_capacity;
};

/**
 * @brief Read-only `soa_device`, built from a `soa_const_view`.
 */
template <concepts::soa_layout layout_t, bool is_in_global_memory = true>
using soa_const_device = soa_device<layout_t, is_in_global_memory, true>;

}  // namespace traccc
