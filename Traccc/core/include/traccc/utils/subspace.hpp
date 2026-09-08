/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"

// Detray include(s).
#include <detray/algebra/common/known_substructure_matrix.hpp>

// System include(s).
#include <array>
#include <cassert>
#include <cstddef>
#include <utility>

namespace traccc {

template <typename algebra_t, detray::dindex_type<algebra_t> kFullSize,
          detray::dindex_type<algebra_t> kSize = 2u>
struct subspace {
  static_assert(1u <= kSize, "Subspace size must be at least 1");
  static_assert(kSize <= kFullSize,
                "Subspace can only be as large as the full space");

 public:
  // Type declarations
  using size_type = detray::dindex_type<algebra_t>;

  // Define the number of bits we use to represent a single index;
  // currently, we hardcode to be three plus one sign bit so that the full
  // size of the subspace is limited to 8. Note that we also reserve one
  // element (all ones) as invalid.
  //
  // TODO: Find a nicer, automated solution to this.
  static constexpr size_type BITS_PER_INDEX = 4;
  static_assert(BITS_PER_INDEX >= 2);

  // Compute the total number of bits necessary to represent a subspace.
  static constexpr size_type TOTAL_BITS = BITS_PER_INDEX * kSize;
  static_assert(TOTAL_BITS <= 64,
                "Subspaces with more than 64 necessary bits are not supported");
  using axes_type = std::conditional_t<TOTAL_BITS <= 32, std::uint_least32_t,
                                       std::uint_least64_t>;

  static constexpr axes_type SIGN_BIT_MASK = 1 << (BITS_PER_INDEX - 1);
  static constexpr axes_type VALUE_BITS_MASK = (1 << (BITS_PER_INDEX - 1)) - 1;
  static constexpr axes_type ELEMENT_BITS_MASK = (1 << BITS_PER_INDEX) - 1;
  static_assert(((1 << (BITS_PER_INDEX - 1)) - 1) >= kFullSize);
  static_assert(ELEMENT_BITS_MASK == (SIGN_BIT_MASK | VALUE_BITS_MASK));

  static constexpr size_type size = kSize;
  static constexpr size_type fullSize = kFullSize;

  /// Construct from a container of axis indices.
  ///
  /// @param indices Unique, ordered indices
  template <typename SIZE_TYPE>
  TRACCC_HOST_DEVICE constexpr subspace(
      const std::array<SIZE_TYPE, kSize>& indices)
      : m_axes(0) {
    for (size_type i = 0u; i < kSize; ++i) {
      assert((indices[i] < kFullSize) and
             "Axis indices must be within the full space");
      set_index(i, static_cast<axes_type>(indices[i]));
    }
  }

  /// Get the projection index for a given element of the subspace.
  TRACCC_HOST_DEVICE
  constexpr size_type get_index(const size_type& i) const {
    const size_type sr = i * BITS_PER_INDEX;
    return (m_axes >> sr) & VALUE_BITS_MASK;
  }

  /// Check whether a given element of the subspace is valid or not.
  ///
  /// Invalid elements are not projected at all and result in an empty
  /// row or column in the projection and expansion matrix.
  TRACCC_HOST_DEVICE
  constexpr bool get_valid(const size_type& i) const {
    const size_type sr = i * BITS_PER_INDEX;
    return ((m_axes >> sr) & VALUE_BITS_MASK) != VALUE_BITS_MASK;
  }

  /// Retrieve the sign of an element of the subspace; if a true value is
  /// returned, the element will be projected in the negative.
  TRACCC_HOST_DEVICE
  constexpr bool get_sign(const size_type& i) const {
    const size_type sr = i * BITS_PER_INDEX;
    return (m_axes >> sr) & SIGN_BIT_MASK;
  }

  /// Set the projection index for a given element.
  TRACCC_HOST_DEVICE
  constexpr void set_index(const size_type& i, const axes_type& j) {
    assert(j == (j & ELEMENT_BITS_MASK));
    const size_type sl = i * BITS_PER_INDEX;
    m_axes = (m_axes & static_cast<axes_type>(~(VALUE_BITS_MASK << sl))) |
             static_cast<axes_type>(j << sl);
  }

  /// Set an element to be invalid so that it is not projected.
  TRACCC_HOST_DEVICE
  constexpr void set_invalid(const size_type& i) {
    // HACK: For reasons not understood by the author, the
    // `VALUE_BITS_MASK` variable cannot be used here.
    set_index(i, ((1 << (BITS_PER_INDEX - 1)) - 1));
  }

  /// Set the sign of an element, where trueish values indicate negative
  /// projection and falseish values indicate positive values.
  TRACCC_HOST_DEVICE
  constexpr void set_sign(const size_type& i, const bool s) {
    const size_type sl = i * BITS_PER_INDEX;
    m_axes = (m_axes & static_cast<axes_type>(~(SIGN_BIT_MASK << sl))) |
             static_cast<axes_type>((s ? SIGN_BIT_MASK : 0) << sl);
  }

  /// The substructure of a projection matrix onto @c D rows, where the
  /// subspace only selects columns with an index below @c MaxIndex.
  template <size_type D, size_type MaxIndex>
  using projection_substructure =
      typename detray::ksm::make_left_columns_substructure<
          D, kFullSize, MaxIndex>::canonical_type;

  /// Projection matrix that maps from the full space into the subspace.
  ///
  /// @tparam D the number of rows to project onto
  /// @tparam MaxIndex the subspace may only select columns with an index
  ///         below @c MaxIndex, so every column from @c MaxIndex on is a
  ///         structural zero of the result. The default promises nothing.
  ///
  /// @note The result is a known-substructure matrix, which a caller reads
  /// through compile-time indices. That is what keeps it in registers: the
  /// column a subspace selects is a runtime value, and writing through it
  /// forces the whole matrix into memory. Call @c to_dense on the result
  /// where a dense matrix is what the caller needs.
  template <size_type D, size_type MaxIndex = kFullSize>
  TRACCC_HOST_DEVICE auto projector() const
      -> detray::ksm::matrix<projection_substructure<D, MaxIndex>,
                             detray::dscalar<algebra_t>> {
    static_assert(D <= kSize,
                  "The dimension of projection should be smaller than "
                  "subspace dimension");
    static_assert(MaxIndex <= kFullSize,
                  "A subspace cannot select a column outside the full space");

    using scalar_t = detray::dscalar<algebra_t>;

    detray::ksm::matrix<projection_substructure<D, MaxIndex>, scalar_t> proj{};

    // A structured matrix is addressable at compile time only, so the
    // selected column cannot index the write. Each free column is compared
    // against it instead.
    const auto set_element = [this, &proj]<std::size_t I, std::size_t J>() {
      assert(!get_valid(I) || get_index(I) < MaxIndex);

      proj.template at<I, J>() =
          (get_index(I) == static_cast<size_type>(J))
              ? (get_sign(I) ? scalar_t{-1} : scalar_t{1})
              : scalar_t{0};
    };

    [&set_element]<std::size_t... Ks>(std::index_sequence<Ks...>) {
      (set_element.template operator()<Ks / MaxIndex, Ks % MaxIndex>(), ...);
    }(std::make_index_sequence<static_cast<std::size_t>(D) * MaxIndex>{});

    return proj;
  }

 private:
  axes_type m_axes;
};

}  // namespace traccc
