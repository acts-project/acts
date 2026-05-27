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

// System include(s).
#include <array>

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
    static_assert(
        TOTAL_BITS <= 64,
        "Subspaces with more than 64 necessary bits are not supported");
    using axes_type = std::conditional_t<TOTAL_BITS <= 32, std::uint_least32_t,
                                         std::uint_least64_t>;

    static constexpr axes_type SIGN_BIT_MASK = 1 << (BITS_PER_INDEX - 1);
    static constexpr axes_type VALUE_BITS_MASK =
        (1 << (BITS_PER_INDEX - 1)) - 1;
    static constexpr axes_type ELEMENT_BITS_MASK = (1 << BITS_PER_INDEX) - 1;
    static_assert(((1 << (BITS_PER_INDEX - 1)) - 1) >= kFullSize);
    static_assert(ELEMENT_BITS_MASK == (SIGN_BIT_MASK | VALUE_BITS_MASK));

    template <size_type ROWS, size_type COLS>
    using matrix_type = detray::dmatrix<algebra_t, ROWS, COLS>;

    using subspace_vector = matrix_type<kSize, 1u>;
    using fullspace_vector = matrix_type<kFullSize, 1u>;
    using projection_matrix = matrix_type<kSize, kFullSize>;
    using expansion_matrix = matrix_type<kFullSize, kSize>;

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

    /// Projection matrix that maps from the full space into the subspace.
    template <size_type D>
    TRACCC_HOST_DEVICE auto projector() const -> matrix_type<D, kFullSize> {
        static_assert(D <= kSize,
                      "The dimension of projection should be smaller than "
                      "subspace dimension");

        auto proj = matrix::zero<matrix_type<D, kFullSize>>();

        for (size_type i = 0u; i < D; ++i) {
            if (get_index(i) < kFullSize) {
                getter::element(proj, i, get_index(i)) =
                    (get_sign(i) ? -1.f : 1.f);
            }
        }
        return proj;
    }

    /// Expansion matrix that maps from the subspace into the full space.
    template <size_type D>
    TRACCC_HOST_DEVICE auto expander() const -> matrix_type<kFullSize, D> {
        static_assert(D <= kSize,
                      "The dimension of projection should be smaller than "
                      "subspace dimension");

        auto expn = matrix::zero<matrix_type<kFullSize, D>>();

        for (size_type i = 0u; i < kSize; ++i) {
            if (get_index(i) < kFullSize) {
                getter::element(expn, get_index(i), i) =
                    (get_sign(i) ? -1.f : 1.f);
            }
        }
        return expn;
    }

    private:
    axes_type m_axes;
};

}  // namespace traccc
