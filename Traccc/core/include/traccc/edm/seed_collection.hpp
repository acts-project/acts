/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/edm/container.hpp>

namespace traccc::edm {

/// Interface for the @c traccc::edm::seed_collection class.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays of the SoA containers, or the variables of the AoS proxies
/// created on top of the SoA containers.
///
template <typename BASE>
class seed : public BASE {

    public:
    /// @name Functions inherited from the base class
    /// @{

    /// Inherit the base class's constructor(s)
    using BASE::BASE;
    /// Inherit the base class's assignment operator(s).
    using BASE::operator=;

    /// @}

    /// @name Seed Information
    /// @{

    /// Index of the bottom spacepoint (non-const)
    ///
    /// @return A (non-const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    auto& bottom_index() { return BASE::template get<0>(); }
    /// Index of the bottom spacepoint (const)
    ///
    /// @return A (const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    const auto& bottom_index() const { return BASE::template get<0>(); }

    /// Index of the middle spacepoint (non-const)
    ///
    /// @return A (non-const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    auto& middle_index() { return BASE::template get<1>(); }
    /// Index of the middle spacepoint (const)
    ///
    /// @return A (const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    const auto& middle_index() const { return BASE::template get<1>(); }

    /// Index of the top spacepoint (non-const)
    ///
    /// @return A (non-const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    auto& top_index() { return BASE::template get<2>(); }
    /// Index of the top spacepoint (const)
    ///
    /// @return A (const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    const auto& top_index() const { return BASE::template get<2>(); }
    /// Quality of the seed (const)
    ///
    /// @return A (const) vector of <tt>float</tt> values
    ///
    TRACCC_HOST_DEVICE
    const auto& quality() const { return BASE::template get<3>(); }
    /// Quality of the seed (non-const)
    ///
    /// @return A (non-const) vector of <tt>float</tt> values
    ///
    TRACCC_HOST_DEVICE
    auto& quality() { return BASE::template get<3>(); }

    /// @}

    /// @name Utility functions
    /// @{

    /// Equality operator
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @param[in] other The object to compare with
    /// @return @c true if the objects are equal, @c false otherwise
    ///
    template <typename T>
    TRACCC_HOST_DEVICE bool operator==(const seed<T>& other) const;

    /// Comparison operator
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @param[in] other The object to compare with
    /// @return A strong ordering object, describing the relation between the
    ///         two objects
    ///
    template <typename T>
    TRACCC_HOST_DEVICE std::strong_ordering operator<=>(
        const seed<T>& other) const;

    /// @}

};  // class seed

/// SoA container describing reconstructed track seeds
using seed_collection =
    vecmem::edm::container<seed, vecmem::edm::type::vector<unsigned int>,
                           vecmem::edm::type::vector<unsigned int>,
                           vecmem::edm::type::vector<unsigned int>,
                           vecmem::edm::type::vector<float> >;

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/seed_collection.ipp"
