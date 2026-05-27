/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"

// VecMem include(s).
#include <vecmem/edm/container.hpp>

// System include(s).
#include <array>
#include <compare>
#include <limits>

namespace traccc::edm {

/// Interface for the @c traccc::edm::spacepoint_collection class.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays of the SoA containers, or the variables of the AoS proxies
/// created on top of the SoA containers.
///
template <typename BASE>
class spacepoint : public BASE {

    public:
    /// @name Functions inherited from the base class
    /// @{

    /// Inherit the base class's constructor(s)
    using BASE::BASE;
    /// Inherit the base class's assignment operator(s).
    using BASE::operator=;

    /// @}

    /// @name Spacepoint Information
    /// @{

    /// Invalid measurement index
    static constexpr unsigned int INVALID_MEASUREMENT_INDEX =
        std::numeric_limits<unsigned int>::max();

    /// The index of the measurement producing this spacepoint (non-const)
    ///
    /// @return A (non-const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    auto& measurement_index_1() { return BASE::template get<0>(); }
    /// The index of the measurement producing this spacepoint (const)
    ///
    /// @return A (const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    const auto& measurement_index_1() const { return BASE::template get<0>(); }

    /// The index of the measurement producing this spacepoint (non-const)
    ///
    /// @return A (non-const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    auto& measurement_index_2() { return BASE::template get<1>(); }
    /// The index of the measurement producing this spacepoint (const)
    ///
    /// @return A (const) vector of <tt>unsigned int</tt> values
    ///
    TRACCC_HOST_DEVICE
    const auto& measurement_index_2() const { return BASE::template get<1>(); }

    /// Global / 3D position of the spacepoint
    ///
    /// @return A (non-const) vector of @c std::array<float,3> values
    ///
    TRACCC_HOST_DEVICE auto& global() { return BASE::template get<2>(); }
    /// Global / 3D position of the spacepoint
    ///
    /// @return A (const) vector of @c std::array<float,3> values
    ///
    TRACCC_HOST_DEVICE const auto& global() const {
        return BASE::template get<2>();
    }

    /// The X position of the spacepoint (non-const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A (non-const) reference to a @c float value
    ///
    TRACCC_HOST_DEVICE
    auto& x();
    /// The X position of the spacepoint (const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A (const) reference to a @c float value
    ///
    TRACCC_HOST_DEVICE
    const auto& x() const;

    /// The Y position of the spacepoint (non-const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A (non-const) reference to a @c float value
    ///
    TRACCC_HOST_DEVICE
    auto& y();
    /// The Y position of the spacepoint (const)
    ///
    /// @return A (const) reference to a @c float value
    ///
    TRACCC_HOST_DEVICE
    const auto& y() const;

    /// The Z position of the spacepoint (non-const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A (non-const) reference to a @c float value
    ///
    TRACCC_HOST_DEVICE
    auto& z();
    /// The Z position of the spacepoint (const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A (const) reference to a @c float value
    ///
    TRACCC_HOST_DEVICE
    const auto& z() const;

    /// The variation on the spacepoint's Z coordinate (non-const)
    ///
    /// @return A (non-const) vector of @c float values
    ///
    TRACCC_HOST_DEVICE auto& z_variance() { return BASE::template get<3>(); }
    /// The variation on the spacepoint's Z coordinate (const)
    ///
    /// @return A (const) vector of @c float values
    ///
    TRACCC_HOST_DEVICE const auto& z_variance() const {
        return BASE::template get<3>();
    }

    /// The radius of the spacepoint in the XY plane (non-const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A @c float value
    ///
    TRACCC_HOST_DEVICE auto radius() const;

    /// The variation on the spacepoint radious (non-const)
    ///
    /// @return A (non-const) vector of @c float values
    ///
    TRACCC_HOST_DEVICE auto& radius_variance() {
        return BASE::template get<4>();
    }
    /// The variation on the spacepoint radious (const)
    ///
    /// @return A (non-const) vector of @c float values
    ///
    TRACCC_HOST_DEVICE const auto& radius_variance() const {
        return BASE::template get<4>();
    }

    /// The azimuthal angle of the spacepoint in the XY plane (non-const)
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return A @c float value
    ///
    TRACCC_HOST_DEVICE auto phi() const;

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
    TRACCC_HOST_DEVICE bool operator==(const spacepoint<T>& other) const;

    /// Comparison operator
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @param[in] other The object to compare with
    /// @return A weak ordering object, describing the relation between the
    ///         two objects
    ///
    template <typename T>
    TRACCC_HOST_DEVICE std::weak_ordering operator<=>(
        const spacepoint<T>& other) const;

    /// @}

};  // class spacepoint

/// SoA container describing reconstructed spacepoints
using spacepoint_collection =
    vecmem::edm::container<spacepoint,
                           // measurement_index_1
                           vecmem::edm::type::vector<unsigned int>,
                           // measurement_index_2
                           vecmem::edm::type::vector<unsigned int>,
                           // global
                           vecmem::edm::type::vector<std::array<float, 3u>>,
                           // z_variance
                           vecmem::edm::type::vector<float>,
                           // radius_variance
                           vecmem::edm::type::vector<float>>;

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/spacepoint_collection.ipp"
