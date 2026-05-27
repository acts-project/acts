/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_parameters.hpp"

// Detray include(s).
#include <detray/definitions/algebra.hpp>

// VecMem include(s).
#include <vecmem/edm/container.hpp>

// System include(s).
#include <concepts>
#include <cstdint>

namespace traccc::edm {

/// Interface for the @c traccc::edm::track_state_collection type.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays of the SoA containers, or the variables of the AoS proxies
/// created on top of the SoA containers.
///
template <typename BASE>
class track_state : public BASE {

    public:
    /// @name Functions inherited from the base class
    /// @{

    /// Inherit the base class's constructor(s)
    using BASE::BASE;
    /// Inherit the base class's assignment operator(s).
    using BASE::operator=;

    /// @}

    /// @name Constants
    /// @{

    /// @c is_hole bit in the state word
    static constexpr std::uint8_t IS_HOLE_MASK = 0x01;
    /// @c is_smoothed bit in the state word
    static constexpr std::uint8_t IS_SMOOTHED_MASK = 0x02;

    /// @}

    /// @name Track State Information
    /// @{

    /// The "state word" of the track (non-const)
    ///
    /// @return A (non-const) vector of unsigned integers
    ///
    TRACCC_HOST_DEVICE
    auto& state() { return BASE::template get<0>(); }
    /// The "state word" of the track (const)
    ///
    /// @return A (const) vector of unsigned integers
    ///
    TRACCC_HOST_DEVICE
    const auto& state() const { return BASE::template get<0>(); }

    /// Chi^2 of the fitered parameters (non-const)
    ///
    /// @return A (non-const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    auto& filtered_chi2() { return BASE::template get<1>(); }
    /// Chi^2 of the filtered parameters (const)
    ///
    /// @return A (const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    const auto& filtered_chi2() const { return BASE::template get<1>(); }

    /// Chi^2 of the smoothed parameters (non-const)
    ///
    /// @return A (non-const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    auto& smoothed_chi2() { return BASE::template get<2>(); }
    /// Chi^2 of the smoothed parameters (const)
    ///
    /// @return A (const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    const auto& smoothed_chi2() const { return BASE::template get<2>(); }

    /// Chi^2 of the backward parameters (non-const)
    ///
    /// @return A (non-const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    auto& backward_chi2() { return BASE::template get<3>(); }
    /// Chi^2 of the backward parameters (const)
    ///
    /// @return A (const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    const auto& backward_chi2() const { return BASE::template get<3>(); }

    /// The filtered parameters of the track on the surface (non-const)
    ///
    /// @return A (non-const) vector of bound track parameters
    ///
    TRACCC_HOST_DEVICE
    auto& filtered_params() { return BASE::template get<4>(); }
    /// The filtered parameters of the track on the surface (const)
    ///
    /// @return A (const) vector of bound track parameters
    ///
    TRACCC_HOST_DEVICE
    const auto& filtered_params() const { return BASE::template get<4>(); }

    /// The smoothed parameters of the track on the surface (non-const)
    ///
    /// @return A (non-const) vector of bound track parameters
    ///
    TRACCC_HOST_DEVICE
    auto& smoothed_params() { return BASE::template get<5>(); }
    /// The smoothed parameters of the track on the surface (const)
    ///
    /// @return A (const) vector of bound track parameters
    ///
    TRACCC_HOST_DEVICE
    const auto& smoothed_params() const { return BASE::template get<5>(); }

    /// The index of the track's measurement on the current surface (non-const)
    ///
    /// @return A (non-const) vector of unsigned integers
    ///
    TRACCC_HOST_DEVICE
    auto& measurement_index() { return BASE::template get<6>(); }
    /// The index of the track's measurement on the current surface (const)
    ///
    /// @return A (const) vector of unsigned integers
    ///
    TRACCC_HOST_DEVICE
    const auto& measurement_index() const { return BASE::template get<6>(); }

    /// @}

    /// @name Utility functions
    /// @{

    /// Check if the track state is on a hole
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return @c true if the track state is on a hole, @c false otherwise
    ///
    TRACCC_HOST_DEVICE
    bool is_hole() const;
    /// Set the track state to be on a hole
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @param value The value to set
    ///
    TRACCC_HOST_DEVICE
    void set_hole(bool value = true);

    /// Check if the track state is smoothed
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @return @c true if the track state is smoothed, @c false otherwise
    ///
    TRACCC_HOST_DEVICE
    bool is_smoothed() const;
    /// Set the track state to be smoothed
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @param value The value to set
    ///
    TRACCC_HOST_DEVICE
    void set_smoothed(bool value = true);

    /// @}

};  // class track_state

/// SoA container describing the fit states of tracks on specific measurements
///
/// @tparam ALGEBRA The algebra type used to describe the tracks
///
template <detray::concepts::algebra ALGEBRA>
using track_state_collection = vecmem::edm::container<
    track_state,
    // state
    vecmem::edm::type::vector<std::uint8_t>,
    // filtered_chi2
    vecmem::edm::type::vector<detray::dscalar<ALGEBRA>>,
    // smoothed_chi2
    vecmem::edm::type::vector<detray::dscalar<ALGEBRA>>,
    // backward_chi2
    vecmem::edm::type::vector<detray::dscalar<ALGEBRA>>,
    // filtered_params
    vecmem::edm::type::vector<bound_track_parameters<ALGEBRA>>,
    // smoothed_params
    vecmem::edm::type::vector<bound_track_parameters<ALGEBRA>>,
    // measurement_index
    vecmem::edm::type::vector<unsigned int>>;

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/track_state_collection.ipp"
