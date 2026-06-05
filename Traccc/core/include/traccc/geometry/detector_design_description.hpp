/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/primitives.hpp"
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/definitions/track_parametrization.hpp"
#include "traccc/utils/subspace.hpp"

// Detray include(s).
#include <detray/geometry/identifier.hpp>

// VecMem include(s).
#include <vecmem/edm/container.hpp>

namespace traccc {

/// Interface for the @c traccc::detector_design_description class.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays defined in @c traccc::detector_design_description.
///
template <typename BASE>
class detector_design_description_interface : public BASE {

    public:
    /// Inherit the base class's constructor(s)
    using BASE::BASE;

    /// @name Detector segmentation information
    /// @{

    /// The identifier of the module's design
    ///
    ///
    /// @return A (const) vector of @c int objects
    ///
    TRACCC_HOST_DEVICE
    auto& design_id() { return BASE::template get<0>(); }
    /// The identifier of the detector module's surface (const)
    ///
    /// Can be used to look up the module in a @c detray::detector object.
    ///
    /// @return A (const) vector of @c detray::geometry::identifier objects
    ///
    TRACCC_HOST_DEVICE
    const auto& design_id() const { return BASE::template get<0>(); }

    /// Reference for local position calculation in X direction const
    ///
    /// The position of a detector element (pixel or strip) is calculated
    /// along the X axis with the formula:
    /// @f$pos_x =  0.5 * (bin_edge_x * index_x + bin_edge_x * (index_x +1)) @f$
    ///
    /// @return A (const) vector of @c traccc::scalar objects
    ///
    TRACCC_HOST_DEVICE
    auto& bin_edges_x() { return BASE::template get<1>(); }
    /// Vector of local cell centres in X direction (const)
    ///
    /// The position of a detector element (pixel or strip) is calculated
    /// along the X axis with the formula:
    /// @f$pos_x = bin_centres_x[cell_index_x]@f$
    ///
    /// @return A (const) vector of @c traccc::scalar objects
    ///
    TRACCC_HOST_DEVICE
    const auto& bin_edges_x() const { return BASE::template get<1>(); }

    /// Reference for local position calculation in Y direction (non-const)
    ///
    /// The position of a detector element (pixel or strip) is calculated
    /// along the Y axis with the formula:
    /// @f$pos_y = bin_centres_y[cell_index_y]@f$
    ///
    /// @return A (non-const) vector of @c traccc::scalar objects
    ///
    TRACCC_HOST_DEVICE
    auto& bin_edges_y() { return BASE::template get<2>(); }
    /// Reference for local position calculation in Y direction (const)
    ///
    /// The position of a detector element (pixel or strip) is calculated
    /// along the Y axis with the formula:
    /// @f$pos_y = reference_y + pitch_y * index_y@f$
    ///
    /// @return A (const) vector of @c traccc::scalar objects
    ///
    TRACCC_HOST_DEVICE
    const auto& bin_edges_y() const { return BASE::template get<2>(); }

    /// The dimensionality (1D/2D) of the detector module (non-const)
    ///
    /// @return A (non-const) vector of @c char objects
    ///
    TRACCC_HOST_DEVICE
    auto& dimensions() { return BASE::template get<3>(); }
    /// The dimensionality (1D/2D) of the detector module (const)
    ///
    /// @return A (const) vector of @c char objects
    ///
    TRACCC_HOST_DEVICE
    const auto& dimensions() const { return BASE::template get<3>(); }

    /// The subspace of measurements on the module
    ///
    /// The "subspace" defines which of the measurement's parameters are
    /// "sensitive", to be used during the track finding/fitting.
    ///
    /// @return A (non-const) vector of @c std::array<uint,2> objects
    ///
    TRACCC_HOST_DEVICE
    auto& subspace() { return BASE::template get<4>(); }

    /// The subspace of measurements on the module
    ///
    /// The "subspace" defines which of the measurement's parameters are
    /// "sensitive", to be used during the track finding/fitting.
    ///
    /// @return A (const) vector of @c std::array<uint,2> objects
    ///
    TRACCC_HOST_DEVICE
    const auto& subspace() const { return BASE::template get<4>(); }

};  // class silicon_detector_description_interface

/// SoA container describing the detector module segmentation information
using detector_design_description = vecmem::edm::container<
    detector_design_description_interface, vecmem::edm::type::vector<int>,
    vecmem::edm::type::jagged_vector<scalar>,
    vecmem::edm::type::jagged_vector<scalar>,
    vecmem::edm::type::vector<unsigned char>,
    vecmem::edm::type::vector<
        std::array<detray::dindex_type<default_algebra>, 2u>>>;

}  // namespace traccc
