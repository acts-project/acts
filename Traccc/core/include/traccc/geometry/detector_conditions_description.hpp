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

/// Interface for the @c traccc::detector_conditions_description class.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays defined in @c traccc::detector_conditions_description.
///
template <typename BASE>
class detector_conditions_description_interface : public BASE {

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
    auto& module_to_design_id() { return BASE::template get<0>(); }
    /// The identifier of the detector module's surface (const)
    ///
    /// Can be used to look up the module in a @c detray::detector object.
    ///
    /// @return A (const) vector of @c detray::geometry::identifier objects
    ///
    TRACCC_HOST_DEVICE
    const auto& module_to_design_id() const { return BASE::template get<0>(); }

    TRACCC_HOST_DEVICE
    auto& geometry_id() { return BASE::template get<1>(); }
    /// The identifier of the detector module's surface (const)
    ///
    /// Can be used to look up the module in a @c detray::detector object.
    ///
    /// @return A (const) vector of @c detray::geometry::identifier objects
    ///
    TRACCC_HOST_DEVICE
    const auto& geometry_id() const { return BASE::template get<1>(); }

    TRACCC_HOST_DEVICE
    auto& acts_geometry_id() { return BASE::template get<2>(); }
    /// The identifier of the detector module's surface (const)
    ///
    /// Can be used to look up the module in a @c detray::detector object.
    ///
    /// @return A (const) vector of @c detray::geometry::identifier objects
    ///
    TRACCC_HOST_DEVICE
    const auto& acts_geometry_id() const { return BASE::template get<2>(); }

    /// The local translation vector to model e.g. Lorentz shifts
    ///
    /// @return A vector by which to translate the measurement in the local
    /// coordinate frame.
    ///
    /// @{
    TRACCC_HOST_DEVICE
    auto& measurement_translation() { return BASE::template get<3>(); }

    TRACCC_HOST_DEVICE
    const auto& measurement_translation() const {
        return BASE::template get<3>();
    }

};  // class silicon_detector_description_interface

/// SoA container describing module to design map and conditions
/// (module specific) data
using detector_conditions_description = vecmem::edm::container<
    detector_conditions_description_interface,
    vecmem::edm::type::vector<unsigned int>,
    vecmem::edm::type::vector<detray::geometry::identifier>,
    vecmem::edm::type::vector<geometry_id>, vecmem::edm::type::vector<vector2>>;

}  // namespace traccc
