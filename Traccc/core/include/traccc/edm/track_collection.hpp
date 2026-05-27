/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/track_constituent_link.hpp"
#include "traccc/edm/track_fit_outcome.hpp"
#include "traccc/edm/track_parameters.hpp"

// Detray include(s).
#include <detray/definitions/algebra.hpp>

// VecMem include(s).
#include <vecmem/edm/container.hpp>

namespace traccc::edm {

/// Interface for the @c traccc::edm::track_collection type.
///
/// It provides the API that users would interact with, while using the
/// columns/arrays of the SoA containers, or the variables of the AoS proxies
/// created on top of the SoA containers.
///
template <typename BASE>
class track : public BASE {

    public:
    /// @name Functions inherited from the base class
    /// @{

    /// Inherit the base class's constructor(s)
    using BASE::BASE;
    /// Inherit the base class's assignment operator(s).
    using BASE::operator=;

    /// @}

    /// @name Track Information
    /// @{

    /// The outcome of the track fit (non-const)
    ///
    /// @return A (non-const) vector of @c traccc::track_fit_outcome
    ///
    TRACCC_HOST_DEVICE
    auto& fit_outcome() { return BASE::template get<0>(); }
    /// The outcome of the track fit (non-const)
    ///
    /// @return A (const) vector of @c traccc::track_fit_outcome
    ///
    TRACCC_HOST_DEVICE
    const auto& fit_outcome() const { return BASE::template get<0>(); }

    /// The parameters of the track (non-const)
    ///
    /// @return A (non-const) vector of bound track parameters
    ///
    TRACCC_HOST_DEVICE
    auto& params() { return BASE::template get<1>(); }
    /// The parameters of the track (const)
    ///
    /// @return A (const) vector of bound track parameters
    ///
    TRACCC_HOST_DEVICE
    const auto& params() const { return BASE::template get<1>(); }

    /// The number of degrees of freedom of the track fit (non-const)
    ///
    /// @return A (non-const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    auto& ndf() { return BASE::template get<2>(); }
    /// The number of degrees of freedom of the track fit (const)
    ///
    /// @return A (const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    const auto& ndf() const { return BASE::template get<2>(); }

    /// The chi square of the track fit (non-const)
    ///
    /// @return A (non-const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    auto& chi2() { return BASE::template get<3>(); }
    /// The chi square of the track fit (const)
    ///
    /// @return A (const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    const auto& chi2() const { return BASE::template get<3>(); }

    /// The p-value of the track fit (non-const)
    ///
    /// @return A (non-const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    auto& pval() { return BASE::template get<4>(); }
    /// The p-value of the track fit (const)
    ///
    /// @return A (const) vector of scalar values
    ///
    TRACCC_HOST_DEVICE
    const auto& pval() const { return BASE::template get<4>(); }

    /// The number of holes in the track pattern (non-const)
    ///
    /// @return A (non-const) vector of unsigned integers
    ///
    TRACCC_HOST_DEVICE
    auto& nholes() { return BASE::template get<5>(); }
    /// The number of holes in the track pattern (const)
    ///
    /// @return A (const) vector of unsigned integers
    ///
    TRACCC_HOST_DEVICE
    const auto& nholes() const { return BASE::template get<5>(); }

    /// Links to the constituents associated with the track (non-const)
    ///
    /// @return A (non-const) jagged vector of links
    ///
    TRACCC_HOST_DEVICE
    auto& constituent_links() { return BASE::template get<6>(); }
    /// Links to the constituents associated with the track (const)
    ///
    /// @return A (const) jagged vector of links
    ///
    TRACCC_HOST_DEVICE
    const auto& constituent_links() const { return BASE::template get<6>(); }

    /// @}

    /// @name Utility functions
    /// @{

    /// Reset the fit quality variables
    TRACCC_HOST_DEVICE
    void reset_quality();

    /// Equality operator
    ///
    /// @note This function must only be used on proxy objects, not on
    ///       containers!
    ///
    /// @param[in] other The object to compare with
    /// @return @c true if the objects are equal, @c false otherwise
    ///
    template <typename T>
    TRACCC_HOST_DEVICE bool operator==(const track<T>& other) const;

    /// @}

};  // class track_fit

/// SoA container describing the fitted tracks
///
/// @tparam ALGEBRA The algebra type used to describe the tracks
///
template <detray::concepts::algebra ALGEBRA>
using track_collection = vecmem::edm::container<
    track,
    // fit_outcome
    vecmem::edm::type::vector<track_fit_outcome>,
    // params
    vecmem::edm::type::vector<bound_track_parameters<ALGEBRA>>,
    // ndf
    vecmem::edm::type::vector<detray::dscalar<ALGEBRA>>,
    // chi2
    vecmem::edm::type::vector<detray::dscalar<ALGEBRA>>,
    // pval
    vecmem::edm::type::vector<detray::dscalar<ALGEBRA>>,
    // nholes
    vecmem::edm::type::vector<unsigned int>,
    // constituent_links
    vecmem::edm::type::jagged_vector<track_constituent_link>>;

}  // namespace traccc::edm

// Include the implementation.
#include "traccc/edm/impl/track_collection.ipp"
