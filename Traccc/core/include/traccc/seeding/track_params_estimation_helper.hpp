/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Library include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/spacepoint_helpers.hpp"
#include "traccc/edm/track_parameters.hpp"

// System include(s).
#include <cassert>
#include <cmath>

namespace traccc {

/// helper functions (for both cpu and gpu) to perform conformal transformation
///
/// @param x is the x value
/// @param y is the y value
/// @return is the conformal transformation result
inline TRACCC_HOST_DEVICE vector2 uv_transform(const scalar& x,
                                               const scalar& y) {
    vector2 uv;
    scalar denominator = x * x + y * y;
    uv[0] = x / denominator;
    uv[1] = y / denominator;
    return uv;
}

/// helper functions (for both cpu and gpu) to calculate bound track parameter
/// at the bottom spacepoint
///
/// @param [out] params the bound track parameter vector to be filled
/// @param [in] measurements is the measurement collection
/// @param [in] spacepoints is the spacepoint collection
/// @param [in] seed is the input seed
/// @param [in] bfield is the magnetic field
///
template <typename T>
inline TRACCC_HOST_DEVICE void seed_to_bound_param_vector(
    bound_track_parameters<>& params,
    const edm::measurement_collection::const_device& measurements,
    const edm::spacepoint_collection::const_device& spacepoints,
    const edm::seed<T>& seed, const vector3& bfield) {

    const edm::spacepoint_collection::const_device::const_proxy_type spB =
        spacepoints.at(seed.bottom_index());
    const edm::spacepoint_collection::const_device::const_proxy_type spM =
        spacepoints.at(seed.middle_index());
    const edm::spacepoint_collection::const_device::const_proxy_type spT =
        spacepoints.at(seed.top_index());

    std::array<vector3, 3> sp_global_positions{
        edm::get_spacepoint_global<default_algebra>(spB),
        edm::get_spacepoint_global<default_algebra>(spM),
        edm::get_spacepoint_global<default_algebra>(spT)};

    // Define a new coordinate frame with its origin at the bottom space
    // point, z axis long the magnetic field direction and y axis
    // perpendicular to vector from the bottom to middle space point.
    // Hence, the projection of the middle space point on the tranverse
    // plane will be located at the x axis of the new frame.
    vector3 relVec = sp_global_positions[1] - sp_global_positions[0];
    vector3 newZAxis = vector::normalize(bfield);
    vector3 newYAxis = vector::normalize(vector::cross(newZAxis, relVec));
    vector3 newXAxis = vector::cross(newYAxis, newZAxis);

    // The center of the new frame is at the bottom space point
    const vector3& translation = sp_global_positions[0];

    transform3 trans(translation, newXAxis, newYAxis, newZAxis);

    // The coordinate of the middle and top space point in the new frame
    const point3 local1 = trans.point_to_local(sp_global_positions[1]);
    const point3 local2 = trans.point_to_local(sp_global_positions[2]);

    // The uv1.y() should be zero
    vector2 uv1 = uv_transform(local1[0], local1[1]);
    vector2 uv2 = uv_transform(local2[0], local2[1]);

    // A,B are slope and intercept of the straight line in the u,v plane
    // connecting the three points
    scalar A = (uv2[1] - uv1[1]) / (uv2[0] - uv1[0]);
    scalar B = uv2[1] - A * uv2[0];

    // Radius (with a sign)
    scalar R = -vector::perp(vector2{1.f, A}) / (2.f * B);
    // The (1/tanTheta) of momentum in the new frame
    scalar invTanTheta =
        local2[2] / (2.f * R * math::asin(vector::perp(local2) / (2.f * R)));

    // The momentum direction in the new frame (the center of the circle
    // has the coordinate (-1.*A/(2*B), 1./(2*B)))
    vector3 transDirection =
        vector3({1.f, A, scalar(vector::perp(vector2{1.f, A})) * invTanTheta});
    // Transform it back to the original frame
    vector3 direction =
        transform3::rotate(trans._data, vector::normalize(transDirection));

    // The estimated phi and theta
    params.set_phi(vector::phi(direction));
    params.set_theta(vector::theta(direction));

    // The measured loc0 and loc1
    assert(spB.measurement_index_2() ==
           edm::spacepoint_collection::device::INVALID_MEASUREMENT_INDEX);
    const edm::measurement meas_for_spB =
        measurements.at(spB.measurement_index_1());
    params.set_surface_link(meas_for_spB.surface_link());
    params.set_bound_local(
        {meas_for_spB.local_position()[0], meas_for_spB.local_position()[1]});

    // The estimated q/pt in [GeV/c]^-1 (note that the pt is the
    // projection of momentum on the transverse plane of the new frame)
    scalar qOverPt = 1.f / (R * vector::norm(bfield));
    // The estimated q/p in [GeV/c]^-1
    params.set_qop(qOverPt / vector::perp(vector2{1.f, invTanTheta}));

    // Make sure the time is a finite value
    assert(std::isfinite(params.time()));
}

}  // namespace traccc
