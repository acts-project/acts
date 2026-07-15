/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/singlet.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

namespace traccc {

inline std::pair<axis2::circular<>, axis2::regular<>> get_axes(
    const spacepoint_grid_config& grid_config, vecmem::memory_resource& mr) {

    detray::dindex phiBins;
    if (grid_config.bFieldInZ == 0) {
        phiBins = 100;
    } else {
        // calculate circle intersections of helix and max detector radius
        scalar minHelixRadius = grid_config.minPt / grid_config.bFieldInZ;

        // sanity check: if yOuter takes the square root of a negative number
        if (minHelixRadius < grid_config.rMax / 2) {
            throw std::domain_error(
                "The value of minHelixRadius cannot be smaller than rMax / 2. "
                "Please "
                "check the configuration of bFieldInZ and minPt");
        }
        scalar maxR2 = grid_config.rMax * grid_config.rMax;
        scalar xOuter = maxR2 / (2 * minHelixRadius);
        scalar yOuter = std::sqrt(maxR2 - xOuter * xOuter);
        scalar outerAngle = std::atan(xOuter / yOuter);

        // intersection of helix and max detector radius minus maximum R
        // distance from middle SP to top SP
        scalar innerAngle = 0;
        scalar rMin = grid_config.rMax;
        if (grid_config.rMax > grid_config.deltaRMax) {
            rMin = grid_config.rMax - grid_config.deltaRMax;
            scalar innerCircleR2 = (grid_config.rMax - grid_config.deltaRMax) *
                                   (grid_config.rMax - grid_config.deltaRMax);
            scalar xInner = innerCircleR2 / (2 * minHelixRadius);
            scalar yInner = std::sqrt(innerCircleR2 - xInner * xInner);
            innerAngle = std::atan(xInner / yInner);
        }

        // evaluating the azimutal deflection including the maximum impact
        // parameter
        scalar deltaAngleWithMaxD0 =
            math::fabs(std::asin(grid_config.impactMax / (rMin)) -
                       std::asin(grid_config.impactMax / grid_config.rMax));

        // evaluating delta Phi based on the inner and outer angle, and the
        // azimutal deflection including the maximum impact parameter Divide by
        // config.phiBinDeflectionCoverage since we combine
        // config.phiBinDeflectionCoverage number of consecutive phi bins in the
        // seed making step. So each individual bin should cover
        // 1/config.phiBinDeflectionCoverage of the maximum expected azimutal
        // deflection
        scalar deltaPhi =
            (outerAngle - innerAngle + deltaAngleWithMaxD0) /
            static_cast<scalar>(grid_config.phiBinDeflectionCoverage);

        // sanity check: if the delta phi is equal to or less than zero, we'll
        // be creating an infinite or a negative number of bins, which would be
        // bad!
        if (deltaPhi <= 0.) {
            throw std::domain_error(
                "Delta phi value is equal to or less than zero, leading to an "
                "impossible number of bins (negative or infinite)");
        }

        // divide 2pi by angle delta to get number of phi-bins
        // size is always 2pi even for regions of interest
        phiBins = static_cast<detray::dindex>(
            std::llround(2 * M_PI / deltaPhi + 0.5));
        // need to scale the number of phi bins accordingly to the number of
        // consecutive phi bins in the seed making step.
        // Each individual bin should be approximately a fraction (depending on
        // this number) of the maximum expected azimutal deflection.
    }

    axis2::circular m_phi_axis{phiBins, grid_config.phiMin, grid_config.phiMax,
                               mr};

    // TODO: can probably be optimized using smaller z bins
    // and returning (multiple) neighbors only in one z-direction for forward
    // seeds
    // FIXME: zBinSize must include scattering

    scalar zBinSize = grid_config.cotThetaMax * grid_config.deltaRMax;
    detray::dindex zBins = std::max(
        static_cast<detray::dindex>(1),
        static_cast<detray::dindex>(
            std::floor((grid_config.zMax - grid_config.zMin) / zBinSize)));

    axis2::regular m_z_axis{zBins, grid_config.zMin, grid_config.zMax, mr};

    return {m_phi_axis, m_z_axis};
}

template <typename T>
inline TRACCC_HOST_DEVICE bool is_valid_sp(const seedfinder_config& config,
                                           const edm::spacepoint<T>& sp) {

    if (sp.z() > config.zMax || sp.z() < config.zMin) {
        return false;
    }
    scalar spPhi = math::atan2(sp.y(), sp.x());
    if (spPhi > config.phiMax || spPhi < config.phiMin) {
        return false;
    }

    return (static_cast<size_t>(vector::perp(vector2{
                sp.x() - config.beamPos[0], sp.y() - config.beamPos[1]})) <
            config.get_num_rbins());
}

}  // namespace traccc
