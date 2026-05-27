/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/doublet.hpp"
#include "traccc/seeding/detail/lin_circle.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/triplet.hpp"
#include "traccc/seeding/triplet_finding_helper.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/memory/memory_resource.hpp>

// System include(s).
#include <cassert>
#include <functional>

namespace traccc::host::details {

/// Triplet finding to search the compatible combintations of two doublets which
/// share same middle spacepoint
struct triplet_finding : public messaging {

    /// Constructor for the triplet finding
    ///
    /// @param finder_config Seed finding configuration parameters
    /// @param filter_config Seed filtering configuration parameters
    /// @param mr The memory resource to use
    ///
    triplet_finding(
        const seedfinder_config& finding_config,
        const seedfilter_config& filter_config, vecmem::memory_resource& mr,
        std::unique_ptr<const Logger> logger = getDummyLogger().clone())
        : messaging(std::move(logger)),
          m_finding_config{finding_config},
          m_filter_config{filter_config},
          m_mr{mr} {}

    /// Callable operator for triplet finding per middle-bottom doublet
    ///
    /// @param spacepoints All the spacepoints in the event
    /// @param sp_grid The spacepoint grid to use
    /// @param mid_bot_doublet is the current middle-bottom doublets
    /// @param mid_bot_lc is transformed coordinate of @c mid_bot_doublet
    /// @param mid_top_doublets is the vector of middle-top doublets which share
    ///                         same middle spacepoint with current
    ///                         middle-bottom doublet
    /// @param mid_top_lcs is transformed coordinates of
    ///                    @c doublets_mid_top
    ///
    /// @return a vector of triplets
    triplet_collection_types::host operator()(
        const edm::spacepoint_collection::const_device& spacepoints,
        const traccc::details::spacepoint_grid_types::host& sp_grid,
        const doublet& mid_bot_doublet, const lin_circle& mid_bot_lc,
        const doublet_collection_types::host& mid_top_doublets,
        const lin_circle_collection_types::host& mid_top_lcs) const {

        // Create the output.
        triplet_collection_types::host result{&(m_mr.get())};

        // Access the middle spacepoint that all the doublets share.
        const edm::spacepoint_collection::const_device::const_proxy_type spM =
            spacepoints.at(sp_grid.bin(
                mid_bot_doublet.sp1.bin_idx)[mid_bot_doublet.sp1.sp_idx]);

        // Calculate quantities that help deciding if two doublets are
        // compatible.
        const scalar iSinTheta2 =
            1.f + mid_bot_lc.cotTheta() * mid_bot_lc.cotTheta();
        scalar scatteringInRegion2 =
            m_finding_config.maxScatteringAngle2 * iSinTheta2;
        scatteringInRegion2 *=
            m_finding_config.sigmaScattering * m_finding_config.sigmaScattering;
        scalar curvature, impact_parameter;

        // Compare this mid-bottom doublet with all the mid-top doublets.
        assert(mid_top_doublets.size() == mid_top_lcs.size());
        for (std::size_t i = 0; i < mid_top_doublets.size(); ++i) {

            const doublet& mid_top_doublet = mid_top_doublets.at(i);
            const lin_circle& mid_top_lc = mid_top_lcs.at(i);

            if (!triplet_finding_helper::isCompatible(
                    spM, mid_bot_lc, mid_top_lc, m_finding_config, iSinTheta2,
                    scatteringInRegion2, curvature, impact_parameter)) {
                continue;
            }

            result.push_back(
                {mid_bot_doublet.sp2,  // bottom
                 mid_bot_doublet.sp1,  // middle
                 mid_top_doublet.sp2,  // top
                 curvature,            // curvature
                 -impact_parameter * m_filter_config.impactWeightFactor,
                 mid_bot_lc.Zo()});
        }

        // Set the triplet weights in a super complicated double loop over the
        // triplets found in the previous step.
        for (std::size_t i = 0; i < result.size(); ++i) {

            triplet& current_triplet = result[i];
            const sp_location& current_spT_idx = current_triplet.sp3;
            const edm::spacepoint_collection::const_device::const_proxy_type
                current_spT = spacepoints.at(sp_grid.bin(
                    current_spT_idx.bin_idx)[current_spT_idx.sp_idx]);
            const scalar currentTop_r = current_spT.radius();

            // if two compatible seeds with high distance in r are found,
            // compatible seeds span 5 layers
            // -> very good seed
            std::vector<scalar> compatibleSeedR;
            scalar lowerLimitCurv = current_triplet.curvature -
                                    m_filter_config.deltaInvHelixDiameter;
            scalar upperLimitCurv = current_triplet.curvature +
                                    m_filter_config.deltaInvHelixDiameter;

            for (std::size_t j = 0; j < result.size(); ++j) {

                if (i == j) {
                    continue;
                }

                const triplet& other_triplet = result[j];
                const sp_location& other_spT_idx = other_triplet.sp3;
                const edm::spacepoint_collection::const_device::const_proxy_type
                    other_spT = spacepoints.at(sp_grid.bin(
                        other_spT_idx.bin_idx)[other_spT_idx.sp_idx]);

                // compared top SP should have at least deltaRMin distance
                const scalar otherTop_r = other_spT.radius();
                const scalar deltaR = currentTop_r - otherTop_r;
                if (std::abs(deltaR) < m_filter_config.deltaRMin) {
                    continue;
                }

                // curvature difference within limits?
                // TODO: how much slower than sorting all vectors by curvature
                // and breaking out of loop? i.e. is vector size large (e.g. in
                // jets?)
                if (other_triplet.curvature < lowerLimitCurv) {
                    continue;
                }
                if (other_triplet.curvature > upperLimitCurv) {
                    continue;
                }

                bool newCompSeed = true;
                for (scalar previousDiameter : compatibleSeedR) {
                    // original ATLAS code uses higher min distance for 2nd
                    // found compatible seed (20mm instead of 5mm) add new
                    // compatible seed only if distance larger than rmin to all
                    // other compatible seeds
                    if (std::abs(previousDiameter - otherTop_r) <
                        m_filter_config.deltaRMin) {
                        newCompSeed = false;
                        break;
                    }
                }

                if (newCompSeed) {
                    compatibleSeedR.push_back(otherTop_r);
                    current_triplet.weight += m_filter_config.compatSeedWeight;
                }

                if (compatibleSeedR.size() >= m_filter_config.compatSeedLimit) {
                    break;
                }
            }
        }

        // Return the reconstructed triplets.
        return result;
    }

    private:
    /// @name Triplet Finding Configuration
    /// @{
    seedfinder_config m_finding_config;
    seedfilter_config m_filter_config;
    /// @}

    /// Memory resource to use
    std::reference_wrapper<vecmem::memory_resource> m_mr;

};  // struct triplet_finding

}  // namespace traccc::host::details
