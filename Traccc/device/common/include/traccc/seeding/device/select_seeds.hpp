/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/device/global_index.hpp"
#include "traccc/edm/device/device_triplet.hpp"
#include "traccc/edm/device/triplet_counter.hpp"

// Project include(s).
#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/seed_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/seeding_config.hpp"
#include "traccc/seeding/detail/spacepoint_grid.hpp"
#include "traccc/seeding/detail/triplet.hpp"

namespace traccc::device {

/// Function used for selecting good triplets to be recorded into seed
/// collection
///
/// @param[in] globalIndex      The index of the current thread
/// @param[in] filter_config    Seed filter config
/// @param[in] spacepoints_view Collection of spacepoints
/// @param[in] sp_grid_view     The spacepoint grid
/// @param[in] spM_tc_view      Collection with the number of triplets per spM
/// @param[in] triplet_view     Collection of triplets
/// @param[in] data     Array for temporary storage of triplets for comparison
/// @param[out] seed_view       Collection of seeds
///
TRACCC_HOST_DEVICE
inline void select_seeds(
    global_index_t globalIndex, const seedfinder_config& finder_config,
    const seedfilter_config& filter_config,
    const edm::spacepoint_collection::const_view& spacepoints_view,
    const traccc::details::spacepoint_grid_types::const_view& sp_grid_view,
    const triplet_counter_spM_collection_types::const_view& spM_tc_view,
    const triplet_counter_collection_types::const_view& tc_view,
    const device_triplet_collection_types::const_view& triplet_view,
    triplet* data, edm::seed_collection::view seed_view);

}  // namespace traccc::device

// Include the implementation.
#include "traccc/seeding/device/impl/select_seeds.ipp"
