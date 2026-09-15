/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s).
#include "traccc/seeding/doublet_finding_helper.hpp"

// VecMem include(s).
#include <vecmem/containers/device_vector.hpp>

namespace traccc::device {

TRACCC_HOST_DEVICE
inline void make_mid_bot_lincircles(
    const global_index_t tid,
    device::device_doublet_collection_types::const_view mb_doublet_view,
    device::doublet_counter_collection_types::const_view doublet_count_view,
    edm::spacepoint_collection::const_view spacepoint_view,
    traccc::details::spacepoint_grid_types::const_view sp_grid_view,
    vecmem::data::vector_view<lin_circle> out_view) {
  const device::device_doublet_collection_types::const_device doublets(
      mb_doublet_view);
  const device::doublet_counter_collection_types::const_device doublet_counts(
      doublet_count_view);
  const edm::spacepoint_collection::const_device spacepoints(spacepoint_view);
  traccc::details::spacepoint_grid_types::const_device sp_grid(sp_grid_view);
  vecmem::device_vector<lin_circle> out(out_view);

  if (tid >= doublets.size()) {
    return;
  }

  const device::device_doublet dub = doublets.at(tid);
  const unsigned int counter_link = dub.counter_link;
  const device::doublet_counter count = doublet_counts.at(counter_link);
  const sp_location spM_loc = count.m_spM;
  const edm::spacepoint_collection::const_device::const_proxy_type spM =
      spacepoints.at(sp_grid.bin(spM_loc.bin_idx)[spM_loc.sp_idx]);
  const sp_location spB_loc = dub.sp2;
  const edm::spacepoint_collection::const_device::const_proxy_type spB =
      spacepoints.at(sp_grid.bin(spB_loc.bin_idx)[spB_loc.sp_idx]);

  out.at(tid) = doublet_finding_helper::transform_coordinates<
      traccc::details::spacepoint_type::bottom>(spM, spB);
}

}  // namespace traccc::device
