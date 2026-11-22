// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/shapes/cuboid3D.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/utils/bounding_volume.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <limits>
#include <vector>

namespace detray::detail {

/// @brief A functor to construct global bounding boxes around masks
///
/// Is called through the mask store visitor. Adds one axis aligned bounding
/// box in global coordinates for every mask in the given index range.
///
/// @tparam algebra_t the algebra type of the detector
template <concepts::algebra algebra_t>
struct bounding_box_creator {
  /// The bounding box type in global coordinates
  using aabb_t = axis_aligned_bounding_volume<cuboid3D, algebra_t>;

  /// Call operator for the mask store visitor
  ///
  /// @tparam mask_group_t the type of mask collection
  /// @tparam idx_range_t the type of the mask index or index range
  ///
  /// @param mask_group the collection of masks of the same type
  /// @param idx_range the index or index range of the masks of one surface
  /// @param envelope the envelope that is added around the mask boundaries
  /// @param trf the placement transform of the surface
  /// @param boxes the bounding boxes in global coordinates (output)
  template <typename mask_group_t, typename idx_range_t>
  DETRAY_HOST_DEVICE inline void operator()(const mask_group_t &mask_group,
                                            const idx_range_t &idx_range,
                                            const dscalar<algebra_t> envelope,
                                            const dtransform3D<algebra_t> &trf,
                                            std::vector<aabb_t> &boxes) const {
    for (const auto &mask : detray::ranges::subrange(mask_group, idx_range)) {
      // Local minimum bounding box
      aabb_t box{mask, boxes.size(), envelope};
      // Bounding box in global coordinates (might no longer be
      // minimum)
      boxes.push_back(box.transform(trf));
    }
  }
};

/// Construct the axis aligned bounding box of a volume in global coordinates
/// from its portals
///
/// @tparam detector_t the type of the detector
///
/// @param det the detector
/// @param vol_idx the index of the volume
/// @param envelope the envelope around the portal masks
///
/// @returns the bounding box of the volume (not necessarily minimal)
template <typename detector_t>
DETRAY_HOST inline auto volume_bounding_box(
    const detector_t &det, const dindex vol_idx,
    const dscalar<typename detector_t::algebra_type> envelope = std::
        numeric_limits<dscalar<typename detector_t::algebra_type>>::epsilon())
    -> axis_aligned_bounding_volume<cuboid3D,
                                    typename detector_t::algebra_type> {
  using algebra_t = typename detector_t::algebra_type;
  using aabb_t = axis_aligned_bounding_volume<cuboid3D, algebra_t>;

  const auto vol = tracking_volume{det, vol_idx};

  std::vector<aabb_t> boxes;
  for (const auto &pt_desc : vol.portals()) {
    det.mask_store().template visit<bounding_box_creator<algebra_t>>(
        pt_desc.mask(), envelope, det.transform_store().at(pt_desc.transform()),
        boxes);
  }

  return aabb_t{boxes, vol_idx, envelope};
}

}  // namespace detray::detail
