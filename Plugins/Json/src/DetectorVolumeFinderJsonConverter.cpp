// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Json/DetectorVolumeFinderJsonConverter.hpp"

#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdaters.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Plugins/Json/UtilitiesJsonConverter.hpp"
#include "Acts/Utilities/GridAxisGenerators.hpp"

#include <array>
#include <memory>
#include <tuple>
#include <vector>

namespace {
/// @brief  The generator struct
struct IndexedVolumesGenerator {
  using value_type = std::size_t;

  /// @brief  Helper function to create and connect the IndexedVolumesImpl
  ///
  /// @tparam grid_type the type of the grid, indicates also the dimension
  ///
  /// @param grid the grid object
  /// @param bv the bin value array
  /// @param transform the transform for the indexed volumes inmplementaiton
  ///
  /// @return a connected DetectorVolumeUpdater object
  template <typename grid_type>
  Acts::Experimental::DetectorVolumeUpdater createUpdater(
      grid_type&& grid,
      const std::array<Acts::BinningValue, grid_type::DIM>& bv,
      const Acts::Transform3& transform) {
    using IndexedDetectorVolumesImpl = Acts::Experimental::IndexedUpdaterImpl<
        grid_type, Acts::Experimental::IndexedDetectorVolumeExtractor,
        Acts::Experimental::DetectorVolumeFiller>;

    auto indexedDetectorVolumeImpl =
        std::make_unique<const IndexedDetectorVolumesImpl>(std::move(grid), bv,
                                                           transform);

    // Create the delegate and connect it
    Acts::Experimental::DetectorVolumeUpdater vFinder;
    vFinder.connect<&IndexedDetectorVolumesImpl::update>(
        std::move(indexedDetectorVolumeImpl));
    return vFinder;
  }
};

}  // namespace

Acts::Experimental::SurfaceCandidatesUpdater
Acts::DetectorVolumeFinderJsonConverter::fromJson(
    const nlohmann::json& jVolumeFinder) {
  // The return object
  auto vFinder = IndexedGridJsonHelper::generateFromJson<
      Experimental::DetectorVolumeUpdater, IndexedVolumesGenerator>(
      jVolumeFinder, "IndexedVolumes");
  if (vFinder.connected()) {
    return vFinder;
  }
  // Return
  return Experimental::tryRootVolumes();
}
