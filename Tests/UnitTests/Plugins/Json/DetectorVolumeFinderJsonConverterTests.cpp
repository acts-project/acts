// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Navigation/DetectorVolumeFinders.hpp"
#include "Acts/Navigation/DetectorVolumeUpdators.hpp"
#include "Acts/Plugins/Json/DetectorVolumeFinderJsonConverter.hpp"

#include <fstream>
#include <memory>
#include <vector>

#include <nlohmann/json.hpp>

BOOST_AUTO_TEST_SUITE(DetectorVolumeFinderJsonConverter)

BOOST_AUTO_TEST_CASE(RzVolumes) {
  std::vector<Acts::ActsScalar> zBoundaries = {-1000., -500, 150.};
  std::vector<Acts::ActsScalar> rBoundaries = {0., 10., 30., 35.};

  using AxesGeneratorType =
      Acts::Experimental::detail::GridAxisGenerators::VarBoundVarBound;

  AxesGeneratorType zrAxes{zBoundaries, rBoundaries};

  // Create the grid with the provided axis generator
  using GridType = typename AxesGeneratorType::template grid_type<std::size_t>;
  GridType grid(zrAxes());

  grid.at(1u) = 0;
  grid.at(2u) = 1;
  grid.at(3u) = 2;
  grid.at(4u) = 3;
  grid.at(5u) = 4;
  grid.at(6u) = 5;

  auto casts = std::array<Acts::BinningValue, 2u>{Acts::binZ, Acts::binR};

  using IndexedDetectorVolumeImpl = Acts::Experimental::IndexedUpdatorImpl<
      GridType, Acts::Experimental::IndexedDetectorVolumeExtractor,
      Acts::Experimental::DetectorVolumeFiller>;

  auto indexedDetectorVolumeImpl =
      std::make_unique<const IndexedDetectorVolumeImpl>(std::move(grid), casts);

  // Return the root volume finder
  Acts::Experimental::DetectorVolumeUpdator rootVolumeFinder;
  rootVolumeFinder.connect<&IndexedDetectorVolumeImpl::update>(
      std::move(indexedDetectorVolumeImpl));

  nlohmann::json rFinderJson =
      Acts::DetectorVolumeFinderJsonConverter::toJson(rootVolumeFinder);

  auto readInRootVolumeFinder =
      Acts::DetectorVolumeFinderJsonConverter::fromJson(rFinderJson);

  BOOST_CHECK(readInRootVolumeFinder.instance() != nullptr);
}

BOOST_AUTO_TEST_SUITE_END()
