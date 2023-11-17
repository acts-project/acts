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

  grid.at(1u) = 0u;
  grid.at(2u) = 1u;
  grid.at(3u) = 2u;
  grid.at(4u) = 3u;
  grid.at(5u) = 4u;
  grid.at(6u) = 5u;

  auto casts = std::array<Acts::BinningValue, 2u>{Acts::binZ, Acts::binR};

  using IndexedDetectorVolumesImpl = Acts::Experimental::IndexedUpdatorImpl<
      GridType, Acts::Experimental::IndexedDetectorVolumeExtractor,
      Acts::Experimental::DetectorVolumeFiller>;

  auto indexedDetectorVolumesImpl =
      std::make_unique<const IndexedDetectorVolumesImpl>(std::move(grid),
                                                         casts);

  // Return the root volume finder
  Acts::Experimental::DetectorVolumeUpdator rootVolumeFinder;
  rootVolumeFinder.connect<&IndexedDetectorVolumesImpl::update>(
      std::move(indexedDetectorVolumesImpl));

  nlohmann::json rFinderJson =
      Acts::DetectorVolumeFinderJsonConverter::toJson(rootVolumeFinder);

  auto readInRootVolumeFinder =
      Acts::DetectorVolumeFinderJsonConverter::fromJson(rFinderJson);

  BOOST_CHECK(readInRootVolumeFinder.instance() != nullptr);

  auto readInIndexedDetectorVolumesImpl =
      dynamic_cast<const IndexedDetectorVolumesImpl*>(
          readInRootVolumeFinder.instance());

  BOOST_CHECK(readInIndexedDetectorVolumesImpl != nullptr);

  const auto& readInGrid = readInIndexedDetectorVolumesImpl->grid;

  BOOST_CHECK_EQUAL(readInGrid.at(1u), 0u);
  BOOST_CHECK_EQUAL(readInGrid.at(2u), 1u);
  BOOST_CHECK_EQUAL(readInGrid.at(3u), 2u);
  BOOST_CHECK_EQUAL(readInGrid.at(4u), 3u);
  BOOST_CHECK_EQUAL(readInGrid.at(5u), 4u);
  BOOST_CHECK_EQUAL(readInGrid.at(6u), 5u);
}

BOOST_AUTO_TEST_SUITE_END()