// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/WhiteBoardUtilities.hpp"
#include "Acts/Utilities/BinUtility.hpp"
#include "Acts/Utilities/Zip.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/Digitization/Smearers.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementReader.hpp"
#include "ActsExamples/Io/Csv/CsvMeasurementWriter.hpp"

#include <fstream>
#include <iostream>
#include <random>

using namespace ActsExamples;
using namespace Acts::Test;

BOOST_AUTO_TEST_CASE(CsvMeasurementRoundTrip) {
  IndexSourceLinkContainer sourceLinksOriginal;
  MeasurementContainer measOriginal;
  ClusterContainer clusterOriginal;
  IndexMultimap<Index> mapOriginal;

  ////////////////////////////
  // Create some dummy data //
  ////////////////////////////
  const std::size_t nMeasurements = 3;
  Acts::GeometryIdentifier someGeoId{298453};

  std::mt19937 gen(23);
  std::uniform_int_distribution<std::uint32_t> disti(1, 10);
  std::uniform_real_distribution<double> distf(0.0, 1.0);

  for (auto i = 0ul; i < nMeasurements; ++i) {
    IndexSourceLink sl(someGeoId, static_cast<Index>(i));
    sourceLinksOriginal.insert(sl);

    Acts::Vector2 p = Acts::Vector2::Random();
    Acts::SquareMatrix2 c = Acts::SquareMatrix2::Random();

    // NOTE this fails:
    // auto m = Acts::makeMeasurement(sl, p, c, eBoundLoc0, eBoundTime)
    // because we don't support non-consecutive parameters here for now
    auto m = Acts::makeMeasurement(Acts::SourceLink{sl}, p, c, Acts::eBoundLoc0,
                                   Acts::eBoundLoc1);

    measOriginal.push_back(m);

    ActsExamples::Cluster cl;

    using Bin2D = ActsFatras::Segmentizer::Bin2D;
    using Seg2D = ActsFatras::Segmentizer::Segment2D;

    // We have two cluster shapes which are displaced randomly
    const auto o = disti(gen);
    cl.channels.emplace_back(
        Bin2D{o + 0, o + 0},
        Seg2D{Acts::Vector2::Random(), Acts::Vector2::Random()}, distf(gen));
    cl.channels.emplace_back(
        Bin2D{o + 0, o + 1},
        Seg2D{Acts::Vector2::Random(), Acts::Vector2::Random()}, distf(gen));
    if (distf(gen) < 0.5) {
      cl.channels.emplace_back(
          Bin2D{o + 0, o + 2},
          Seg2D{Acts::Vector2::Random(), Acts::Vector2::Random()}, distf(gen));
      cl.sizeLoc0 = 1;
      cl.sizeLoc1 = 3;
    } else {
      cl.channels.emplace_back(
          Bin2D{o + 1, o + 0},
          Seg2D{Acts::Vector2::Random(), Acts::Vector2::Random()}, distf(gen));
      cl.sizeLoc0 = 2;
      cl.sizeLoc1 = 2;
    }

    clusterOriginal.push_back(cl);

    // Just generate some random hitid
    mapOriginal.insert(std::pair<Index, Index>{i, disti(gen)});
  }

  ///////////
  // Write //
  ///////////
  CsvMeasurementWriter::Config writerConfig;
  writerConfig.inputClusters = "clusters";
  writerConfig.inputMeasurements = "meas";
  writerConfig.inputMeasurementSimHitsMap = "map";
  writerConfig.outputDir = "";

  CsvMeasurementWriter writer(writerConfig, Acts::Logging::WARNING);

  auto writeTool =
      GenericReadWriteTool<>()
          .add(writerConfig.inputMeasurements, measOriginal)
          .add(writerConfig.inputClusters, clusterOriginal)
          .add(writerConfig.inputMeasurementSimHitsMap, mapOriginal);

  writeTool.write(writer);

  //////////////////
  // Write & Read //
  //////////////////
  CsvMeasurementReader::Config readerConfig;
  readerConfig.inputDir = writerConfig.outputDir;
  readerConfig.outputMeasurements = writerConfig.inputMeasurements;
  readerConfig.outputMeasurementSimHitsMap =
      writerConfig.inputMeasurementSimHitsMap;
  readerConfig.outputClusters = writerConfig.inputClusters;
  readerConfig.outputSourceLinks = "sourcelinks";

  CsvMeasurementReader reader(readerConfig, Acts::Logging::WARNING);

  auto readTool =
      writeTool.add(readerConfig.outputSourceLinks, sourceLinksOriginal);

  const auto [measRead, clusterRead, mapRead, sourceLinksRead] =
      readTool.read(reader);

  ///////////
  // Check //
  ///////////
  auto checkMeasurementClose = [](const auto &m1, const auto &m2) {
    constexpr auto SizeA = std::decay_t<decltype(m1)>::size();
    constexpr auto SizeB = std::decay_t<decltype(m2)>::size();
    if constexpr (SizeA == SizeB) {
      CHECK_CLOSE_REL(m1.parameters(), m2.parameters(), 1e-4);
    }
  };

  static_assert(
      std::is_same_v<std::decay_t<decltype(measRead)>, decltype(measOriginal)>);
  BOOST_REQUIRE(measRead.size() == measOriginal.size());
  for (const auto &[a, b] : Acts::zip(measRead, measOriginal)) {
    std::visit(checkMeasurementClose, a, b);
  }

  static_assert(std::is_same_v<std::decay_t<decltype(clusterRead)>,
                               decltype(clusterOriginal)>);
  BOOST_REQUIRE(clusterRead.size() == clusterOriginal.size());
  for (auto [a, b] : Acts::zip(clusterRead, clusterOriginal)) {
    BOOST_REQUIRE(a.sizeLoc0 == b.sizeLoc0);
    BOOST_REQUIRE(a.sizeLoc1 == b.sizeLoc1);

    for (const auto &ca : a.channels) {
      auto match = [&](const auto &cb) {
        return ca.bin == cb.bin &&
               std::abs(ca.activation - cb.activation) < 1.e-4;
      };

      BOOST_CHECK(std::any_of(b.channels.begin(), b.channels.end(), match));
    }
  }

  static_assert(
      std::is_same_v<std::decay_t<decltype(mapRead)>, decltype(mapOriginal)>);
  BOOST_REQUIRE(mapRead.size() == mapOriginal.size());
  for (const auto &[a, b] : Acts::zip(mapRead, mapOriginal)) {
    BOOST_REQUIRE(a == b);
  }

  static_assert(std::is_same_v<std::decay_t<decltype(sourceLinksRead)>,
                               decltype(sourceLinksOriginal)>);
  BOOST_REQUIRE(sourceLinksRead.size() == sourceLinksOriginal.size());
  for (const auto &[a, b] : Acts::zip(sourceLinksRead, sourceLinksOriginal)) {
    BOOST_REQUIRE(a.geometryId() == b.geometryId());
    BOOST_REQUIRE(a.index() == b.index());
  }
}
