// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/Utilities/CudaHoughTransformUtils.hpp"

#include "Acts/Seeding/HoughTransformUtils.hpp"

#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <vector>

namespace ActsTests {

namespace CudaHT = ActsExamples::CudaHoughTransformUtils;

BOOST_AUTO_TEST_SUITE(CudaHoughTransformUtilsSuite)

namespace {

struct DriftCircleInput {
  double y;
  double z;
  double r;
  double uncert;
};

std::vector<DriftCircleInput> driftCircleInputs() {
  constexpr double uncert = 0.3;

  return {
      {-427.981, -225.541, 14.5202, uncert},
      {-412.964, -199.530, 1.66237, uncert},
      {-427.981, -173.519, 12.3176, uncert},
      {-427.981, 173.519, 1.5412, uncert},
      {-442.999, 199.530, 12.3937, uncert},
      {-427.981, 225.541, 3.77967, uncert},
  };
}

ActsExamples::CudaMuonSpacePointContainer makeBatchedDriftCircleContainer(
    std::size_t nBuckets) {
  const std::vector<DriftCircleInput> driftCircles = driftCircleInputs();
  const std::size_t hitsPerBucket = driftCircles.size();

  ActsExamples::CudaMuonSpacePointContainer container{nBuckets * hitsPerBucket};

  for (std::size_t bucket = 0; bucket < nBuckets; ++bucket) {
    const double bucketYOffset = 10.0 * static_cast<double>(bucket);

    const std::size_t start = bucket * hitsPerBucket;
    const std::size_t end = start + hitsPerBucket;

    for (std::size_t local = 0; local < hitsPerBucket; ++local) {
      const std::size_t index = start + local;
      const auto& dc = driftCircles[local];

      container.setGeometryId(index, index);
      container.setId(index, 0u);

      container.defineCoordinates(index,
                                  Acts::Vector3{0.0, dc.y + bucketYOffset,
                                                dc.z},
                                  Acts::Vector3{1.0, 0.0, 0.0},
                                  Acts::Vector3{0.0, 1.0, 0.0});

      container.setRadius(index, dc.r);
      container.setTime(index, 0.0);
      container.setCovariance(index, dc.uncert * dc.uncert,
                              dc.uncert * dc.uncert, 0.0);

      container.setLogicalLayer(index, static_cast<std::uint32_t>(local));
    }

    container.addBucket(start, end);
  }

  return container;
}

// Utility to save data to CSNV  for pyhton visualization
void writeDriftCircleCsv(const std::filesystem::path& path) {
  const std::vector<DriftCircleInput> driftCircles = driftCircleInputs();

  std::ofstream out{path};
  BOOST_REQUIRE_MESSAGE(out, "Failed to open " << path.string());

  out << std::setprecision(17);
  out << "hitIndex,x,y,z,r,uncert,layer\n";

  for (std::size_t i = 0; i < driftCircles.size(); ++i) {
    const auto& dc = driftCircles[i];

    out << i << ","
        << 0.0 << ","
        << dc.y << ","
        << dc.z << ","
        << dc.r << ","
        << dc.uncert << ","
        << i << "\n";
  }
}

// Utility to save data to CSNV  for pyhton visualization
void writeHoughHistogramCsv(
    const std::filesystem::path& path,
    const CudaHT::CudaHoughPlaneBatch& plane,
    const Acts::HoughTransformUtils::HoughAxisRanges& axisRanges) {
  std::ofstream out{path};
  BOOST_REQUIRE_MESSAGE(out, "Failed to open " << path.string());

  out << std::setprecision(17);
  out << "xBin,yBin,tanTheta,interceptY,nHits,nLayers,layerMask\n";

  // Only one bucket that is 0th
  std::uint32_t bucketId = 0;

  for (std::size_t yBin = 0; yBin < plane.nBinsY(); ++yBin) {
    for (std::size_t xBin = 0; xBin < plane.nBinsX(); ++xBin) {
      const double tanTheta = Acts::HoughTransformUtils::binCenter(
          axisRanges.xMin, axisRanges.xMax, plane.nBinsX(), xBin);
      const double interceptY = Acts::HoughTransformUtils::binCenter(
          axisRanges.yMin, axisRanges.yMax, plane.nBinsY(), yBin);

      out << xBin << ","
          << yBin << ","
          << tanTheta << ","
          << interceptY << ","
          << plane.nHits(bucketId, xBin, yBin) << ","
          << plane.nLayers(bucketId, xBin, yBin) << ","
          << static_cast<unsigned long long>(plane.layerMask(bucketId, xBin, yBin))
          << "\n";
    }
  }
}

}  // namespace

// This checks the batched cell model.
//
// We fill a cell in bucket 1 with four hit contributions:
//   layers: 1, 2, 2, 4
//
// nHits counts all contributions and becomes 4.
// nLayers counts unique layers and becomes 3.
BOOST_AUTO_TEST_CASE(cuda_hough_batch_bit_mask_layer_counting) {
  CudaHT::CudaHoughPlaneBatch batch{{4, 4}, 3};

  batch.fillBin(1, 1, 2, 1, 1.0f);
  batch.fillBin(1, 1, 2, 2, 1.0f);
  batch.fillBin(1, 1, 2, 2, 1.0f);
  batch.fillBin(1, 1, 2, 4, 1.0f);

  BOOST_CHECK_EQUAL(batch.nHits(1, 1, 2), 4.0f);
  BOOST_CHECK_EQUAL(batch.nLayers(1, 1, 2), 3.0f);

  BOOST_CHECK(batch.hasLayer(1, 1, 2, 1));
  BOOST_CHECK(batch.hasLayer(1, 1, 2, 2));
  BOOST_CHECK(batch.hasLayer(1, 1, 2, 4));
  BOOST_CHECK(!batch.hasLayer(1, 1, 2, 3));

  BOOST_CHECK_EQUAL(batch.nHits(0, 1, 2), 0.0f);
  BOOST_CHECK_EQUAL(batch.nHits(2, 1, 2), 0.0f);
}

// This checks that the batched CUDA MDT eta fill agrees with the batched CPU
// reference implementation for all buckets.
BOOST_AUTO_TEST_CASE(cuda_hough_batch_eta_mdt_fill_matches_host_reference) {
  constexpr std::size_t nBuckets = 10;

  auto spacePointsForHost = makeBatchedDriftCircleContainer(nBuckets);
  auto spacePointsForCuda = makeBatchedDriftCircleContainer(nBuckets);

  Acts::HoughTransformUtils::HoughAxisRanges axisRanges{
      -0.20, 0.20, -650.0, -150.0};

  CudaHT::CudaHoughPlaneBatch hostBatch{{40, 40}, nBuckets};
  CudaHT::CudaHoughPlaneBatch cudaBatch{{40, 40}, nBuckets};

  hostBatch.fillEtaDriftCirclesHost(spacePointsForHost, axisRanges);

  cudaBatch.fillEtaDriftCirclesOnDevice(spacePointsForCuda, axisRanges,
                                        3.0, 1.0, 1.0f, 128);
  cudaBatch.moveToHost();

  for (std::size_t bucket = 0; bucket < nBuckets; ++bucket) {
    BOOST_CHECK_EQUAL(cudaBatch.maxHits(bucket), hostBatch.maxHits(bucket));
    BOOST_CHECK_EQUAL(cudaBatch.maxLayers(bucket), hostBatch.maxLayers(bucket));

    for (std::size_t x = 0; x < hostBatch.nBinsX(); ++x) {
      for (std::size_t y = 0; y < hostBatch.nBinsY(); ++y) {
        BOOST_CHECK_EQUAL(cudaBatch.nHits(bucket, x, y),
                          hostBatch.nHits(bucket, x, y));
        BOOST_CHECK_EQUAL(cudaBatch.nLayers(bucket, x, y),
                          hostBatch.nLayers(bucket, x, y));
        BOOST_CHECK_EQUAL(cudaBatch.layerMask(bucket, x, y),
                          hostBatch.layerMask(bucket, x, y));
      }
    }
  }
}


// This is the visual test, that exports cvs files for python script.
//
// The known simulated line from the original ACTS test is approximately:
//   tan(theta)  = -0.0401472 / 0.994974
//   interceptY  = -422.612
BOOST_AUTO_TEST_CASE(cuda_hough_eta_drift_circle_csv_visual_example) {
  auto spacePoints = makeBatchedDriftCircleContainer(1);

  const double expectedTanTheta = -0.0401472 / 0.994974;
  const double expectedInterceptY = -422.612;

  // Only one bucket that is 0th
  std::uint32_t bucketId = 0;

  // The axis has to be relatively small, as plane is very sparse.
  Acts::HoughTransformUtils::HoughAxisRanges axisRanges{
      -0.1, 0.0, -500.0, -400.0};

  CudaHT::CudaHoughPlaneBatch plane{{15, 15}, 1};

  plane.fillEtaDriftCirclesOnDevice(spacePoints, axisRanges);
  plane.moveToHost();

  const auto [xBin, yBin] = plane.locMaxHits(bucketId);

  const double foundTanTheta = Acts::HoughTransformUtils::binCenter(
      axisRanges.xMin, axisRanges.xMax, plane.nBinsX(), xBin);
  const double foundInterceptY = Acts::HoughTransformUtils::binCenter(
      axisRanges.yMin, axisRanges.yMax, plane.nBinsY(), yBin);

  const std::filesystem::path outDir = std::filesystem::current_path();
  const std::filesystem::path hitsCsv = outDir / "cuda_hough_visual_hits.csv";
  const std::filesystem::path histCsv =
      outDir / "cuda_hough_visual_histogram.csv";

  writeDriftCircleCsv(hitsCsv);
  writeHoughHistogramCsv(histCsv, plane, axisRanges);

  BOOST_TEST_MESSAGE("Wrote Hough visual debug hits to: "
                     << hitsCsv.string());
  BOOST_TEST_MESSAGE("Wrote Hough visual debug histogram to: "
                     << histCsv.string());

  BOOST_TEST_MESSAGE("Maximum bin: x=" << xBin << ", y=" << yBin);
  BOOST_TEST_MESSAGE("Maximum parameters: tanTheta=" << foundTanTheta
                                                     << ", interceptY="
                                                     << foundInterceptY);
  BOOST_TEST_MESSAGE("Max hits: " << plane.maxHits(bucketId));
  BOOST_TEST_MESSAGE("Max layers: " << plane.maxLayers(0));

  BOOST_CHECK_CLOSE(foundTanTheta, expectedTanTheta, 50.0);
  BOOST_CHECK_CLOSE(foundInterceptY, expectedInterceptY, 10.0);

  BOOST_CHECK_GE(plane.maxHits(bucketId), 3.0f);
  BOOST_CHECK_GE(plane.maxLayers(bucketId), 3.0f);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
