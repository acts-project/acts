// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Seeding/HoughTransformUtils.hpp"
#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/Root/RootCudaMuonSpacePointReader.hpp"
#include "ActsExamples/Utilities/CudaHoughTransformUtils.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
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

      container.defineCoordinates(
          index, Acts::Vector3{0.0, dc.y + bucketYOffset, dc.z},
          Acts::Vector3{1.0, 0.0, 0.0}, Acts::Vector3{0.0, 1.0, 0.0});

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

// Utility to save data to CSNV  for python visualization
void writeDriftCircleCsv(const std::filesystem::path& path) {
  const std::vector<DriftCircleInput> driftCircles = driftCircleInputs();

  std::ofstream out{path};
  BOOST_REQUIRE_MESSAGE(out, "Failed to open " << path.string());

  out << std::setprecision(17);
  out << "hitIndex,x,y,z,r,uncert,layer\n";

  for (std::size_t i = 0; i < driftCircles.size(); ++i) {
    const auto& dc = driftCircles[i];

    out << i << "," << 0.0 << "," << dc.y << "," << dc.z << "," << dc.r << ","
        << dc.uncert << "," << i << "\n";
  }
}

// Utility to save data to CSNV  for python visualization
void writeHoughHistogramCsv(
    const std::filesystem::path& path, const CudaHT::CudaHoughPlaneBatch& plane,
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

      out << xBin << "," << yBin << "," << tanTheta << "," << interceptY << ","
          << plane.nHits(bucketId, xBin, yBin) << ","
          << plane.nLayers(bucketId, xBin, yBin) << ","
          << static_cast<unsigned long long>(
                 plane.layerMask(bucketId, xBin, yBin))
          << "\n";
    }
  }
}

std::uint32_t rawMuonIdLayer(std::uint32_t rawId) {
  static constexpr std::uint32_t fourBit = 0xFu;
  static constexpr std::uint32_t layerShift = 17u;

  return (rawId >> layerShift) & fourBit;
}

void writeFirstBucketHitsCsv(
    const std::filesystem::path& path,
    const ActsExamples::CudaMuonSpacePointContainer& container) {
  std::ofstream out{path};
  BOOST_REQUIRE_MESSAGE(out, "Failed to open " << path.string());

  BOOST_REQUIRE_GT(container.bucketCount(), 0u);

  const std::size_t bucketId = 0;
  const std::size_t start = container.bucketStart(bucketId);
  const std::size_t end = container.bucketEnd(bucketId);

  out << std::setprecision(17);
  out << "hitIndex,x,y,z,r,uncert,layer\n";

  for (std::size_t i = start; i < end; ++i) {
    auto sp = container[i];

    const Acts::Vector3& pos = sp->localPosition();
    const std::array<double, 3>& cov = sp->covariance();

    const double uncert = std::sqrt(std::max(cov[1], 0.0));

    out << (i - start) << "," << pos.x() << "," << pos.y() << "," << pos.z()
        << "," << sp->driftRadius() << "," << uncert << ","
        << rawMuonIdLayer(container.muonId(i)) << "\n";
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

  Acts::HoughTransformUtils::HoughAxisRanges axisRanges{-0.20, 0.20, -650.0,
                                                        -150.0};

  CudaHT::CudaHoughPlaneBatch hostBatch{{40, 40}, nBuckets};
  CudaHT::CudaHoughPlaneBatch cudaBatch{{40, 40}, nBuckets};

  hostBatch.fillEtaDriftCirclesHost(spacePointsForHost, axisRanges);

  cudaBatch.fillEtaDriftCirclesOnDevice(spacePointsForCuda, axisRanges, 3.0,
                                        1.0, 1.0f, 128);
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
  Acts::HoughTransformUtils::HoughAxisRanges axisRanges{-0.1, 0.0, -500.0,
                                                        -400.0};

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

  BOOST_TEST_MESSAGE("Wrote Hough visual debug hits to: " << hitsCsv.string());
  BOOST_TEST_MESSAGE(
      "Wrote Hough visual debug histogram to: " << histCsv.string());

  BOOST_TEST_MESSAGE("Maximum bin: x=" << xBin << ", y=" << yBin);
  BOOST_TEST_MESSAGE("Maximum parameters: tanTheta="
                     << foundTanTheta << ", interceptY=" << foundInterceptY);
  BOOST_TEST_MESSAGE("Max hits: " << plane.maxHits(bucketId));
  BOOST_TEST_MESSAGE("Max layers: " << plane.maxLayers(0));

  BOOST_CHECK_CLOSE(foundTanTheta, expectedTanTheta, 50.0);
  BOOST_CHECK_CLOSE(foundInterceptY, expectedInterceptY, 10.0);

  BOOST_CHECK_GE(plane.maxHits(bucketId), 3.0f);
  BOOST_CHECK_GE(plane.maxLayers(bucketId), 3.0f);
}

// Visual test with realistic data
// Bucket 0 exported for visualization
BOOST_AUTO_TEST_CASE(cuda_hough_eta_mdt_root_first_event_whole_event_sanity) {
  const char* envPath = std::getenv("ACTS_CUDA_MUON_SP_ROOT");

  const std::filesystem::path filePath = std::filesystem::path{
      "/data/mgawlas/ACTS/acts/Tests/Data/ParticleGun_MU0.root"};

  if (!std::filesystem::exists(filePath)) {
    BOOST_TEST_MESSAGE(
        "Skipping ROOT CUDA Hough test. File does not exist: " << filePath);
    return;
  }

  ActsExamples::RootCudaMuonSpacePointReader::Config cfg{};
  cfg.filePath = filePath.string();
  cfg.treeName = "MuonSpacePoints";
  cfg.outputSpacePoints = "CudaMuonSpacePoints";

  ActsExamples::RootCudaMuonSpacePointReader reader{cfg, Acts::Logging::INFO};

  const auto [firstEvent, lastEvent] = reader.availableEvents();

  BOOST_CHECK_EQUAL(firstEvent, 0u);
  BOOST_REQUIRE_GT(lastEvent, 0u);

  ActsExamples::WhiteBoard eventStore;
  ActsExamples::AlgorithmContext context{0, 0, eventStore, 0};

  BOOST_CHECK(reader.read(context) == ActsExamples::ProcessCode::SUCCESS);
  BOOST_REQUIRE(eventStore.exists(cfg.outputSpacePoints));

  ActsExamples::ReadDataHandle<ActsExamples::CudaMuonSpacePointContainer>
      outputHandle{&reader, "OutputSpacePoints"};
  outputHandle.initialize(cfg.outputSpacePoints);

  const auto& constContainer = outputHandle(context);

  BOOST_REQUIRE_GT(constContainer.size(), 0u);
  BOOST_REQUIRE_GT(constContainer.bucketCount(), 0u);

  // The CUDA fill moves the container to device, so it needs a non-const
  // object.
  auto& spacePoints =
      const_cast<ActsExamples::CudaMuonSpacePointContainer&>(constContainer);

  for (std::size_t bucket = 0; bucket < spacePoints.bucketCount(); ++bucket) {
    BOOST_CHECK_LE(spacePoints.bucketStart(bucket),
                   spacePoints.bucketEnd(bucket));
    BOOST_CHECK_LE(spacePoints.bucketEnd(bucket), spacePoints.size());
  }

  const Acts::HoughTransformUtils::HoughAxisRanges axisRanges{-1.00, 1.20,
                                                              -550.0, 450.0};

  CudaHT::CudaHoughPlaneBatch plane{{15, 15}, spacePoints.bucketCount()};

  plane.fillEtaDriftCirclesOnDevice(spacePoints, axisRanges);
  plane.moveToHost();

  std::size_t bucketsWithHits = 0;
  std::size_t totalNonEmptyCells = 0;

  std::size_t bestBucket = 0;
  CudaHT::YieldType bestMaxHits = 0.0f;
  CudaHT::YieldType bestMaxLayers = 0.0f;

  for (std::size_t bucket = 0; bucket < plane.nBuckets(); ++bucket) {
    const auto maxHits = plane.maxHits(bucket);
    const auto maxLayers = plane.maxLayers(bucket);

    BOOST_CHECK(std::isfinite(maxHits));
    BOOST_CHECK(std::isfinite(maxLayers));
    BOOST_CHECK_LE(maxLayers, maxHits);

    const auto nonEmpty = plane.nonEmptyBins(bucket);
    totalNonEmptyCells += nonEmpty.size();

    if (maxHits > 0.0f) {
      ++bucketsWithHits;
    }

    if (maxHits > bestMaxHits) {
      bestMaxHits = maxHits;
      bestMaxLayers = maxLayers;
      bestBucket = bucket;
    }
  }

  BOOST_CHECK_GT(bucketsWithHits, 0u);
  BOOST_CHECK_GT(totalNonEmptyCells, 0u);
  BOOST_CHECK_GT(bestMaxHits, 0.0f);

  const auto [bestXBin, bestYBin] = plane.locMaxHits(bestBucket);

  const double bestTanTheta = Acts::HoughTransformUtils::binCenter(
      axisRanges.xMin, axisRanges.xMax, plane.nBinsX(), bestXBin);
  const double bestInterceptY = Acts::HoughTransformUtils::binCenter(
      axisRanges.yMin, axisRanges.yMax, plane.nBinsY(), bestYBin);

  BOOST_TEST_MESSAGE("Events in file: " << lastEvent);
  BOOST_TEST_MESSAGE("Space points in first event: " << spacePoints.size());
  BOOST_TEST_MESSAGE("Buckets in first event: " << spacePoints.bucketCount());
  BOOST_TEST_MESSAGE("Hough axis ranges: tanTheta=["
                     << axisRanges.xMin << ", " << axisRanges.xMax
                     << "], interceptY=[" << axisRanges.yMin << ", "
                     << axisRanges.yMax << "]");
  BOOST_TEST_MESSAGE("Buckets with non-empty Hough response: "
                     << bucketsWithHits << " / " << plane.nBuckets());
  BOOST_TEST_MESSAGE("Total non-empty cells: " << totalNonEmptyCells);
  BOOST_TEST_MESSAGE("Best bucket: " << bestBucket);
  BOOST_TEST_MESSAGE("Best maximum bin: x=" << bestXBin << ", y=" << bestYBin);
  BOOST_TEST_MESSAGE("Best maximum parameters: tanTheta="
                     << bestTanTheta << ", interceptY=" << bestInterceptY);
  BOOST_TEST_MESSAGE("Best max hits: " << bestMaxHits);
  BOOST_TEST_MESSAGE("Best max layers: " << bestMaxLayers);

  // Dump bucket 0 for visualization.
  const std::filesystem::path outDir = std::filesystem::current_path();
  const std::filesystem::path hitsCsv =
      outDir / "cuda_hough_root_first_bucket_hits.csv";
  const std::filesystem::path histCsv =
      outDir / "cuda_hough_root_first_bucket_histogram.csv";

  writeFirstBucketHitsCsv(hitsCsv, spacePoints);
  writeHoughHistogramCsv(histCsv, plane, axisRanges);

  BOOST_TEST_MESSAGE("Wrote first ROOT bucket hits to: " << hitsCsv.string());
  BOOST_TEST_MESSAGE(
      "Wrote first ROOT bucket Hough histogram to: " << histCsv.string());

  const auto [bucket0XBin, bucket0YBin] = plane.locMaxHits(0);

  BOOST_TEST_MESSAGE("Bucket 0 max bin: x=" << bucket0XBin
                                            << ", y=" << bucket0YBin);
  BOOST_TEST_MESSAGE("Bucket 0 max hits: " << plane.maxHits(0));
  BOOST_TEST_MESSAGE("Bucket 0 max layers: " << plane.maxLayers(0));
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
