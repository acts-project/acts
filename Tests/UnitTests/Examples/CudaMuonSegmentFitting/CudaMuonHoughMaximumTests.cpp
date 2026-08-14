// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/EventData/CudaMuonHoughMaximum.hpp"
#include "ActsExamples/Utilities/CudaUtilities.hpp"

#include <array>
#include <cstdint>
#include <stdexcept>

#include <cuda_runtime.h>

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(CudaMuonHoughMaximumHostConstruction) {
  using MaximumBatch = ActsExamples::CudaHoughMaximumBatch<2u>;

  MaximumBatch maxima{3u};

  BOOST_CHECK_EQUAL(maxima.nBuckets(), 3u);
  BOOST_CHECK_EQUAL(maxima.capacityPerBucket(), 2u);
  BOOST_CHECK_EQUAL(maxima.totalCapacity(), 6u);
  BOOST_CHECK_EQUAL(maxima.totalAssociatedHits(), 0u);

  BOOST_CHECK(!maxima.empty());
  BOOST_CHECK(!maxima.isOnDevice());
  BOOST_CHECK(!maxima.hasAssociationStorage());

  BOOST_CHECK_EQUAL(maxima.nMaxima(0u), 0u);
  BOOST_CHECK_EQUAL(maxima.nMaxima(1u), 0u);
  BOOST_CHECK_EQUAL(maxima.nMaxima(2u), 0u);

  BOOST_CHECK_THROW(maxima.nAssociatedHits(0u, 0u), std::logic_error);
  BOOST_CHECK_THROW(maxima.associatedHitIndices(0u, 0u), std::logic_error);
}

BOOST_AUTO_TEST_CASE(CudaMuonHoughMaximumExactAssociationStorage) {
  int deviceCount = 0;

  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess || deviceCount == 0) {
    BOOST_TEST_MESSAGE("No CUDA device found, skipping CUDA runtime test");
    return;
  }

  using MaximumBatch = ActsExamples::CudaHoughMaximumBatch<2u>;

  MaximumBatch maxima{3u};
  ActsExamples::CudaStream stream;
  maxima.moveToDevice(stream.get());

  BOOST_REQUIRE(maxima.isOnDevice());

  auto device = maxima.deviceArrays();

  BOOST_REQUIRE(device.nMaxima != nullptr);
  BOOST_REQUIRE(device.nAssociatedHits != nullptr);
  BOOST_CHECK(device.associatedHitOffsets == nullptr);
  BOOST_CHECK(device.associatedHitIndices == nullptr);
  BOOST_CHECK_EQUAL(device.totalAssociatedHits, 0u);

  // Bucket 0 has two maxima, bucket 1 has one, bucket 2 has none.
  const std::array<std::uint32_t, 3u> nMaxima{2u, 1u, 0u};

  // Counts are stored bucket by bucket:
  // bucket 0 -> 3, 1
  // bucket 1 -> 2, unused
  // bucket 2 -> unused, unused
  const std::array<std::uint32_t, 6u> nAssociatedHits{3u, 1u, 2u, 0u, 0u, 0u};

  const std::array<ActsExamples::CoordType, 6u> tanBeta{42.0, 43.0, 44.0,
                                                        0.0,  0.0,  0.0};

  BOOST_REQUIRE_EQUAL(cudaMemcpy(device.nMaxima, nMaxima.data(),
                                 sizeof(nMaxima), cudaMemcpyHostToDevice),
                      cudaSuccess);

  BOOST_REQUIRE_EQUAL(
      cudaMemcpy(device.nAssociatedHits, nAssociatedHits.data(),
                 sizeof(nAssociatedHits), cudaMemcpyHostToDevice),
      cudaSuccess);

  BOOST_REQUIRE_EQUAL(cudaMemcpy(device.tanBeta, tanBeta.data(),
                                 sizeof(tanBeta), cudaMemcpyHostToDevice),
                      cudaSuccess);

  // Only nMaxima and nAssociatedHits are copied here.
  maxima.copyAssociationMetadataToHost(stream.get());

  BOOST_CHECK_EQUAL(maxima.nMaxima(0u), 2u);
  BOOST_CHECK_EQUAL(maxima.nMaxima(1u), 1u);
  BOOST_CHECK_EQUAL(maxima.nMaxima(2u), 0u);

  BOOST_CHECK_EQUAL(maxima.nAssociatedHits(0u, 0u), 3u);
  BOOST_CHECK_EQUAL(maxima.nAssociatedHits(0u, 1u), 1u);
  BOOST_CHECK_EQUAL(maxima.nAssociatedHits(1u, 0u), 2u);

  // tanBeta was deliberately not copied by the metadata-only transfer.
  BOOST_CHECK_EQUAL(maxima.tanBeta(0u, 0u), 0.0);

  maxima.allocateAssociationStorage(stream.get());

  BOOST_REQUIRE(maxima.hasAssociationStorage());
  BOOST_CHECK_EQUAL(maxima.totalAssociatedHits(), 6u);

  device = maxima.deviceArrays();

  BOOST_REQUIRE(device.associatedHitOffsets != nullptr);
  BOOST_REQUIRE(device.associatedHitIndices != nullptr);
  BOOST_CHECK_EQUAL(device.totalAssociatedHits, 6u);

  std::array<std::uint32_t, 7u> offsets{};

  BOOST_REQUIRE_EQUAL(cudaMemcpy(offsets.data(), device.associatedHitOffsets,
                                 sizeof(offsets), cudaMemcpyDeviceToHost),
                      cudaSuccess);

  const std::array<std::uint32_t, 7u> expectedOffsets{0u, 3u, 4u, 6u,
                                                      6u, 6u, 6u};

  BOOST_CHECK_EQUAL_COLLECTIONS(offsets.begin(), offsets.end(),
                                expectedOffsets.begin(), expectedOffsets.end());

  const std::array<std::uint32_t, 6u> associatedIndices{10u, 11u, 12u,
                                                        20u, 30u, 31u};

  BOOST_REQUIRE_EQUAL(
      cudaMemcpy(device.associatedHitIndices, associatedIndices.data(),
                 sizeof(associatedIndices), cudaMemcpyHostToDevice),
      cudaSuccess);

  maxima.copyAssociatedHitIndicesToHost(stream.get());

  const std::array<std::uint32_t, 3u> expectedFirst{10u, 11u, 12u};
  const std::array<std::uint32_t, 1u> expectedSecond{20u};
  const std::array<std::uint32_t, 2u> expectedThird{30u, 31u};

  const auto first = maxima.associatedHitIndices(0u, 0u);
  const auto second = maxima.associatedHitIndices(0u, 1u);
  const auto third = maxima.associatedHitIndices(1u, 0u);

  BOOST_CHECK_EQUAL_COLLECTIONS(first.begin(), first.end(),
                                expectedFirst.begin(), expectedFirst.end());

  BOOST_CHECK_EQUAL_COLLECTIONS(second.begin(), second.end(),
                                expectedSecond.begin(), expectedSecond.end());

  BOOST_CHECK_EQUAL_COLLECTIONS(third.begin(), third.end(),
                                expectedThird.begin(), expectedThird.end());

  // The ordinary maximum values are copied only by moveToHost(stream).
  maxima.moveToHost(stream.get());

  BOOST_CHECK_EQUAL(maxima.tanBeta(0u, 0u), 42.0);
  BOOST_CHECK_EQUAL(maxima.tanBeta(0u, 1u), 43.0);
  BOOST_CHECK_EQUAL(maxima.tanBeta(1u, 0u), 44.0);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
