// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/Io/Root/RootCudaMuonSpacePointReader.hpp"

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/UnitVectors.hpp"

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <cmath>
#include <filesystem>
#include <memory>
#include <vector>

using namespace Acts;

BOOST_AUTO_TEST_CASE(CudaMuonSpacePointLoadRootFileToHost) {
  using namespace Acts::UnitLiterals;

  const std::filesystem::path fullPath{
      "/data/mgawlas/ACTS/acts/Tests/Data/ParticleGun_MU0.root"};

  std::filesystem::path filePath{};

  for (std::filesystem::path path = std::filesystem::current_path();
       !path.empty(); path = path.parent_path()) {
    const auto candidate = fullPath;
    std::cout << "Candidate path: " << candidate << std::endl;

    BOOST_TEST_MESSAGE("Candidate path: " << candidate);

    if (std::filesystem::exists(candidate)) {
      filePath = candidate;
      break;
    }

    if (path == path.root_path()) {
      break;
    }
  }

  BOOST_REQUIRE_MESSAGE(!filePath.empty(),
                        "Could not find ROOT input file: "
                            << fullPath.string());

  std::unique_ptr<TFile> file{TFile::Open(filePath.c_str(), "READ")};

  BOOST_REQUIRE(file != nullptr);
  BOOST_REQUIRE(!file->IsZombie());

  TTreeReader reader{"MuonSpacePoints", file.get()};

  BOOST_REQUIRE(reader.GetTree() != nullptr);

  TTreeReaderValue<std::uint32_t> eventId{reader, "event_id"};

  TTreeReaderValue<std::vector<Acts::GeometryIdentifier::Value>> geometryId{
      reader, "spacePoint_geometryId"};
  TTreeReaderValue<std::vector<std::uint16_t>> bucketId{
      reader, "spacePoint_bucketId"};
  TTreeReaderValue<std::vector<std::uint32_t>> muonId{
      reader, "spacePoint_muonId"};

  TTreeReaderValue<std::vector<float>> localPositionX{
      reader, "spacePoint_localPosX"};
  TTreeReaderValue<std::vector<float>> localPositionY{
      reader, "spacePoint_localPosY"};
  TTreeReaderValue<std::vector<float>> localPositionZ{
      reader, "spacePoint_localPosZ"};

  TTreeReaderValue<std::vector<float>> sensorDirectionTheta{
      reader, "spacePoint_sensorDirTheta"};
  TTreeReaderValue<std::vector<float>> sensorDirectionPhi{
      reader, "spacePoint_sensorDirPhi"};

  TTreeReaderValue<std::vector<float>> toNextSensorTheta{
      reader, "spacePoint_toNextDirTheta"};
  TTreeReaderValue<std::vector<float>> toNextSensorPhi{
      reader, "spacePoint_toNextDirPhi"};

  TTreeReaderValue<std::vector<float>> covLoc0{reader, "spacePoint_covLoc0"};
  TTreeReaderValue<std::vector<float>> covLoc1{reader, "spacePoint_covLoc1"};
  TTreeReaderValue<std::vector<float>> covT{reader, "spacePoint_covT"};

  TTreeReaderValue<std::vector<float>> driftR{
      reader, "spacePoint_driftRadius"};
  TTreeReaderValue<std::vector<float>> time{reader, "spacePoint_time"};

  const auto nEvents = static_cast<std::size_t>(reader.GetEntries());
  BOOST_REQUIRE_GT(nEvents, 0u);

  BOOST_TEST_MESSAGE("ROOT file: " << filePath.string());
  BOOST_TEST_MESSAGE("Number of events: " << nEvents);

  BOOST_REQUIRE(reader.Next());

  const std::size_t nSpacePoints = bucketId->size();

  BOOST_TEST_MESSAGE("First event id: " << *eventId);
  BOOST_TEST_MESSAGE("Space points in first event: " << nSpacePoints);

  BOOST_REQUIRE_GT(nSpacePoints, 0u);

  ActsExamples::CudaMuonSpacePointContainer container{nSpacePoints};

  std::size_t bucketStart = 0;
  std::uint16_t currentBucket = 0;
  bool hasOpenBucket = false;

  for (std::size_t spIdx = 0; spIdx < nSpacePoints; ++spIdx) {
    const auto currentBucketId = bucketId->at(spIdx);

    if (!hasOpenBucket) {
      bucketStart = spIdx;
      currentBucket = currentBucketId;
      hasOpenBucket = true;
    } else if (currentBucketId != currentBucket) {
      container.addBucket(bucketStart, spIdx);
      bucketStart = spIdx;
      currentBucket = currentBucketId;
    }

    container.setGeometryId(spIdx, geometryId->at(spIdx));
    container.setId(spIdx, muonId->at(spIdx));

    const Vector3 position{localPositionX->at(spIdx),
                           localPositionY->at(spIdx),
                           localPositionZ->at(spIdx)};

    const Vector3 sensorDir{makeDirectionFromPhiTheta<double>(
        sensorDirectionPhi->at(spIdx) * 1._degree,
        sensorDirectionTheta->at(spIdx) * 1._degree)};

    const Vector3 toNext{makeDirectionFromPhiTheta<double>(
        toNextSensorPhi->at(spIdx) * 1._degree,
        toNextSensorTheta->at(spIdx) * 1._degree)};

    container.defineCoordinates(spIdx, position, sensorDir, toNext);

    container.setRadius(spIdx, driftR->at(spIdx));
    container.setTime(spIdx, time->at(spIdx));
    container.setCovariance(spIdx, covLoc0->at(spIdx), covLoc1->at(spIdx),
                            covT->at(spIdx));
  }

  if (hasOpenBucket) {
    container.addBucket(bucketStart, nSpacePoints);
  }

  BOOST_CHECK_EQUAL(container.size(), nSpacePoints);
  BOOST_CHECK_GT(container.bucketCount(), 0u);
  BOOST_CHECK_EQUAL(container.bucketStart(0), 0u);
  BOOST_CHECK_LE(container.bucketEnd(0), container.size());

  auto firstSpacePoint = container[0];

  const Vector3& firstPosition = firstSpacePoint->localPosition();
  const std::array<double, 3>& firstCovariance = firstSpacePoint->covariance();

  BOOST_TEST_MESSAGE("Buckets in first event: " << container.bucketCount());
  BOOST_TEST_MESSAGE("First bucket range: [" << container.bucketStart(0) << ", "
                                             << container.bucketEnd(0) << ")");
  BOOST_TEST_MESSAGE("First space point position: (" << firstPosition.x()
                                                     << ", "
                                                     << firstPosition.y()
                                                     << ", "
                                                     << firstPosition.z()
                                                     << ")");
  BOOST_TEST_MESSAGE("First space point drift radius: "
                     << firstSpacePoint->driftRadius());
  BOOST_TEST_MESSAGE("First space point time: " << firstSpacePoint->time());
  BOOST_TEST_MESSAGE("First space point covariance: (" << firstCovariance[0]
                                                       << ", "
                                                       << firstCovariance[1]
                                                       << ", "
                                                       << firstCovariance[2]
                                                       << ")");

  BOOST_CHECK(std::isfinite(firstPosition.x()));
  BOOST_CHECK(std::isfinite(firstPosition.y()));
  BOOST_CHECK(std::isfinite(firstPosition.z()));
}
