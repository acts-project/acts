// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <cuda_runtime.h>

using namespace Acts;

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(CudaMuonSpacePointHostAccess) {
  static_assert(Acts::Experimental::CompositeSpacePoint<
                ActsExamples::CudaMuonSpacePointProxy>);
  static_assert(Acts::Experimental::CompositeSpacePointPtr<
                ActsExamples::CudaMuonSpacePointPtr>);
  static_assert(Acts::Experimental::CompositeSpacePointContainer<
                ActsExamples::CudaMuonSpacePointContainer>);

  ActsExamples::CudaMuonSpacePointContainer container{2};

  ActsExamples::MuonSpacePoint::MuonId muonId{};
  muonId.setChamber(ActsExamples::MuonSpacePoint::MuonId::StationName::BIS,
                    ActsExamples::MuonSpacePoint::MuonId::DetSide::A, 1,
                    ActsExamples::MuonSpacePoint::MuonId::TechField::Mdt);
  muonId.setLayAndCh(1, 1);
  muonId.setCoordFlags(false, true, true);

  container.setGeometryId(0, GeometryIdentifier{42}.value());
  container.setId(0, muonId.toInt());

  const Vector3 position{1.0, 2.0, 3.0};
  const Vector3 sensorDirection{1.0, 0.0, 0.0};
  const Vector3 toNextSensor{0.0, 1.0, 0.0};

  container.defineCoordinates(0, position, sensorDirection, toNextSensor);
  container.setRadius(0, 4.0);
  container.setTime(0, 5.0);
  container.setCovariance(0, 6.0, 7.0, 8.0);

  container.addBucket(0, 1);
  container.addBucket(1, 2);

  BOOST_CHECK_EQUAL(container.size(), 2u);
  BOOST_CHECK_EQUAL(container.bucketCount(), 2u);
  BOOST_CHECK_EQUAL(container.bucketStart(0), 0u);
  BOOST_CHECK_EQUAL(container.bucketEnd(0), 1u);
  BOOST_CHECK_EQUAL(container.bucketStart(1), 1u);
  BOOST_CHECK_EQUAL(container.bucketEnd(1), 2u);

  auto spacePoint = container[0];

  const Vector3& loadedPosition = spacePoint->localPosition();
  BOOST_CHECK_EQUAL(loadedPosition.x(), 1.0);
  BOOST_CHECK_EQUAL(loadedPosition.y(), 2.0);
  BOOST_CHECK_EQUAL(loadedPosition.z(), 3.0);

  const Vector3& loadedNormal = spacePoint->planeNormal();
  BOOST_CHECK_EQUAL(loadedNormal.x(), 0.0);
  BOOST_CHECK_EQUAL(loadedNormal.y(), 0.0);
  BOOST_CHECK_EQUAL(loadedNormal.z(), 1.0);

  const std::array<double, 3>& covariance = spacePoint->covariance();
  BOOST_CHECK_EQUAL(covariance[0], 6.0);
  BOOST_CHECK_EQUAL(covariance[1], 7.0);
  BOOST_CHECK_EQUAL(covariance[2], 8.0);

  BOOST_CHECK_EQUAL(spacePoint->geometryId().value(),
                    GeometryIdentifier{42}.value());
  BOOST_CHECK(spacePoint->isStraw());
  BOOST_CHECK(spacePoint->hasTime());
  BOOST_CHECK(!spacePoint->measuresLoc1());
  BOOST_CHECK(spacePoint->measuresLoc0());
  BOOST_CHECK_EQUAL(spacePoint->driftRadius(), 4.0);
  BOOST_CHECK_EQUAL(spacePoint->time(), 5.0);
}

BOOST_AUTO_TEST_CASE(CudaMuonSpacePointDeviceTransfer) {
  int deviceCount = 0;

  if (cudaGetDeviceCount(&deviceCount) != cudaSuccess || deviceCount == 0) {
    BOOST_TEST_MESSAGE("No CUDA device found, skipping CUDA runtime test");
    return;
  }

  ActsExamples::CudaMuonSpacePointContainer container{1};

  ActsExamples::MuonSpacePoint::MuonId muonId{};
  muonId.setChamber(ActsExamples::MuonSpacePoint::MuonId::StationName::BIS,
                    ActsExamples::MuonSpacePoint::MuonId::DetSide::A, 1,
                    ActsExamples::MuonSpacePoint::MuonId::TechField::Mdt);
  muonId.setLayAndCh(1, 1);
  muonId.setCoordFlags(false, true, true);

  container.setGeometryId(0, GeometryIdentifier{42}.value());
  container.setId(0, muonId.toInt());
  container.defineCoordinates(0, Vector3{1.0, 2.0, 3.0}, Vector3{1.0, 0.0, 0.0},
                              Vector3{0.0, 1.0, 0.0});
  container.setRadius(0, 4.0);
  container.setTime(0, 5.0);
  container.setCovariance(0, 6.0, 7.0, 8.0);
  container.addBucket(0, 1);

  container.moveToDevice();

  BOOST_CHECK(container.isOnDevice());
  BOOST_CHECK(container.deviceArrays().localPositionX != nullptr);
  BOOST_CHECK(container.deviceArrays().localPositionY != nullptr);
  BOOST_CHECK(container.deviceArrays().localPositionZ != nullptr);
  BOOST_CHECK(container.deviceArrays().bucketStart != nullptr);
  BOOST_CHECK(container.deviceArrays().bucketEnd != nullptr);

  // Invalidate the host memory
  container.defineCoordinates(0, Vector3{0.0, 0.0, 0.0}, Vector3{0.0, 1.0, 0.0},
                              Vector3{0.0, 0.0, 1.0});

  container.moveToHost();

  auto spacePoint = container[0];
  const Vector3& position = spacePoint->localPosition();

  BOOST_CHECK_EQUAL(position.x(), 1.0);
  BOOST_CHECK_EQUAL(position.y(), 2.0);
  BOOST_CHECK_EQUAL(position.z(), 3.0);
}

BOOST_AUTO_TEST_CASE(CudaMuonSpacePointConstructFromMuonSpacePointContainer) {
  ActsExamples::MuonSpacePoint::MuonId muonId{};
  muonId.setChamber(ActsExamples::MuonSpacePoint::MuonId::StationName::BIS,
                    ActsExamples::MuonSpacePoint::MuonId::DetSide::A, 1,
                    ActsExamples::MuonSpacePoint::MuonId::TechField::Mdt);
  muonId.setLayAndCh(2, 17);
  muonId.setCoordFlags(true, false, true);

  ActsExamples::MuonSpacePoint spacePoint{};
  spacePoint.setGeometryId(GeometryIdentifier{42});
  spacePoint.setId(muonId);
  spacePoint.defineCoordinates(Vector3{1.0, 2.0, 3.0}, Vector3{1.0, 0.0, 0.0},
                               Vector3{0.0, 1.0, 0.0});
  spacePoint.setRadius(4.0);
  spacePoint.setTime(5.0);
  spacePoint.setCovariance(6.0, 7.0, 8.0);

  ActsExamples::MuonSpacePointContainer input{};
  input.emplace_back();
  input.back().push_back(std::move(spacePoint));

  ActsExamples::CudaMuonSpacePointContainer container{input};

  BOOST_CHECK_EQUAL(container.size(), 1u);
  BOOST_CHECK_EQUAL(container.bucketCount(), 1u);
  BOOST_CHECK_EQUAL(container.bucketStart(0), 0u);
  BOOST_CHECK_EQUAL(container.bucketEnd(0), 1u);
  BOOST_CHECK(!container.isOnDevice());

  auto converted = container[0];

  BOOST_CHECK_EQUAL(converted->geometryId().value(),
                    GeometryIdentifier{42}.value());
  BOOST_CHECK_EQUAL(converted->id().toInt(), muonId.toInt());

  BOOST_CHECK_EQUAL(converted->localPosition().x(), 1.0);
  BOOST_CHECK_EQUAL(converted->localPosition().y(), 2.0);
  BOOST_CHECK_EQUAL(converted->localPosition().z(), 3.0);

  BOOST_CHECK_EQUAL(converted->sensorDirection().x(), 1.0);
  BOOST_CHECK_EQUAL(converted->sensorDirection().y(), 0.0);
  BOOST_CHECK_EQUAL(converted->sensorDirection().z(), 0.0);

  BOOST_CHECK_EQUAL(converted->toNextSensor().x(), 0.0);
  BOOST_CHECK_EQUAL(converted->toNextSensor().y(), 1.0);
  BOOST_CHECK_EQUAL(converted->toNextSensor().z(), 0.0);

  BOOST_CHECK_EQUAL(converted->planeNormal().x(), 0.0);
  BOOST_CHECK_EQUAL(converted->planeNormal().y(), 0.0);
  BOOST_CHECK_EQUAL(converted->planeNormal().z(), 1.0);

  BOOST_CHECK_EQUAL(converted->driftRadius(), 4.0);
  BOOST_CHECK_EQUAL(converted->time(), 5.0);

  BOOST_CHECK_EQUAL(converted->covariance()[0], 6.0);
  BOOST_CHECK_EQUAL(converted->covariance()[1], 7.0);
  BOOST_CHECK_EQUAL(converted->covariance()[2], 8.0);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
