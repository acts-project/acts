// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "ActsExamples/EventData/CudaMuonSpacePoint.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Io/Root/RootCudaMuonSpacePointReader.hpp"
#include "ActsExamples/Io/Root/RootCudaMuonSpacePointWriter.hpp"

#include <cmath>
#include <filesystem>
#include <memory>
#include <vector>

#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

using namespace Acts;

namespace {

ActsExamples::CudaMuonSpacePointContainer makeTestContainer() {
  ActsExamples::CudaMuonSpacePointContainer container{2};

  ActsExamples::MuonSpacePoint::MuonId firstId{};
  firstId.setChamber(ActsExamples::MuonSpacePoint::MuonId::StationName::BIS,
                     ActsExamples::MuonSpacePoint::MuonId::DetSide::A, 1,
                     ActsExamples::MuonSpacePoint::MuonId::TechField::Mdt);
  firstId.setLayAndCh(1, 1);
  firstId.setCoordFlags(false, true, true);

  ActsExamples::MuonSpacePoint::MuonId secondId{};
  secondId.setChamber(ActsExamples::MuonSpacePoint::MuonId::StationName::BIL,
                      ActsExamples::MuonSpacePoint::MuonId::DetSide::C, 2,
                      ActsExamples::MuonSpacePoint::MuonId::TechField::Rpc);
  secondId.setLayAndCh(2, 3);
  secondId.setCoordFlags(true, false, false);

  container.setGeometryId(0, Acts::GeometryIdentifier{42}.value());
  container.setId(0, firstId.toInt());
  container.defineCoordinates(0, Acts::Vector3{1.0, 2.0, 3.0},
                              Acts::Vector3{1.0, 0.0, 0.0},
                              Acts::Vector3{0.0, 1.0, 0.0});
  container.setRadius(0, 4.0);
  container.setTime(0, 5.0);
  container.setCovariance(0, 6.0, 7.0, 8.0);

  container.setGeometryId(1, Acts::GeometryIdentifier{43}.value());
  container.setId(1, secondId.toInt());
  container.defineCoordinates(1, Acts::Vector3{10.0, 20.0, 30.0},
                              Acts::Vector3{0.0, 1.0, 0.0},
                              Acts::Vector3{0.0, 0.0, 1.0});
  container.setRadius(1, 40.0);
  container.setTime(1, 50.0);
  container.setCovariance(1, 60.0, 70.0, 80.0);

  container.addBucket(0, 1);
  container.addBucket(1, 2);

  return container;
}

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_CASE(RootCudaMuonSpacePointWriterWritesAdHocContainer) {
  const auto filePath = std::filesystem::temp_directory_path() /
                        "ActsRootCudaMuonSpacePointWriterAdHoc.root";

  std::filesystem::remove(filePath);

  ActsExamples::RootCudaMuonSpacePointWriter::Config cfg{};
  cfg.inputSpacePoints = "CudaMuonSpacePoints";
  cfg.filePath = filePath.string();
  cfg.fileMode = "RECREATE";
  cfg.treeName = "muonSpacePoints";

  ActsExamples::RootCudaMuonSpacePointWriter writer{cfg, Acts::Logging::INFO};

  ActsExamples::WhiteBoard eventStore;
  ActsExamples::AlgorithmContext context{0, 0, eventStore, 0};

  ActsExamples::WriteDataHandle<ActsExamples::CudaMuonSpacePointContainer>
      inputHandle{&writer, "InputSpacePoints"};
  inputHandle.initialize(cfg.inputSpacePoints);
  inputHandle(context, makeTestContainer());

  BOOST_CHECK(writer.write(context) == ActsExamples::ProcessCode::SUCCESS);
  BOOST_CHECK(writer.finalize() == ActsExamples::ProcessCode::SUCCESS);

  BOOST_REQUIRE(std::filesystem::exists(filePath));

  std::unique_ptr<TFile> file{TFile::Open(filePath.c_str(), "READ")};

  BOOST_REQUIRE(file != nullptr);
  BOOST_REQUIRE(!file->IsZombie());

  TTreeReader reader{cfg.treeName.c_str(), file.get()};

  BOOST_REQUIRE(reader.GetTree() != nullptr);
  BOOST_REQUIRE_EQUAL(reader.GetEntries(), 1);

  TTreeReaderValue<std::uint32_t> eventId{reader, "event_id"};
  TTreeReaderValue<std::vector<std::uint16_t>> bucketId{reader,
                                                        "spacePoint_bucketId"};
  TTreeReaderValue<std::vector<Acts::GeometryIdentifier::Value>> geometryId{
      reader, "spacePoint_geometryId"};
  TTreeReaderValue<std::vector<std::uint32_t>> muonId{reader,
                                                      "spacePoint_muonId"};
  TTreeReaderValue<std::vector<float>> localPositionX{reader,
                                                      "spacePoint_localPosX"};
  TTreeReaderValue<std::vector<float>> localPositionY{reader,
                                                      "spacePoint_localPosY"};
  TTreeReaderValue<std::vector<float>> localPositionZ{reader,
                                                      "spacePoint_localPosZ"};
  TTreeReaderValue<std::vector<float>> driftRadius{reader,
                                                   "spacePoint_driftRadius"};
  TTreeReaderValue<std::vector<float>> time{reader, "spacePoint_time"};

  BOOST_REQUIRE(reader.Next());

  BOOST_CHECK_EQUAL(*eventId, 0u);

  BOOST_REQUIRE_EQUAL(bucketId->size(), 2u);
  BOOST_REQUIRE_EQUAL(geometryId->size(), 2u);
  BOOST_REQUIRE_EQUAL(muonId->size(), 2u);
  BOOST_REQUIRE_EQUAL(localPositionX->size(), 2u);
  BOOST_REQUIRE_EQUAL(localPositionY->size(), 2u);
  BOOST_REQUIRE_EQUAL(localPositionZ->size(), 2u);

  BOOST_CHECK_EQUAL(bucketId->at(0), 0u);
  BOOST_CHECK_EQUAL(bucketId->at(1), 1u);

  BOOST_CHECK_EQUAL(geometryId->at(0), Acts::GeometryIdentifier{42}.value());
  BOOST_CHECK_EQUAL(geometryId->at(1), Acts::GeometryIdentifier{43}.value());

  BOOST_CHECK_EQUAL(localPositionX->at(0), 1.0f);
  BOOST_CHECK_EQUAL(localPositionY->at(0), 2.0f);
  BOOST_CHECK_EQUAL(localPositionZ->at(0), 3.0f);

  BOOST_CHECK_EQUAL(localPositionX->at(1), 10.0f);
  BOOST_CHECK_EQUAL(localPositionY->at(1), 20.0f);
  BOOST_CHECK_EQUAL(localPositionZ->at(1), 30.0f);

  BOOST_CHECK_EQUAL(driftRadius->at(0), 4.0f);
  BOOST_CHECK_EQUAL(time->at(0), 5.0f);

  std::filesystem::remove(filePath);
}

/*
/// Test is problematic as it relies on existence of ParticleGun_MU0.root file
BOOST_AUTO_TEST_CASE(RootCudaMuonSpacePointReaderReadsFirstEvent) {
  const std::filesystem::path filePath{
      "/data/mgawlas/ACTS/acts/Tests/Data/ParticleGun_MU0.root"};

  BOOST_REQUIRE(std::filesystem::exists(filePath));

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

  const auto& container = outputHandle(context);

  BOOST_REQUIRE_GT(container.size(), 0u);
  BOOST_REQUIRE_GT(container.bucketCount(), 0u);

  BOOST_CHECK_EQUAL(container.bucketStart(0), 0u);
  BOOST_CHECK_LE(container.bucketEnd(0), container.size());

  auto firstSpacePoint = container[0];

  const Acts::Vector3& position = firstSpacePoint->localPosition();
  const std::array<double, 3>& covariance = firstSpacePoint->covariance();

  BOOST_TEST_MESSAGE("Events in file: " << lastEvent);
  BOOST_TEST_MESSAGE("Space points in first event: " << container.size());
  BOOST_TEST_MESSAGE("Buckets in first event: " << container.bucketCount());
  BOOST_TEST_MESSAGE("First bucket range: [" << container.bucketStart(0)
                                             << ", " << container.bucketEnd(0)
                                             << ")");
  BOOST_TEST_MESSAGE("First position: (" << position.x() << ", "
                                         << position.y() << ", "
                                         << position.z() << ")");
  BOOST_TEST_MESSAGE("First covariance: (" << covariance[0] << ", "
                                           << covariance[1] << ", "
                                           << covariance[2] << ")");

  BOOST_CHECK(std::isfinite(position.x()));
  BOOST_CHECK(std::isfinite(position.y()));
  BOOST_CHECK(std::isfinite(position.z()));
  BOOST_CHECK(std::isfinite(firstSpacePoint->driftRadius()));
  BOOST_CHECK(std::isfinite(firstSpacePoint->time()));
}
*/

BOOST_AUTO_TEST_CASE(RootCudaMuonSpacePointWriterReaderRoundTrip) {
  const auto filePath = std::filesystem::temp_directory_path() /
                        "ActsRootCudaMuonSpacePointWriterReaderRoundTrip.root";

  std::filesystem::remove(filePath);

  {
    ActsExamples::RootCudaMuonSpacePointWriter::Config writerCfg{};
    writerCfg.inputSpacePoints = "CudaMuonSpacePoints";
    writerCfg.filePath = filePath.string();
    writerCfg.fileMode = "RECREATE";
    writerCfg.treeName = "MuonSpacePoints";

    ActsExamples::RootCudaMuonSpacePointWriter writer{writerCfg,
                                                      Acts::Logging::INFO};

    ActsExamples::WhiteBoard eventStore;
    ActsExamples::AlgorithmContext context{0, 0, eventStore, 0};

    ActsExamples::WriteDataHandle<ActsExamples::CudaMuonSpacePointContainer>
        inputHandle{&writer, "InputSpacePoints"};
    inputHandle.initialize(writerCfg.inputSpacePoints);
    inputHandle(context, makeTestContainer());

    BOOST_CHECK(writer.write(context) == ActsExamples::ProcessCode::SUCCESS);
    BOOST_CHECK(writer.finalize() == ActsExamples::ProcessCode::SUCCESS);
  }

  BOOST_REQUIRE(std::filesystem::exists(filePath));

  ActsExamples::RootCudaMuonSpacePointReader::Config readerCfg{};
  readerCfg.filePath = filePath.string();
  readerCfg.treeName = "MuonSpacePoints";
  readerCfg.outputSpacePoints = "ReadCudaMuonSpacePoints";

  ActsExamples::RootCudaMuonSpacePointReader reader{readerCfg,
                                                    Acts::Logging::INFO};

  const auto [firstEvent, lastEvent] = reader.availableEvents();

  BOOST_CHECK_EQUAL(firstEvent, 0u);
  BOOST_REQUIRE_EQUAL(lastEvent, 1u);

  ActsExamples::WhiteBoard eventStore;
  ActsExamples::AlgorithmContext context{0, 0, eventStore, 0};

  BOOST_CHECK(reader.read(context) == ActsExamples::ProcessCode::SUCCESS);
  BOOST_REQUIRE(eventStore.exists(readerCfg.outputSpacePoints));

  ActsExamples::ReadDataHandle<ActsExamples::CudaMuonSpacePointContainer>
      outputHandle{&reader, "OutputSpacePoints"};
  outputHandle.initialize(readerCfg.outputSpacePoints);

  const auto& container = outputHandle(context);

  BOOST_REQUIRE_EQUAL(container.size(), 2u);
  BOOST_REQUIRE_EQUAL(container.bucketCount(), 2u);

  BOOST_CHECK_EQUAL(container.bucketStart(0), 0u);
  BOOST_CHECK_EQUAL(container.bucketEnd(0), 1u);
  BOOST_CHECK_EQUAL(container.bucketStart(1), 1u);
  BOOST_CHECK_EQUAL(container.bucketEnd(1), 2u);

  auto firstSpacePoint = container[0];
  const Acts::Vector3& firstPosition = firstSpacePoint->localPosition();

  BOOST_CHECK_EQUAL(firstSpacePoint->geometryId().value(),
                    Acts::GeometryIdentifier{42}.value());

  BOOST_CHECK_EQUAL(firstPosition.x(), 1.0);
  BOOST_CHECK_EQUAL(firstPosition.y(), 2.0);
  BOOST_CHECK_EQUAL(firstPosition.z(), 3.0);

  BOOST_CHECK_EQUAL(firstSpacePoint->driftRadius(), 4.0);
  BOOST_CHECK_EQUAL(firstSpacePoint->time(), 5.0);

  const std::array<double, 3>& firstCovariance = firstSpacePoint->covariance();

  BOOST_CHECK_EQUAL(firstCovariance[0], 6.0);
  BOOST_CHECK_EQUAL(firstCovariance[1], 7.0);
  BOOST_CHECK_EQUAL(firstCovariance[2], 8.0);

  BOOST_CHECK(firstSpacePoint->isStraw());
  BOOST_CHECK(firstSpacePoint->hasTime());
  BOOST_CHECK(firstSpacePoint->measuresLoc0());
  BOOST_CHECK(!firstSpacePoint->measuresLoc1());

  auto secondSpacePoint = container[1];
  const Acts::Vector3& secondPosition = secondSpacePoint->localPosition();

  BOOST_CHECK_EQUAL(secondSpacePoint->geometryId().value(),
                    Acts::GeometryIdentifier{43}.value());

  BOOST_CHECK_EQUAL(secondPosition.x(), 10.0);
  BOOST_CHECK_EQUAL(secondPosition.y(), 20.0);
  BOOST_CHECK_EQUAL(secondPosition.z(), 30.0);

  std::filesystem::remove(filePath);
}

}  // namespace ActsTests
