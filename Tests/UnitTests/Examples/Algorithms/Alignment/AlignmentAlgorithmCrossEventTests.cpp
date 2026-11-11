// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Alignment.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsAlignment/Kernel/AlignmentError.hpp"
#include "ActsExamples/Alignment/AlignmentAlgorithm.hpp"
#include "ActsExamples/DetectorCommons/AlignmentContext.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <array>
#include <filesystem>
#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

using namespace Acts;
using namespace ActsExamples;

namespace ActsTests {
namespace {

class RecordingAlignment final : public AlignmentAlgorithm::AlignmentFunction {
 public:
  mutable std::vector<std::vector<IndexSourceLink>> sourceLinks;
  mutable TrackParametersContainer parameters;
  mutable std::size_t calls = 0;
  bool fail = false;
  bool markNotConverged = false;
  bool readFirstMeasurement = false;
  mutable double firstParameter = 0.;
  mutable double firstVariance = 0.;
  std::unordered_map<const Surface*, std::size_t> indexedSurfaces;

  AlignmentAlgorithm::AlignmentResult operator()(
      const std::vector<std::vector<IndexSourceLink>>& links,
      const TrackParametersContainer& initialParameters,
      const ActsAlignment::AlignmentOptions<
          AlignmentAlgorithm::TrackFitterOptions>& options) const override {
    sourceLinks = links;
    parameters = initialParameters;
    ++calls;
    if (fail) {
      return AlignmentAlgorithm::AlignmentResult::failure(
          ActsAlignment::AlignmentError::AlignmentParametersUpdateFailure);
    }

    if (readFirstMeasurement && !links.empty() && !links.front().empty()) {
      Acts::VectorMultiTrajectory trajectory;
      auto state = trajectory.makeTrackState(Acts::TrackStatePropMask::All);
      options.fitOptions.extensions.calibrator(
          GeometryContext::dangerouslyDefaultConstruct(), CalibrationContext{},
          SourceLink{links.front().front()}, state);
      firstParameter = state.template calibrated<1>()[0];
      firstVariance = state.template calibratedCovariance<1>()(0, 0);
    }

    ActsAlignment::AlignmentResult output;
    output.deltaAlignmentParameters =
        DynamicVector::Zero(static_cast<Eigen::Index>(eAlignmentSize));
    output.deltaAlignmentParameters(0) = 0.25;
    output.alignmentCovariance =
        DynamicMatrix::Identity(static_cast<Eigen::Index>(eAlignmentSize),
                                static_cast<Eigen::Index>(eAlignmentSize));
    output.idxedAlignSurfaces = indexedSurfaces;
    output.deltaChi2 = 1.5;
    if (markNotConverged) {
      output.result = Acts::Result<void>::failure(
          ActsAlignment::AlignmentError::ConvergeFailure);
    }
    return Acts::Result<ActsAlignment::AlignmentResult>::success(
        std::move(output));
  }
};

std::shared_ptr<const TrackingGeometry> makeGeometry() {
  auto bounds = std::make_shared<CylinderVolumeBounds>(0., 20., 50.);
  auto volume =
      std::make_shared<TrackingVolume>(Transform3::Identity(), bounds, "World");
  return std::make_shared<TrackingGeometry>(std::move(volume));
}

TrackParameters trackAt(double loc0) {
  auto perigee = Surface::makeShared<PerigeeSurface>(Vector3::Zero());
  BoundVector values = BoundVector::Zero();
  values[eBoundLoc0] = loc0;
  values[eBoundTheta] = 1.;
  values[eBoundQOverP] = 1.;
  return TrackParameters(std::move(perigee), values, std::nullopt,
                         ParticleHypothesis::pion());
}

MeasurementContainer measurementsOf(GeometryIdentifier geoId,
                                    std::size_t count) {
  MeasurementContainer container;
  for (std::size_t i = 0; i < count; ++i) {
    container.makeMeasurement(2, geoId);
  }
  return container;
}

MeasurementContainer measurementOf(GeometryIdentifier geoId, double value,
                                   double variance) {
  MeasurementContainer container;
  auto measurement = container.makeMeasurement(1, geoId);
  measurement.setSubspaceIndices(
      std::array<std::uint8_t, 1>{static_cast<std::uint8_t>(eBoundLoc0)});
  measurement.parameters()(0) = value;
  measurement.covariance()(0, 0) = variance;
  return container;
}

ProtoTrackContainer protoTracks(std::initializer_list<std::size_t> indices) {
  ProtoTrackContainer tracks;
  ProtoTrack track;
  for (std::size_t index : indices) {
    track.push_back(index);
  }
  tracks.push_back(std::move(track));
  return tracks;
}

ProtoTrackContainer protoTrackList(
    std::initializer_list<std::initializer_list<std::size_t>> tracks) {
  ProtoTrackContainer container;
  for (const auto& indices : tracks) {
    ProtoTrack track;
    for (std::size_t index : indices) {
      track.push_back(index);
    }
    container.push_back(std::move(track));
  }
  return container;
}

// The handles register their address on the algorithm and are not
// unregistered on destruction, so they have to outlive every use of the
// algorithm, including writeHandles().
struct EventHandles {
  WriteDataHandle<MeasurementContainer> measurements;
  WriteDataHandle<ProtoTrackContainer> protoTracks;
  WriteDataHandle<TrackParametersContainer> parameters;

  explicit EventHandles(AlignmentAlgorithm& algorithm)
      : measurements{&algorithm, "TestMeasurements"},
        protoTracks{&algorithm, "TestProto"},
        parameters{&algorithm, "TestParameters"} {
    measurements.initialize("measurements");
    protoTracks.initialize("proto");
    parameters.initialize("parameters");
  }

  EventHandles(const EventHandles&) = delete;
  EventHandles(EventHandles&&) = delete;
};

void storeEvent(EventHandles& handles, const AlgorithmContext& ctx,
                MeasurementContainer measurements, ProtoTrackContainer tracks,
                TrackParametersContainer parameters) {
  handles.measurements(ctx, std::move(measurements));
  handles.protoTracks(ctx, std::move(tracks));
  handles.parameters(ctx, std::move(parameters));
}

AlignmentAlgorithm::Config configFor(
    const std::shared_ptr<RecordingAlignment>& align,
    const std::shared_ptr<const TrackingGeometry>& geometry,
    const std::filesystem::path& directory) {
  AlignmentAlgorithm::Config cfg;
  cfg.inputMeasurements = "measurements";
  cfg.inputProtoTracks = "proto";
  cfg.inputInitialTrackParameters = "parameters";
  cfg.outputAlignmentFile = directory / "aligned_result.txt";
  cfg.outputAlignmentIndexFile = directory / "index_map.txt";
  cfg.trackingGeometry = geometry;
  cfg.align = align;
  cfg.maxNumTracks = 2;
  return cfg;
}

}  // namespace

BOOST_AUTO_TEST_SUITE(AlignmentAlgorithmCrossEventSuite)

BOOST_AUTO_TEST_CASE(CollectsUntilTheTrackLimitAndWritesMp2) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(10., 20.));
  auto geoId = GeometryIdentifier()
                   .withVolume(1)
                   .withBoundary(4)
                   .withLayer(2)
                   .withApproach(5)
                   .withSensitive(3);
  surface->assignGeometryId(geoId);
  align->indexedSurfaces.emplace(surface.get(), 0);

  const auto directory =
      std::filesystem::temp_directory_path() / "acts-alignment-cross-event";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard firstStore;
  AlgorithmContext first{0, 0, firstStore, 0};
  storeEvent(handles, first, measurementsOf(geoId, 1), protoTracks({0}),
             {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(first) == ProcessCode::SUCCESS);

  WhiteBoard secondStore;
  AlgorithmContext second{0, 1, secondStore, 0};
  storeEvent(handles, second, measurementsOf(geoId, 2), protoTracks({1}),
             {trackAt(20.)});
  BOOST_CHECK(algorithm.execute(second) == ProcessCode::SUCCESS);

  WhiteBoard thirdStore;
  AlgorithmContext third{0, 2, thirdStore, 0};
  storeEvent(handles, third, measurementsOf(geoId, 1), protoTracks({0}),
             {trackAt(30.)});
  BOOST_CHECK(algorithm.execute(third) == ProcessCode::SUCCESS);

  BOOST_CHECK(algorithm.finalize() == ProcessCode::SUCCESS);
  BOOST_CHECK_EQUAL(align->calls, 1u);
  BOOST_REQUIRE_EQUAL(align->sourceLinks.size(), 2u);
  BOOST_REQUIRE_EQUAL(align->sourceLinks[0].size(), 1u);
  BOOST_REQUIRE_EQUAL(align->sourceLinks[1].size(), 1u);
  // Event 1 keeps only the measurement its track uses, so that hit is the
  // second measurement in the merged container.
  BOOST_CHECK_EQUAL(align->sourceLinks[0][0].index(), 0u);
  BOOST_CHECK_EQUAL(align->sourceLinks[1][0].index(), 1u);
  BOOST_CHECK_EQUAL(align->parameters[0].parameters()[eBoundLoc0], 10.);
  BOOST_CHECK_EQUAL(align->parameters[1].parameters()[eBoundLoc0], 20.);

  std::ifstream result(directory / "aligned_result.txt");
  BOOST_REQUIRE(result);
  std::string line;
  BOOST_REQUIRE(std::getline(result, line));
  BOOST_CHECK(line.find("99") != std::string::npos);

  std::ifstream index(directory / "index_map.txt");
  BOOST_REQUIRE(index);
  std::string header;
  std::string firstLabel;
  BOOST_REQUIRE(std::getline(index, header));
  BOOST_REQUIRE(std::getline(index, firstLabel));
  std::istringstream labelStream(firstLabel);
  std::string label;
  std::string surfaceIndexText;
  std::string dof;
  std::string name;
  std::string geoIdText;
  labelStream >> label >> surfaceIndexText >> dof >> name >> geoIdText;
  BOOST_CHECK_EQUAL(name, "dx");
  BOOST_CHECK_EQUAL(geoIdText, "vol=1|bnd=4|lay=2|apr=5|sen=3");

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(TruncatesTracksInsideOneEvent) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-truncate";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard store;
  AlgorithmContext ctx{0, 0, store, 0};
  // maxNumTracks is 2. The third track, and the measurement only it uses,
  // are dropped. Kept hits are copied in track order, so proto indices 2
  // and 0 become collected indices 0 and 1.
  storeEvent(handles, ctx, measurementsOf(GeometryIdentifier(), 3),
             protoTrackList({{2}, {0}, {1}}),
             {trackAt(10.), trackAt(20.), trackAt(30.)});
  BOOST_CHECK(algorithm.execute(ctx) == ProcessCode::SUCCESS);
  BOOST_CHECK(algorithm.finalize() == ProcessCode::SUCCESS);

  BOOST_CHECK_EQUAL(align->calls, 1u);
  BOOST_REQUIRE_EQUAL(align->sourceLinks.size(), 2u);
  BOOST_REQUIRE_EQUAL(align->sourceLinks[0].size(), 1u);
  BOOST_REQUIRE_EQUAL(align->sourceLinks[1].size(), 1u);
  BOOST_CHECK_EQUAL(align->sourceLinks[0][0].index(), 0u);
  BOOST_CHECK_EQUAL(align->sourceLinks[1][0].index(), 1u);
  BOOST_CHECK_EQUAL(align->parameters[0].parameters()[eBoundLoc0], 10.);
  BOOST_CHECK_EQUAL(align->parameters[1].parameters()[eBoundLoc0], 20.);

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(AbortsWhenTheResultFileCannotBeOpened) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-unwritable";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  auto cfg = configFor(align, geometry, directory);
  const auto resultPath = directory / "missing" / "aligned_result.txt";
  cfg.outputAlignmentFile = resultPath;

  AlignmentAlgorithm algorithm{
      std::move(cfg),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard store;
  AlgorithmContext ctx{0, 0, store, 0};
  storeEvent(handles, ctx, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(ctx) == ProcessCode::SUCCESS);
  BOOST_CHECK(algorithm.finalize() == ProcessCode::ABORT);
  BOOST_CHECK_EQUAL(align->calls, 1u);
  BOOST_CHECK(!std::filesystem::exists(resultPath));

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(KeepsLowerEventsWhenALaterEventArrivesFirst) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-order";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  auto cfg = configFor(align, geometry, directory);
  cfg.maxNumTracks = 3;
  AlignmentAlgorithm algorithm{
      std::move(cfg),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard laterStore;
  AlgorithmContext later{0, 2, laterStore, 0};
  storeEvent(handles, later, measurementsOf(GeometryIdentifier(), 2),
             protoTrackList({{0}, {1}}), {trackAt(30.), trackAt(31.)});
  BOOST_CHECK(algorithm.execute(later) == ProcessCode::SUCCESS);

  WhiteBoard earlierStore;
  AlgorithmContext earlier{0, 0, earlierStore, 0};
  storeEvent(handles, earlier, measurementsOf(GeometryIdentifier(), 2),
             protoTrackList({{0}, {1}}), {trackAt(10.), trackAt(11.)});
  BOOST_CHECK(algorithm.execute(earlier) == ProcessCode::SUCCESS);

  BOOST_CHECK(algorithm.finalize() == ProcessCode::SUCCESS);
  BOOST_REQUIRE_EQUAL(align->parameters.size(), 3u);
  BOOST_CHECK_EQUAL(align->parameters[0].parameters()[eBoundLoc0], 10.);
  BOOST_CHECK_EQUAL(align->parameters[1].parameters()[eBoundLoc0], 11.);
  BOOST_CHECK_EQUAL(align->parameters[2].parameters()[eBoundLoc0], 30.);

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(PreservesMeasurementValuesInEventOrder) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  align->readFirstMeasurement = true;
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-values";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};
  const auto geoId =
      GeometryIdentifier().withVolume(1).withLayer(2).withSensitive(3);

  WhiteBoard laterStore;
  AlgorithmContext later{0, 1, laterStore, 0};
  storeEvent(handles, later, measurementOf(geoId, 8.5, 9.), protoTracks({0}),
             {trackAt(20.)});
  BOOST_CHECK(algorithm.execute(later) == ProcessCode::SUCCESS);

  WhiteBoard earlierStore;
  AlgorithmContext earlier{0, 0, earlierStore, 0};
  storeEvent(handles, earlier, measurementOf(geoId, 1.25, 4.), protoTracks({0}),
             {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(earlier) == ProcessCode::SUCCESS);

  BOOST_CHECK(algorithm.finalize() == ProcessCode::SUCCESS);
  BOOST_CHECK_EQUAL(align->parameters[0].parameters()[eBoundLoc0], 10.);
  BOOST_CHECK_EQUAL(align->parameters[1].parameters()[eBoundLoc0], 20.);
  BOOST_CHECK_EQUAL(align->firstParameter, 1.25);
  BOOST_CHECK_EQUAL(align->firstVariance, 4.);

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(AbortsOnAMeasurementIndexPastTheContainer) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-index";
  std::filesystem::remove_all(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard store;
  AlgorithmContext ctx{0, 0, store, 0};
  storeEvent(handles, ctx, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({5}), {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(ctx) == ProcessCode::ABORT);
  BOOST_CHECK_EQUAL(align->calls, 0u);
}

BOOST_AUTO_TEST_CASE(AbortsWhenTheAlignmentStoreChanges) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-store";
  std::filesystem::remove_all(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  GeoIdAlignmentStore storeA{
      std::unordered_map<GeometryIdentifier, Transform3>{}};
  GeoIdAlignmentStore storeB{
      std::unordered_map<GeometryIdentifier, Transform3>{}};

  WhiteBoard firstStore;
  AlgorithmContext first{0, 0, firstStore, 0};
  first.recoGeoContext = AlignmentContext{&storeA};
  storeEvent(handles, first, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(first) == ProcessCode::SUCCESS);

  WhiteBoard secondStore;
  AlgorithmContext second{0, 1, secondStore, 0};
  second.recoGeoContext = AlignmentContext{&storeB};
  storeEvent(handles, second, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(20.)});
  BOOST_CHECK(algorithm.execute(second) == ProcessCode::ABORT);
  BOOST_CHECK(algorithm.finalize() == ProcessCode::ABORT);
  BOOST_CHECK_EQUAL(align->calls, 0u);
}

BOOST_AUTO_TEST_CASE(LeavesTheResultFileUntouchedWhenTheFitFails) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  align->fail = true;
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-keep";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);
  {
    std::ofstream existing(directory / "aligned_result.txt");
    existing << "keep-me\n";
  }

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard store;
  AlgorithmContext ctx{0, 0, store, 0};
  storeEvent(handles, ctx, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(ctx) == ProcessCode::SUCCESS);
  BOOST_CHECK(algorithm.finalize() == ProcessCode::ABORT);
  BOOST_CHECK_EQUAL(align->calls, 1u);

  std::ifstream result(directory / "aligned_result.txt");
  std::string line;
  BOOST_REQUIRE(std::getline(result, line));
  BOOST_CHECK_EQUAL(line, "keep-me");
  BOOST_CHECK(!std::filesystem::exists(directory / "index_map.txt"));

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(WritesTheResultWhenAlignmentDidNotConverge) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  align->markNotConverged = true;
  auto surface = Surface::makeShared<PlaneSurface>(
      Transform3::Identity(), std::make_shared<RectangleBounds>(1., 1.));
  align->indexedSurfaces.emplace(surface.get(), 0);
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-converge";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  WhiteBoard store;
  AlgorithmContext ctx{0, 0, store, 0};
  storeEvent(handles, ctx, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(ctx) == ProcessCode::SUCCESS);
  BOOST_CHECK(algorithm.finalize() == ProcessCode::SUCCESS);
  BOOST_CHECK_EQUAL(align->calls, 1u);

  std::ifstream result(directory / "aligned_result.txt");
  std::string line;
  BOOST_REQUIRE(std::getline(result, line));
  BOOST_CHECK(line.find("99") != std::string::npos);

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_CASE(InitializeDropsCollectedTracksAndTheStoreCheck) {
  auto geometry = makeGeometry();
  auto align = std::make_shared<RecordingAlignment>();
  const auto directory = std::filesystem::temp_directory_path() /
                         "acts-alignment-cross-event-initialize";
  std::filesystem::remove_all(directory);
  std::filesystem::create_directories(directory);

  AlignmentAlgorithm algorithm{
      configFor(align, geometry, directory),
      getDefaultLogger("AlignmentAlgorithmCrossEvent", Logging::INFO)};
  EventHandles handles{algorithm};

  GeoIdAlignmentStore storeA{
      std::unordered_map<GeometryIdentifier, Transform3>{}};
  GeoIdAlignmentStore storeB{
      std::unordered_map<GeometryIdentifier, Transform3>{}};

  WhiteBoard firstStore;
  AlgorithmContext first{0, 0, firstStore, 0};
  first.recoGeoContext = AlignmentContext{&storeA};
  storeEvent(handles, first, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(10.)});
  BOOST_CHECK(algorithm.execute(first) == ProcessCode::SUCCESS);
  BOOST_CHECK(algorithm.initialize() == ProcessCode::SUCCESS);

  WhiteBoard secondStore;
  AlgorithmContext second{0, 1, secondStore, 0};
  second.recoGeoContext = AlignmentContext{&storeB};
  storeEvent(handles, second, measurementsOf(GeometryIdentifier(), 1),
             protoTracks({0}), {trackAt(20.)});
  BOOST_CHECK(algorithm.execute(second) == ProcessCode::SUCCESS);
  BOOST_CHECK(algorithm.finalize() == ProcessCode::SUCCESS);
  BOOST_CHECK_EQUAL(align->calls, 1u);
  BOOST_REQUIRE_EQUAL(align->parameters.size(), 1u);
  BOOST_CHECK_EQUAL(align->parameters[0].parameters()[eBoundLoc0], 20.);

  std::filesystem::remove_all(directory);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
