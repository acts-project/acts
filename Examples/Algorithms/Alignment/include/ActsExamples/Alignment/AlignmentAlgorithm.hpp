// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/TrackFitting/KalmanFitter.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsAlignment/Kernel/Alignment.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"

#include <filesystem>
#include <map>
#include <memory>
#include <mutex>
#include <optional>
#include <vector>

namespace Acts {
class MagneticFieldProvider;
}

namespace ActsExamples {

class IAlignmentStore;

class AlignmentGroup {
 public:
  AlignmentGroup(const std::string& name,
                 const std::vector<Acts::GeometryIdentifier>& geoIds)
      : m_name(name), m_map(constructHierarchyMap(geoIds)) {}

  // Access the name of the group
  std::string getNameOfGroup() const { return m_name; }

  // Useful for testing
  bool has(Acts::GeometryIdentifier geoId) {
    auto it = m_map.find(geoId);
    return (it == m_map.end()) ? false : *it;
  }

 private:
  std::string m_name;  //  storing the name in the class
  Acts::GeometryHierarchyMap<bool> m_map;

  Acts::GeometryHierarchyMap<bool> constructHierarchyMap(
      const std::vector<Acts::GeometryIdentifier>& geoIds) {
    std::vector<Acts::GeometryHierarchyMap<bool>::InputElement> ies;
    for (const auto& geoId : geoIds) {
      ies.emplace_back(geoId, true);
    }
    return Acts::GeometryHierarchyMap<bool>(ies);
  }
};

class AlignmentAlgorithm final : public IAlgorithm {
 public:
  using AlignmentResult = Acts::Result<ActsAlignment::AlignmentResult>;
  /// Alignment function that takes sets of input measurements, initial
  /// trackstate and alignment options and returns some alignment-specific
  /// result.
  using TrackFitterOptions =
      Acts::KalmanFitterOptions<Acts::VectorMultiTrajectory>;

  /// Alignment function that takes the above parameters and runs alignment
  /// @note This is separated into a virtual interface to keep compilation units
  /// small
  class AlignmentFunction {
   public:
    virtual ~AlignmentFunction() = default;
    virtual AlignmentResult operator()(
        const std::vector<std::vector<IndexSourceLink>>&,
        const TrackParametersContainer&,
        const ActsAlignment::AlignmentOptions<TrackFitterOptions>&) const = 0;
  };

  /// Create the alignment function implementation.
  ///
  /// The magnetic field is intentionally given by-value since the variant
  /// contains shared_ptr anyway.
  static std::shared_ptr<AlignmentFunction> makeAlignmentFunction(
      std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
      std::shared_ptr<const Acts::MagneticFieldProvider> magneticField,
      Acts::Logging::Level logLevel = Acts::Logging::INFO);

  struct Config {
    /// Input measurements collection.
    std::string inputMeasurements;
    /// Input proto tracks collection, i.e. groups of hit indices.
    std::string inputProtoTracks;
    /// Input initial track parameter estimates for each proto track.
    std::string inputInitialTrackParameters;
    /// Millepede-style 6-column alignment result. Unset skips the file.
    /// Written after the fit returns successfully.
    std::optional<std::filesystem::path> outputAlignmentFile;
    /// Surface-index map for the labels in @c outputAlignmentFile. Unset
    /// skips the file. Written after the fit returns successfully.
    std::optional<std::filesystem::path> outputAlignmentIndexFile;
    /// Geometry used to resolve source-link surfaces. Pass the same geometry
    /// as @c makeAlignmentFunction.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
    /// Type erased fitter function.
    std::shared_ptr<AlignmentFunction> align;
    /// The aligned transform updater
    ActsAlignment::AlignedTransformUpdater alignedTransformUpdater;
    /// The surfaces (with detector elements) to be aligned
    std::vector<Acts::SurfacePlacementBase*> alignedDetElements;
    /// The alignment mask at each iteration
    std::map<unsigned int, std::bitset<6>> iterationState;
    /// Cutoff value for average chi2/ndf
    double chi2ONdfCutOff = 0.10;
    /// Cutoff value for delta of average chi2/ndf within a couple of iterations
    std::pair<std::size_t, double> deltaChi2ONdfCutOff = {10, 0.00001};
    /// Maximum number of iterations
    std::size_t maxNumIterations = 100;
    /// Track cap for the fit. Non-positive disables it. The kept sample is the
    /// leading tracks in event-number order; within one event, its leading
    /// tracks. Events are held until @c finalize(), including when they run
    /// in parallel.
    int maxNumTracks = -1;
    std::vector<AlignmentGroup> m_groups;
  };

  /// Constructor of the alignment algorithm
  ///
  /// @param cfg is the config struct to configure the algorithm
  /// @param logger is the logging instance
  explicit AlignmentAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  /// Copy this event's measurements and track parameters.
  /// The fit runs once in @c finalize().
  ///
  /// @param ctx is the algorithm context that holds event-wise information
  /// @return a process code to steer the algorithm flow
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Fit the kept sample and write the alignment result.
  ProcessCode finalize() override;

  /// Drop tracks and the alignment-store check from a previous run.
  ProcessCode initialize() override;

 private:
  struct CollectedChunk {
    Acts::GeometryContext geoContext{
        Acts::GeometryContext::dangerouslyDefaultConstruct()};
    Acts::MagneticFieldContext magFieldContext;
    Acts::CalibrationContext calibContext;
    MeasurementContainer measurements;
    std::vector<std::vector<IndexSourceLink>> sourceLinks;
    TrackParametersContainer parameters;
  };

  /// Apply the track cap and merge. Contexts come from the lowest event.
  CollectedChunk takeCollectedSample();

  Config m_cfg;

  ReadDataHandle<MeasurementContainer> m_inputMeasurements{this,
                                                           "InputMeasurements"};
  ReadDataHandle<TrackParametersContainer> m_inputInitialTrackParameters{
      this, "InputInitialTrackParameters"};
  ReadDataHandle<ProtoTrackContainer> m_inputProtoTracks{this,
                                                         "InputProtoTracks"};

  // execute() may run concurrently. The mutex guards the chunks and the
  // alignment-store check. initialize() clears both. finalize() runs after
  // event processing.
  mutable std::mutex m_collectionMutex;
  mutable std::map<std::size_t, CollectedChunk> m_collectedChunks;
  mutable const IAlignmentStore* m_alignmentStore = nullptr;
  mutable bool m_alignmentStoreSeen = false;
  mutable bool m_alignmentStoreMismatch = false;
};

}  // namespace ActsExamples
