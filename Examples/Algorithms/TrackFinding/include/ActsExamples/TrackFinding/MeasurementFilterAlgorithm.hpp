// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

/// Records which measurements were consumed by a tracking pass and produces a
/// MeasurementSubset for the next pass.
///
/// All source-link indices in the output subset are in original-container
/// space, so the MeasurementParticlesMap from digitization is valid for truth
/// matching in every pass without any remapping.
///
/// Typical two-pass sequence:
///
///   Digitization -> MeasurementSubset("measurement_subset", all measurements)
///
///   Pass 1:  SpacePointMaker + TrackFindingAlgorithm
///                (inputMeasurementSubset = "measurement_subset")
///            -> MeasurementFilterAlgorithm(
///                   inputTracks             = "tracks",
///                   inputMeasurementSubset  = "measurement_subset",
///                   outputMeasurementSubset = "lrt_measurement_subset")
///
///   Pass 2:  SpacePointMaker + TrackFindingAlgorithm
///                (inputMeasurementSubset = "lrt_measurement_subset")
class MeasurementFilterAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input track collection (output of CKF or ambiguity resolution).
    std::string inputTracks;
    /// Input MeasurementSubset (initial full subset from digitization for pass
    /// 1; filtered subset from the previous MeasurementFilterAlgorithm for
    /// subsequent passes).
    std::string inputMeasurementSubset;
    /// Output MeasurementSubset for the next pass (indices in
    /// original-container space).
    std::string outputMeasurementSubset;
    /// If true, outlier track states are also excluded from the output subset.
    bool includeOutliers = false;
  };

  explicit MeasurementFilterAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<MeasurementSubset> m_inputMeasurementSubset{
      this, "InputMeasurementSubset"};
  WriteDataHandle<MeasurementSubset> m_outputMeasurementSubset{
      this, "OutputMeasurementSubset"};
};

}  // namespace ActsExamples
