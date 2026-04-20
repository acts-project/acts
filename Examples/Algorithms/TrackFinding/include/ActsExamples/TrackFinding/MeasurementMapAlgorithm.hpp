// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/EventData/MeasurementMap.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <string>

namespace ActsExamples {

/// Records which measurements were used by a tracking pass and accumulates
/// them into a MeasurementMap that downstream passes can use to skip those
/// measurements.
///
/// All indices stored in the output MeasurementMap are relative to the
/// original (unfiltered) MeasurementContainer produced by digitization.
/// If the tracks were found using a filtered container, supply
/// inputIndexRemapping so indices are translated back to original space.
///
/// Typical multi-pass sequence:
///   Pass N:  TrackFindingAlgorithm
///            -> MeasurementMapAlgorithm(inputMeasurementMap=passN-1Map,
///                                       inputIndexRemapping=passN-1Remap,
///                                       outputMeasurementMap=passNMap)
///   Pass N+1: MeasurementFilterAlgorithm(inputMeasurementMap=passNMap,
///                                        outputIndexRemapping=passNRemap)
///             -> TrackFindingAlgorithm
///             -> MeasurementMapAlgorithm(inputMeasurementMap=passNMap,
///                                        inputIndexRemapping=passNRemap, ...)
class MeasurementMapAlgorithm final : public IAlgorithm {
 public:
  struct Config {
    /// Input track collection (output of CKF or ambiguity resolution).
    std::string inputTracks;
    /// Optional input MeasurementMap from the previous tracking pass.
    /// When provided the output map is the union of this map and the
    /// measurements used in the current pass.
    std::string inputMeasurementMap;
    /// Optional index remapping produced by MeasurementFilterAlgorithm.
    /// Element i contains the original-container index of measurement i in
    /// the filtered container that was used for this tracking pass.
    /// Required when tracks were found using a filtered MeasurementContainer.
    std::string inputIndexRemapping;
    /// Output MeasurementMap for this pass (indices in original-container
    /// space).
    std::string outputMeasurementMap;
  };

  explicit MeasurementMapAlgorithm(
      Config cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ProcessCode execute(const AlgorithmContext& ctx) const final;

  const Config& config() const { return m_cfg; }

 private:
  Config m_cfg;

  ReadDataHandle<ConstTrackContainer> m_inputTracks{this, "InputTracks"};
  ReadDataHandle<MeasurementMap> m_inputMeasurementMap{this,
                                                       "InputMeasurementMap"};
  ReadDataHandle<MeasurementIndexRemapping> m_inputIndexRemapping{
      this, "InputIndexRemapping"};
  WriteDataHandle<MeasurementMap> m_outputMeasurementMap{this,
                                                         "OutputMeasurementMap"};
};

}  // namespace ActsExamples
