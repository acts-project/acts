// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsPlugins/EDM4hep/PodioUtil.hpp"

#include <optional>
#include <unordered_map>

namespace ActsExamples {

/// ConversionHelper implementation for DD4hep-based detectors
///
/// This class provides the mapping between ACTS surfaces/sourcelinks
/// and their DD4hep CellID identifiers.
class DD4hepPodioConversionHelper
    : public ActsPlugins::PodioUtil::ConversionHelper {
 public:
  /// Constructor
  /// @param detector The DD4hep detector reference
  /// @param measurements The measurement container
  explicit DD4hepPodioConversionHelper(
      const DD4hepDetector& detector, const MeasurementContainer& measurements);

  std::optional<ActsPlugins::PodioUtil::Identifier> surfaceToIdentifier(
      const Acts::Surface& surface) const override;

  const Acts::Surface* identifierToSurface(
      ActsPlugins::PodioUtil::Identifier identifier) const override;

  Acts::SourceLink identifierToSourceLink(
      ActsPlugins::PodioUtil::Identifier identifier) const override;

  ActsPlugins::PodioUtil::Identifier sourceLinkToIdentifier(
      const Acts::SourceLink& sourceLink) const override;

 private:
  const DD4hepDetector* m_detector;
  const MeasurementContainer* m_measurements;
};

}  // namespace ActsExamples
