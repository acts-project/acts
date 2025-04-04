// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

namespace ActsExamples {

class InternallyAlignedDetectorElement;
class InternalAlignmentDecorator;

class AlignedDetector : public GenericDetector {
 public:
  struct Config : public GenericDetector::Config {
    /// Seed for the decorator random numbers.
    unsigned int seed = 1324354657;
    /// Size of a valid IOV.
    unsigned int iovSize = 100;
    /// Span until garbage collection is active.
    unsigned int flushSize = 200;
    /// Run the garbage collection?
    bool doGarbageCollection = true;
    /// Sigma of the in-plane misalignment
    double sigmaInPlane = 100 * Acts::UnitConstants::um;
    /// Sigma of the out-of-plane misalignment
    double sigmaOutPlane = 50 * Acts::UnitConstants::um;
    /// Sigma of the in-plane rotation misalignment
    double sigmaInRot = 20 * 0.001;  // millirad
    /// Sigma of the out-of-plane rotation misalignment
    double sigmaOutRot = 0;
    /// Keep the first iov batch nominal.
    bool firstIovNominal = false;
    /// Log level for the decorator
    Acts::Logging::Level decoratorLogLevel = Acts::Logging::INFO;

    enum class Mode { Internal, External };
    Mode mode = Mode::Internal;
  };

  explicit AlignedDetector(const Config& cfg);

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
