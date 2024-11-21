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
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include <cstddef>
#include <memory>

namespace ActsExamples {
class InternallyAlignedDetectorElement;
class InternalAlignmentDecorator;

class AlignedDetector : public DetectorBase {
 public:
  struct Config : public GenericDetector::Config {
    /// Seed for the decorator random numbers.
    std::size_t seed = 1324354657;
    /// Size of a valid IOV.
    std::size_t iovSize = 100;
    /// Span until garbage collection is active.
    std::size_t flushSize = 200;
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

    std::shared_ptr<const Acts::IMaterialDecorator> materialDecorator;
  };

  explicit AlignedDetector(const Config& cfg);

  Gen1GeometryHolder buildGen1Geometry() override;

 private:
  Config m_cfg;
};

}  // namespace ActsExamples
