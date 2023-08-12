// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"

#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
#include <random>  // --> random distributions
#include <algorithm> //  --> std::find

namespace Acts {
class TrackingGeometry;
class IMaterialDecorator;
}  // namespace Acts

namespace ActsExamples {
class IContextDecorator;
namespace Generic {
class GenericDetectorElement;
}  // namespace Generic
}  // namespace ActsExamples

namespace ActsExamples::Contextual {
class InternallyAlignedDetectorElement;
class InternalAlignmentDecorator;

class AlignedDetector {
 public:
  using ContextDecorators =
      std::vector<std::shared_ptr<ActsExamples::IContextDecorator>>;
  using TrackingGeometryPtr = std::shared_ptr<const Acts::TrackingGeometry>;

  struct Config : public GenericDetector::Config {
    /// Seed for the decorator random numbers.
    size_t seed = 1324354657;
    /// Size of a valid IOV.
    size_t iovSize = 100;
    /// Span until garbage collection is active.
    size_t flushSize = 200;
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

  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(
      const Config& cfg,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);

  std::vector<std::vector<std::shared_ptr<Generic::GenericDetectorElement>>>&
  detectorStore() {
    return m_detectorStore;
  }

 private:
  /// The Store of the detector elements (lifetime: job)
  std::vector<std::vector<std::shared_ptr<Generic::GenericDetectorElement>>>
      m_detectorStore;
  
    // helper function
  std::tuple<double, double, double, double> generateCorrelatedMisalignments(
      double sigmaInPlane, double sigmaOutPlane, double sigmaInRot,
      double sigmaOutRot, std::mt19937& rng);

 public:
  std::pair<TrackingGeometryPtr, ContextDecorators> finalize(
      const Config& cfg,
      std::shared_ptr<const Acts::IMaterialDecorator> mdecorator);

};

}  // namespace ActsExamples::Contextual


//TODO: add implementation of generateCorrelatedMisalignments --> outside the class ?
std::tuple<double, double, double, double>
ActsExamples::Contextual::AlignedDetector::generateCorrelatedMisalignments(
    double sigmaInPlane, double sigmaOutPlane, double sigmaInRot,
    double sigmaOutRot, std::mt19937& rng) {
  std::normal_distribution<double> inPlaneDistribution(0.0, sigmaInPlane);
  std::normal_distribution<double> outPlaneDistribution(0.0, sigmaOutPlane);
  std::normal_distribution<double> inRotDistribution(0.0, sigmaInRot);
  std::normal_distribution<double> outRotDistribution(0.0, sigmaOutRot);

  double inPlaneMisalignment = inPlaneDistribution(rng);
  double outPlaneMisalignment = outPlaneDistribution(rng);
  double inRotMisalignment = inRotDistribution(rng);
  double outRotMisalignment = outRotDistribution(rng);

  return std::make_tuple(inPlaneMisalignment, outPlaneMisalignment,
                         inRotMisalignment, outRotMisalignment);
}

// Add the implementation of finalize outside the class
std::pair<ActsExamples::Contextual::AlignedDetector::TrackingGeometryPtr,
          ActsExamples::Contextual::AlignedDetector::ContextDecorators>
ActsExamples::Contextual::AlignedDetector::finalize(
    const Config& cfg,
    std::shared_ptr<const Acts::IMaterialDecorator> mdecorator) {


  // local random number generator
  std::mt19937 rng(cfg.seed);

  // For each detector element apply misalignment
  for (auto& detectorElements : m_detectorStore) {
    // correlated misalignment values
    auto misalignments = generateCorrelatedMisalignments(
        cfg.sigmaInPlane, cfg.sigmaOutPlane, cfg.sigmaInRot, cfg.sigmaOutRot,
        rng);

    // 
    for (auto& detectorElement : detectorElements) {
      detectorElement->applyMisalignment(std::get<0>(misalignments),
                                         std::get<1>(misalignments),
                                         std::get<2>(misalignments),
                                         std::get<3>(misalignments));
    }
  }

}
