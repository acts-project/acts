// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryHierarchyMap.hpp"
#include "ActsExamples/Digitization/DigitizationConfig.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/Framework/BareAlgorithm.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsFatras/Digitization/Channelizer.hpp"
#include "ActsFatras/Digitization/PlanarSurfaceDrift.hpp"
#include "ActsFatras/Digitization/PlanarSurfaceMask.hpp"

#include <memory>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

namespace Acts {
class Surface;
class TrackingGeometry;
}  // namespace Acts

namespace ActsExamples {

/// Algorithm that turns simulated hits into measurements by truth smearing.
class DigitizationAlgorithm final : public BareAlgorithm {
 public:
  struct Config {
    /// Input collection of simulated hits.
    std::string inputSimHits;
    /// Output source links collection.
    std::string outputSourceLinks;
    /// Output measurements collection.
    std::string outputMeasurements;
    /// Output cluster collection.
    std::string outputClusters;
    /// Output collection to map measured hits to contributing particles.
    std::string outputMeasurementParticlesMap;
    /// Output collection to map measured hits to simulated hits.
    std::string outputMeasurementSimHitsMap;
    /// Tracking geometry required to access global-to-local transforms.
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
    /// Random numbers tool.
    std::shared_ptr<const RandomNumbers> randomNumbers = nullptr;
    /// The digitizers per GeometryIdentifier
    Acts::GeometryHierarchyMap<DigitizationConfig> digitizationConfigs;
  };

  /// Construct the smearing algorithm.
  ///
  /// @param cfg is the algorithm configuration
  /// @param lvl is the logging level
  DigitizationAlgorithm(Config cfg, Acts::Logging::Level lvl);

  /// Build measurement from simulation hits at input.
  ///
  /// @param ctx is the algorithm context with event information
  /// @return a process code indication success or failure
  ProcessCode execute(const AlgorithmContext& ctx) const final override;

 private:
  /// Nested struct for digitized parameters
  struct DigitizedParameters {
    std::vector<Acts::BoundIndices> indices = {};
    std::vector<Acts::ActsScalar> values = {};
    std::vector<Acts::ActsScalar> variances = {};

    Cluster cluster;
  };

  /// Helper method for the geometric channelizing part
  ///
  /// @param geoCfg is the geometric digitization configuration
  /// @param hit the Simultated hit
  /// @param surface the Surface on which this is supposed to happen
  /// @param gctx the Geometry context
  /// @param rng the Random number engine for the drift smearing
  ///
  /// @return the list of channels
  std::vector<ActsFatras::Channelizer::ChannelSegment> channelizing(
      const GeometricDigitizationConfig& geoCfg, const SimHit& hit,
      const Acts::Surface& surface, const Acts::GeometryContext& gctx,
      RandomEngine& rng) const;

  /// Helper method for creating digitized parameters from clusters
  ///
  /// @todo ADD random smearing
  /// @param geoCfg is the geometric digitization configuration
  /// @param channels are the input channels
  /// @param rng the Random number engine for the charge generation smearing
  ///
  /// @return the list of digitized parameters
  DigitizedParameters localParameters(
      const GeometricDigitizationConfig& geoCfg,
      const std::vector<ActsFatras::Channelizer::ChannelSegment>& channels,
      RandomEngine& rng) const;

  /// Helper method for created a measurement from digitized parameters
  ///
  /// @param dParams The digitized parameters of variable size
  /// @param isl The indexed source link for the measurement
  ///
  /// @return a variant measurement
  Measurement createMeasurement(const DigitizedParameters& dParams,
                                const IndexSourceLink& isl) const
      noexcept(false);

  /// Nested smearer struct that holds geometric digitizer and smearing
  /// Support up to 4 dimensions.
  template <size_t kSmearDIM>
  struct CombinedDigitizer {
    GeometricDigitizationConfig geometric;
    ActsFatras::BoundParametersSmearer<RandomEngine, kSmearDIM> smearing;
  };

  // Support max 4 digitization dimensions - either digital or smeared
  using Digitizer = std::variant<CombinedDigitizer<0>, CombinedDigitizer<1>,
                                 CombinedDigitizer<2>, CombinedDigitizer<3>,
                                 CombinedDigitizer<4>>;

  /// Configuration of the Algorithm
  Config m_cfg;
  /// Digitizers within geometry hierarchy
  Acts::GeometryHierarchyMap<Digitizer> m_digitizers;
  /// Geometric digtizers
  ActsFatras::PlanarSurfaceDrift m_surfaceDrift;
  ActsFatras::PlanarSurfaceMask m_surfaceMask;
  ActsFatras::Channelizer m_channelizer;

  /// Contruct the constituents of a measurement.
  ///
  /// @tparam kMeasDIM the full dimension of the measurement
  ///
  /// @param dParams the struct of arrays of parameters to be created
  ///
  /// @return a tuple of constituents for a measurement
  template <size_t kMeasDIM>
  std::tuple<std::array<Acts::BoundIndices, kMeasDIM>,
             Acts::ActsVector<kMeasDIM>, Acts::ActsSymMatrix<kMeasDIM>>
  measurementConstituents(const DigitizedParameters& dParams) const {
    std::array<Acts::BoundIndices, kMeasDIM> indices;
    Acts::ActsVector<kMeasDIM> par;
    Acts::ActsSymMatrix<kMeasDIM> cov =
        Acts::ActsSymMatrix<kMeasDIM>::Identity();
    for (Eigen::Index ei = 0; ei < static_cast<Eigen::Index>(kMeasDIM); ++ei) {
      indices[ei] = dParams.indices[ei];
      par[ei] = dParams.values[ei];
      cov(ei, ei) = dParams.variances[ei];
    }
    return {indices, par, cov};
  }

  /// Construct a fixed-size smearer from a configuration.
  ///
  /// It's templated on the smearing dimention given by @tparam kSmearDIM
  ///
  /// @param cfg Is the digitization configuration input
  ///
  /// @return a variant of a Digitizer
  template <size_t kSmearDIM>
  static Digitizer makeDigitizer(const DigitizationConfig& cfg) {
    CombinedDigitizer<kSmearDIM> impl;
    // Copy the geometric configuration
    impl.geometric = cfg.geometricDigiConfig;
    // Prepare the smearing configuration
    for (int i = 0; i < static_cast<int>(kSmearDIM); ++i) {
      impl.smearing.indices[i] = cfg.smearingDigiConfig.at(i).index;
      impl.smearing.smearFunctions[i] =
          cfg.smearingDigiConfig.at(i).smearFunction;
    }
    return impl;
  }
};

}  // namespace ActsExamples
