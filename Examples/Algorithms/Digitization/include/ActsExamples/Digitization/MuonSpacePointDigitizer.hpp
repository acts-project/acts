// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

namespace ActsExamples {

/// Algorithm that turns simulated hits into measurements by truth smearing.
class MuonSpacePointDigitizer final : public IAlgorithm {
 public:
  struct Config {
    /// @brief Name of the input simulated hits collection
    std::string inputSimHits{"simhits"};
    /// @brief Name of the input simulated particles collection
    std::string inputParticles{"particles_simulated"};
    /// @brief Name of the output space points collection
    std::string outputSpacePoints{"MuonSpacePoints"};
    /// @brief Random number generator service
    std::shared_ptr<const RandomNumbers> randomNumbers{};
    /// @brief Pointer to the muon calibrator to fetch the smearing constants
    std::shared_ptr<const MuonSpacePointCalibrator> calibrator{};
    /// @brief Pointer to the tracking geometry to fetch the surfaces
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry{};
    /// @brief Add strip times
    bool digitizeTime{false};
    /// @brief Visualize the digitization
    bool dumpVisualization{true};
    /// @brief Applied dead time between two consecutive straw hits
    double strawDeadTime{1. * Acts::UnitConstants::ms};
    /// @brief Applied dead time between two consecutive rpc hits
    double rpcDeadTime{50. * Acts::UnitConstants::ns};
  };
  /// @brief Constructor
  MuonSpacePointDigitizer(const Config& cfg, Acts::Logging::Level lvl);

  /// @brief Initialize the digitizer
  ProcessCode initialize() override;
  /// @brief Execute the digitization
  ProcessCode execute(const AlgorithmContext& ctx) const override;
  /// @brief  Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief Returns the reference to the configured calibrator
  const MuonSpacePointCalibrator& calibrator() const;
  /// @brief Returns the reference to the passed tracking geometry
  const Acts::TrackingGeometry& trackingGeometry() const;
  /// @brief  Returns the transformation from the local hit frame into the
  ///         chamber's surface frame
  /// @param gctx Geometry context to access the local -> global transform of the surface
  /// @param hitId Geometry identifier of the hit of interest
  Acts::Transform3 toSpacePointFrame(
      const Acts::GeometryContext& gctx,
      const Acts::GeometryIdentifier& hitId) const;
  /// @brief Visualizes the digitized space point bucket and plots the
  ///        measurements on a y-z planes. Truth trajectories of muons
  ///        are also depicted on the canvas
  /// @param ctx: Context to fetch the sim particle container
  /// @param gctx: Geometry context to construct the toSpacePoint transforms
  /// @param bucket: All space points in the chamber to draw.
  void visualizeBucket(const AlgorithmContext& ctx,
                       const Acts::GeometryContext& gctx,
                       const MuonSpacePointBucket& bucket) const;
  /// @brief Determines the axis ranges in the z-y plane to draw
  ///        the measurements on the canvas
  /// @param bucket: List of all measurements in the chamber
  Acts::RangeXD<2, double> canvasRanges(
      const MuonSpacePointBucket& bucket) const;
  /// @brief Returns whether a surface should be drawn on the canvas
  /// @param gctx: Geometry context to construct the toSpacePoint transforms
  /// @param surface: Reference to the surface of question to depict
  /// @param canvasBoundaries: Visible section on the canvas
  bool isSurfaceToDraw(const Acts::GeometryContext& gctx,
                       const Acts::Surface& surface,
                       const Acts::RangeXD<2, double>& canvasBoundaries) const;
  /// @brief Configuration of the digitizer
  Config m_cfg;
  /// @brief Data handle for the input simulated hits
  ReadDataHandle<SimHitContainer> m_inputSimHits{this, "InputSimHits"};
  /// @brief Data handle for the input simulated particles
  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "G4Particles"};
  /// @brief Data handle for the output space points
  WriteDataHandle<MuonSpacePointContainer> m_outputSpacePoints{this,
                                                               "SpacePoints"};
};
}  // namespace ActsExamples
