// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/MaterialInteraction.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/EventData/PropagationSummary.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/EventData/SimParticleDecay.hpp"
#include "ActsExamples/Framework/DataHandle.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Geant4/Geant4ConstructionOptions.hpp"
#include "ActsExamples/Geant4/SensitiveSurfaceMapper.hpp"

#include <cstddef>
#include <functional>
#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

class G4RunManager;
class G4VUserPrimaryGeneratorAction;
class G4VUserDetectorConstruction;
class G4UserRunAction;
class G4UserEventAction;
class G4UserTrackingAction;
class G4UserSteppingAction;
class G4VUserPhysicsList;

namespace Acts {
class Volume;
}  // namespace Acts

namespace ActsExamples {
struct Geant4Handle;

namespace Geant4 {
class EventContext;
struct EventResult;
class SimulationSession;
class NativeSimulationSession;
}  // namespace Geant4

/// Abstracts common Geant4 Acts algorithm behaviour.
class Geant4SimulationBase : public IAlgorithm {
 public:
  /// Nested configuration struct for the Geant4 simulation
  struct Config {
    // Name of the input particle collection
    std::string inputParticles = "";

    /// Random number service.
    std::shared_ptr<const RandomNumbers> randomNumbers;

    /// Geant4 construction options.
    Geant4ConstructionOptions constructionOptions;

    /// Detector instance to access Geant4 geometry construction.
    std::shared_ptr<Detector> detector;

    /// Optional Geant4 handle for sequential reuse with the same physics list.
    /// Each new configuration replaces the previous geometry and actions.
    /// Finish the active run before constructing a new configuration.
    std::shared_ptr<Geant4Handle> geant4Handle;

    /// Cap the largest acceptable step in the propagator
    /// if inf use Geant4 default (100m)
    double propagatorLargestAcceptableStep =
        std::numeric_limits<double>::infinity();  // mm
  };

  Geant4SimulationBase(const Config& cfg, const std::string& name,
                       std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~Geant4SimulationBase() override;

  /// Open one Geant4 run until finalize(). Geometry, materials, regions and
  /// cuts must remain fixed. Serial transport requires the initializing thread.
  ProcessCode initialize() final;

  /// Close the Geant4 run and invoke the end-of-run action once.
  ProcessCode finalize() final;

  /// Readonly access to the configuration
  virtual const Config& config() const = 0;

  std::shared_ptr<Geant4Handle> geant4Handle() const;

 protected:
  void commonInitialization();

  G4RunManager& runManager() const;

  std::shared_ptr<Geant4::EventContext> eventContext() const;

  /// Transport the framework input and return collections owned by this event.
  Geant4::EventResult simulateEvent(
      const AlgorithmContext& ctx,
      const SimParticleDecayGraph* decays = nullptr) const;

  int m_geant4Level{};

  std::unique_ptr<Geant4::SimulationSession> m_session;
  std::unique_ptr<Geant4::NativeSimulationSession> m_nativeSession;

  /// Detector construction object.
  /// G4RunManager will take care of deletion
  G4VUserDetectorConstruction* m_detectorConstruction{};

  ReadDataHandle<SimParticleContainer> m_inputParticles{this, "InputParticles"};
};

/// Algorithm to run Geant4 simulation in the ActsExamples framework
class Geant4Simulation final : public Geant4SimulationBase {
 public:
  using SensitiveCandidates = std::function<std::vector<const Acts::Surface*>(
      const Acts::GeometryContext&, const Acts::Vector3&)>;

  /// Create a stepping action for one simulation instance. During a stepping
  /// callback the context gives read-only access to the active event.
  /// Event references are valid only during that event on the session thread.
  using SteppingActionFactory =
      std::function<std::unique_ptr<G4UserSteppingAction>(
          std::shared_ptr<const Geant4::EventContext>)>;

  struct Config : public Geant4SimulationBase::Config {
    /// Number of native Geant4 transport workers. Zero selects the serial
    /// session.
    unsigned int numGeant4Workers = 0;

    /// Optional prescribed decays keyed by the barcodes of the input particles.
    std::string inputDecayGraph;

    /// Name of the output collection : hits
    std::string outputSimHits = "simhits";

    /// Name of the output collection : simulated particles
    std::string outputParticles = "particles_simulated";

    /// Name of the output collection : propagation records (debugging)
    std::string outputPropagationSummaries = "propagation_summaries";

    /// The ACTS sensitive surfaces in a mapper, used for hit creation
    std::shared_ptr<const Geant4::SensitiveSurfaceMapper>
        sensitiveSurfaceMapper = nullptr;

    /// The ACTS Magnetic field provider
    std::shared_ptr<const Acts::MagneticFieldProvider> magneticField = nullptr;

    /// If a physics list has to be instantiated this one is chosen.
    std::string physicsList = "FTFP_BERT";

    std::shared_ptr<const Acts::Volume> killVolume;
    double killAfterTime = std::numeric_limits<double>::infinity();
    bool killSecondaries = false;

    bool recordHitsOfCharged = true;

    bool recordHitsOfNeutrals = false;

    bool recordHitsOfPrimaries = true;

    bool recordHitsOfSecondaries = true;

    /// Prescribed parents are always retained to preserve the decay ancestry.
    bool keepParticlesWithoutHits = true;

    bool recordPropagationSummaries = false;

    /// Invoked during initialization, once per native worker or serial session.
    /// Serial retries may invoke the factories again. Each factory creates
    /// an independent action. Actions run in configuration
    /// order after the built-in particle-kill and hit actions.
    std::vector<SteppingActionFactory> steppingActionFactories;
  };

  /// Simulation constructor
  ///
  /// @param config is the configuration struct
  /// @param level is the logging level to be used
  explicit Geant4Simulation(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~Geant4Simulation() override;

  /// Algorithm execute method, called once per event with context
  ///
  /// @param ctx the AlgorithmContext for this event
  ProcessCode execute(const AlgorithmContext& ctx) const final;

  /// Readonly access to the configuration
  const Config& config() const final { return m_cfg; }

 private:
  Config m_cfg;

  WriteDataHandle<SimParticleContainer> m_outputParticles{this,
                                                          "OutputParticles"};
  ReadDataHandle<SimParticleDecayGraph> m_inputDecayGraph{this,
                                                          "InputDecayGraph"};

  WriteDataHandle<SimHitContainer> m_outputSimHits{this, "OutputSimHIts"};

  WriteDataHandle<PropagationSummaries> m_outputPropagationSummaries{
      this, "OutputPropagationSummaries"};
};

class Geant4MaterialRecording final : public Geant4SimulationBase {
 public:
  struct Config : public Geant4SimulationBase::Config {
    /// Name of the output collection: material tracks
    std::string outputMaterialTracks = "material_tracks";

    /// Materials to exclude from the recording.
    std::vector<std::string> excludeMaterials = {"Air", "Vacuum"};
    bool recordElementFractions = false;
  };

  /// Material recording constructor
  ///
  /// @param config is the configuration struct
  /// @param level is the logging level to be used
  explicit Geant4MaterialRecording(
      const Config& cfg, std::unique_ptr<const Acts::Logger> logger = nullptr);

  ~Geant4MaterialRecording() override;

  /// Algorithm execute method, called once per event with context
  ///
  /// @param ctx the AlgorithmContext for this event
  ProcessCode execute(const AlgorithmContext& ctx) const override;

  /// Readonly access to the configuration
  const Config& config() const final { return m_cfg; }

 private:
  Config m_cfg;

  WriteDataHandle<std::unordered_map<std::size_t, Acts::RecordedMaterialTrack>>
      m_outputMaterialTracks{this, "OutputMaterialTracks"};
};

}  // namespace ActsExamples
