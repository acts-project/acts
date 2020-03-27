// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "FatrasSimulationBase.hpp"

#include <boost/program_options.hpp>

#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Fatras/FatrasAlgorithm.hpp"
#include "ACTFW/Fatras/FatrasOptions.hpp"
#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/Framework/Sequencer.hpp"
#include "ACTFW/Generators/FlattenEvent.hpp"
#include "ACTFW/Generators/ParticleSelector.hpp"
#include "ACTFW/Io/Csv/CsvParticleWriter.hpp"
#include "ACTFW/Io/Root/RootParticleWriter.hpp"
#include "ACTFW/Io/Root/RootSimHitWriter.hpp"
#include "ACTFW/Options/CommonOptions.hpp"
#include "ACTFW/Plugins/BField/BFieldOptions.hpp"
#include "ACTFW/Plugins/BField/ScalableBField.hpp"
#include "ACTFW/Utilities/Paths.hpp"
#include "Acts/Geometry/GeometryID.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/InterpolatedBFieldMap.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsFatras/Kernel/PhysicsList.hpp"
#include "ActsFatras/Kernel/Process.hpp"
#include "ActsFatras/Kernel/Simulator.hpp"
#include "ActsFatras/Physics/EnergyLoss/BetheBloch.hpp"
#include "ActsFatras/Physics/EnergyLoss/BetheHeitler.hpp"
#include "ActsFatras/Physics/Scattering/Highland.hpp"
#include "ActsFatras/Selectors/ChargeSelectors.hpp"
#include "ActsFatras/Selectors/KinematicCasts.hpp"
#include "ActsFatras/Selectors/SelectorHelpers.hpp"

namespace {

/// Simple struct to select surfaces where hits should be generated.
struct HitSurfaceSelector {
  bool sensitive = false;
  bool material = false;
  bool passive = false;

  /// Check if the surface should be used.
  bool operator()(const Acts::Surface& surface) const {
    if (sensitive and surface.associatedDetectorElement()) {
      return true;
    }
    if (material and surface.surfaceMaterial()) {
      return true;
    }
    if (passive) {
      return true;
    }
    return false;
  }
};

/// @brief Simulation setup for the FatrasAlgorithm
///
/// @tparam bfield_t Type of the bfield for the simulation to be set up
///
/// @param fieldMap The field map for the simulation setup
/// @param sequencer The framework sequencer
/// @param variables The boost variable map to resolve
/// @param tGeometry The TrackingGeometry for the tracking setup
/// @param barcodesSvc The barcode service to be used for the simulation
/// @param randomNumberSvc The random number service to be used for the
/// simulation
template <typename magnetic_field_t>
void setupSimulationAlgorithms(
    const FW::Options::Variables& variables, FW::Sequencer& sequencer,
    std::shared_ptr<const FW::RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry,
    magnetic_field_t&& magneticField) {
  // Read the log level
  Acts::Logging::Level logLevel = FW::Options::readLogLevel(variables);

  // Convert generated events to selected particles
  auto select = FW::ParticleSelector::readConfig(variables);
  select.inputEvent = "event_generated";
  select.outputEvent = "event_selected";
  sequencer.addAlgorithm(
      std::make_shared<FW::ParticleSelector>(select, logLevel));
  FW::FlattenEvent::Config flatten;
  flatten.inputEvent = select.outputEvent;
  flatten.outputParticles = "particles_selected";
  sequencer.addAlgorithm(std::make_shared<FW::FlattenEvent>(flatten, logLevel));

  // setup propagator-related types
  // use the default navigation
  using Navigator = Acts::Navigator;
  // propagate charged particles numerically in the given magnetic field
  using ChargedStepper = Acts::EigenStepper<magnetic_field_t>;
  using ChargedPropagator = Acts::Propagator<ChargedStepper, Navigator>;
  // propagate neutral particles with just straight lines
  using NeutralStepper = Acts::StraightLineStepper;
  using NeutralPropagator = Acts::Propagator<NeutralStepper, Navigator>;

  // setup simulator-related types
  using MinP = ActsFatras::Min<ActsFatras::Casts::P>;
  // charged particles w/ standard em physics list and selectable hits
  using ChargedSelector =
      ActsFatras::CombineAnd<ActsFatras::ChargedSelector, MinP>;
  using ChargedSimulator = ActsFatras::ParticleSimulator<
      ChargedPropagator, ActsFatras::ChargedElectroMagneticPhysicsList,
      HitSurfaceSelector>;
  // neutral particles w/o physics and no hits
  using NeutralSelector =
      ActsFatras::CombineAnd<ActsFatras::NeutralSelector, MinP>;
  using NeutralSimulator = ActsFatras::ParticleSimulator<
      NeutralPropagator, ActsFatras::PhysicsList<>, ActsFatras::NoSurface>;
  // full simulator type for charged and neutrals
  using Simulator = ActsFatras::Simulator<ChargedSelector, ChargedSimulator,
                                          NeutralSelector, NeutralSimulator>;
  // final algorihm type
  using SimulationAlgorithm = FW::FatrasAlgorithm<Simulator>;

  // construct the simulator
  Navigator navigator(trackingGeometry);
  // construct the charged simulator
  ChargedStepper chargedStepper(std::move(magneticField));
  ChargedPropagator chargedPropagator(std::move(chargedStepper), navigator);
  ChargedSimulator chargedSimulator(std::move(chargedPropagator), logLevel);
  // construct the neutral simulator
  NeutralStepper neutralStepper;
  NeutralPropagator neutralPropagator(std::move(neutralStepper), navigator);
  NeutralSimulator neutralSimulator(std::move(neutralPropagator), logLevel);
  // construct the combined simulator
  Simulator simulator(std::move(chargedSimulator), std::move(neutralSimulator));

  // construct/add the simulation algorithm
  auto fatras = FW::Options::readFatrasConfig(variables, std::move(simulator));
  fatras.inputParticles = flatten.outputParticles;
  fatras.outputParticlesInitial = "particles_initial";
  fatras.outputParticlesFinal = "particles_final";
  fatras.outputHits = "hits";
  fatras.randomNumbers = randomNumbers;
  sequencer.addAlgorithm(
      std::make_shared<SimulationAlgorithm>(fatras, logLevel));

  // Output directory
  const auto outputDir = FW::ensureWritableDirectory(
      variables["output-dir"].template as<std::string>());

  // Write simulation information as CSV files
  if (variables["output-csv"].template as<bool>()) {
    FW::CsvParticleWriter::Config writeInitial;
    writeInitial.inputParticles = fatras.outputParticlesInitial;
    writeInitial.outputDir = outputDir;
    writeInitial.outputStem = fatras.outputParticlesInitial;
    sequencer.addWriter(
        std::make_shared<FW::CsvParticleWriter>(writeInitial, logLevel));
    FW::CsvParticleWriter::Config writeFinal;
    writeFinal.inputParticles = fatras.outputParticlesFinal;
    writeFinal.outputDir = outputDir;
    writeFinal.outputStem = fatras.outputParticlesFinal;
    sequencer.addWriter(
        std::make_shared<FW::CsvParticleWriter>(writeFinal, logLevel));
  }

  // Write simulation information as ROOT files
  if (variables["output-root"].template as<bool>()) {
    // write initial simulated particles
    FW::RootParticleWriter::Config writeInitial;
    writeInitial.inputParticles = fatras.outputParticlesInitial;
    writeInitial.filePath =
        FW::joinPaths(outputDir, fatras.outputParticlesInitial + ".root");
    sequencer.addWriter(
        std::make_shared<FW::RootParticleWriter>(writeInitial, logLevel));

    // write final simulated particles
    FW::RootParticleWriter::Config writeFinal;
    writeFinal.inputParticles = fatras.outputParticlesFinal;
    writeFinal.filePath =
        FW::joinPaths(outputDir, fatras.outputParticlesFinal + ".root");
    sequencer.addWriter(
        std::make_shared<FW::RootParticleWriter>(writeFinal, logLevel));

    // write simulated hits
    FW::RootSimHitWriter::Config writeHits;
    writeHits.inputSimulatedHits = fatras.outputHits;
    writeHits.filePath = FW::joinPaths(outputDir, fatras.outputHits + ".root");
    sequencer.addWriter(
        std::make_shared<FW::RootSimHitWriter>(writeHits, logLevel));
  }
}

}  // namespace

void FW::setupSimulation(
    const FW::Options::Variables& variables, FW::Sequencer& sequencer,
    std::shared_ptr<const RandomNumbers> randomNumbers,
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry) {
  auto magneticFieldVariant = FW::Options::readBField(variables);
  std::visit(
      [&](auto&& inputField) {
        using magnetic_field_t =
            typename std::decay_t<decltype(inputField)>::element_type;
        Acts::SharedBField<magnetic_field_t> magneticField(inputField);
        setupSimulationAlgorithms(variables, sequencer, randomNumbers,
                                  trackingGeometry, std::move(magneticField));
      },
      magneticFieldVariant);
}
