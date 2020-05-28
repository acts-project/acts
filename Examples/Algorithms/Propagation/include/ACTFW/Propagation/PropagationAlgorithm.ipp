// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename propagator_t>
std::optional<Acts::BoundSymMatrix>
PropagationAlgorithm<propagator_t>::generateCovariance(
    FW::RandomEngine& rnd, std::normal_distribution<double>& gauss) const {
  if (m_cfg.covarianceTransport) {
    // We start from the correlation matrix
    Acts::BoundSymMatrix newCov(m_cfg.correlations);
    // Then we draw errors according to the error values
    Acts::BoundVector covs_smeared = m_cfg.covariances;
    for (size_t k = 0; k < size_t(covs_smeared.size()); ++k) {
      covs_smeared[k] *= gauss(rnd);
    }
    // and apply a double loop
    for (size_t i = 0; i < size_t(newCov.rows()); ++i) {
      for (size_t j = 0; j < size_t(newCov.cols()); ++j) {
        (newCov)(i, j) *= covs_smeared[i];
        (newCov)(i, j) *= covs_smeared[j];
      }
    }
    return newCov;
  }
  return std::nullopt;
}

template <typename propagator_t>
PropagationAlgorithm<propagator_t>::PropagationAlgorithm(
    const PropagationAlgorithm<propagator_t>::Config& cfg,
    Acts::Logging::Level loglevel)
    : BareAlgorithm("PropagationAlgorithm", loglevel), m_cfg(cfg) {}

/// Templated execute test method for
/// charged and netural particles
/// @param [in] context is the contextual data of this event
/// @param [in] startParameters the start parameters
/// @param [in] pathLength the maximal path length to go
template <typename propagator_t>
template <typename parameters_t>
PropagationOutput PropagationAlgorithm<propagator_t>::executeTest(
    const AlgorithmContext& context, const parameters_t& startParameters,
    double pathLength, FW::RandomEngine& rnd,
    std::normal_distribution<double>& gauss) const {
  ACTS_DEBUG("Test propagation/extrapolation starts");

  PropagationOutput pOutput;

  // This is the outside in mode
  if (m_cfg.mode == 0) {
    // The step length logger for testing & end of world aborter
    using MaterialInteractor = Acts::MaterialInteractor;
    using SteppingLogger = Acts::detail::SteppingLogger;
    using DebugOutput = Acts::DebugOutputActor;
    using EndOfWorld = Acts::EndOfWorldReached;
    using Generator = Acts::MeasurementGenerator;
    using Approacher = Acts::MeasurementApproacher;

    // Action list and abort list
    using ActionList = Acts::ActionList<Generator, SteppingLogger,
                                        MaterialInteractor, DebugOutput>;
    using AbortList = Acts::AbortList<EndOfWorld>;
    using PropagatorOptions =
        Acts::DenseStepperPropagatorOptions<ActionList, AbortList>;

    PropagatorOptions options(context.geoContext, context.magFieldContext);
    options.pathLimit = pathLength;
    options.debug = m_cfg.debugOutput;

    // Activate loop protection at some pt value
    options.loopProtection =
        (Acts::VectorHelpers::perp(startParameters.momentum()) <
         m_cfg.ptLoopers);

    // Set the radii for the generation
    auto& mGenerator = options.actionList.get<Generator>();
    mGenerator.radii = {200, 205, 210, 215, 220, 225, 230, 235, 240, 245};

    // Switch the material interaction on/off & eventually into logging mode
    auto& mInteractor = options.actionList.get<MaterialInteractor>();
    mInteractor.multipleScattering = m_cfg.multipleScattering;
    mInteractor.energyLoss = m_cfg.energyLoss;
    mInteractor.recordInteractions = m_cfg.recordMaterialInteractions;

    // Set a maximum step size
    options.maxStepSize = m_cfg.maxStepSize;

    // Propagate using the propagator
    const auto& result =
        m_cfg.propagator.propagate(startParameters, options).value();
    auto steppingResults = result.template get<SteppingLogger::result_type>();

    // Get the output generation result & create smeared space points
    std::vector<Acts::Vector3D> smearedSpacePoints;
    double smearT = 2.;
    double smearL = 5.;
    auto generatorResult = result.template get<Generator::result_type>();
    ACTS_DEBUG("Generator created " << generatorResult.spacepoints.size()
                                    << " measurements.");

    smearedSpacePoints.reserve(generatorResult.spacepoints.size());
    for (auto& pos : generatorResult.spacepoints) {
      Acts::Vector3D spos(pos.x() + smearT * gauss(rnd),
                          pos.y() + smearT * gauss(rnd),
                          pos.z() + smearL * gauss(rnd));
      ACTS_VERBOSE(" --> Measurement at radius = "
                   << Acts::VectorHelpers::perp(pos) << " created.");
      smearedSpacePoints.push_back(spos);
    }

    // Set the stepping result
    pOutput.first = std::move(steppingResults.steps);
    // Also set the material recording result - if configured
    if (m_cfg.recordMaterialInteractions) {
      auto materialResult =
          result.template get<MaterialInteractor::result_type>();
      pOutput.second = std::move(materialResult);
    }

    // screen output if requested
    if (m_cfg.debugOutput) {
      auto& debugResult = result.template get<DebugOutput::result_type>();
      ACTS_VERBOSE(debugResult.debugString);
    }

    if (not generatorResult.spacepoints.empty()) {
      // Action list and abort list
      using ActionRecList = Acts::ActionList<Approacher, SteppingLogger,
                                             MaterialInteractor, DebugOutput>;
      using PropagatorRecOptions =
          Acts::DenseStepperPropagatorOptions<ActionRecList, AbortList>;

      PropagatorRecOptions recOptions(context.geoContext,
                                      context.magFieldContext);
      recOptions.pathLimit = pathLength;
      recOptions.debug = m_cfg.debugOutput;

      // Activate loop protection at some pt value
      recOptions.loopProtection =
          (Acts::VectorHelpers::perp(startParameters.momentum()) <
           m_cfg.ptLoopers);

      // Set the radii for the generation
      auto& mApproacher = recOptions.actionList.get<Approacher>();
      mApproacher.measurements = smearedSpacePoints;

      // Switch the material interaction on/off & eventually into logging mode
      auto& mRecInteractor = recOptions.actionList.get<MaterialInteractor>();
      mRecInteractor.multipleScattering = m_cfg.multipleScattering;
      mRecInteractor.energyLoss = m_cfg.energyLoss;
      mRecInteractor.recordInteractions = m_cfg.recordMaterialInteractions;

      // Set a maximum step size
      recOptions.maxStepSize = m_cfg.maxStepSize;

      // Propagate using the propagator - for the reconstruction step
      const auto& recResult =
          m_cfg.propagator.propagate(startParameters, recOptions).value();

      auto approacherResult = recResult.template get<Approacher::result_type>();
      ACTS_DEBUG("Approacher found " << approacherResult.approaches.size()
                                     << " measurements.");
    
      if (approacherResult.approaches.size() ==
          mApproacher.measurements.size()) {
        for (size_t is = 0; is < approacherResult.approaches.size(); ++is) {
          auto approach = approacherResult.approaches[is];
          auto measurement = mApproacher.measurements[is];
          std::cout << approach.x() << "," << approach.y() << ","
                    << approach.z() << "," << measurement.x() << "," << measurement.y()
                    << "," << measurement.z() << ","
                    << (approach - measurement).norm() << std::endl;
        }
      }
    }
  }
  return pOutput;
}

template <typename propagator_t>
ProcessCode PropagationAlgorithm<propagator_t>::execute(
    const AlgorithmContext& context) const {
  // Create a random number generator
  FW::RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);

  // Standard gaussian distribution for covarianmces
  std::normal_distribution<double> gauss(0., 1.);

  // Setup random number distributions for some quantities
  std::uniform_real_distribution<double> phiDist(m_cfg.phiRange.first,
                                                 m_cfg.phiRange.second);
  std::uniform_real_distribution<double> etaDist(m_cfg.etaRange.first,
                                                 m_cfg.etaRange.second);
  std::uniform_real_distribution<double> ptDist(m_cfg.ptRange.first,
                                                m_cfg.ptRange.second);
  std::uniform_real_distribution<double> qDist(0., 1.);

  std::shared_ptr<const Acts::PerigeeSurface> surface =
      Acts::Surface::makeShared<Acts::PerigeeSurface>(
          Acts::Vector3D(0., 0., 0.));

  // Output : the propagation steps
  std::vector<std::vector<Acts::detail::Step>> propagationSteps;
  propagationSteps.reserve(m_cfg.ntests);

  // Output (optional): the recorded material
  std::vector<RecordedMaterialTrack> recordedMaterial;
  if (m_cfg.recordMaterialInteractions) {
    recordedMaterial.reserve(m_cfg.ntests);
  }

  // loop over number of particles
  for (size_t it = 0; it < m_cfg.ntests; ++it) {
    /// get the d0 and z0
    double d0 = m_cfg.d0Sigma * gauss(rng);
    double z0 = m_cfg.z0Sigma * gauss(rng);
    double phi = phiDist(rng);
    double eta = etaDist(rng);
    double theta = 2 * atan(exp(-eta));
    double pt = ptDist(rng);
    double p = pt / sin(theta);
    double charge = qDist(rng) > 0.5 ? 1. : -1.;
    double qop = charge / p;
    double t = m_cfg.tSigma * gauss(rng);
    // parameters
    Acts::BoundVector pars;
    pars << d0, z0, phi, theta, qop, t;
    // some screen output

    Acts::Vector3D sPosition(0., 0., 0.);
    Acts::Vector3D sMomentum(0., 0., 0.);

    // The covariance generation
    auto cov = generateCovariance(rng, gauss);

    // execute the test for charged particles
    PropagationOutput pOutput;
    if (charge) {
      // charged extrapolation - with hit recording
      Acts::BoundParameters startParameters(context.geoContext, std::move(cov),
                                            std::move(pars), surface);
      sPosition = startParameters.position();
      sMomentum = startParameters.momentum();
      pOutput = executeTest<Acts::TrackParameters>(context, startParameters,
                                                   1000., rng, gauss);
    } else {
      // execute the test for neeutral particles
      Acts::NeutralBoundParameters neutralParameters(
          context.geoContext, std::move(cov), std::move(pars), surface);
      sPosition = neutralParameters.position();
      sMomentum = neutralParameters.momentum();
      pOutput = executeTest<Acts::NeutralParameters>(context, neutralParameters,
                                                     1000., rng, gauss);
    }
    // Record the propagator steps
    propagationSteps.push_back(std::move(pOutput.first));
    if (m_cfg.recordMaterialInteractions &&
        pOutput.second.materialInteractions.size()) {
      // Create a recorded material track
      RecordedMaterialTrack rmTrack;
      // Start position
      rmTrack.first.first = std::move(sPosition);
      // Start momentum
      rmTrack.first.second = std::move(sMomentum);
      // The material
      rmTrack.second = std::move(pOutput.second);
      // push it it
      recordedMaterial.push_back(std::move(rmTrack));
    }
  }

  // Write the propagation step data to the event store
  context.eventStore.add(m_cfg.propagationStepCollection,
                         std::move(propagationSteps));

  // Write the recorded material to the event store
  if (m_cfg.recordMaterialInteractions) {
    context.eventStore.add(m_cfg.propagationMaterialCollection,
                           std::move(recordedMaterial));
  }

  return ProcessCode::SUCCESS;
}
