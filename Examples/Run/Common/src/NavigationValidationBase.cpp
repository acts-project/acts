// This file is part of the Acts project.
//
// Copyright (C) 2019-2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/SharedBField.hpp"
#include "Acts/Propagator/FSMNavigator.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "ActsExamples/Detector/IBaseDetector.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/MagneticField/MagneticFieldOptions.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Options/ParticleGunOptions.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Propagation/PropagatorInterface.hpp"
#include "ActsExamples/Utilities/Paths.hpp"

#include <memory>

#include <boost/program_options.hpp>

using namespace ActsExamples;
using namespace Acts::UnitLiterals;

#define CHECK(x)                                 \
  do {                                           \
    if (auto c = x; c != ProcessCode::SUCCESS) { \
      return c;                                  \
    }                                            \
  } while (false)

namespace {

/// This is a simple cache struct to mimic a Stepper
struct PseudoStepper {
  // comply with concept
  using Jacobian = Acts::BoundMatrix;
  using Covariance = Acts::BoundSymMatrix;
  using BoundState = std::tuple<Acts::BoundTrackParameters, Jacobian, double>;
  using CurvilinearState =
      std::tuple<Acts::CurvilinearTrackParameters, Jacobian, double>;
  using BField = int;

  template <typename, typename>
  using return_parameter_type = void;

  /// This is a simple cache struct to mimic the
  /// Stepper cache in the propagation
  struct State {
    /// Position
    Acts::Vector4 pos4 = Acts::Vector4(0., 0., 0., 0.);

    /// Direction
    Acts::Vector3 dir = Acts::Vector3(1., 0., 0.);

    /// Momentum
    double p;

    /// Charge
    double q;

    /// the navigation direction
    Acts::NavigationDirection navDir = Acts::NavigationDirection::Forward;

    // accummulated path length cache
    double pathAccumulated = 0.;

    // adaptive sep size of the runge-kutta integration
    Acts::ConstrainedStep stepSize = Acts::ConstrainedStep(100_cm);

    // Previous step size for overstep estimation (ignored here)
    double previousStepSize = 0.;

    /// The tolerance for the stepping
    double tolerance = Acts::s_onSurfaceTolerance;

    Acts::GeometryContext geoContext = Acts::GeometryContext();
  };

  double getStepSize(const State& state,
                     Acts::ConstrainedStep::Type stype) const {
    return state.stepSize.value(stype);
  }

  void releaseStepSize(State& state) const {
    state.stepSize.release(Acts::ConstrainedStep::actor);
  }

  void step(State& sstate, double fraction = 1) {
    // update the cache position
    double ssize = sstate.stepSize * fraction;
    Acts::Vector4 prev = sstate.pos4;
    sstate.pos4[Acts::ePos0] += ssize * sstate.dir[Acts::eMom0];
    sstate.pos4[Acts::ePos1] += ssize * sstate.dir[Acts::eMom1];
    sstate.pos4[Acts::ePos2] += ssize * sstate.dir[Acts::eMom2];

    std::cout << "PseudoStepper: Performing step with size: " << ssize
              << " along [" << sstate.dir.transpose() << "]: " << std::endl;

    auto rz = [](const Acts::Vector4& v) -> std::string {
      return std::to_string(Acts::VectorHelpers::perp(v)) + "," +
             std::to_string(v[Acts::eFreePos2]);
    };
    std::cout << "               [" << prev.transpose();
    std::cout << "] -> [" << sstate.pos4.transpose() << "]" << std::endl;

    std::cout << "               [" << rz(prev);
    std::cout << "] -> [" << rz(sstate.pos4.transpose()) << "]" << std::endl;
    // create navigation parameters
    return;
  }

  /// State resetter
  void resetState(State& /*unused*/, const Acts::BoundVector& /*unused*/,
                  const Acts::BoundSymMatrix& /*unused*/,
                  const Acts::Surface& /*unused*/,
                  const Acts::NavigationDirection /*unused*/,
                  const double /*unused*/) const {}

  /// Global particle position accessor
  Acts::Vector3 position(const State& state) const {
    return state.pos4.segment<3>(Acts::ePos0);
  }

  /// Time access
  double time(const State& state) const { return state.pos4[Acts::eTime]; }

  /// Momentum direction accessor
  Acts::Vector3 direction(const State& state) const { return state.dir; }

  /// Momentum accessor
  double momentum(const State& state) const { return state.p; }

  /// Charge access
  double charge(const State& state) const { return state.q; }

  /// Overstep limit access
  double overstepLimit(const State& /*state*/) const {
    return Acts::s_onSurfaceTolerance;
  }

  Acts::Intersection3D::Status updateSurfaceStatus(
      State& state, const Acts::Surface& surface,
      const Acts::BoundaryCheck& bcheck, Acts::LoggerWrapper logger) const {
    return Acts::detail::updateSingleSurfaceStatus<PseudoStepper>(
        *this, state, surface, bcheck, logger);
  }

  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      bool release = true) const {
    Acts::detail::updateSingleStepSize<PseudoStepper>(state, oIntersection,
                                                      release);
  }

  void setStepSize(
      State& state, double stepSize,
      Acts::ConstrainedStep::Type stype = Acts::ConstrainedStep::actor,
      bool release = true) const {
    state.previousStepSize = state.stepSize;
    state.stepSize.update(stepSize, stype, release);
  }

  std::string outputStepSize(const State& state) const {
    return state.stepSize.toString();
  }

  Acts::Result<BoundState> boundState(
      State& state, const Acts::Surface& surface, bool /*unused*/,
      const Acts::FreeToBoundCorrection& /*unused*/
  ) const {
    auto bound = Acts::BoundTrackParameters::create(
        surface.getSharedPtr(), state.geoContext, state.pos4, state.dir,
        state.p, state.q);
    if (!bound.ok()) {
      return bound.error();
    }
    BoundState bState{std::move(*bound), Jacobian::Identity(),
                      state.pathAccumulated};
    return bState;
  }

  CurvilinearState curvilinearState(State& state, bool /*unused*/
  ) const {
    Acts::CurvilinearTrackParameters parameters(state.pos4, state.dir, state.p,
                                                state.q);
    // Create the bound state
    CurvilinearState curvState{std::move(parameters), Jacobian::Identity(),
                               state.pathAccumulated};
    return curvState;
  }

  void update(State& /*state*/, const Acts::FreeVector& /*freePars*/,
              const Acts::BoundVector& /*boundPars*/, const Covariance& /*cov*/,
              const Acts::Surface& /*surface*/) const {}

  void update(State& /*state*/, const Acts::Vector3& /*uposition*/,
              const Acts::Vector3& /*udirection*/, double /*up*/,
              double /*time*/) const {}

  void transportCovarianceToCurvilinear(State& /*state*/) const {}

  void transportCovarianceToBound(
      State& /*unused*/, const Acts::Surface& /*surface*/,
      const Acts::FreeToBoundCorrection& /*freeToBoundCorrection*/ =
          Acts::FreeToBoundCorrection{false}) const {}

  Acts::Result<Acts::Vector3> getField(State& /*state*/,
                                       const Acts::Vector3& /*pos*/) const {
    // get the field from the cell
    return Acts::Result<Acts::Vector3>::success({0., 0., 0.});
  }
};

static_assert(Acts::StepperConcept<PseudoStepper>,
              "Dummy stepper does not fulfill concept");

/// This is a simple cache struct to mimic the
/// Propagator cache
template <typename navigator_t>
struct PropagatorState {
  /// emulate the options template
  struct Options {
    /// Debug output
    /// the string where debug messages are stored (optionally)
    bool debug = false;
    std::string debugString = "";
    /// buffer & formatting for consistent output
    size_t debugPfxWidth = 30;
    size_t debugMsgWidth = 50;

    Acts::LoggerWrapper logger{Acts::getDummyLogger()};
  };

  /// Navigation cache: the start surface
  const Acts::Surface* startSurface = nullptr;

  /// Navigation cache: the current surface
  const Acts::Surface* currentSurface = nullptr;

  /// Navigation cache: the target surface
  const Acts::Surface* targetSurface = nullptr;
  bool targetReached = false;

  /// Give some options
  Options options;

  /// The Stepper state - internal statew of the Stepper
  PseudoStepper::State stepping;

  /// Navigation state - internal state of the Navigator
  typename navigator_t::State navigation;

  // The context cache for this propagation
  Acts::GeometryContext geoContext = Acts::GeometryContext();
};

class NavigationValidationAlgorithm : public BareAlgorithm {
 public:
  struct Config {
    std::string inputParticles{"particles"};
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry;
  };
  NavigationValidationAlgorithm(const Config& config,
                                Acts::Logging::Level level)
      : BareAlgorithm{"NavValAlg", level}, m_cfg{config} {}

  ProcessCode execute(const AlgorithmContext& context) const override {
    const auto& particles =
        context.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);

    PseudoStepper stepper;
    Acts::Navigator navA{{m_cfg.trackingGeometry}};
    Acts::FSMNavigator navB{
        {m_cfg.trackingGeometry},
        Acts::getDefaultLogger("FSMNavigator", Acts::Logging::Level::VERBOSE)};
    auto llogger =
        Acts::getDefaultLogger("Navigator", Acts::Logging::Level::VERBOSE);

    ACTS_VERBOSE("Have " << particles.size() << " particles");
    for (const auto& particle : particles) {
      Acts::CurvilinearTrackParameters params{
          particle.fourPosition(), particle.unitDirection(),
          particle.absoluteMomentum(), particle.charge()};
      const auto* startSurface = &params.referenceSurface();

      ACTS_VERBOSE("- particle: " << particle.fourPosition().transpose() << " "
                                  << particle.unitDirection().transpose());

      PropagatorState<Acts::Navigator> stateA;
      stateA.navigation.startSurface = startSurface;
      stateA.options.logger = Acts::LoggerWrapper{*llogger};
      stateA.stepping.pos4 = particle.fourPosition();
      stateA.stepping.dir = particle.unitDirection();

      PropagatorState<Acts::FSMNavigator> stateB;
      stateB.navigation.startSurface = startSurface;
      stateB.options.logger = Acts::LoggerWrapper{*llogger};
      stateB.stepping.pos4 = particle.fourPosition();
      stateB.stepping.dir = particle.unitDirection();

      auto fmtObj = [](const auto* o) -> std::string {
        if (o == nullptr) {
          return "0x0";
        } else {
          std::stringstream sstr;
          sstr << o->geometryId();
          return sstr.str();
        }
      };

      size_t s = 1;

      auto checkConsistency = [&]() {
        ACTS_INFO("Checking navigator consistency at step " << s);
        auto* sA = stateA.navigation.currentSurface;
        auto* sB = stateB.navigation.currentSurface;

        if (sA != sB) {
          ACTS_ERROR("Navigation currentSurface inconsistent: "
                     << fmtObj(sA) << " <-> " << fmtObj(sB));
          return ProcessCode::ABORT;
        }
        auto* vA = stateA.navigation.currentVolume;
        auto* vB = stateB.navigation.currentVolume;

        if (vA != vB) {
          ACTS_ERROR("Navigation currentVolume inconsistent: "
                     << fmtObj(vA) << " <-> " << fmtObj(vB));
          return ProcessCode::ABORT;
        }

        if (stateA.stepping.pos4 != stateB.stepping.pos4) {
          ACTS_ERROR("Navigation inconsistent pos4: "
                     << stateA.stepping.pos4.transpose());
          ACTS_ERROR("                          vs. "
                     << stateB.stepping.pos4.transpose());
          return ProcessCode::ABORT;
        }
        if (stateA.stepping.dir != stateB.stepping.dir) {
          ACTS_ERROR("Navigation inconsistent dir: "
                     << stateA.stepping.dir.transpose());
          ACTS_ERROR("                          vs. "
                     << stateB.stepping.dir.transpose());
          return ProcessCode::ABORT;
        }
        ACTS_INFO("Navigators consistent at step " << s);
        s++;
        return ProcessCode::SUCCESS;
      };

      // auto step = [&](const auto& nav, auto& state, double f) {
      //   ACTS_VERBOSE("----> STATUS CALL");
      //   nav.status(state, stepper);
      //   ACTS_VERBOSE("----> TARGET CALL");
      //   nav.target(state, stepper);

      //   if (state.navigation.navigationBreak) {
      //     return;
      //   }

      //   ACTS_VERBOSE("----> Step with fraction " << f);
      //   stepper.step(state.stepping, f);
      // };

      size_t i = 0;
      while (!(stateA.navigation.navigationBreak ||
               stateB.navigation.navigationBreak)) {
        ACTS_VERBOSE("-------- Navigation iteration begin: #" << i
                                                              << " --------");
        ACTS_VERBOSE("----> STATUS CALL");
        ACTS_VERBOSE("------> Navigator A");
        navA.status(stateA, stepper);
        ACTS_VERBOSE("------> Navigator B");
        navB.status(stateB, stepper);

        CHECK(checkConsistency());

        if (stateA.navigation.navigationBreak ||
            stateB.navigation.navigationBreak) {
          break;
        }

        ACTS_VERBOSE("----> TARGET CALL");
        ACTS_VERBOSE("------> Navigator A");
        navA.target(stateA, stepper);
        ACTS_VERBOSE("------> Navigator B");
        navB.target(stateB, stepper);

        CHECK(checkConsistency());

        ACTS_VERBOSE("----> Step with fraction 0.5");
        ACTS_VERBOSE("------> Navigator A");
        stepper.step(stateA.stepping, 0.5);
        ACTS_VERBOSE("------> Navigator B");
        stepper.step(stateB.stepping, 0.5);

        CHECK(checkConsistency());

        ACTS_VERBOSE("----> STATUS CALL");
        ACTS_VERBOSE("------> Navigator A");
        navA.status(stateA, stepper);
        ACTS_VERBOSE("------> Navigator B");
        navB.status(stateB, stepper);

        CHECK(checkConsistency());

        ACTS_VERBOSE("----> TARGET CALL");
        ACTS_VERBOSE("------> Navigator A");
        navA.target(stateA, stepper);
        ACTS_VERBOSE("------> Navigator B");
        navB.target(stateB, stepper);

        CHECK(checkConsistency());

        stepper.step(stateA.stepping, 1.0);
        stepper.step(stateB.stepping, 1.0);

        CHECK(checkConsistency());

        ACTS_VERBOSE("-------- Navigation iteration end: #" << i
                                                            << " --------");
        i++;
      }

      if (stateA.navigation.navigationBreak !=
          stateB.navigation.navigationBreak) {
        ACTS_ERROR("Navigation breaks inconsistent");
        return ProcessCode::ABORT;
      }
    }

    return ProcessCode::SUCCESS;
  }

 private:
  Config m_cfg;
};
}  // namespace

int navigationValidation(int argc, char* argv[],
                         ActsExamples::IBaseDetector& detector) {
  // Setup and parse options
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addGeometryOptions(desc);
  ActsExamples::Options::addMaterialOptions(desc);
  ActsExamples::Options::addMagneticFieldOptions(desc);
  ActsExamples::Options::addRandomNumbersOptions(desc);
  ActsExamples::Options::addParticleGunOptions(desc);
  ActsExamples::Options::addOutputOptions(
      desc, ActsExamples::OutputFormat::Root | ActsExamples::OutputFormat::Obj);

  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = ActsExamples::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }
  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // The geometry, material and decoration
  auto geometry = ActsExamples::Geometry::build(vm, detector);
  auto tGeometry = geometry.first;
  auto contextDecorators = geometry.second;
  // Add the decorator to the sequencer
  for (auto cdr : contextDecorators) {
    sequencer.addContextDecorator(cdr);
  }

  // Create the random number engine
  auto randomNumberSvcCfg = ActsExamples::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<ActsExamples::RandomNumbers>(randomNumberSvcCfg);

  EventGenerator::Config evgen = Options::readParticleGunOptions(vm);
  evgen.outputParticles = "particles";
  evgen.randomNumbers = randomNumberSvc;
  sequencer.addReader(std::make_shared<EventGenerator>(evgen, logLevel));

  // Create BField service
  ActsExamples::Options::setupMagneticFieldServices(vm, sequencer);
  auto bField = ActsExamples::Options::readMagneticField(vm);

  NavigationValidationAlgorithm::Config valCfg;
  valCfg.trackingGeometry = tGeometry;
  sequencer.addAlgorithm(std::make_shared<NavigationValidationAlgorithm>(
      valCfg, Acts::Logging::VERBOSE));

  // Check what output exists, if none exists, the SteppingLogger
  // will switch to sterile.
  // bool rootOutput = vm["output-root"].template as<bool>();
  // bool objOutput = vm["output-obj"].template as<bool>();

  // ---------------------------------------------------------------------------------
  // Output directory
  // std::string outputDir = vm["output-dir"].template as<std::string>();
  // auto psCollection = vm["prop-step-collection"].as<std::string>();

  // if (rootOutput) {
  //   // Write the propagation steps as ROOT TTree
  //   ActsExamples::RootPropagationStepsWriter::Config pstepWriterRootConfig;
  //   pstepWriterRootConfig.collection = psCollection;
  //   pstepWriterRootConfig.filePath =
  //       ActsExamples::joinPaths(outputDir, psCollection + ".root");
  //   sequencer.addWriter(
  //       std::make_shared<ActsExamples::RootPropagationStepsWriter>(
  //           pstepWriterRootConfig));
  // }

  // if (objOutput) {
  //   using PropagationSteps = Acts::detail::Step;
  //   using ObjPropagationStepsWriter =
  //       ActsExamples::ObjPropagationStepsWriter<PropagationSteps>;

  //   // Write the propagation steps as Obj TTree
  //   ObjPropagationStepsWriter::Config pstepWriterObjConfig;
  //   pstepWriterObjConfig.collection = psCollection;
  //   pstepWriterObjConfig.outputDir = outputDir;
  //   sequencer.addWriter(
  //       std::make_shared<ObjPropagationStepsWriter>(pstepWriterObjConfig));
  // }

  return sequencer.run();
}
