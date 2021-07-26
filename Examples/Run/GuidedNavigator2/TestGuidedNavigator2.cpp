#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/detail/SteppingLogger.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsExamples/TelescopeDetector/BuildTelescopeDetector.hpp"

#include <vector>

#include "Guiders.hpp"
#include "GuidedNavigator2.hpp"

using namespace Acts::UnitLiterals;

const Acts::GeometryContext gctx;
const Acts::MagneticFieldContext mctx;

int main() {
  // Magnetic Field
  auto magField =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3(0.0, 0.0, 0.0));

  // Detector
  const typename ActsExamples::Telescope::TelescopeDetectorElement::ContextType
      detectorContext;
  std::vector<
      std::shared_ptr<ActsExamples::Telescope::TelescopeDetectorElement>>
      detectorElementStorage;
  const std::vector<double> distances = {100_mm, 200_mm, 300_mm,
                                         400_mm, 500_mm, 600_mm};
  const std::array<double, 2> offsets = {0.0_mm, 0.0_mm};
  const std::array<double, 2> bounds = {100._mm, 100._mm};
  const double thickness = 10._mm;
  const auto type = ActsExamples::Telescope::TelescopeSurfaceType::Plane;
  const auto detectorDirection = Acts::BinningValue::binX;

  auto detector = std::shared_ptr(ActsExamples::Telescope::buildDetector(
      detectorContext, detectorElementStorage, distances, offsets, bounds,
      thickness, type, detectorDirection));

  // Make a start surface
  auto startBounds = std::make_shared<Acts::RectangleBounds>(10, 10);
  Acts::Transform3 trafo =
      Acts::Transform3::Identity() *
      Eigen::AngleAxisd(0.5 * M_PI, Eigen::Vector3d::UnitY());
  auto startSurface =
      Acts::Surface::makeShared<Acts::PlaneSurface>(trafo, startBounds);

  // Find Target surface
  const Acts::Surface *targetSurface = nullptr;
  detector->visitSurfaces([&](auto surface) {
    if (surface->center(gctx)[0] == 600._mm) {
      targetSurface = surface;
    }
  });

  throw_assert(targetSurface, "must have a target surface");

  // Make Surface sequence
  std::vector<const Acts::Surface *> surfaceSequence;
  detector->visitSurfaces(
      [&](auto surface) { surfaceSequence.push_back(surface); });

  std::sort(surfaceSequence.begin(), surfaceSequence.end(),
            [](auto s1, auto s2) {
              return s1->center(gctx)[0] < s2->center(gctx)[0];
            });

  std::cout << "surfaceSequence.size() = " << surfaceSequence.size() << "\n";
  throw_assert(surfaceSequence.size() > 0, "need some surfaceSequence");
  // Logger
  const auto logger = Acts::getDefaultLogger("Single", Acts::Logging::VERBOSE);

  // Stepper
  using Stepper = Acts::EigenStepper<>;

  // Start parameters
  const double l0{0.}, l1{0.}, theta{0.5 * M_PI}, phi{0.}, p{50._GeV}, q{-1.},
      t{0.};
  Acts::BoundVector pars;
  pars << l0, l1, phi, theta, q / p, t;

  Acts::BoundTrackParameters startPars(startSurface, pars, std::nullopt);

  // Store SteppingLogger results
  std::array<std::vector<Acts::detail::Step>, 3> r;

  ///////////////////////////////
  // Run with normal Navigator
  ///////////////////////////////
  {
    Acts::Navigator::Config cfg;
    cfg.trackingGeometry = detector;
    Acts::Navigator navigator(cfg);

    using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
    Propagator prop(Stepper(magField), navigator);

    using Actors = Acts::ActionList<Acts::detail::SteppingLogger>;
    using Aborters = Acts::AbortList<Acts::EndOfWorldReached>;

    Acts::PropagatorOptions<Actors, Aborters> propOptions(
        gctx, mctx, Acts::LoggerWrapper(*logger));

    auto res = prop.propagate(startPars, propOptions);

    r[0] = (*res).get<Acts::detail::SteppingLogger::result_type>().steps;
  }

  std::cout << "===========================================================\n";

  /////////////////////////////
  // Run with Guided - Direct
  /////////////////////////////
  {
    using Navigator = Acts::GuidedNavigator2<Acts::DirectGuider>;
    Navigator navigator;

    using Propagator = Acts::Propagator<Stepper, Navigator>;
    Propagator prop(Stepper(magField), navigator);

    using Actors =
        Acts::ActionList<Navigator::Initializer, Acts::detail::SteppingLogger>;
    using Aborters = Acts::AbortList<>;

    Acts::PropagatorOptions<Actors, Aborters> propOptions(
        gctx, mctx, Acts::LoggerWrapper(*logger));

    auto &dInitializer =
        propOptions.actionList.template get<Navigator::Initializer>();
    dInitializer.navSurfaces = surfaceSequence;

    auto res = prop.propagate(startPars, *targetSurface, propOptions);

    r[1] = (*res).get<Acts::detail::SteppingLogger::result_type>().steps;
  }

  std::cout << "==========================================================\n";

  ////////////////////////////////////
  // Run with Guided - TrialAndError
  ////////////////////////////////////
  {
    using Navigator = Acts::GuidedNavigator2<Acts::TrialAndErrorGuider>;
    Navigator navigator;
    navigator.guider.possibleSurfaces = surfaceSequence;

    using Propagator = Acts::Propagator<Stepper, Navigator>;
    Propagator prop(Stepper(magField), navigator);

    using Actors =
        Acts::ActionList<Navigator::Initializer, Acts::detail::SteppingLogger>;
    using Aborters = Acts::AbortList<>;

    Acts::PropagatorOptions<Actors, Aborters> propOptions(
        gctx, mctx, Acts::LoggerWrapper(*logger));

    auto &dInitializer =
        propOptions.actionList.template get<Navigator::Initializer>();
    dInitializer.navSurfaces = surfaceSequence;

    auto res = prop.propagate(startPars, *targetSurface, propOptions);

    r[2] = (*res).get<Acts::detail::SteppingLogger::result_type>().steps;
  }
  
  std::cout << "\nStandard navigator\t\tDirect-Guided navigator\t\tTrial&Error navigator\n";
  std::cout << "------------------------------------------------------------------------\n";
  
  for(const auto &step : r[0])
  {
    if( !step.surface )
        continue;
    
    std::cout << step.surface->geometryId() << "\t\t";
      
    const auto inDirect = std::find_if(r[1].begin(), r[1].end(), [&](const auto &s){ return s.surface && step.surface->geometryId() == s.surface->geometryId(); });
    const auto inTrialError = std::find_if(r[2].begin(), r[2].end(), [&](const auto &s){ return s.surface && step.surface->geometryId() == s.surface->geometryId(); });
    
    if( inDirect != r[1].end() && inDirect->surface)
        std::cout << inDirect->surface->geometryId() << "\t\t";
    
    if( inDirect != r[2].end() && inTrialError->surface)
        std::cout << inTrialError->surface->geometryId();
    
    std::cout << "\n";
  }
}
