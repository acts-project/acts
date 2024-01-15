// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/old/interface.hpp>
#include <boost/test/tools/output_test_stream.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>

namespace bdata = boost::unit_test::data;
using namespace Acts::UnitLiterals;

namespace Acts {
namespace Test {

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;

inline Material makeLiquidArgon() {
  return Material::fromMassDensity(140.034_mm, 857.064_mm, 39.948, 18,
                                   1.396 * 1_g / 1_cm3);
}

inline std::tuple<std::shared_ptr<const Acts::TrackingGeometry>,
                  std::vector<const Acts::Surface*>>
makeDetector() {
  CuboidVolumeBuilder::Config conf;
  conf.position = {0., 0., 0.};
  conf.length = {4_m, 2_m, 2_m};

  {
    CuboidVolumeBuilder::VolumeConfig start;
    start.position = {-1_m, 0, 0};
    start.length = {1_m, 1_m, 1_m};

    conf.volumeCfg.push_back(start);
  }

  {
    CuboidVolumeBuilder::VolumeConfig dense;
    dense.position = {0, 0, 0};
    dense.length = {1_m, 1_m, 1_m};
    dense.volumeMaterial =
        std::make_shared<const HomogeneousVolumeMaterial>(makeLiquidArgon());

    conf.volumeCfg.push_back(dense);
  }

  {
    CuboidVolumeBuilder::SurfaceConfig surface;
    surface.position = {1.5_m, 0, 0};
    surface.rotation =
        Eigen::AngleAxisd(M_PI / 2, Eigen::Vector3d::UnitY()).matrix();
    surface.rBounds = std::make_shared<RectangleBounds>(1_m, 1_m);

    CuboidVolumeBuilder::LayerConfig layer;
    layer.surfaceCfg.push_back(surface);

    CuboidVolumeBuilder::VolumeConfig end;
    end.position = {1_m, 0, 0};
    end.length = {1_m, 1_m, 1_m};
    end.layerCfg.push_back(layer);

    conf.volumeCfg.push_back(end);
  }

  CuboidVolumeBuilder cvb(conf);

  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  auto detector = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);

  std::vector<const Acts::Surface*> surfaces;
  detector->visitSurfaces(
      [&](const Acts::Surface* surface) { surfaces.push_back(surface); });

  return {std::move(detector), std::move(surfaces)};
}

auto makePropagator(std::shared_ptr<const Acts::TrackingGeometry> detector,
                    std::shared_ptr<const Acts::MagneticFieldProvider> bfield) {
  using Stepper = Acts::EigenStepper<
      Acts::StepperExtensionList<Acts::DefaultExtension,
                                 Acts::DenseEnvironmentExtension>,
      Acts::detail::HighestValidAuctioneer>;
  using Navigator = Acts::Navigator;
  using Propagator = Acts::Propagator<Stepper, Navigator>;
  using Actors = ActionList<>;
  using Aborters = AbortList<EndOfWorldReached>;
  using Options = DenseStepperPropagatorOptions<Actors, Aborters>;

  Navigator navigator({detector, true, true, false},
                      getDefaultLogger("nav", Logging::INFO));
  Stepper stepper(bfield);
  Propagator propagator(std::move(stepper), std::move(navigator),
                        getDefaultLogger("prop", Logging::INFO));

  return propagator;
}

BOOST_DATA_TEST_CASE(dense_propagator_test,
                     bdata::make({1_GeV, 10_GeV, 100_GeV}), p) {
  const double q = 1;

  auto bfield =
      std::make_shared<Acts::ConstantBField>(Acts::Vector3{0., 0., 0.});
  auto [detector, surfaces] = makeDetector();

  auto propagator = makePropagator(detector, bfield);

  auto particleHypothesis = ParticleHypothesis::muon();
  double qOverP = particleHypothesis.qOverP(p, q);
  CurvilinearTrackParameters startParams(
      Vector4(-1.5_m, 0, 0, 0), Vector3(1, 0, 0), qOverP,
      BoundVector::Constant(1e-16).asDiagonal(), particleHypothesis);

  DenseStepperPropagatorOptions<> options(geoCtx, magCtx);
  options.maxStepSize = 1_m;
  options.maxSteps = 10000;

  const Acts::Surface& target = *surfaces.back();

  auto result = propagator.propagate(startParams, target, options);

  BOOST_CHECK(result.ok());
  CHECK_CLOSE_REL(3_m, result->pathLength, 1e-6);
  BOOST_CHECK(result->endParameters);

  BoundTrackParameters endParams = result->endParameters.value();

  BOOST_CHECK(endParams.covariance());
  CHECK_CLOSE_ABS(startParams.position(geoCtx) + Vector3(3_m, 0, 0),
                  endParams.position(geoCtx), 1e-6);
  CHECK_CLOSE_ABS(startParams.direction(), endParams.direction(), 1e-6);

  const auto& cov = endParams.covariance().value();

  double endP = endParams.absoluteMomentum();
  double endVarX = cov(eBoundLoc0, eBoundLoc0);
  double endVarY = cov(eBoundLoc1, eBoundLoc1);
  double endVarQOverP = cov(eBoundQOverP, eBoundQOverP);
  double endVarP =
      std::pow(q / std::pow(endParams.qOverP(), 2), 2) * endVarQOverP;

  std::cout << "input p = " << p << std::endl;
  std::cout << "output p = " << endP << std::endl;
  std::cout << "output std x = " << std::sqrt(endVarX) << std::endl;
  std::cout << "output std y = " << std::sqrt(endVarY) << std::endl;
  std::cout << "output std q/p = " << std::sqrt(endVarQOverP) << std::endl;
  std::cout << "output std p = " << std::sqrt(endVarP) << std::endl;

  float theta0 = computeMultipleScatteringTheta0(
      MaterialSlab(makeLiquidArgon(), 1_m), particleHypothesis.absolutePdg(),
      particleHypothesis.mass(), qOverP, particleHypothesis.absoluteCharge());

  std::cout << "theta0 = " << theta0 << std::endl;
}

}  // namespace Test
}  // namespace Acts
