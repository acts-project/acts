// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"

#include <limits>

#include "PropagationDatasets.hpp"
#include "PropagationTests.hpp"

namespace {

namespace ds = ActsTests::PropagationDatasets;
using namespace Acts::UnitLiterals;

using MagneticField = Acts::ConstantBField;
using Stepper = Acts::EigenStepper<
    MagneticField, Acts::StepperExtensionList<Acts::DenseEnvironmentExtension>>;
using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;

constexpr auto epsPos = 1_um;
constexpr auto epsDir = 0.125_mrad;
constexpr auto epsMom = 1_keV;

const Acts::GeometryContext geoCtx;
const Acts::MagneticFieldContext magCtx;

inline Propagator makePropagator(double Bz) {
  using namespace Acts;

  // avoid rebuilding the tracking geometry for every propagator
  static std::shared_ptr<const Acts::TrackingGeometry> detector;
  if (not detector) {
    CuboidVolumeBuilder::VolumeConfig vConf;
    vConf.position = {1.5_m, 0., 0.};
    vConf.length = {3_m, 1_m, 1_m};
    vConf.volumeMaterial = std::make_shared<const HomogeneousVolumeMaterial>(
        Acts::Test::makeBeryllium());
    CuboidVolumeBuilder::Config conf;
    conf.volumeCfg.push_back(vConf);
    conf.position = {1.5_m, 0., 0.};
    conf.length = {3_m, 1_m, 1_m};
    CuboidVolumeBuilder cvb(conf);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto&) {
          return cvb.trackingVolume(context, inner, nullptr);
        });
    detector = TrackingGeometryBuilder(tgbCfg).trackingGeometry(geoCtx);
  }

  MagneticField magField(Acts::Vector3D(0.0, 0.0, Bz));
  Stepper stepper(std::move(magField));
  return Propagator(std::move(stepper), Acts::Navigator(detector));
}

}  // namespace

BOOST_AUTO_TEST_SUITE(PropagationDenseConstant)

BOOST_DATA_TEST_CASE(ForwardBackward,
                     ds::phi* ds::theta* ds::absMomentum* ds::chargeNonZero*
                         ds::pathLength* ds::magneticField,
                     phi, theta, p, q, s, b) {
  const auto propagator = makePropagator(b);
  // phi is ill-defined in forward/backward tracks
  const auto phiFixed = ((0 < theta) and (theta < M_PI)) ? phi : 0.0;
  const auto initial = makeParametersCurvilinear(phiFixed, theta, p, q);
  runForwardBackwardTest<Acts::DenseStepperPropagatorOptions>(
      propagator, geoCtx, magCtx, initial, s, epsPos, epsDir, epsMom);
}

BOOST_AUTO_TEST_SUITE_END()
