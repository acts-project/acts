// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/NeutralTrackParameters.hpp"
#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/CylinderVolumeBuilder.hpp"
#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialGridHelper.hpp"
#include "Acts/Material/MaterialMapUtils.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Material/ProtoVolumeMaterial.hpp"
#include "Acts/Material/VolumeMaterialMapper.hpp"
#include "Acts/Propagator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tests/CommonHelpers/PredefinedMaterials.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

#include <limits>
#include <random>
#include <vector>

namespace Acts {

/// @brief Collector of material and position along propagation
struct MaterialCollector {
  struct this_result {
    std::vector<Material> matTrue;
    std::vector<Vector3> position;
  };
  using result_type = this_result;

  template <typename propagator_state_t, typename stepper_t,
            typename navigator_t>
  void operator()(propagator_state_t& state, const stepper_t& stepper,
                  const navigator_t& navigator, result_type& result,
                  const Logger& /*logger*/) const {
    if (navigator.currentVolume(state.navigation) != nullptr) {
      auto position = stepper.position(state.stepping);
      result.matTrue.push_back(
          (navigator.currentVolume(state.navigation)->volumeMaterial() !=
           nullptr)
              ? navigator.currentVolume(state.navigation)
                    ->volumeMaterial()
                    ->material(position)
              : Material());

      result.position.push_back(position);
    }
  }
};

namespace Test {

/// Test the filling and conversion
BOOST_AUTO_TEST_CASE(SurfaceMaterialMapper_tests) {
  using namespace Acts::UnitLiterals;

  BinUtility bu1(4, 0_m, 1_m, open, binX);
  bu1 += BinUtility(2, -0.5_m, 0.5_m, open, binY);
  bu1 += BinUtility(2, -0.5_m, 0.5_m, open, binZ);

  BinUtility bu2(4, 1_m, 2_m, open, binX);
  bu2 += BinUtility(2, -0.5_m, 0.5_m, open, binY);
  bu2 += BinUtility(2, -0.5_m, 0.5_m, open, binZ);

  BinUtility bu3(4, 2_m, 3_m, open, binX);
  bu3 += BinUtility(2, -0.5_m, 0.5_m, open, binY);
  bu3 += BinUtility(2, -0.5_m, 0.5_m, open, binZ);

  // Build a vacuum volume
  CuboidVolumeBuilder::VolumeConfig vCfg1;
  vCfg1.position = Vector3(0.5_m, 0., 0.);
  vCfg1.length = Vector3(1_m, 1_m, 1_m);
  vCfg1.name = "Vacuum volume";
  vCfg1.volumeMaterial = std::make_shared<const ProtoVolumeMaterial>(bu1);

  // Build a material volume
  CuboidVolumeBuilder::VolumeConfig vCfg2;
  vCfg2.position = Vector3(1.5_m, 0., 0.);
  vCfg2.length = Vector3(1_m, 1_m, 1_m);
  vCfg2.name = "First material volume";
  vCfg2.volumeMaterial = std::make_shared<const ProtoVolumeMaterial>(bu2);

  // Build another material volume with different material
  CuboidVolumeBuilder::VolumeConfig vCfg3;
  vCfg3.position = Vector3(2.5_m, 0., 0.);
  vCfg3.length = Vector3(1_m, 1_m, 1_m);
  vCfg3.name = "Second material volume";
  vCfg3.volumeMaterial = std::make_shared<const ProtoVolumeMaterial>(bu3);

  // Configure world
  CuboidVolumeBuilder::Config cfg;
  cfg.position = Vector3(1.5_m, 0., 0.);
  cfg.length = Vector3(3_m, 1_m, 1_m);
  cfg.volumeCfg = {vCfg1, vCfg2, vCfg3};

  GeometryContext gc;

  // Build a detector
  CuboidVolumeBuilder cvb(cfg);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::shared_ptr<const TrackingGeometry> tGeometry = tgb.trackingGeometry(gc);

  /// We need a Navigator, Stepper to build a Propagator
  Navigator navigator({tGeometry});
  StraightLineStepper stepper;
  VolumeMaterialMapper::StraightLinePropagator propagator(stepper,
                                                          std::move(navigator));

  /// The config object
  Acts::VolumeMaterialMapper::Config vmmConfig;
  Acts::VolumeMaterialMapper vmMapper(
      vmmConfig, std::move(propagator),
      getDefaultLogger("VolumeMaterialMapper", Logging::VERBOSE));

  /// Create some contexts
  GeometryContext gCtx;
  MagneticFieldContext mfCtx;

  /// Now create the mapper state
  auto mState = vmMapper.createState(gCtx, mfCtx, *tGeometry);

  /// Test if this is not null
  BOOST_CHECK_EQUAL(mState.materialBin.size(), 3u);
}

/// @brief Test case for comparison between the mapped material and the
/// associated material by propagation
BOOST_AUTO_TEST_CASE(VolumeMaterialMapper_comparison_tests) {
  using namespace Acts::UnitLiterals;

  // Build a vacuum volume
  CuboidVolumeBuilder::VolumeConfig vCfg1;
  vCfg1.position = Vector3(0.5_m, 0., 0.);
  vCfg1.length = Vector3(1_m, 1_m, 1_m);
  vCfg1.name = "Vacuum volume";
  vCfg1.volumeMaterial =
      std::make_shared<const HomogeneousVolumeMaterial>(Material());

  // Build a material volume
  CuboidVolumeBuilder::VolumeConfig vCfg2;
  vCfg2.position = Vector3(1.5_m, 0., 0.);
  vCfg2.length = Vector3(1_m, 1_m, 1_m);
  vCfg2.name = "First material volume";
  vCfg2.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(makeSilicon());

  // Build another material volume with different material
  CuboidVolumeBuilder::VolumeConfig vCfg3;
  vCfg3.position = Vector3(2.5_m, 0., 0.);
  vCfg3.length = Vector3(1_m, 1_m, 1_m);
  vCfg3.name = "Second material volume";
  vCfg3.volumeMaterial =
      std::make_shared<const HomogeneousVolumeMaterial>(Material());

  // Configure world
  CuboidVolumeBuilder::Config cfg;
  cfg.position = Vector3(1.5_m, 0., 0.);
  cfg.length = Vector3(3_m, 1_m, 1_m);
  cfg.volumeCfg = {vCfg1, vCfg2, vCfg3};

  GeometryContext gc;

  // Build a detector
  CuboidVolumeBuilder cvb(cfg);
  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& context, const auto& inner, const auto&) {
        return cvb.trackingVolume(context, inner, nullptr);
      });
  TrackingGeometryBuilder tgb(tgbCfg);
  std::unique_ptr<const TrackingGeometry> detector = tgb.trackingGeometry(gc);

  // Set up the grid axes
  Acts::MaterialGridAxisData xAxis{0_m, 3_m, 7};
  Acts::MaterialGridAxisData yAxis{-0.5_m, 0.5_m, 7};
  Acts::MaterialGridAxisData zAxis{-0.5_m, 0.5_m, 7};

  // Set up a random engine for sampling material
  std::random_device rd;
  std::mt19937 gen(42);
  std::uniform_real_distribution<> disX(0., 3_m);
  std::uniform_real_distribution<> disYZ(-0.5_m, 0.5_m);

  // Sample the Material in the detector
  RecordedMaterialVolumePoint matRecord;
  for (unsigned int i = 0; i < 1e4; i++) {
    Vector3 pos(disX(gen), disYZ(gen), disYZ(gen));
    std::vector<Vector3> volPos;
    volPos.push_back(pos);
    Material tv =
        (detector->lowestTrackingVolume(gc, pos)->volumeMaterial() != nullptr)
            ? (detector->lowestTrackingVolume(gc, pos)->volumeMaterial())
                  ->material(pos)
            : Material();
    MaterialSlab matProp(tv, 1);
    matRecord.push_back(std::make_pair(matProp, volPos));
  }

  // Build the material grid
  Grid3D Grid = createGrid(xAxis, yAxis, zAxis);
  std::function<Vector3(Vector3)> transfoGlobalToLocal =
      [](Vector3 pos) -> Vector3 {
    return {pos.x(), pos.y(), pos.z()};
  };

  // Walk over each properties
  for (const auto& rm : matRecord) {
    // Walk over each point associated with the properties
    for (const auto& point : rm.second) {
      // Search for fitting grid point and accumulate
      Acts::Grid3D::index_t index =
          Grid.localBinsFromLowerLeftEdge(transfoGlobalToLocal(point));
      Grid.atLocalBins(index).accumulate(rm.first);
    }
  }

  MaterialGrid3D matGrid = mapMaterialPoints(Grid);

  // Construct a simple propagation through the detector
  StraightLineStepper sls;
  Navigator::Config navCfg;
  navCfg.trackingGeometry = std::move(detector);
  Navigator nav(navCfg);
  Propagator<StraightLineStepper, Navigator> prop(sls, nav);

  // Set some start parameters
  Vector4 pos4(0., 0., 0., 42_ns);
  Vector3 dir(1., 0., 0.);
  NeutralCurvilinearTrackParameters sctp(pos4, dir, 1 / 1_GeV);

  MagneticFieldContext mc;
  // Launch propagation and gather result
  PropagatorOptions<ActionList<MaterialCollector>, AbortList<EndOfWorldReached>>
      po(gc, mc);
  po.maxStepSize = 1._mm;
  po.maxSteps = 1e6;

  const auto& result = prop.propagate(sctp, po).value();
  const MaterialCollector::this_result& stepResult =
      result.get<typename MaterialCollector::result_type>();

  // Collect the material as given by the grid and test it
  std::vector<Material> matvector;
  double gridX0 = 0., gridL0 = 0., trueX0 = 0., trueL0 = 0.;
  for (unsigned int i = 0; i < stepResult.position.size(); i++) {
    matvector.push_back(matGrid.atPosition(stepResult.position[i]));
    gridX0 += 1 / matvector[i].X0();
    gridL0 += 1 / matvector[i].L0();
    trueX0 += 1 / stepResult.matTrue[i].X0();
    trueL0 += 1 / stepResult.matTrue[i].L0();
  }
  CHECK_CLOSE_REL(gridX0, trueX0, 1e-1);
  CHECK_CLOSE_REL(gridL0, trueL0, 1e-1);
}

}  // namespace Test
}  // namespace Acts
