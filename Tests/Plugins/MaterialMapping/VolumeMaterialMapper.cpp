// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// clang-format off
#define BOOST_TEST_MODULE material utils test
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
// clang-format on

#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Utilities/MaterialMapUtils.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include <limits>
#include <random>
#include <vector>
#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/EventData/SingleCurvilinearTrackParameters.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Plugins/MaterialMapping/VolumeMaterialMapper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/Propagator/detail/StandardAborters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"
#include "Acts/Tools/CuboidVolumeBuilder.hpp"
#include "Acts/Tools/TrackingGeometryBuilder.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/detail/Grid.hpp"

namespace bdata = boost::unit_test::data;

using Acts::VectorHelpers::perp;

namespace Acts {

using namespace detail;

namespace Test {

  using RecordedMaterial
      = std::vector<std::pair<Acts::Material, Acts::Vector3D>>;
  using EAxis = Acts::detail::EquidistantAxis;
  using Grid2D
      = Acts::detail::Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis>;
  using Grid3D = Acts::detail::
      Grid<Acts::AccumulatedVolumeMaterial, EAxis, EAxis, EAxis>;
  using MaterialGrid2D = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis>;
  using MaterialGrid3D
      = Acts::detail::Grid<Acts::ActsVectorF<5>, EAxis, EAxis, EAxis>;

  /// @brief This function assigns material to the 2D bin number that represents
  /// the local index of the first axis to the material point.
  ///
  /// @param [in] matPos Position of the material
  /// @param [in] grid Grid that is used for the look-up
  ///
  /// @return Local grid point with the closest distance to @p matPos along the
  /// first axis
  Grid2D::index_t
  mapToBin2D(const Vector3D& matPos, const Grid2D& grid)
  {
    double dist  = std::numeric_limits<double>::max();
    size_t index = 0;
    // Loop through all elements in the first axis
    for (size_t i = 0; i < grid.getNBins()[0]; i++) {
      // Search the closest distance - elements are ordered
      if (std::abs(grid.getUpperRightBinEdge({{i, 0}})[0] - matPos.x())
          < dist) {
        // Store distance and index
        dist  = std::abs(grid.getUpperRightBinEdge({{i, 0}})[0] - matPos.x());
        index = i;
      } else {  // Break if distance becomes larger
        break;
      }
    }
    return {{index, 0}};
  }

  /// @brief This function assigns material to the 3D bin number that represents
  /// the local index of the first axis to the material point.
  ///
  /// @param [in] matPos Position of the material
  /// @param [in] grid Grid that is used for the look-up
  ///
  /// @return Local grid point with the closest distance to @p matPos along the
  /// first axis
  Grid3D::index_t
  mapToBin3D(const Vector3D& matPos, const Grid3D& grid)
  {
    double dist  = std::numeric_limits<double>::max();
    size_t index = 0;
    // Loop through all elements in the first axis
    for (size_t i = 0; i < grid.getNBins()[0]; i++) {
      // Search the closest distance - elements are ordered
      if (std::abs(grid.getUpperRightBinEdge({{i, 0, 0}})[0] - matPos.x())
          < dist) {
        // Store distance and index
        dist = std::abs(grid.getUpperRightBinEdge({{i, 0, 0}})[0] - matPos.x());
        index = i;
      } else {  // Break if distance becomes larger
        break;
      }
    }
    return {{index, 0, 0}};
  }

  /// @brief Rough searcher for closest point
  ///
  /// @param [in] matPos Position of the material
  /// @param [in] grid Grid that is used for the look-up
  ///
  /// @return Local grid point with the closest distance to @p matPos
  Grid3D::index_t
  mapMaterial3D(const Vector3D& matPos, const Grid3D& grid)
  {
    double dist   = std::numeric_limits<double>::max();
    size_t indexX = 0, indexY = 0, indexZ = 0;
    // Loop through all elements
    for (size_t i = 0; i < grid.getNBins()[0]; i++) {
      for (size_t j = 0; j < grid.getNBins()[1]; j++) {
        for (size_t k = 0; k < grid.getNBins()[2]; k++) {
          // Search the closest distance - elements are ordered
          double dX = grid.getUpperRightBinEdge({{i, j, k}})[0] - matPos.x();
          double dY = grid.getUpperRightBinEdge({{i, j, k}})[1] - matPos.y();
          double dZ = grid.getUpperRightBinEdge({{i, j, k}})[2] - matPos.z();

          if (std::sqrt(dX * dX + dY * dY + dZ * dZ) < dist) {
            // Store distance and index
            dist   = std::sqrt(dX * dX + dY * dY + dZ * dZ);
            indexX = i;
            indexY = j;
            indexZ = k;
          } else {  // Break if distance becomes larger
            break;
          }
        }
      }
    }
    return {{indexX, indexY, indexZ}};
  }

  /// @brief Collector of material and position along propagation
  struct MaterialCollector
  {
    struct this_result
    {
      std::vector<Material> matTrue;
      std::vector<Vector3D> position;
    };
    using result_type = this_result;

    template <typename propagator_state_t, typename stepper_t>
    void
    operator()(propagator_state_t& state,
               const stepper_t&    stepper,
               result_type&        result) const
    {
      if (state.navigation.currentVolume != nullptr) {
        result.matTrue.push_back(
            (state.navigation.currentVolume->material() != nullptr)
                ? *state.navigation.currentVolume->material()
                : Material());

        result.position.push_back(stepper.position(state.stepping));
      }
    }
  };

  /// @brief Various test cases of the VolumeMaterialMapper functions
  BOOST_AUTO_TEST_CASE(VolumeMaterialMapper_tests)
  {

    // Define some axes and grid points
    std::array<double, 3> axis1 = {0., 1., 2.};
    std::array<double, 3> axis2 = {2., 4., 3.};
    std::array<double, 3> axis3 = {5., 6., 3.};

    // Make some materials
    std::vector<std::pair<Material, Vector3D>> matRecord;
    Material mat1(1., 2., 3., 4., 5.);
    Material mat2(6., 7., 8., 9., 10.);

    Material       vacuum;
    ActsVectorF<5> matMix;
    matMix << 3.5, 4.5, 5.5, 6.5, 7.5;

    //
    // Test the production chain in 2D
    //
    matRecord.clear();
    matRecord.push_back(std::make_pair(mat1, Vector3D(0., 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.4, 0., 0.)));
    matRecord.push_back(std::make_pair(mat2, Vector3D(0.6, 0., 0.)));

    MaterialGrid2D mgrid2d
        = createMaterialGrid(axis1, axis2, matRecord, mapToBin2D);

    // Test sizes
    BOOST_CHECK_EQUAL(mgrid2d.size(), (axis1[2] + 2) * (axis2[2] + 2));
    for (size_t index = 0; index < mgrid2d.size(); index++) {
      // Check the contained data
      if (index == 0) {
        BOOST_CHECK_EQUAL(mgrid2d.at(index), matMix);
        continue;
      }
      if (index == 5) {
        BOOST_CHECK_EQUAL(mgrid2d.at(index), mat2.classificationNumbers());
        continue;
      }
      BOOST_CHECK_EQUAL(mgrid2d.at(index), vacuum.classificationNumbers());
    }

    //
    // Test the production chain in 3D
    //
    MaterialGrid3D mgrid3d
        = createMaterialGrid(axis1, axis2, axis3, matRecord, mapToBin3D);

    // Test sizes
    BOOST_CHECK_EQUAL(mgrid3d.size(),
                      (axis1[2] + 2) * (axis2[2] + 2) * (axis3[2] + 2));
    for (size_t index = 0; index < mgrid3d.size(); index++) {
      // Check the contained data
      if (index == 0) {
        BOOST_CHECK_EQUAL(mgrid3d.at(index), matMix);
        continue;
      }
      if (index == 25) {
        BOOST_CHECK_EQUAL(mgrid3d.at(index), mat2.classificationNumbers());
        continue;
      }
      BOOST_CHECK_EQUAL(mgrid3d.at(index), vacuum.classificationNumbers());
    }
  }

  /// @brief Test case for comparison between the mapped material and the
  /// associated material by propagation
  BOOST_AUTO_TEST_CASE(VolumeMaterialMapper_comparison_tests)
  {
    // Build a vacuum volume
    CuboidVolumeBuilder::VolumeConfig vCfg1;
    vCfg1.position = Vector3D(0.5 * units::_m, 0., 0.);
    vCfg1.length   = Vector3D(1. * units::_m, 1. * units::_m, 1. * units::_m);
    vCfg1.name     = "Vacuum volume";
    vCfg1.material
        = std::make_shared<const Material>(352.8, 407., 9.012, 4., 1.848e-3);

    // Build a material volume
    CuboidVolumeBuilder::VolumeConfig vCfg2;
    vCfg2.position = Vector3D(1.5 * units::_m, 0., 0.);
    vCfg2.length   = Vector3D(1. * units::_m, 1. * units::_m, 1. * units::_m);
    vCfg2.name     = "First material volume";
    vCfg2.material
        = std::make_shared<const Material>(95.7, 465.2, 28.03, 14., 2.32e-3);

    // Build another material volume with different material
    CuboidVolumeBuilder::VolumeConfig vCfg3;
    vCfg3.position = Vector3D(2.5 * units::_m, 0., 0.);
    vCfg3.length   = Vector3D(1. * units::_m, 1. * units::_m, 1. * units::_m);
    vCfg3.name     = "Second material volume";
    vCfg3.material
        = std::make_shared<const Material>(352.8, 407., 9.012, 4., 1.848e-3);

    // Configure world
    CuboidVolumeBuilder::Config cfg;
    cfg.position  = Vector3D(1.5 * units::_m, 0., 0.);
    cfg.length    = Vector3D(3. * units::_m, 1. * units::_m, 1. * units::_m);
    cfg.volumeCfg = {vCfg1, vCfg2, vCfg3};

    GeometryContext gc;

    // Build a detector
    CuboidVolumeBuilder             cvb(cfg);
    TrackingGeometryBuilder::Config tgbCfg;
    tgbCfg.trackingVolumeBuilders.push_back(
        [=](const auto& context, const auto& inner, const auto&) {
          return cvb.trackingVolume(context, inner, nullptr);
        });
    TrackingGeometryBuilder                 tgb(tgbCfg);
    std::unique_ptr<const TrackingGeometry> detector = tgb.trackingGeometry(gc);

    // Set up the grid axes
    std::array<double, 3> xAxis{0. * units::_m, 3. * units::_m, 7};
    std::array<double, 3> yAxis{-0.5 * units::_m, 0.5 * units::_m, 7};
    std::array<double, 3> zAxis{-0.5 * units::_m, 0.5 * units::_m, 7};

    // Set up a random engine for sampling material
    std::random_device               rd;
    std::mt19937                     gen(rd());
    std::uniform_real_distribution<> disX(0., 3. * units::_m);
    std::uniform_real_distribution<> disYZ(-0.5 * units::_m, 0.5 * units::_m);

    // Sample the Material in the detector
    RecordedMaterial matRecord;
    for (unsigned int i = 0; i < 1e4; i++) {
      Vector3D pos(disX(gen), disYZ(gen), disYZ(gen));
      Material tv
          = (detector->lowestTrackingVolume(gc, pos)->material() != nullptr)
          ? *(detector->lowestTrackingVolume(gc, pos)->material())
          : Material();
      matRecord.push_back(std::make_pair(tv, pos));
    }

    // Build the material grid
    MaterialGrid3D grid
        = createMaterialGrid(xAxis, yAxis, zAxis, matRecord, mapMaterial3D);

    // Construct a simple propagation through the detector
    StraightLineStepper sls;
    Navigator           nav(std::move(detector));
    Propagator<StraightLineStepper, Navigator> prop(sls, nav);

    // Set some start parameters
    Vector3D pos(0., 0., 0.);
    Vector3D mom(1. * units::_GeV, 0., 0.);
    SingleCurvilinearTrackParameters<NeutralPolicy> sctp(nullptr, pos, mom);

    MagneticFieldContext mc;

    // Launch propagation and gather result
    PropagatorOptions<ActionList<MaterialCollector>,
                      AbortList<detail::EndOfWorldReached>>
        po(gc, mc);
    po.maxStepSize     = 1. * units::_mm;
    po.maxSteps        = 1e6;
    const auto& result = prop.propagate(sctp, po).value();
    const MaterialCollector::this_result& stepResult
        = result.get<typename MaterialCollector::result_type>();

    // Collect the material as given by the grid and test it
    std::vector<Material> matGrid;
    double                gridX0 = 0., gridL0 = 0., trueX0 = 0., trueL0 = 0.;
    for (unsigned int i = 0; i < stepResult.position.size(); i++) {
      matGrid.push_back(grid.at(stepResult.position[i]));
      gridX0 += matGrid[i].X0();
      gridL0 += matGrid[i].L0();
      trueX0 += stepResult.matTrue[i].X0();
      trueL0 += stepResult.matTrue[i].L0();
    }
    CHECK_CLOSE_REL(gridX0, trueX0, 1e-1);
    CHECK_CLOSE_REL(gridL0, trueL0, 1e-1);
  }
}  // namespace Test
}  // namespace Acts