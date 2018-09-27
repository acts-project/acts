// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///  Boost include(s)
#define BOOST_TEST_MODULE StepActor Tests

#include <boost/test/included/unit_test.hpp>
// leave blank line

#include "Acts/Detector/TrackingGeometry.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/CuboidVolumeBounds.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Detector/TrackingVolume.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Extrapolator/StepActor.hpp"
#include "Acts/Extrapolator/Navigator.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Tools/LayerArrayCreator.hpp"
#include "Acts/Layers/PlaneLayer.hpp"

namespace Acts {
namespace Test {

  ///
  /// @brief Aborter for the case that a particle leaves the detector
  ///
  struct EndOfWorld
  {
    /// @brief Constructor
    EndOfWorld() = default;

    /// @brief Main call operator for the abort operation
    ///
    /// @tparam propagator_state_t State of the propagator
    /// @param [in] state State of the propagation
    /// @return Boolean statement if the particle is still in the detector
    template <typename propagator_state_t>
    bool
    operator()(propagator_state_t& state) const
    {
      if (std::abs(state.stepping.position().x()) > 0.5 * units::_m
          || std::abs(state.stepping.position().y()) > 0.5 * units::_m
          || std::abs(state.stepping.position().z()) > 2. * units::_m)
        return true;
      return false;
    }
  };

struct StepCollector
{

  struct this_result
  {
    std::vector<Vector3D> position;
    std::vector<Vector3D> momentum;
    std::vector<ActsSymMatrixD<5>> cov;
  };

  using result_type = this_result;

  template <typename propagator_state_t>
  void
  operator()(propagator_state_t& state, result_type& result) const
  {
	  std::cout << state.options.debugString << std::endl;
	  result.position.push_back(state.stepping.position());
	  result.momentum.push_back(state.stepping.momentum());
	  result.cov.push_back(state.stepping.cov);
  }
};

std::shared_ptr<TrackingGeometry>
buildVacDetector()
{
	// Build vacuum block
	Transform3D trafoVac(Transform3D::Identity());
    trafoVac.translation() = Vector3D(0., 0., 1. * units::_m);
    
    std::shared_ptr<const PlanarBounds> rBounds(new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    LayerPtr dummyLayer = PlaneLayer::create(std::make_shared<const Transform3D>(trafoVac),
											 rBounds);
         
	LayerArrayCreator layArrCreator(getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    std::unique_ptr<const LayerArray> layArr(
        layArrCreator.layerArray({dummyLayer},
                                 0.,
                                 2. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binZ));
    
    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(0.5 * units::_m, 0.5 * units::_m, 1. * units::_m);
    
    //~ auto trackingVac = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVac), boundsVol, nullptr, "Vacuum");
    auto trackingVac = TrackingVolume::create(std::make_shared<const Transform3D>(trafoVac), boundsVol, nullptr, std::move(layArr), {}, {}, {}, "Vacuum");
    

    return std::shared_ptr<TrackingGeometry>(new TrackingGeometry(trackingVac));
}

std::shared_ptr<TrackingGeometry>
buildMatDetector()
{                                 
	// Build material block
    Transform3D trafoMat(Transform3D::Identity());
    trafoMat.translation() = Vector3D(0., 0., 1. * units::_m);
    
    std::shared_ptr<const PlanarBounds> rBounds(new RectangleBounds(0.5 * units::_m, 0.5 * units::_m));
    LayerPtr dummyLayer = PlaneLayer::create(std::make_shared<const Transform3D>(trafoMat),
											 rBounds);
											 
    LayerArrayCreator layArrCreator(getDefaultLogger("LayerArrayCreator", Logging::VERBOSE));
    std::unique_ptr<const LayerArray> layArr(
        layArrCreator.layerArray({dummyLayer},
                                 0.,
                                 2. * units::_m,
                                 BinningType::arbitrary,
                                 BinningValue::binZ));
    
    auto boundsVol = std::make_shared<const CuboidVolumeBounds>(0.5 * units::_m, 0.5 * units::_m, 1. * units::_m);
     
    std::shared_ptr<const Material> mat(new Material(352.8, 407., 9.012, 4., 1.848e-3));
    
    auto trackingMat = TrackingVolume::create(std::make_shared<const Transform3D>(trafoMat), boundsVol, mat, std::move(layArr), {}, {}, {}, "Material");

	return std::shared_ptr<TrackingGeometry>(new TrackingGeometry(trackingMat));
}

BOOST_AUTO_TEST_CASE(step_actor_test)
{
	{
	std::shared_ptr<TrackingGeometry> vacuum = buildVacDetector();

	// Build navigators
    Navigator naviVac(vacuum);
    naviVac.resolvePassive   = true;
    naviVac.resolveMaterial  = true;
    naviVac.resolveSensitive = true;

    // Set initial parameters for the particle track
    ActsSymMatrixD<5> cov;
    cov << 1. * units::_mm, 0., 0., 0., 0., 
		   0., 1. * units::_mm, 0., 0., 0.,
		   0., 0., 1., 0., 0., 
		   0., 0., 0., 1., 0., 
		   0., 0., 0., 0., 1.;
    auto     covPtr = std::make_unique<const ActsSymMatrixD<5>>(cov);
    Vector3D startParams(0., 0., 0.), startMom(0., 0., 1. * units::_GeV);
    SingleCurvilinearTrackParameters<ChargedPolicy> sbtp(
        std::move(covPtr), startParams, startMom, 1.);

    // Create action list for surface collection
    ActionList<StepCollector, StepActor> aList;
	AbortList<EndOfWorld> abortList;
	
    // Set options for propagator
	Propagator<EigenStepper<ConstantBField>, Navigator>::
        Options<ActionList<StepCollector, StepActor>, AbortList<EndOfWorld>> propOpts;
    propOpts.actionList     = aList;
    propOpts.stopConditions = abortList;
    propOpts.maxSteps       = 1e6;
    propOpts.maxStepSize = 2. * units::_m;
    propOpts.debug = true;

    // Re-configure propagation with B-field
    ConstantBField               bField(Vector3D(0., 0., 0.));
    EigenStepper<ConstantBField> es(bField);
    Propagator<EigenStepper<ConstantBField>, Navigator> prop(es, naviVac);
    
    const auto& result = prop.propagate(sbtp, propOpts);
    const StepCollector::this_result& stepResult
        = result.get<typename StepCollector::result_type>();

    for(const auto& pos : stepResult.position)
    {
		BOOST_TEST(pos.x() == 0.);
		BOOST_TEST(pos.y() == 0.);
    }
    for(const auto& mom : stepResult.momentum)
    {
		BOOST_TEST(mom.x() == 0.);
		BOOST_TEST(mom.y() == 0.);
		BOOST_TEST(mom.z() == 1. * units::_GeV);
	}
	for(const auto& c : stepResult.cov)
	{
		BOOST_TEST(c == ActsSymMatrixD<5>::Identity());
	}
}
	//~ {
	//~ std::shared_ptr<TrackingGeometry> material = buildMatDetector();
	
    //~ Navigator naviMat(material);
    //~ naviMat.resolvePassive   = true;
    //~ naviMat.resolveMaterial  = true;
    //~ naviMat.resolveSensitive = true;
    
    
    //~ ConstantBField               bField(Vector3D(0., 0.5 * units::_T, 0.));
    //~ EigenStepper<ConstantBField> es(bField);
    //~ Propagator<EigenStepper<ConstantBField>, Navigator> propB(es, naviVac);
	//~ }
}
 
}  // namespace Test
}  // namespace Acts
