// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Propagator/StraightLineStepper.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/NeutralParameters.hpp"
#include "Acts/Tests/CommonHelpers/FloatComparisons.hpp"

namespace tt = boost::test_tools;

namespace Acts {
namespace Test {

using Covariance = BoundSymMatrix;
using Jacobian = BoundMatrix;

struct PropState
{
	PropState(StraightLineStepper::State sState) : stepping(sState) {}
	
	StraightLineStepper::State stepping;
	
	struct {
		double mass = 42.;
	} options;
};

/// These tests are aiming to test whether the state setup is working properly
BOOST_AUTO_TEST_CASE(straight_line_stepper_state_test) {
	// Set up some variables
  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  
  Vector3D pos(1., 2., 3.);
  Vector3D mom(4., 5., 6.);
  double time = 7.;
  double charge = -1.;
  
  // Test charged parameters without covariance matrix
  CurvilinearParameters cp(std::nullopt, pos, mom, charge, time);
	StraightLineStepper::State slsState(tgContext, mfContext, cp, ndir, stepSize, tolerance);
	
	// Test the result & compare with the input/test for reasonable members
	BOOST_TEST(slsState.jacToGlobal == BoundToFreeMatrix::Zero());
	BOOST_TEST(slsState.jacTransport == FreeMatrix::Identity());
	BOOST_TEST(slsState.derivative == FreeVector::Zero());
	BOOST_TEST(!slsState.covTransport);
	BOOST_TEST(slsState.cov == Covariance::Zero());
	BOOST_TEST(slsState.pos == pos);
	BOOST_TEST(slsState.dir == mom.normalized());
	BOOST_TEST(slsState.p == mom.norm());
	BOOST_TEST(slsState.q == charge);
	BOOST_TEST(slsState.t == time);
	BOOST_TEST(slsState.navDir == ndir);
	BOOST_TEST(slsState.pathAccumulated == 0.);
	BOOST_TEST(slsState.stepSize == ndir * stepSize);
	BOOST_TEST(slsState.previousStepSize == 0.);
	BOOST_TEST(slsState.tolerance == tolerance);
  
    // Test without charge and covariance matrix
    NeutralCurvilinearParameters ncp(std::nullopt, pos, mom, time);
    slsState = StraightLineStepper::State(tgContext, mfContext, ncp, ndir, stepSize, tolerance);
    BOOST_TEST(slsState.q == 0.);
    
    // Test with covariance matrix
    Covariance cov = 8. * Covariance::Identity();
    ncp = NeutralCurvilinearParameters(cov, pos, mom, time);
    slsState = StraightLineStepper::State(tgContext, mfContext, ncp, ndir, stepSize, tolerance);
    BOOST_TEST(slsState.jacToGlobal != BoundToFreeMatrix::Zero());
	BOOST_TEST(slsState.covTransport);
    BOOST_TEST(slsState.cov == cov);
}

/// These tests are aiming to test the functions of the StraightLineStepper
/// The numerical correctness of the stepper is tested in the integration tests
BOOST_AUTO_TEST_CASE(straight_line_stepper_test) {
	
	  GeometryContext tgContext = GeometryContext();
  MagneticFieldContext mfContext = MagneticFieldContext();
  NavigationDirection ndir = backward;
  double stepSize = 123.;
  double tolerance = 234.;
  
  Vector3D pos(1., 2., 3.);
  Vector3D mom(4., 5., 6.);
  double time = 7.;
  double charge = -1.;
    Covariance cov = 8. * Covariance::Identity();
  
  CurvilinearParameters cp(cov, pos, mom, charge, time);
	StraightLineStepper::State slsState(tgContext, mfContext, cp, ndir, stepSize, tolerance);

	StraightLineStepper sls;

	BOOST_TEST(sls.position(slsState) == slsState.pos);
	BOOST_TEST(sls.direction(slsState) == slsState.dir);
	BOOST_TEST(sls.momentum(slsState) == slsState.p);
	BOOST_TEST(sls.charge(slsState) == slsState.q);
	BOOST_TEST(sls.time(slsState) == slsState.t);

    BOOST_TEST(sls.overstepLimit(slsState) == tolerance);
    
  //~ double overstepLimit(const State& /*state*/) const {
    //~ return s_onSurfaceTolerance;
  //~ }
  

	const std::string originalStepSize = slsState.stepSize.toString();
	
	sls.setStepSize(slsState, 1337.);
	BOOST_TEST(slsState.previousStepSize == ndir * stepSize);
	BOOST_TEST(slsState.stepSize == 1337.);
	
	sls.releaseStepSize(slsState);
	BOOST_TEST(slsState.stepSize == -123.);
	BOOST_TEST(sls.outputStepSize(slsState) == originalStepSize);

                        
	auto curvState = sls.curvilinearState(slsState, true);
	auto curvPars = std::get<0>(curvState);
	CHECK_CLOSE_ABS(curvPars.position(), cp.position(), 1e-6);
	CHECK_CLOSE_ABS(curvPars.momentum(), cp.momentum(), 1e-6);
	CHECK_CLOSE_ABS(curvPars.charge(), cp.charge(), 1e-6);
	CHECK_CLOSE_ABS(curvPars.time(), cp.time(), 1e-6);
	BOOST_TEST(curvPars.covariance().has_value());
	BOOST_TEST(*curvPars.covariance() != cov);
	CHECK_CLOSE_COVARIANCE(std::get<1>(curvState), BoundMatrix(BoundMatrix::Identity()), 1e-6);
	CHECK_CLOSE_ABS(std::get<2>(curvState), 0., 1e-6);
	



	Vector3D newPos(2., 4., 8.);
	Vector3D newMom(3., 9., 27.);
	double newTime(321.);
	sls.update(slsState, newPos, newMom.normalized(), newMom.norm(), newTime);
	BOOST_TEST(slsState.pos == newPos);
	BOOST_TEST(slsState.dir == newMom.normalized());
	BOOST_TEST(slsState.p == newMom.norm());
	BOOST_TEST(slsState.q == charge);
	BOOST_TEST(slsState.t == newTime);
	
	slsState.cov = cov;
    sls.covarianceTransport(slsState, true);
    BOOST_TEST(slsState.cov != cov);
	BOOST_TEST(slsState.jacToGlobal != BoundToFreeMatrix::Zero());
	BOOST_TEST(slsState.jacTransport == FreeMatrix::Identity());
	BOOST_TEST(slsState.derivative == FreeVector::Zero());

	

	slsState.cov = cov;
	PropState ps(slsState);

	ps.stepping.covTransport = false;	
	double h = sls.step(ps).value();
	BOOST_TEST(ps.stepping.stepSize == ndir * stepSize);
	BOOST_TEST(ps.stepping.stepSize == h);
	CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
	BOOST_TEST(ps.stepping.pos.norm() > newPos.norm());
	BOOST_TEST(ps.stepping.dir == newMom.normalized());
	BOOST_TEST(ps.stepping.q == charge);
	BOOST_TEST(ps.stepping.t < newTime);
	BOOST_TEST(ps.stepping.derivative == FreeVector::Zero());
	BOOST_TEST(ps.stepping.jacTransport == FreeMatrix::Identity());
		
	ps.stepping.covTransport = true;
	double h2 = sls.step(ps).value();
	BOOST_TEST(ps.stepping.stepSize == ndir * stepSize);
	BOOST_TEST(h2 == h);
	CHECK_CLOSE_COVARIANCE(ps.stepping.cov, cov, 1e-6);
	BOOST_TEST(ps.stepping.pos.norm() > newPos.norm());
	BOOST_TEST(ps.stepping.dir == newMom.normalized());
	BOOST_TEST(ps.stepping.q == charge);
	BOOST_TEST(ps.stepping.t < newTime);
	BOOST_TEST(ps.stepping.derivative != FreeVector::Zero());
	BOOST_TEST(ps.stepping.jacTransport != FreeMatrix::Identity());
	
	
	
  //~ Intersection::Status updateSurfaceStatus(State& state, const Surface& surface,
                                           //~ const BoundaryCheck& bcheck) const {
    //~ return detail::updateSingleSurfaceStatus<StraightLineStepper>(
        //~ *this, state, surface, bcheck);
  //~ }
  //~ template <typename object_intersection_t>
  //~ void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      //~ bool release = true) const {
    //~ detail::updateSingleStepSize<StraightLineStepper>(state, oIntersection,
                                                      //~ release);
  //~ }
  
  
  
  //~ BoundState boundState(State& state, const Surface& surface,
                        //~ bool reinitialize) const;
                        
                        
                          //~ void update(State& state, const BoundParameters& pars) const;
                          
                          
                          
  //~ void covarianceTransport(State& state, const Surface& surface,
                           //~ bool reinitialize = false) const;


}
}  // namespace Test
}  // namespace Acts