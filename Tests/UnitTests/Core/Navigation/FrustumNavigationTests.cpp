// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "ActsExamples/MuonSpectrometerMockupDetector/GeoMuonMockupExperiment.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelMuonMockupBuilder.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"

using namespace Acts;
using namespace ActsExamples;
using namespace Acts::GeoModel;

GeometryContext gContext;

BOOST_AUTO_TEST_SUITE(Experimental)
auto logger = getDefaultLogger("FrustumNavigationTests", Logging::VERBOSE);

// This tests the frustum navigation policy for gen3 geometry interface
BOOST_AUTO_TEST_CASE(Frustum_NavigationPolicy) {
	//create mockup detector
	GeoMuonMockupExperiment::Config detCfg;
	GeoMuonMockupExperiment mockDet(detCfg,"GeoMockUpMS");
	Acts::GeoModelTree gmTree=mockDet.constructMS();
	BOOST_CHECK(gmTree.publisher->getPublishedVol("Muon").size()>0);

	GeoModelDetector::Config gmdConfig;
	gmdConfig.geoModelTree=gmTree;
	GeoModelDetector muonDetector(gmdConfig);

	GeoModelDetectorObjectFactory::Config gmdObjFactoryCfg;
	gmdObjFactoryCfg.nameList={"RpcGasGap","MDTDriftGas"};
	gmdObjFactoryCfg.convertSubVolumes=true;
	gmdObjFactoryCfg.convertBox="MDT";
	GeoModelDetectorObjectFactory gmdObjFactory(gmdObjFactoryCfg);
	GeoModelDetectorObjectFactory::Options options;
	options.queries={"Muon"};
	GeoModelDetectorObjectFactory::Cache cache;
	gmdObjFactory.construct(cache,gContext,gmTree,options);

	GeoModelMuonMockupBuilder::Config builderConfig;
	builderConfig.stationNames = {"Inner", "Middle", "Outer"};
	builderConfig.volumeBoxFPVs = cache.boundingBoxes;
	GeoModelMuonMockupBuilder trackingGeometryBuilder(builderConfig,"GeoModelMuonMockupBuilder");

	field=ConstantBField(Vector3(0, 0, 0 * u.T));

	std::shared_ptr<const TrackingGeometry> trackingGeometry = detector.buildTrackingGeometry(gContext, trackingGeometryBuilder);

	// check the navigation policy
  	NavigationStream main;
  	AppendOnlyNavigationStream stream{main};
  	Vector3 startPos = {0., 0., 0.};
  	Vector3 startDir = {0., 1., 0.};
  	NavigationArguments args{startPos, startDir};
}

BOOST_AUTO_TEST_SUITE_END()
