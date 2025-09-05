// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/MuonSpectrometerMockupDetector/GeoMuonMockupExperiment.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelDetector.hpp"
#include "ActsExamples/GeoModelDetector/GeoModelMuonMockupBuilder.hpp"
#include "Acts/Plugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "Acts/Plugins/GeoModel/GeoModelTree.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Navigation/FrustumNavigationPolicy.hpp"

using namespace Acts;
using namespace ActsExamples;
using namespace Acts::GeoModel;
using namespace Acts::Experimental;

GeometryContext gContext;

auto logger = getDefaultLogger("FrustumNavigationTests", Logging::VERBOSE);

// This tests the frustum navigation policy for gen3 geometry interface
int main() {
	//create mockup detector
	GeoMuonMockupExperiment::Config detCfg;
	detCfg.nEtaStations=4;
	detCfg.nSectors=12;
	GeoMuonMockupExperiment mockDet(detCfg,getDefaultLogger("FrustumNavigationTestsGMM", Logging::VERBOSE));
	Acts::GeoModelTree gmTree=mockDet.constructMS();

	GeoModelDetector::Config gmdConfig;
	gmdConfig.geoModelTree=gmTree;
	GeoModelDetector muonDetector(gmdConfig);

	GeoModelDetectorObjectFactory::Config gmdObjFactoryCfg;
	gmdObjFactoryCfg.nameList={"RpcGasGap","MDTDriftGas"};
	gmdObjFactoryCfg.convertSubVolumes=true;
	gmdObjFactoryCfg.convertBox={"MDT"};
	GeoModelDetectorObjectFactory gmdObjFactory(gmdObjFactoryCfg,getDefaultLogger("FrustumNavigationTestsFactory", Logging::VERBOSE));
	GeoModelDetectorObjectFactory::Options options;
	options.queries={"Muon"};
	GeoModelDetectorObjectFactory::Cache cache;
	gmdObjFactory.construct(cache,gContext,gmTree,options);

	GeoModelMuonMockupBuilder::Config builderConfig;
	builderConfig.stationNames = {"Inner", "Middle", "Outer"};
	builderConfig.volumeBoxFPVs = cache.volumeBoxFPVs;
	GeoModelMuonMockupBuilder trackingGeometryBuilder(builderConfig,getDefaultLogger("FrustumNavigationTestsBuilder", Logging::VERBOSE));

	ConstantBField field=ConstantBField(Vector3(0, 0, 0 ) * UnitConstants::T);

	std::shared_ptr<const TrackingGeometry> trackingGeometry = muonDetector.buildTrackingGeometry(gContext, trackingGeometryBuilder);

	// check the navigation policy
  	NavigationStream main;
  	AppendOnlyNavigationStream stream{main};
  	Vector3 startPos = {0., 0., 0.};
  	Vector3 startDir = {0., 1., 0.};
  	NavigationArguments args{startPos, startDir};

	FrustumNavigationPolicy::Config frustumConfig;
	FrustumNavigationPolicy frustumNav(gContext,*trackingGeometry->highestTrackingVolume(),*logger,frustumConfig);
	NavigationDelegate delegate;
  	frustumNav.connect(delegate);
	delegate(args,stream, *logger);
	return 0;
}

