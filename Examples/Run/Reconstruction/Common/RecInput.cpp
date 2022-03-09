// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Detector/IBaseDetector.hpp"
#ifdef ACTS_PLUGIN_ONNX
#include "Acts/Plugins/Onnx/MLTrackClassifier.hpp"
#endif
#include "ActsExamples/Digitization/DigitizationOptions.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"
#include "ActsExamples/Io/Json/JsonDigitizationConfig.hpp"
#include "ActsExamples/Io/Performance/CKFPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/SeedingPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFinderPerformanceWriter.hpp"
#include "ActsExamples/Io/Performance/TrackFitterPerformanceWriter.hpp"
#include "ActsExamples/Io/Root/RootTrajectoryStatesWriter.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TrackFinding/SeedingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/SpacePointMaker.hpp"
#include "ActsExamples/TrackFinding/TrackFindingAlgorithm.hpp"
#include "ActsExamples/TrackFinding/TrackFindingOptions.hpp"
#include "ActsExamples/TrackFitting/SurfaceSortingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingAlgorithm.hpp"
#include "ActsExamples/TrackFitting/TrackFittingOptions.hpp"
#include "ActsExamples/TruthTracking/TruthTrackFinder.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include "RecInput.hpp"

ActsExamples::CsvSimHitReader::Config setupSimHitReading(
																												 const ActsExamples::Options::Variables& vars,
																												 ActsExamples::Sequencer& sequencer) {
	using namespace ActsExamples;
	
	// Read some standard options
	auto logLevel = Options::readLogLevel(vars);
	
	// Read truth hits from CSV files
	auto simHitReaderCfg = Options::readCsvSimHitReaderConfig(vars);
	simHitReaderCfg.inputStem = "hits";
	simHitReaderCfg.outputSimHits = "hits";
	sequencer.addReader(
											std::make_shared<CsvSimHitReader>(simHitReaderCfg, logLevel));
	
	return simHitReaderCfg;
}

ActsExamples::CsvParticleReader::Config setupParticleReading(
																														 const ActsExamples::Options::Variables& vars,
																														 ActsExamples::Sequencer& sequencer) {
	using namespace ActsExamples;
	
	// Read some standard options
	auto logLevel = Options::readLogLevel(vars);
	
	// Read particles (initial states) and clusters from CSV files
	auto particleReader = Options::readCsvParticleReaderConfig(vars);
	particleReader.inputStem = "particles_initial";
	particleReader.outputParticles = "particles_initial";
	sequencer.addReader(
											std::make_shared<CsvParticleReader>(particleReader, logLevel));
	
	return particleReader;
}

ActsExamples::CsvSpacePointReader::Config setupSpacePointReading(
																																 const ActsExamples::Options::Variables& vars,
																																 ActsExamples::Sequencer& sequencer,
																																 const std::string& inputCollectionName) {
	using namespace ActsExamples;
	
	// Read some standard options
	auto logLevel = Options::readLogLevel(vars);
	
	// Read truth hits from CSV files
	auto spacePointReaderCfg = Options::readCsvSpacePointReaderConfig(vars);
	spacePointReaderCfg.inputStem = "spacepoints";
	spacePointReaderCfg.inputCollection = inputCollectionName;
	sequencer.addReader(
											std::make_shared<CsvSpacePointReader>(spacePointReaderCfg, logLevel));
	
	return spacePointReaderCfg;
}
