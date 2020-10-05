// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string>
#include <fstream>
#include <chrono>
#include <boost/program_options.hpp>
#include "ActsExamples/DD4hepDetector/DD4hepDetectorOptions.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"
#include "ActsExamples/Geant4/InteractionProcessRecording.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/DD4hepG4/DD4hepToG4Svc.hpp"
#include "ActsExamples/Plugins/Root/RootMaterialTrackWriter.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"

#include "ActsExamples/Plugins/Geant4/OREventAction.hpp"
#include "ActsExamples/Framework/WriterT.hpp"
#include "ActsExamples/EventData/Barcode.hpp"

#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/Geometry/CommonGeometry.hpp"


namespace po = boost::program_options;

using PartRec = FW::Geant4::Collection;
//~ using PartRec = std::map<int, std::vector<FW::Geant4::ParticleRecord>>;
//~ using PartRec = 	std::map<int, std::vector<FW::Data::SimVertex<FW::Data::SimParticle>>;
using WriteIt = FW::WriterT<PartRec>;

namespace FW
{
class ParticleRecordWriting : public WriteIt
{
public:

	struct Config
    {
      std::string collection;              ///< particle collection to write
      std::shared_ptr<BarcodeSvc>
             barcodeSvc;          ///< the barcode service to decode (optional)
    };
    
    /// Constructor
    ///
    /// @param cfg Configuration struct
    /// @param level Message level declaration
    ParticleRecordWriting(const Config&        cfg,
                       Acts::Logging::Level level = Acts::Logging::INFO) : WriteIt(cfg.collection, "ParticleRecordWriter", level), m_cfg(cfg) 
    {
		using namespace std::chrono;
		milliseconds ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
		ofs.open("geantOut" + std::to_string(ms.count()) + ".txt");		
	}

    /// Virtual destructor
    ~ParticleRecordWriting() override 
    {
		ofs.close();
	}

    /// End-of-run hook
    ProcessCode
    endRun() final override 
    {
		  return ProcessCode::SUCCESS;
	} // TODO: Maybe write the particles here to file

  protected:
    /// @brief Write method called by the base class
    /// @param [in] context is the algorithm context for event information
    /// @param [in] vertices is the process vertex collection for the
    /// particles to be attached
    ProcessCode
    writeT(const AlgorithmContext&       context,
           const PartRec& collection) final override
	{
		int eventNr = context.eventNumber;
		
		std::cout << "Writing event " << eventNr << std::endl;
		ofs << eventNr << " " << collection.pdg << " " << collection.momentum << " " << collection.phi << " " << -log(tan(collection.theta * 0.5)) << "\n";
		
		std::set<std::string> procsElse;
		for(const auto& part : collection.particles)
		{
			//~ check the number of charged particles in the event, the type of particles, their energy spectrum, the angular distribution			
			if(part.second.back().volume == "No volume")
			{			
				writeToFile(part.second.back());
				
				std::set<std::string> procs;
				for(const auto& p : part.second)
				{
					procs.insert(p.process);
				}
				for(auto it = procs.cbegin(); it != procs.cend(); it++)
					ofs << " " << *it;
				ofs << "\n";
			}
			else
			{
				for(const auto& p : part.second)
				{
					procsElse.insert(p.process);
				}
			}
		}
		if(procsElse.empty())
			ofs << "#\n";
		else
		{
			ofs << "999 ";
			for(auto it = procsElse.cbegin(); it != procsElse.cend(); it++)
				ofs << *it << " ";
			ofs << "\n#\n";
		}
		
		return ProcessCode::SUCCESS;
	}

  private:
    Config     m_cfg;         ///< The config class
    std::ofstream ofs;
    
    void
    writeToFile(const FW::Geant4::ParticleRecord& p)
    {
		const double r = sqrt(p.momentum[0] * p.momentum[0] + p.momentum[1] * p.momentum[1] + p.momentum[2] * p.momentum[2]);
		Acts::Vector3D dir;
		dir << p.momentum[0] / r, p.momentum[1] / r, p.momentum[2] / r;
		
		const double theta = acos(dir.z());
		const double eta = -log(tan(theta * 0.5));
		
		const double phi = atan2(dir.y(), dir.x());
		
		ofs << p.pdg << " " << p.charge << " " << r << " " << phi << " " << eta;
	}
};
}

void
g4SequencerBuild(boost::program_options::variables_map& vm)
{
  FW::Sequencer g4Sequencer(FW::Options::readSequencerConfig(vm));
  
  using namespace std::chrono;
  milliseconds ms = duration_cast< milliseconds >(system_clock::now().time_since_epoch());
  int    randomSeed1 = ms.count();
  int    randomSeed2 = ms.count() + 10;
  
  Acts::GeometryContext geoContext;

  // DETECTOR:
  // --------------------------------------------------------------------------------
  // DD4Hep detector definition
  // read the detector config & dd4hep detector
  auto dd4HepDetectorConfig
      = FW::Options::readDD4hepConfig<po::variables_map>(vm);
  auto geometrySvc = std::make_shared<FW::DD4hep::DD4hepGeometryService>(
      dd4HepDetectorConfig);
  std::shared_ptr<const Acts::TrackingGeometry> tGeometry
      = geometrySvc->trackingGeometry(geoContext);

  // DD4Hep to Geant4 conversion
  //
  FW::DD4hepG4::DD4hepToG4Svc::Config dgConfig("DD4hepToG4",
                                               Acts::Logging::INFO);
  dgConfig.dd4hepService = geometrySvc;
  auto dd4hepToG4Svc = std::make_shared<FW::DD4hepG4::DD4hepToG4Svc>(dgConfig);

  // --------------------------------------------------------------------------------
  // Geant4 JOB:
  // --------------------------------------------------------------------------------
  // set up the writer for
  // ---------------------------------------------------------------------------------

  // set up the algorithm writing out the material map
  FW::InteractionProcessRecording::Config g4rConfig = FW::Options::readInteractionProcessRecordingConfig(vm);
  g4rConfig.geant4Service  = dd4hepToG4Svc;
  g4rConfig.seed1          = randomSeed1;
  g4rConfig.seed2          = randomSeed2;
  g4rConfig.particleCollection = "geant-outcome-tracks";
  
  // create the geant4 algorithm
  auto g4rAlgorithm
      = std::make_shared<FW::InteractionProcessRecording>(g4rConfig, Acts::Logging::INFO);

  // Output directory
  std::string particleCollection = g4rConfig.particleCollection;

  // Write the propagation steps as ROOT TTree
  FW::ParticleRecordWriting::Config config;
  config.collection  = particleCollection;
  auto writer
	= std::make_shared<FW::ParticleRecordWriting>(
		config);
  g4Sequencer.addWriter(writer);

  // Append the algorithm and run
  g4Sequencer.addAlgorithm(g4rAlgorithm);
  g4Sequencer.run();
}

int
main(int argc, char* argv[])
{
	using po::value;

  DD4hepDetector detector;

  // Declare the supported program options.
  // Setup and parse options
  auto desc = FW::Options::makeDefaultOptions();
  FW::Options::addSequencerOptions(desc);
  FW::Options::addOutputOptions(desc);
  FW::Options::addGeometryOptions(desc);
  FW::Options::addParticleGunOptions(desc); // TODO: Replace whatever is given here and in outcomerecoptions
  FW::Options::addBFieldOptions(desc);
  FW::Options::addFatrasOptions(desc);
  FW::Options::addRandomNumbersOptions(desc);
  FW::Options::addMaterialOptions(desc);
  FW::Options::addDigitizationOptions(desc);
  desc.add_options()("evg-input-type",
                     value<std::string>()->default_value("gun"),
                     "Type of evgen input 'gun', 'pythia8'");
                     
  // Add specific options for this geometry
  detector.addOptions(desc);
  auto vm = FW::Options::parse(desc, argc, argv);
  if (vm.empty()) {
    return EXIT_FAILURE;
  }

std::cout << "Building g4 sequencer" << std::endl;
  g4SequencerBuild(vm);
std::cout << "Done" << std::endl;
}