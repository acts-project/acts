//basijng on seeding ortho and rosie test 
#include "ActsExamples/TrackFinding/SeedingFTFAlgorithm.hpp"

#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iostream>
#include <map>
#include <fstream>
#include <vector> 
#include <sstream>

using namespace std;

//constructor: 
ActsExamples::SeedingFTFAlgorithm::SeedingFTFAlgorithm(
    ActsExamples::SeedingFTFAlgorithm::Config cfg, 
    Acts::Logging::Level lvl) 
    : ActsExamples::IAlgorithm("SeedingAlgorithm", lvl), 
      m_cfg(std::move(cfg)) {
    //fill config struct
    m_cfg.layerMappingFile = m_cfg.layerMappingFile ; 
    //now interal units error on SeedFilter so using on all 3 

    m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
 
    //got error "SeedFInderOrthoConfig not in ACTS internal units": 
    m_cfg.seedFinderConfig =
    m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();
    //because this is ortho type uses function in this class 

    m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);


  for (const auto &spName : m_cfg.inputSpacePoints) {
    if (spName.empty()) {
      throw std::invalid_argument("Invalid space point input collection");
    }

    auto &handle = m_inputSpacePoints.emplace_back(
        std::make_unique<ReadDataHandle<SimSpacePointContainer>>(
            this,
            "InputSpacePoints#" + std::to_string(m_inputSpacePoints.size())));
    handle->initialize(spName);
  }


    //throw statements for differnet cases 
    m_outputSeeds.initialize(m_cfg.outputSeeds);


    //filling the memeber of class (not const) from member of const 
    //at the ned of both 
    m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(m_cfg.seedFilterConfig));
    
    //this is now calling to a core algorithm 
    m_seedFinder = Acts::SeedFinderFTF<SimSpacePoint>(m_cfg.seedFinderConfig);
      } //this is not FTF type because it is a meber of the algs config, which is of type FTF cofig  

//exectute: 

ActsExamples::ProcessCode ActsExamples::SeedingFTFAlgorithm::execute(
    const AlgorithmContext& ctx) const {

  //defining things need for seeds function, need to go through to understand
  //try changing to type originally needed for loadsp?  
  std::vector<const SimSpacePoint *> spacePoints;

//for loop filling space points, seeing if can do without 
  for (const auto &isp : m_inputSpacePoints) {
    for (const auto &spacePoint : (*isp)(ctx)) {
      spacePoints.push_back(&spacePoint);
    }
  } 

  Acts::SeedFinderFTF<SimSpacePoint> finder(m_cfg.seedFinderConfig); 
  //this is now calling on a core algorithm 

  //need this function as create_coords is needed for seeds 
  std::function<std::pair<Acts::Vector3, Acts::Vector2>(
      const SimSpacePoint *sp)>
      create_coordinates = [](const SimSpacePoint *sp) {
        Acts::Vector3 position(sp->x(), sp->y(), sp->z());
        Acts::Vector2 variance(sp->varianceR(), sp->varianceZ());
        return std::make_pair(position, variance);
      };  


      //output of function needed for seed

  //no output to this function, input is vector of simspacepoints  





//create map from csv 
  map<std::pair<int, int>,int> ACTS_FTF;
  //0 in this file refers to no FTF ID 
  std::ifstream data(m_cfg.layerMappingFile);
  std::string line;
  std::vector<std::vector<std::string> > parsedCsv;
  while(std::getline(data,line))
  {
      std::stringstream lineStream(line);
      std::string cell;
      std::vector<std::string> parsedRow;
      while(std::getline(lineStream,cell,','))
      {
          parsedRow.push_back(cell); 
      }

      parsedCsv.push_back(parsedRow);
            
  }

  for (auto i : parsedCsv){
        
      int FTF = stoi(i[0]); 
      int ACTS_vol = stoi(i[1]); 
      int ACTS_lay = stoi(i[2]);
      ACTS_FTF.insert({{ACTS_vol,ACTS_lay},FTF}) ; 
  }

  
  //make FTF_SP vector, keeping same const pointer type as before 
  //reserve space here 
  std::vector<FTF_SP<SimSpacePoint>> FTF_spacePoints; //vector to the objects not pointers 
  //defined input to LSP std::vector<const FTF_SP<external_spacepoint_t>*>
  FTF_spacePoints.reserve(m_inputSpacePoints.size()); //not sure if this is enough, each one has several sp 


   //loop over space points, call on map 
  for (const auto &isp : m_inputSpacePoints) {
    for (const auto &spacePoint : (*isp)(ctx)) {

      auto source_link = spacePoint.sourceLinks() ; 
      //warning if source link empty 
      if (source_link.empty()){
        //warning in officaial acts format 
        ACTS_WARNING("warning source link vector is empty");
        continue; 
      }


      int ACTS_vol_id = source_link.front().geometryId().volume() ;
      int ACTS_lay_id = source_link.front().geometryId().layer() ; 
      auto key = std::make_pair(ACTS_vol_id,ACTS_lay_id) ; 
      auto Find = ACTS_FTF.find(key) ;
        
      //dont want strips or HGTD, at the moment never called so commented out 
      // if (ACTS_vol_id == 2 or  ACTS_vol_id == 22 or ACTS_vol_id == 23 or ACTS_vol_id == 24){ 
      //   continue; 
      // }

      //warning if key not in map 
      if (Find == ACTS_FTF.end()){
        ACTS_WARNING("Key not found in FTF map for volume id: " << ACTS_vol_id << " and layer id: " << ACTS_lay_id  );
        continue; 
      } 

      //now should be pixel with FTF ID: 
      int FTF_id = Find->second ;

      //backup warning shouldnt need as no 0 in csv 
      if (FTF_id == 0) {
        ACTS_WARNING("No assigned FTF ID for key for volume id: " << ACTS_vol_id << " and layer id: " << ACTS_lay_id) ;
      }
 

      //fill FTF vector with current sapce point and ID 
      //when it fills spacePoints it uses & before   
      FTF_spacePoints.emplace_back(&spacePoint,FTF_id); 
      // FTF_spacePoints.emplace_back(const FTF_SP<SimSpacePoint>(&spacePoint,FTF_id));
    }
  }
  ACTS_VERBOSE("Space points successfully assigned FTF ID") ; 




  finder.loadSpacePoints(FTF_spacePoints); 
  //this createseeds function is defined in SeedFInderFTF in core 
  SimSeedContainer seeds = finder.createSeeds(m_cfg.seedFinderOptions,
                                              spacePoints, create_coordinates);

                             

  // //debug statement
  
  // ctx.eventStore.add(m_cfg.outputSeeds, std::move(seeds));
  // ctx.eventStore.add(m_cfg.outputProtoTracks, std::move(protoTracks));
  m_outputSeeds(ctx, std::move(seeds));



  return ActsExamples::ProcessCode::SUCCESS;
}

//defenitions of print functions, leaving for now 
