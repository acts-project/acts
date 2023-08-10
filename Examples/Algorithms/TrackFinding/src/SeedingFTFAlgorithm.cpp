//basing on seeding ortho 
#include "ActsExamples/TrackFinding/SeedingFTFAlgorithm.hpp"

#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SimSeed.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
//might need for layer numbering funciton 
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "ActsExamples/EventData/Measurement.hpp"

#include <iostream>
#include <map>
#include <fstream>
#include <vector> 
#include <sstream>
#include <random>


using namespace std;

template class Acts::TrigFTF_GNN_Layer<ActsExamples::SimSpacePoint> ;
template class Acts::TrigFTF_GNN_Geometry<ActsExamples::SimSpacePoint> ;
template class Acts::TrigFTF_GNN_Node<ActsExamples::SimSpacePoint> ;
template class Acts::TrigFTF_GNN_EtaBin<ActsExamples::SimSpacePoint> ;
template class Acts::FTF_SP<ActsExamples::SimSpacePoint> ;
template class Acts::TrigFTF_GNN_DataStorage<ActsExamples::SimSpacePoint> ;
template class Acts::TrigFTF_GNN_Edge<ActsExamples::SimSpacePoint> ;

//constructor: 
ActsExamples::SeedingFTFAlgorithm::SeedingFTFAlgorithm(
    ActsExamples::SeedingFTFAlgorithm::Config cfg, 
    Acts::Logging::Level lvl) 
    : ActsExamples::IAlgorithm("SeedingAlgorithm", lvl), 
      m_cfg(std::move(cfg)) {

    // std::cout << "in initialise 1" << std::endl ;     
    //fill config struct
    m_cfg.layerMappingFile = m_cfg.layerMappingFile ; 
    // std::cout << "in initialise 2 " << std::endl ;     

    m_cfg.seedFilterConfig = m_cfg.seedFilterConfig.toInternalUnits();
    // std::cout << "in initialise 3 " << std::endl ;     

    m_cfg.seedFinderConfig =
    m_cfg.seedFinderConfig.toInternalUnits().calculateDerivedQuantities();
    // std::cout << "in initialise 4 " << std::endl ;     

    m_cfg.seedFinderOptions =
      m_cfg.seedFinderOptions.toInternalUnits().calculateDerivedQuantities(
          m_cfg.seedFinderConfig);
    // std::cout << "in initialise 5 " << std::endl ;     
    
    
    m_cfg.ACTS_FTF_Map = Make_ACTS_FTF_Map() ;
    // std::cout << "in initialise 6 " << std::endl ;     


    m_cfg.seedFinderConfig.input_vector = LayerNumbering() ; 
    // std::cout << "in initialise 7 " << std::endl ;     


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
    // std::cout << "in initialise 8 " << std::endl ;     

    m_outputSeeds.initialize(m_cfg.outputSeeds);
    // std::cout << "in initialise 9 " << std::endl ;     


    m_inputSourceLinks.initialize(m_cfg.inputSourceLinks);
    // std::cout << "in initialise 9 " << std::endl ;     


    m_cfg.seedFinderConfig.seedFilter =
      std::make_unique<Acts::SeedFilter<SimSpacePoint>>(
          Acts::SeedFilter<SimSpacePoint>(m_cfg.seedFilterConfig));
    // std::cout << "in initialise 10 " << std::endl ;     

    
    // m_seedFinder = Acts::SeedFinderFTF<SimSpacePoint>(m_cfg.seedFinderConfig);
    // std::cout << "in initialise 11 " << std::endl ;     

      } //this is not FTF config type because it is a meber of the algs config, which is of type FTF cofig  




//exectute: 
ActsExamples::ProcessCode ActsExamples::SeedingFTFAlgorithm::execute(
    const AlgorithmContext& ctx) const {

  // std::cout << "in execute 1 " << std::endl ;     
    
  std::vector<Acts::FTF_SP<SimSpacePoint>> FTF_spacePoints = Make_FTF_spacePoints(ctx, m_cfg.ACTS_FTF_Map);


  // //loop over FTF_SP vector to plot the space points 
  // fstream fout;
  // fout.open("FTF_SP_output.csv", ios::out | ios::app);
  // for (auto sp : FTF_spacePoints){
  //   fout << sp.FTF_ID << ", "
  //     << sp.SP->z() << ", "
  //     << sp.SP->r() 
  //     << "\n";
  // } 
  // std::cout << "in execute 2 " << std::endl ;     


  for (auto sp : FTF_spacePoints){
    ACTS_DEBUG("FTF space points: " << " FTF_id: "<< sp.FTF_ID << " z: " 
        << sp.SP->z() << " r: " << sp.SP->r() << " ACTS volume:  " 
        << sp.SP->sourceLinks().front().geometryId().volume() << "\n");
  }   
  
  // std::cout << "in execute 3 " << std::endl ;     

  
  //this is now calling on a core algorithm
  Acts::SeedFinderFTF<SimSpacePoint> finder(m_cfg.seedFinderConfig);
  // std::cout << "in execute 4 " << std::endl ;     

  //need this function as create_coords is needed for seeds 
  std::function<std::pair<Acts::Vector3, Acts::Vector2>(
      const SimSpacePoint *sp)>
      create_coordinates = [](const SimSpacePoint *sp) {
        Acts::Vector3 position(sp->x(), sp->y(), sp->z());
        Acts::Vector2 variance(sp->varianceR(), sp->varianceZ());
        return std::make_pair(position, variance);
      };  
      //output of function needed for seed
  // std::cout << "in execute 5 " << std::endl ;     

  //crashing here 
  finder.loadSpacePoints(FTF_spacePoints); 
  // std::cout << "in execute 6 " << std::endl ;     


  //still to develop, input will be FTF_spacePoints 
  SimSeedContainer seeds = finder.createSeeds(m_cfg.seedFinderOptions,
                                              FTF_spacePoints, create_coordinates);

  // std::cout << "in execute 7 " << std::endl ;     
                             

  m_outputSeeds(ctx, std::move(seeds));


  return ActsExamples::ProcessCode::SUCCESS;
}


std::map<std::pair<int, int>,std::pair<int, int>> ActsExamples::SeedingFTFAlgorithm::Make_ACTS_FTF_Map() const {

 //create map from csv 
  // map<std::pair<int, int>,int> ACTS_FTF;
  map<std::pair<int, int>,std::pair<int, int>> ACTS_FTF;
  std::ifstream data(m_cfg.layerMappingFile);  //0 in this file refers to no FTF ID 
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
  // file in format ACTS_vol,ACTS_lay,ACTS_mod,FTF_id
  for (auto i : parsedCsv){
        
      int ACTS_vol = stoi(i[0]); 
      int ACTS_lay = stoi(i[1]);
      int ACTS_mod = stoi(i[2]);
      int FTF = stoi(i[5]); 
      int eta_mod = stoi(i[6]);
      int ACTS_joint = ACTS_vol*100 + ACTS_lay ;
      // ACTS_FTF.insert({{ACTS_joint,ACTS_mod},FTF}) ; 
      ACTS_FTF.insert({{ACTS_joint,ACTS_mod},{FTF,eta_mod}}) ; 
      //here think should have pair of new number, vol*100+layer, mod which is 0 or number 
  }

  return ACTS_FTF ; 

}

std::vector<Acts::FTF_SP<ActsExamples::SimSpacePoint>> ActsExamples::SeedingFTFAlgorithm::Make_FTF_spacePoints(const AlgorithmContext& ctx, std::map<std::pair<int, int>,std::pair<int, int>> map) const {
  
    //create space point vectors 
  std::vector<const ActsExamples::SimSpacePoint *> spacePoints;
  std::vector<Acts::FTF_SP<ActsExamples::SimSpacePoint>> FTF_spacePoints ; 
  FTF_spacePoints.reserve(m_inputSpacePoints.size()); //not sure if this is enough, each one has several sp 

 //for loop filling space
  for (const auto &isp : m_inputSpacePoints) {
    for (const auto &spacePoint : (*isp)(ctx)) {
      //fill originial space point vector 
      spacePoints.push_back(&spacePoint);

      //FTF space point vector 
       //loop over space points, call on map 
      auto source_link = spacePoint.sourceLinks() ; 
      //warning if source link empty 
      if (source_link.empty()){
        //warning in officaial acts format 
        ACTS_WARNING("warning source link vector is empty");
        continue; 
      }

      int ACTS_vol_id = source_link.front().geometryId().volume() ;
      int ACTS_lay_id = source_link.front().geometryId().layer() ; 
      int ACTS_mod_id = source_link.front().geometryId().sensitive() ; 

      //dont want strips or HGTD, at the moment never called so commented out 
      if (ACTS_vol_id == 2 or  ACTS_vol_id == 22 or ACTS_vol_id == 23 or ACTS_vol_id == 24){ 
        continue; 
      }
      
      // Search for vol, lay 0n if doesnt esist (end) then search for full thing 
      //vol*100+lay as first number in pair then 0 or mod id 
      auto ACTS_joint_id = ACTS_vol_id*100 + ACTS_lay_id ;
      auto key = std::make_pair(ACTS_joint_id,0) ; //here the key needs to be pair of(vol*100+lay, 0) 
      auto Find = map.find(key) ;

      if (Find == map.end()){ //if end then make new key of (vol*100+lay, modid)
        key = std::make_pair(ACTS_joint_id,ACTS_mod_id) ; //mod ID 
        Find = map.find(key) ;
      } 
        
      //warning if key not in map 
      if (Find == map.end()){
        ACTS_WARNING("Key not found in FTF map for volume id: " << ACTS_vol_id << " and layer id: " << ACTS_lay_id  );
        continue; 
      } 

      //now should be pixel with FTF ID: 
      // int FTF_id = Find->second ; 
      int FTF_id = Find->second.first ; //new map the item is a pair so want first from it 


      //backup warning shouldnt need as no 0 in csv 
      if (FTF_id == 0) {
        ACTS_WARNING("No assigned FTF ID for key for volume id: " << ACTS_vol_id << " and layer id: " << ACTS_lay_id) ;
      }

      //need to have combined id here 
      //access eta mod from map 
      int eta_mod = Find->second.second ;

      int combined_id = FTF_id*1000 + eta_mod ; 


      //fill FTF vector with current sapce point and ID 
      FTF_spacePoints.emplace_back(&spacePoint,FTF_id,combined_id); 
    }
  }
  ACTS_VERBOSE("Space points successfully assigned FTF ID") ; 



  return FTF_spacePoints ; //not sure if I also need to return normal space point vector 
}


std::vector<Acts::TrigInDetSiLayer> ActsExamples::SeedingFTFAlgorithm::LayerNumbering() const {
 
  std::vector<Acts::TrigInDetSiLayer> input_vector ;
  std::vector<size_t> count_vector ;

  // //make a map of combined ID, vol, lay, z, r, mod id 
  // std::map<int,std::vector<int>> filling_map ;
  // fstream fout;
  // fout.open("Module_loop_sensitives.csv", ios::out | ios::app);  

  for (Acts::GeometryIdentifier geoId : m_cfg.geometrySelection) { //think this is the volume loop 

    m_cfg.trackingGeometry->visitSurfaces([this, &input_vector, &count_vector](const Acts::Surface* surface) {

        Acts::GeometryIdentifier geoid = surface->geometryId() ; 
        auto ACTS_vol_id = geoid.volume() ;
        auto ACTS_lay_id = geoid.layer() ;
        auto mod_id = geoid.sensitive() ; 
        auto bounds_vect = surface->bounds().values() ; 
        auto center = surface->center(Acts::GeometryContext()) ; //vector of 3 values //this needs geo context 
    

        // std::cout << "seeing what mod info is " << mod_id << std::endl ; 
        //make bounds global 
        Acts::Vector3 globalFakeMom(1, 1, 1); //fake mom, idea from hough 
        Acts::Vector2 min_bound_local = Acts::Vector2(bounds_vect[0],bounds_vect[1]) ; 
        Acts::Vector2 max_bound_local = Acts::Vector2(bounds_vect[2],bounds_vect[3]) ; 
        Acts::Vector3 min_bound_global = surface->localToGlobal(Acts::GeometryContext(),min_bound_local, globalFakeMom) ; 
        Acts::Vector3 max_bound_global = surface->localToGlobal(Acts::GeometryContext(),max_bound_local, globalFakeMom) ; 
        
        if(min_bound_global(0)>max_bound_global(0)){ //checking that not wrong way round 
          min_bound_global.swap(max_bound_global) ; 
        }

        float rc=0.0;
        float minBound = 100000.0;
        float maxBound =-100000.0;



        //convert to FTF ID 
        auto ACTS_joint_id = ACTS_vol_id*100 + ACTS_lay_id ;
        auto key = std::make_pair(ACTS_joint_id,0) ; //here the key needs to be pair of(vol*100+lay, 0) 
        auto Find = m_cfg.ACTS_FTF_Map.find(key) ;
        // int FTF_id = Find->second ;
        int FTF_id = Find->second.first ; //new map, item is pair want first 
        if (Find == m_cfg.ACTS_FTF_Map.end()){ //if end then make new key of (vol*100+lay, modid)
          key = std::make_pair(ACTS_joint_id,mod_id) ; //mod ID 
          Find = m_cfg.ACTS_FTF_Map.find(key) ; 
          FTF_id = Find->second.first; 
          // std::cout << "keys not for mod 0 " << ACTS_joint_id << " mod id: " << mod_id << " FTF: " << FTF_id<< std::endl ;
               
        } 

        short barrel_ec = 0; //a variable that says if barrrel, 0 = barrel, need to find what the n values are 
        // int eta_mod = 0 ; //could set eta mod to map memeber here, 80s have 0 in map too 
        int eta_mod = Find->second.second ;

        //assign barrel_ec depending on FTF_layer 
        if(79 < FTF_id && FTF_id < 85 ){ //80s, barrel 
          barrel_ec  = 0 ; 
          // eta_mod = 0 ;
        }
        else if(89 < FTF_id && FTF_id < 99) {//90s positive 
          barrel_ec  = 2 ;
          // eta_mod = Find->second->second ; //access second member of item 
        }
        else { //70s negative 
          barrel_ec  = -2 ;
          // eta_mod = Find->second->second ;
          // std::cout << "checking some eta mod values are correct " << eta_mod << std::endl ;

        }



        if(barrel_ec == 0){
            rc = sqrt(center(0)*center(0)+center(1)*center(1)); //barrel center in r 
            //bounds of z 
            if(min_bound_global(2) < minBound) minBound = min_bound_global(2);
            if(max_bound_global(2) > maxBound) maxBound = max_bound_global(2);
        }
        else{
            rc = center(2); //check this is right, not barrel center in Z 
            //bounds of r 
            float min = sqrt(min_bound_global(0)*min_bound_global(0)+min_bound_global(1)*min_bound_global(1)) ;
            float max = sqrt(max_bound_global(0)*max_bound_global(0)+max_bound_global(1)*max_bound_global(1)) ;
            if( min < minBound) minBound = min;
            if( max > maxBound) maxBound = max ;
        }


        int combined_id = FTF_id*1000 + eta_mod ; 
        auto current_index = find_if(input_vector.begin(), input_vector.end(),[combined_id](auto n) { return n.m_subdet == combined_id; }); 
        if (current_index != input_vector.end()) { //not end so does exist 
        //   //add to old, need index. FTF ID and type stay same! 
          size_t index = std::distance(input_vector.begin(), current_index);
          input_vector[index].m_refCoord += rc ; 
          input_vector[index].m_minBound += minBound ; 
          input_vector[index].m_maxBound += maxBound ; 
          count_vector[index] += 1 ; //increase count at the index 

        } 
        else { //end so doesnt exists 
          //make new if one with FTF ID doesnt exist: 
          Acts::TrigInDetSiLayer new_FTF_ID(combined_id, barrel_ec, rc, minBound, maxBound) ; 
          input_vector.push_back(new_FTF_ID) ; 
          count_vector.push_back(1) ; //so the element exists and not divinding by 0 

        }

        
        // //temporary making map 
        // int ACTS_IDs = ACTS_vol_id*100 + ACTS_lay_id; 
        // if (filling_map.find(ACTS_IDs) != filling_map.end()) { //not end so does exist 
        //   filling_map[ACTS_IDs][2] += center(2) ; //z 
        //   filling_map[ACTS_IDs][3] += sqrt(center(0)*center(0)+center(1)*center(1)) ; //r 
        //   filling_map[ACTS_IDs][4] += 1 ; //increase count

        // } 
        // else { //end so doesnt exists 
        //   filling_map[ACTS_IDs] = {ACTS_vol_id, ACTS_lay_id,center(2),sqrt(center(0)*center(0)+center(1)*center(1)) ,1 } ;
        // }

        //print to csv for each module 

        // fout << ACTS_vol_id << ", " //vol
        // << ACTS_lay_id << ", " //lay 
        // << mod_id << ", " //module 
        // << center(2)  << ", " //z 
        // << sqrt(center(0)*center(0)+center(1)*center(1))  //r 
        // << "\n";




    });


  } 

  // fstream fout;
  // fout.open("Module_loop.csv", ios::out | ios::app);
  // for (auto  i : filling_map){
  //         fout << i.second[0] << ", " //vol
  //           << i.second[1] << ", " //lay 
  //           << i.second[2]/i.second[4]  << ", " //z 
  //           << i.second[3]/i.second[4]  //r 
  //           << "\n";
  // }   
  
  for (int i = 0; i < input_vector.size(); i++){ 
    input_vector[i].m_refCoord = input_vector[i].m_refCoord/count_vector[i] ; 
    //do I need to also divide bounds?  
    // std::cout << "checking vector is correct" << input_vector[i].m_subdet << std::endl ;  

  }  
  // std::cout << "checking vector is filled" << input_vector.size() << count_vector.size();  

  return input_vector ; 

  
} 

