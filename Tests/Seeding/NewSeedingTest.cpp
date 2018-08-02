#include "Acts/Seeding/New_Seedmaker.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"

#include "SpacePoint.hpp"
#include "ATLASCuts.hpp"

#include <boost/type_erasure/any_cast.hpp>
 
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>


std::vector<const SpacePoint*> readFile(std::string filename){

  std::string line;
  int   layer;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()){
    while(!spFile.eof())
    {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r, covr, covz;
      if (linetype == "lxyz"){
        ss >> layer >> x >> y >> z >> covr >> covz;
        r = std::sqrt(x*x+y*y);
        SpacePoint * sp = new SpacePoint{x,y,z,r,layer,covr,covz};
   //     if(r < 200.){
   //       sp->setClusterList(1,0);
   //     }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

int main(){
  std::vector<const SpacePoint*> spVec = readFile("sp.txt");
  std::cout << "size of read SP: " <<spVec.size() << std::endl;

  Acts::SeedmakerConfig config;
// silicon detector max
  config.rMax = 600.;
  config.deltaRMin = 5.;
  config.deltaRMax = 270.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  //2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 400.;

  config.beamPos={-.5,-.5};

  auto bottomBinFinder = std::make_shared<Acts::BinFinder>(Acts::BinFinder());
  auto topBinFinder = std::make_shared<Acts::BinFinder>(Acts::BinFinder());
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts atlasCuts = Acts::ATLASCuts();
  config.seedFilter = std::make_unique<Acts::SeedFilter>(Acts::SeedFilter(sfconf, &atlasCuts));
  Acts::New_Seedmaker a(config);

  // covariance tool, sets covariances per spacepoint as required
  std::function<Acts::Vector2D(const SpacePoint*,float,float,float)> ct = [=]
    (const SpacePoint* sp,float zAlign,float rAlign,float sigma=1)
    -> Acts::Vector2D
    {
      return {sp->covz,sp->covr};
    };

  std::shared_ptr<Acts::SeedmakerState> state = a.initState(spVec, ct, bottomBinFinder, topBinFinder);
  auto start = std::chrono::system_clock::now();
  for(Acts::SeedingStateIterator it = state->begin(); !(it == state->end()); ++it){
    a.createSeedsForSP(it, state);
  }
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "time to create seeds: " << elapsed_seconds.count() << std::endl;
  std::cout << "Seeds created: "<<state->outputQueue.size() << std::endl;
  while(!(state->outputQueue.empty())){
    std::unique_ptr<const Acts::InternalSeed> seed = std::move(state->outputQueue.front());
    state->outputQueue.pop();
    const Acts::InternalSpacePoint* spC = seed->spacepoint0();
    const SpacePoint* sp = spVec[spC->spIndex()];
//    std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    spC = seed->spacepoint1();
    sp = spVec[spC->spIndex()];
//    std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    spC = seed->spacepoint2();
    sp = spVec[spC->spIndex()];
//    std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
//    std::cout << std::endl;
  }
  return 0;
}
