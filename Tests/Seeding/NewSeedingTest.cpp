#include "Acts/Seeding/New_Seedmaker.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SpacePointConcept.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/SPForSeed.hpp"

#include "SpacePoint.hpp"

#include <boost/type_erasure/any_cast.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>


std::vector<const SpacePoint*> readFile(std::string filename){

  std::string line;
  int   layer, ns, hitid;
  float r, phi, z;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()){
    while(!spFile.eof())
    {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y;
      if (linetype == "lxyz"){
        ss >> layer >> x >> y >> z;
        r = std::sqrt(x*x+y*y);
        SpacePoint * sp = new SpacePoint{x,y,z,r,layer,0.003,0.003};
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
  config.maxSeedsPerSpM = 10;
  //2.7 eta
  config.cotThetaMax = 7.40627;

  config.minPt = 400.;

  config.bottomBinFinder = std::make_unique<Acts::BinFinder>(Acts::BinFinder());
  config.topBinFinder = std::make_unique<Acts::BinFinder>(Acts::BinFinder());
  Acts::SeedFilterConfig sfconf;
  config.seedFilter = std::make_unique<Acts::SeedFilter>(Acts::SeedFilter(sfconf));
  Acts::New_Seedmaker a(config);

  std::function<Acts::Vector2D(const SpacePoint*,float,float,float)> ct = [=]
          (const SpacePoint* sp,float zAlign,float rAlign,float sigma=1)
          -> Acts::Vector2D
          { Acts::Vector2D cov;
            cov[0] = ((*sp).covr + rAlign*rAlign) * sigma;
            cov[1] = ((*sp).covz + zAlign*zAlign) * sigma;
            return cov;
          };

  std::shared_ptr<Acts::SeedmakerState> state = a.initState(spVec, ct);
  a.createSeeds(state);
  while(!(state->outputQueue.empty())){
    std::shared_ptr<Acts::InternalSeed> seed = state->outputQueue.front();
    state->outputQueue.pop();
    std::shared_ptr<Acts::SPForSeed> spC = seed->spacepoint0();
    const SpacePoint* sp = spVec[spC->spIndex()];
    std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    spC = seed->spacepoint1();
    sp = spVec[spC->spIndex()];
    std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    spC = seed->spacepoint2();
    sp = spVec[spC->spIndex()];
    std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    std::cout << std::endl;
  }
  return 0;
}
