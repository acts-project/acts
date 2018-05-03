#include "ACTS/Seeding/New_Seedmaker.hpp"
#include "ACTS/Seeding/BinFinder.hpp"
#include "ACTS/Seeding/SeedFilter.hpp"

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

namespace InDetDD{
  struct DetectorElement{
    bool m_isBarrel = true;
    bool m_isPixel = true;
    bool isBarrel(){
      return m_isBarrel;
    }
    bool isPixel(){
      return m_isPixel;
    }
  };
}

namespace InDet{
  struct SiCluster{
    std::shared_ptr<InDetDD::DetectorElement> de;
  };
}

struct SpacePoint{
  float x;
  float y;
  float z;
  float r;
  std::pair<shared_ptr<InDet::SiCluster>,shared_ptr<InDet::SiCluster> > m_clusterList(
                       std::make_shared<InDet::SiCluster>(InDet::SiCluster()),
                       std::make_shared<InDet::SiCluster>(InDet::SiCluster()));
  void setClusterList(int first, int second) {
    if(second==0){
      m_clusterList.second.reset();
    }
  }
  const std::pair<InDet::SiCluster,InDet::SiCluster> clusterList() const {
    return m_clusterList;
  }
  int surface;
};


std::vector<SpacePoint*> readFile(std::string filename){

  std::string line;
  int   layer, ns, hitid;
  float r, phi, z;
  std::vector<SpacePoint*> readSP;

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
        SpacePoint* sp = new SpacePoint(x,y,z,r,layer);
        if(r < 200.){
          sp->setClusterList(1,0);
        }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

int main(){
  Acts::Seeding::New_Seedmaker a;
  std::vector<SpacePoint*> spVec = readFile("sp.txt");
  std::cout << "size of read SP: " <<spVec.size() << std::endl;

  std::shared_ptr<Acts::Seeding::Config> config = std::make_shared<Acts::Seeding::Config>();
// silicon detector max
  config->rMax = 600.;
  config->deltaRMin = 5.;
  config->deltaRMax = 270.;
  config->collisionRegionMin = -250.;
  config->collisionRegionMax = 250.;
  config->zMin = -2800.;
  config->zMax = 2800.;
  config->maxSeedsPerSpM = 10;
  //2.7 eta
  config->cotThetaMax = 7.40627;

//  this is wrong on so many levels
  config->bottomBinFinder = std::make_unique<Acts::Seeding::BinFinder>(Acts::Seeding::BinFinder());
  config->topBinFinder = std::make_unique<Acts::Seeding::BinFinder>(Acts::Seeding::BinFinder());
  config->seedFilter = std::make_unique<Acts::Seeding::SeedFilter>(Acts::Seeding::SeedFilter(-100.));
  std::shared_ptr<Acts::Seeding::Cache> cache = a.initialize(config);
  a.newEvent(spVec,cache,config);
  a.production3Sp(cache,config);
  for(int i = 0; i < cache->outputSeeds.size(); i++){
    for (const SpacePoint* sp : cache->outputSeeds.at(i)->spacePoints()){
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    }
    std::cout << std::endl;
  }
  return 0;
}
