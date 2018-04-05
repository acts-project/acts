#include "ACTS/Seeding/New_Seedmaker.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>

//struct SpacePoint{
//  float x;
//  float y;
//  float z;
//  float r;
//  std::pair<int,int> m_clusterList = std::pair<int,int>(1,1);
//  void setClusterList(int first, int second) {
//    m_clusterList = std::pair<int,int>(first,second);
//  }
//  const std::pair<int,int> clusterList() const {
//    return m_clusterList;
//  }
//  int surface;
//
//};

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
        SpacePoint* sp = new SpacePoint();
        sp->x()=x;
        sp->y()=y;
        sp->z()=z;
//        sp->x=r; //straight line for debug purposes
        sp->r()=r;
        sp->surface = layer;
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

  std::unique_ptr<Acts::Seeding::Config> config = make_unique(Acts::Seeding::Config);
  Acts::Seeding::NewSeedmaker::Cache cache = a.initialize(
                        std::unique_ptr<Acts::Seeding::Config > config)
  a.newEvent(spVec,cache,config);
  a.production3Sp(cache,config);
  for(auto seed : cache->outputSeeds){
    for (const SpacePoint* const sp : seed->spacePoints()){
      std::cout << sp->surface << " (" << sp->x << ", " << sp->y << ", " << sp->z << ") ";
    }
    std::cout << std::endl;
    seed = a.next();
  }

  a.newEvent(1,spVec.begin(),spVec.end());
  return 0;
}
