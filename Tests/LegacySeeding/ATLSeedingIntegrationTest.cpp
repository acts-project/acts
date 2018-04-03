#include "ACTS/Seeding/ATL_Seedmaker.hpp"
#include <fstream>
#include <sstream>

#define BOOST_TEST_MODULE SeedmakerIntegrationTest
#include <boost/test/included/unit_test.hpp>

struct SpacePoint{
  float x;
  float y;
  float z;
  float r;
  float covr=0.03;
  float covz=0.03;
  std::pair<int,int> m_clusterList = std::pair<int,int>(1,1);
  void setClusterList(int first, int second) {
    m_clusterList = std::pair<int,int>(first,second);
  }
  const std::pair<int,int> clusterList() const {
    return m_clusterList;
  }
  int surface;

};

std::vector<SpacePoint*> readFile(std::string filename){

  std::string line;
  int   layer;//, ns, hitid; are never used
  std::vector<SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()){
    while(!spFile.eof())
    {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r;
      if (linetype == "lxyz"){
        ss >> layer >> x >> y >> z;
        r = std::sqrt(x*x+y*y);
        SpacePoint* sp = new SpacePoint();
        sp->x=x;
        sp->y=y;
        sp->z=z;
//        sp->x=r; //straight line for debug purposes
        sp->r=r;
        sp->surface = layer;
        if(r < 200){
          sp->setClusterList(1,0);
        }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

int runSeeding(std::vector<SpacePoint*> spVec ){

  Acts::ATL_Seedmaker<SpacePoint> seedMaker;
  seedMaker.newEvent(0,spVec.begin(),spVec.end());
  seedMaker.find3Sp();
  const Acts::Seed<SpacePoint>* seed = seedMaker.next();
  int numSeeds = 0;
  while(seed != 0){
    numSeeds++;
    seed = seedMaker.next();
  }
  for(auto sp : spVec){
    delete sp;
  }
  return numSeeds;
}

BOOST_AUTO_TEST_CASE(number_of_seeds_correct_){
  std::vector<SpacePoint*> sp = readFile("/afs/cern.ch/user/r/rlangenb/workspace/acts2/acts/Tests/LegacySeeding/singleTrack.txt");
  // compare with result obtained from original ATLAS code
  BOOST_CHECK(runSeeding(sp)==8);
}
