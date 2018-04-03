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
  return numSeeds;
}

BOOST_AUTO_TEST_CASE(number_of_seeds_correct_){
  std::vector<SpacePoint*> spVec;
  std::vector<int> layerVec{1,2,2,3,4,11,13,14};
  std::vector<float> xVec{-33.3403 ,-48.2369 ,-49.4129 ,-88.8567 ,-122.5566, -283.169, -412.277, -462.5564};
  std::vector<float> yVec{2.7288 ,4.5193 ,4.6755 ,11.1935,18.7696,83.1666,179.1006,232.9765};
  std::vector<float> zVec{-74.5553 , -91.9763, -93.3541, -139.779, -179.889, -381.403, -568.641, -654.2494};
  for (unsigned int i = 0; i< layerVec.size(); i++){
    SpacePoint * sp = new SpacePoint();
    sp->surface = layerVec.at(i);
    sp->x = xVec.at(i);
    sp->y = yVec.at(i);
    sp->z = zVec.at(i);
    sp->r = std::sqrt(sp->x*sp->x+sp->y*sp->y);
    if(sp->r < 200.){
      sp->setClusterList(1,0);
    }
    spVec.push_back(sp);
  }

  // compare with result obtained from original ATLAS code
  BOOST_CHECK(runSeeding(spVec)==8);
  for(auto sp : spVec){
    delete sp;
  }
}
