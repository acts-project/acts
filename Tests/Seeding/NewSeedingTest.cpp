#include "ACTS/Seeding/New_Seedmaker.hpp"
#include "ACTS/Seeding/BinFinder.hpp"
#include "ACTS/Seeding/SeedFilter.hpp"
#include "ACTS/Seeding/SpacePointConcept.hpp"

#include "TestQualityTool.hpp"
#include "TestCovarianceTool.hpp"
#include "SpacePoint.hpp"

#include <boost/type_erasure/any_cast.hpp>

#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>


std::vector<const Acts::concept::AnySpacePoint<>*> readFile(std::string filename){

  std::string line;
  int   layer, ns, hitid;
  float r, phi, z;
  std::vector<const Acts::concept::AnySpacePoint<>*> readSP;

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
        Acts::concept::AnySpacePoint<>* sp = new Acts::concept::AnySpacePoint<>(SpacePoint{x,y,z,r,layer,0.003,0.003});
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
  Acts::Seeding::New_Seedmaker a;
  std::vector<const Acts::concept::AnySpacePoint<>*> spVec = readFile("sp.txt");
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
  Acts::Seeding::SeedFilterConfig sfconf;
  auto qTool = std::make_shared<Acts::Seeding::QualityTool>(Acts::Seeding::QualityTool());
  config->seedFilter = std::make_unique<Acts::Seeding::SeedFilter>(Acts::Seeding::SeedFilter(sfconf,qTool));
  config->covarianceTool = std::make_unique<Acts::Seeding::CovarianceTool>(Acts::Seeding::CovarianceTool());
  std::shared_ptr<Acts::Seeding::Cache> cache = a.initialize(config);
  a.newEvent(spVec,cache,config);
  auto outputSeeds = a.production3Sp(cache,config);
  for(int i = 0; i < outputSeeds.size(); i++){
    for (auto spC : outputSeeds.at(i)->spacePoints()){
      const SpacePoint* sp = boost::type_erasure::any_cast<const SpacePoint *>(spC);
      std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", " << sp->z() << ") ";
    }
    std::cout << std::endl;
  }
  return 0;
}
