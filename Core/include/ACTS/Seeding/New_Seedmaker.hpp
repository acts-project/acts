#pragma once

#include <list>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <utility>
#include <array>
#include <memory>

#include "ACTS/Seeding/SPForSeed.hpp"
#include "ACTS/Seeding/InternalSeed.hpp"
#include "ACTS/Seeding/SeedmakerConfig.hpp"
#include "ACTS/Seeding/SeedmakerCache.hpp"
#include "ACTS/Seeding/SpacePointConcept.hpp"

namespace Acts {
  namespace Seeding{
  class New_Seedmaker 
    {
      ///////////////////////////////////////////////////////////////////
      // Public methods:
      ///////////////////////////////////////////////////////////////////
      
    public:
    
      New_Seedmaker () = default;
      ~New_Seedmaker() = default;

      ///////////////////////////////////////////////////////////////////
      // Methods to initialize tool for new event or region
      ///////////////////////////////////////////////////////////////////
      std::shared_ptr<Acts::Seeding::Cache>
      initialize(std::shared_ptr<Acts::Seeding::Config> config);

      void
      newEvent
       (std::vector<const Acts::concept::AnySpacePoint<>*> spVec,
        std::shared_ptr<Acts::Seeding::Cache> cache,
        std::shared_ptr<Acts::Seeding::Config> config);

      
      std::vector<std::shared_ptr<Seed> > 
      production3Sp (std::shared_ptr<Acts::Seeding::Cache> cache,
                     std::shared_ptr<Acts::Seeding::Config> config);

      
    protected:
              /**    @name Disallow default instantiation, copy, assignment */
  //@{
  New_Seedmaker(const New_Seedmaker&) = delete;
  New_Seedmaker &operator=(const New_Seedmaker&) = delete;
  //@}

      ///////////////////////////////////////////////////////////////////
      // Protected methods
      ///////////////////////////////////////////////////////////////////


      std::vector<std::shared_ptr<Seed> > 
      production3Sp (std::vector<std::shared_ptr<SPForSeed > > currentBin,
                          std::set<size_t > bottomBins,
                          std::set<size_t > topBins,
                          std::shared_ptr<Acts::Seeding::Cache> cache,
                          std::shared_ptr<Acts::Seeding::Config> config)   ;

      void transformCoordinates (std::vector<std::shared_ptr<SPForSeed> >& vec,
                                 std::shared_ptr<SPForSeed> spM,
                                 bool bottom,
                                 std::vector<LinCircle>& linCircle);
   };

  } // end of Seeding namespace
} // end of Acts namespace

///////////////////////////////////////////////////////////////////
// Object-function for curvature seeds comparison
///////////////////////////////////////////////////////////////////

class comCurvature  {
public:
  
  bool operator ()
  (const std::pair<float,Acts::Seeding::SPForSeed*>& i1, 
   const std::pair<float,Acts::Seeding::SPForSeed*>& i2)
  {
    return i1.first < i2.first;
  }
};

class comR {
public:
  
  bool operator()(const Acts::concept::AnySpacePoint<>* i1, const Acts::concept::AnySpacePoint<>* i2){
    return i1->r() < i2->r();
  }
};
#include "ACTS/Seeding/New_Seedmaker.ipp"
