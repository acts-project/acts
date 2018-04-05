#ifndef New_Seedmaker_H
#define New_Seedmaker_H

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
#include "ACTS/Seeding/SpacePoint.hpp"

/// ist das kunst oder kann das weg?
#include <iostream>

namespace Acts {
  namespace Seeding{
  class New_Seedmaker 
    {
      ///////////////////////////////////////////////////////////////////
      // Public methods:
      ///////////////////////////////////////////////////////////////////
      
    public:
    


      
      ///////////////////////////////////////////////////////////////////
      // Standard tool methods
      ///////////////////////////////////////////////////////////////////
      
      New_Seedmaker ();
      virtual ~New_Seedmaker();

      ///////////////////////////////////////////////////////////////////
      // Methods to initialize tool for new event or region
      ///////////////////////////////////////////////////////////////////
      std::unique_ptr<Acts::Seeding::Cache>
      initialize(std::unique_ptr<Acts::Seeding::Config> config);

      void
      newEvent
       (std::vector<SpacePoint*> spVec,
        std::unique_ptr<Acts::Seeding::Cache>& cache,
        std::unique_ptr<Acts::Seeding::Config>& config);

      
    protected:
              /**    @name Disallow default instantiation, copy, assignment */
  //@{
  New_Seedmaker(const New_Seedmaker&) = delete;
  New_Seedmaker &operator=(const New_Seedmaker&) = delete;
  //@}

      ///////////////////////////////////////////////////////////////////
      // Protected methods
      ///////////////////////////////////////////////////////////////////

      void production3Sp (std::unique_ptr<Acts::Seeding::Cache>& cache,
                          std::unique_ptr<Acts::Seeding::Config>& config);

      void production3Sp (std::vector<std::shared_ptr<SPForSeed > > currentBin,
                          std::set<size_t > bottomBins,
                          std::set<size_t > topBins,
                          std::unique_ptr<Acts::Seeding::Cache>& cache,
                          std::unique_ptr<Acts::Seeding::Config>& config)   ;

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
  
  bool operator()(SpacePoint* i1, SpacePoint* i2){
    return i1->r() < i2->r();
  }
};
#include "ACTS/Seeding/New_Seedmaker.ipp"

#endif // New_Seedmaker_H
