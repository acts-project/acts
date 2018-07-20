#pragma once

#include <list>

namespace Acts {
  template<typename SpacePoint>
  class Seed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    Seed();
    Seed
      (const SpacePoint*&, const SpacePoint*&, const SpacePoint*&, const double&);
    Seed(const Seed&) = default;
    Seed& operator = (const Seed&);
    virtual ~Seed();
    
    std::list<const SpacePoint*> m_spacepoints;
    double                       m_zvertex    ;  
  };

  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////
 

  
  
  ///////////////////////////////////////////////////////////////////////////////
  // Constructors
  ///////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint> 
  Seed<SpacePoint>::Seed(){}

  template <typename SpacePoint>
  Seed<SpacePoint>::Seed(const SpacePoint*& b, const SpacePoint*& m, const SpacePoint*& u, const double& vertex){
    m_zvertex = vertex;
    m_spacepoints.push_back(b);
    m_spacepoints.push_back(m);
    m_spacepoints.push_back(u);

  }

  template <typename SpacePoint>
  Seed<SpacePoint>::~Seed(){
  }
} // end of Acts namespace
