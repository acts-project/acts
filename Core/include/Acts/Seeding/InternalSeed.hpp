#pragma once

#include "Acts/Seeding/SPForSeed.hpp"
#include "Acts/Seeding/Seed.hpp"

#include <memory>

namespace Acts {
  class InternalSeed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    InternalSeed();
    InternalSeed(std::shared_ptr<SPForSeed>,std::shared_ptr<SPForSeed>,std::shared_ptr<SPForSeed>,float);
    InternalSeed(const InternalSeed&);
    virtual ~InternalSeed();
    InternalSeed& operator  = (const InternalSeed&);

    std::shared_ptr<SPForSeed> spacepoint0() {return m_s0;}
    std::shared_ptr<SPForSeed> spacepoint1() {return m_s1;}
    std::shared_ptr<SPForSeed> spacepoint2() {return m_s2;}
    float                      z()           {return m_z ;}

  protected:
    
    std::shared_ptr<SPForSeed> m_s0;
    std::shared_ptr<SPForSeed> m_s1;
    std::shared_ptr<SPForSeed> m_s2;
    float                      m_z ;
  };


  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////

  inline InternalSeed::InternalSeed ()
    {
      m_s0 = 0 ;
      m_s1 = 0 ;
      m_s2 = 0 ;
      m_z  = 0.;
    }

  inline InternalSeed& InternalSeed::operator = 
    (const InternalSeed& sp) 
    {
      if(&sp!=this) {

    m_z   = sp.m_z ;
    m_s0  = sp.m_s0;
    m_s1  = sp.m_s1;
    m_s2  = sp.m_s2;
      }
      return(*this);
    }
 
  inline InternalSeed::InternalSeed
    (std::shared_ptr<SPForSeed> s0,std::shared_ptr<SPForSeed> s1,std::shared_ptr<SPForSeed> s2,float z)
    {
      m_z   = z ;
      m_s0  = s0;
      m_s1  = s1;
      m_s2  = s2;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Copy constructor
  /////////////////////////////////////////////////////////////////////////////////

  inline InternalSeed::InternalSeed (const InternalSeed& sp): m_s0(sp.m_s0),m_s1(sp.m_s1),m_s2(sp.m_s2)
    {
      *this = sp;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Destructor
  /////////////////////////////////////////////////////////////////////////////////

  inline InternalSeed::~InternalSeed() 
  {
  }

} // end of Acts namespace

