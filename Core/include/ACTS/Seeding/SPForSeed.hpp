#pragma once

#include <cmath>
#include <limits>
#include <array>

#include "ACTS/Seeding/SpacePointConcept.hpp" 

namespace Acts {
namespace Seeding{
  class SPForSeed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    SPForSeed()                                             ;
    SPForSeed(Acts::concept::AnySpacePoint<>* const&, const std::array<float,2>);
    SPForSeed(const SPForSeed&)                             ;
    ~SPForSeed() = default                                  ;
    SPForSeed& operator  = (const SPForSeed&)               ;

    void set(Acts::concept::AnySpacePoint<>* const&, const std::array<float,2>)   ;
    void setQuality(float)                                  ;

    const Acts::concept::AnySpacePoint<>* spacepoint                            ; 
    const float&          x() const {return m_x;}
    const float&          y() const {return m_y;}
    const float&          z() const {return m_z;}
    const float&     radius() const {return m_r;}
          float         phi() const {return atan2(m_y,m_x);}
    const float&       covr() const {return m_covr;}
    const float&       covz() const {return m_covz;}

  protected:
    
    float m_x   ; // x-coordinate in beam system coordinates  
    float m_y   ; // y-coordinate in beam system coordinates
    float m_z   ; // z-coordinate in beam system coordinetes
    float m_r   ; // radius       in beam system coordinates
    float m_covr; //
    float m_covz; //

  };


  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////

  
  inline SPForSeed::SPForSeed ()
    {
      spacepoint = 0;
      m_x     = 0.;
      m_y     = 0.;
      m_z     = 0.;
      m_r     = 0.;
      m_covr  = 0.;
      m_covz  = 0.;
   }

  
  inline SPForSeed& SPForSeed::operator = 
    (const SPForSeed& sp) 
    {
      if(&sp!=this) {
    spacepoint  = sp.spacepoint;
    m_x         = sp.m_x       ;
    m_y         = sp.m_y       ;
    m_z         = sp.m_z       ;
    m_r         = sp.m_r       ;
    m_covr      = sp.m_covr    ;
    m_covz      = sp.m_covz    ;
      }
      return(*this);
    }
 
  
  inline SPForSeed::SPForSeed
    (Acts::concept::AnySpacePoint<>* const& sp, const std::array<float,2> r) 
    {
      set(sp,r); 
    }


  /////////////////////////////////////////////////////////////////////////////////
  // Copy constructor
  /////////////////////////////////////////////////////////////////////////////////

  
  inline SPForSeed::SPForSeed (const SPForSeed& sp)
    {
      *this = sp;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Set
  /////////////////////////////////////////////////////////////////////////////////

  
  inline void SPForSeed::set
    (Acts::concept::AnySpacePoint<>* const& sp,const std::array<float,2> r)
    {
      spacepoint = sp  ;
      m_x        = sp->x()+r[0];
      m_y        = sp->y()+r[1];
      m_z        = sp->z();
      m_r        =sqrt(m_x*m_x+m_y*m_y);
    } 
 
} // end of namespace Seeding
} // end of namespace Acts
