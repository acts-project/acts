#ifndef SPForSeed_h
#define SPForSeed_h

#include "ACTS/Seeding/SpacePoint.hpp"

#include <cmath>
#include <limits>
#include <array>


//COLLECTION OF MAGIC NUMBERS IN HERE:
//
//
//m_q        = 100000.;
//
//float cov = wid*wid*.08333;
//
//Pixel Case
//    if(de->isBarrel()) {m_covz = 9.*cov; m_covr = .06;} 
//    else               {m_covr = 9.*cov; m_covz = .06;}
//
//SCT Case
//    if(de->isBarrel()) {m_covz = 8.*f22; m_covr = .1;}
//    else               {m_covr = 8.*f22; m_covz = .1;}
// currently "calibrationTool" function doesn't set either

namespace Acts {
namespace Seeding{
  class SPForSeed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    SPForSeed();
    SPForSeed(SpacePoint* const&, const std::array<float,2>)             ;
    SPForSeed(const SPForSeed&)                             ;
    virtual ~SPForSeed()                                    ;
    SPForSeed& operator  = (const SPForSeed&)               ;

    void set(SpacePoint* const&, const std::array<float,2>)   ;
    void setQuality(float)                                  ;

    const SpacePoint* spacepoint                            ; 
    const float&          x() const {return m_x;}
    const float&          y() const {return m_y;}
    const float&          z() const {return m_z;}
    const float&     radius() const {return m_r;}
          float         phi() const {return atan2(m_y,m_x);}
    const float&       covr() const {return m_covr;}
    const float&       covz() const {return m_covz;}
    const float&    quality() const {return m_q ;}

    const int&      surface() const {return m_surface;}

  protected:
    
    float m_x   ; // x-coordinate in beam system coordinates  
    float m_y   ; // y-coordinate in beam system coordinates
    float m_z   ; // z-coordinate in beam system coordinetes
    float m_r   ; // radius       in beam system coordinates
    float m_covr; //
    float m_covz; //
    float m_q   ;

    int   m_surface; // surface identifier
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
      m_q     = 0.;
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
    m_q         = sp.m_q       ;
      }
      return(*this);
    }
 
  
  inline SPForSeed::SPForSeed
    (SpacePoint* const& sp, const std::array<float,2> r) 
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
    (SpacePoint* const& sp,const std::array<float,2> r)
    {
      spacepoint = sp  ;
      m_x        = sp->x()+r[0];
      m_y        = sp->y()+r[1];
      m_z        = sp->z();
      m_r        =sqrt(m_x*m_x+m_y*m_y);
      m_surface  = sp->surface;
      // setting quality low
      m_q        = std::numeric_limits<float>::min();
    } 
    
  
  inline void  SPForSeed::setQuality(float q)
    {
      if(q <= m_q) m_q = q;
    }
 
} // end of namespace Seeding
} // end of namespace Acts

#endif  // SPForSeed_h
