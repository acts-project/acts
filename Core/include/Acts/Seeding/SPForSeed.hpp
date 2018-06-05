// This file is part of the ACTS project.
//
// Copyright (C) 2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include <limits>
#include <array>

#include "Acts/Seeding/SpacePointConcept.hpp" 

namespace Acts {
  class SPForSeed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    SPForSeed() = delete;
    SPForSeed(const Acts::concept::AnySpacePoint<>*,
              const std::array<float,2> offsetXY,
              const std::array<float,2> cov)                ;
    SPForSeed(const SPForSeed&)                             ;
    ~SPForSeed() = default                                  ;
    SPForSeed& operator  = (const SPForSeed&)               ;

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

  inline SPForSeed& SPForSeed::operator = 
    (const SPForSeed& sp) 
    {
      spacepoint  = sp.spacepoint;
      m_x         = sp.m_x       ;
      m_y         = sp.m_y       ;
      m_z         = sp.m_z       ;
      m_r         = sp.m_r       ;
      m_covr      = sp.m_covr    ;
      m_covz      = sp.m_covz    ;
      return(*this);
    }
 
  inline SPForSeed::SPForSeed
    (const Acts::concept::AnySpacePoint<>* sp, const std::array<float,2> offsetXY, const std::array<float,2> cov) 
    {
      spacepoint = sp  ;
      m_x        = sp->x()+offsetXY[0];
      m_y        = sp->y()+offsetXY[1];
      m_z        = sp->z();
      m_r        = sqrt(m_x*m_x+m_y*m_y);
      m_covr     = cov[0];
      m_covz     = cov[1];
    }


  /////////////////////////////////////////////////////////////////////////////////
  // Copy constructor
  /////////////////////////////////////////////////////////////////////////////////

  
  inline SPForSeed::SPForSeed (const SPForSeed& sp)
    {
      *this = sp;
    }

} // end of namespace Acts
