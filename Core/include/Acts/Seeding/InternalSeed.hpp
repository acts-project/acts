// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"

#include <memory>

namespace Acts {
  template<typename SpacePoint>
  class InternalSeed {
    
    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////
    
  public:
    
    InternalSeed();
    InternalSeed(const InternalSpacePoint<SpacePoint>*,const InternalSpacePoint<SpacePoint>*,const InternalSpacePoint<SpacePoint>*,float);
    InternalSeed(const InternalSeed&);
    virtual ~InternalSeed();
    InternalSeed& operator  = (const InternalSeed&);

    const std::array<const InternalSpacePoint<SpacePoint> *, 3> sp;
    float                      z()          const {return m_z ;}

  protected:
    
    float                      m_z ;
  };


  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::InternalSeed ()
    {
      sp={0,0,0};
      m_z  = 0.;
    }

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>& InternalSeed<SpacePoint>::operator = 
    (const InternalSeed<SpacePoint>& seed) 
    {
      if(&seed!=this) {

    m_z   = seed.m_z ;
    sp    = seed.sp;
      }
      return(*this);
    }
 
  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::InternalSeed
    (const InternalSpacePoint<SpacePoint>* s0,const InternalSpacePoint<SpacePoint>* s1,const InternalSpacePoint<SpacePoint>* s2,float z):sp({s0,s1,s2})
    {
      m_z   = z ;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Copy constructor
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::InternalSeed (const InternalSeed<SpacePoint>& seed): sp(seed.sp)
    {
      *this = seed;
    }

  /////////////////////////////////////////////////////////////////////////////////
  // Destructor
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::~InternalSeed() 
  {
  }

} // end of Acts namespace

