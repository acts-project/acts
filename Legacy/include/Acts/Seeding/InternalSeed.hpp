// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// InternalSeed.hpp Acts project
///////////////////////////////////////////////////////////////////

#pragma once
#include "Acts/Seeding/SPForSeed.hpp"
#include "Acts/Seeding/Seed.hpp"

namespace Acts {
namespace Seeding {
  template <typename SpacePoint>
  class InternalSeed
  {

    /////////////////////////////////////////////////////////////////////////////////
    // Public methods:
    /////////////////////////////////////////////////////////////////////////////////

  public:
    InternalSeed();
    InternalSeed(SPForSeed<SpacePoint>*&,
                 SPForSeed<SpacePoint>*&,
                 SPForSeed<SpacePoint>*&,
                 float);
    InternalSeed(const InternalSeed<SpacePoint>&);
    virtual ~InternalSeed();
    InternalSeed<SpacePoint>&
    operator=(const InternalSeed<SpacePoint>&);

    SPForSeed<SpacePoint>*
    spacepoint0()
    {
      return m_s0;
    }
    SPForSeed<SpacePoint>*
    spacepoint1()
    {
      return m_s1;
    }
    SPForSeed<SpacePoint>*
    spacepoint2()
    {
      return m_s2;
    }
    const float&
    z() const
    {
      return m_z;
    }
    const float&
    quality() const
    {
      return m_q;
    }

    void
    set(SPForSeed<SpacePoint>*&,
        SPForSeed<SpacePoint>*&,
        SPForSeed<SpacePoint>*&,
        float);

    bool
    setQuality(float);

    bool
    set3(Acts::Seeding::Seed<SpacePoint>&);

  protected:
    SPForSeed<SpacePoint>* m_s0;
    SPForSeed<SpacePoint>* m_s1;
    SPForSeed<SpacePoint>* m_s2;
    float                  m_z;
    float                  m_q;
  };

  /////////////////////////////////////////////////////////////////////////////////
  // Inline methods
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::InternalSeed()
  {
    m_s0 = 0;
    m_s1 = 0;
    m_s2 = 0;
    m_z  = 0.;
    m_q  = 0.;
  }

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>&
  InternalSeed<SpacePoint>::operator=(const InternalSeed& sp)
  {
    if (&sp != this) {

      m_z  = sp.m_z;
      m_q  = sp.m_q;
      m_s0 = sp.m_s0;
      m_s1 = sp.m_s1;
      m_s2 = sp.m_s2;
    }
    return (*this);
  }

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::InternalSeed(SPForSeed<SpacePoint>*& s0,
                                                SPForSeed<SpacePoint>*& s1,
                                                SPForSeed<SpacePoint>*& s2,
                                                float                   z)
  {
    set(s0, s1, s2, z);
    m_q = 0.;
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Copy constructor
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::InternalSeed(const InternalSeed& sp)
    : m_s0(sp.m_s0), m_s1(sp.m_s1), m_s2(sp.m_s2)
  {
    *this = sp;
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Destructor
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline InternalSeed<SpacePoint>::~InternalSeed()
  {
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Set
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline void
  InternalSeed<SpacePoint>::set(SPForSeed<SpacePoint>*& s0,
                                SPForSeed<SpacePoint>*& s1,
                                SPForSeed<SpacePoint>*& s2,
                                float                   z)
  {
    m_z  = z;
    m_s0 = s0;
    m_s1 = s1;
    m_s2 = s2;
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Set three space points seed
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline bool
  InternalSeed<SpacePoint>::set3(Acts::Seeding::Seed<SpacePoint>& s)
  {

    bool pixb = !m_s0->spacepoint->clusterList().second;
    bool pixt = !m_s2->spacepoint->clusterList().second;

    if (pixb != pixt) {
      if (m_q > m_s0->quality() && m_q > m_s1->quality()
          && m_q > m_s2->quality()) {
        return false;
      }
    }

    m_s0->setQuality(m_q);
    m_s1->setQuality(m_q);
    m_s2->setQuality(m_q);

    s.erase();
    s.add(m_s0->spacepoint);
    s.add(m_s1->spacepoint);
    s.add(m_s2->spacepoint);
    s.setZVertex(double(m_z));
    return true;
  }

  /////////////////////////////////////////////////////////////////////////////////
  // Set quality in pro seed
  /////////////////////////////////////////////////////////////////////////////////

  template <typename SpacePoint>
  inline bool
  InternalSeed<SpacePoint>::setQuality(float q)
  {
    m_q       = q;
    bool pixb = !m_s0->spacepoint->clusterList().second;
    bool pixt = !m_s2->spacepoint->clusterList().second;
    if (pixb == pixt) {
      m_s0->setQuality(q);
      m_s1->setQuality(q);
      m_s2->setQuality(q);
      return true;
    }
    if (q < m_s0->quality() || q < m_s1->quality() || q < m_s2->quality()) {
      return true;
    }
    return false;
  }

}  // end of Seeding namespace
}  // end of Acts namespace