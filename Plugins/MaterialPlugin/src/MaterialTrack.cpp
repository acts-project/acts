// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialTrack.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialPlugins/MaterialTrack.hpp"

Acts::MaterialTrack::MaterialTrack(const MaterialStep::Position& startPos,
                                   double                        theta,
                                   double                        phi,
                                   std::vector<MaterialStep>     materialSteps,
                                   double                        tX0,
                                   double                        tL0)
  : m_startPosition(startPos)
  , m_theta(theta)
  , m_phi(phi)
  , m_tX0(tX0)
  , m_tL0(tL0)
  , m_materialSteps(materialSteps)
{
}

Acts::MaterialTrack::MaterialTrack(const MaterialTrack& mtrecord)
  : m_startPosition(mtrecord.m_startPosition)
  , m_theta(mtrecord.m_theta)
  , m_phi(mtrecord.m_phi)
  , m_tX0(mtrecord.m_tX0)
  , m_tL0(mtrecord.m_tL0)
  , m_materialSteps(mtrecord.m_materialSteps)
{
}

Acts::MaterialTrack&
Acts::MaterialTrack::operator=(const MaterialTrack& mtrecord)
{
  if (this != &mtrecord) {
    m_startPosition = mtrecord.m_startPosition;
    m_theta         = mtrecord.m_theta;
    m_phi           = mtrecord.m_phi;
    m_materialSteps = mtrecord.m_materialSteps;
    m_tX0           = mtrecord.m_tX0;
    m_tL0           = mtrecord.m_tL0;
  }
  return (*this);
}
