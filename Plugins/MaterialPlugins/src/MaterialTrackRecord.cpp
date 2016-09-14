// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialTrackRecord.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/MaterialTrackRecord.hpp"

Acts::MaterialTrackRecord::MaterialTrackRecord(
    const MaterialStep::Position& startPos,
    double                        theta,
    double                        phi,
    std::vector<MaterialStep>     materialSteps)
  : m_startPosition(startPos)
  , m_theta(theta)
  , m_phi(phi)
  , m_materialSteps(materialSteps)
{
}

Acts::MaterialTrackRecord::MaterialTrackRecord(
    const MaterialTrackRecord& mtrecord)
  : m_startPosition(mtrecord.m_startPosition)
  , m_theta(mtrecord.m_theta)
  , m_phi(mtrecord.m_phi)
  , m_materialSteps(mtrecord.m_materialSteps)
{
}

const Acts::MaterialTrackRecord*
Acts::MaterialTrackRecord::clone() const
{
  return (new MaterialTrackRecord(*this));
}

Acts::MaterialTrackRecord&
Acts::MaterialTrackRecord::operator=(const MaterialTrackRecord& mtrecord)
{
  if (this != &mtrecord) {
    m_startPosition = mtrecord.m_startPosition;
    m_theta         = mtrecord.m_theta;
    m_phi           = mtrecord.m_phi;
    m_materialSteps = mtrecord.m_materialSteps;
  }
  return (*this);
}
