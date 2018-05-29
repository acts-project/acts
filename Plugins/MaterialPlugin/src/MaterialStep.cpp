// This file is part of the Acts project.
//
// Copyright (C) 2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStep.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialPlugins/MaterialStep.hpp"

Acts::MaterialStep::MaterialStep()
  : m_position(Position()), m_material(), m_geoID(0)
{
}

Acts::MaterialStep::MaterialStep(const Acts::MaterialProperties& mat,
                                 const Position&                 pos,
                                 uint64_t                        geoID)
  : m_position(pos), m_material(mat), m_geoID(geoID)
{
}

Acts::MaterialStep::MaterialStep(const MaterialStep& mstep)
  : m_position(mstep.m_position), m_material(mstep.m_material)
{
}

Acts::MaterialStep&
Acts::MaterialStep::operator=(const Acts::MaterialStep& mstep)
{
  if (this != &mstep) {
    m_position = mstep.m_position;
    m_material = mstep.m_material;
    m_geoID    = mstep.m_geoID;
  }
  return (*this);
}
