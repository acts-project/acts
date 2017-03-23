// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// MaterialStep.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Plugins/MaterialPlugins/MaterialStep.hpp"

Acts::MaterialStep::MaterialStep() : m_position(Position()), m_material()
{
}

Acts::MaterialStep::MaterialStep(const Acts::MaterialProperties& mat,
                                 const Position&                 pos)
  : m_position(pos), m_material(mat)
{
}

Acts::MaterialStep::MaterialStep(const MaterialStep& mstep)
  : m_position(mstep.m_position), m_material(mstep.m_material)
{
}

Acts::MaterialStep*
Acts::MaterialStep::clone() const
{
  return (new MaterialStep(*this));
}

Acts::MaterialStep&
Acts::MaterialStep::operator=(const Acts::MaterialStep& mstep)
{
  if (this != &mstep) {
    m_position = mstep.m_position;
    m_material = mstep.m_material;
  }
  return (*this);
}
