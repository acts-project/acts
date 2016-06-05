// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// VolumeExcluder.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Volumes/VolumeExcluder.hpp"

// Default constructor
Acts::VolumeExcluder::VolumeExcluder() : m_vol(0)
{
}

// constructor with volume
Acts::VolumeExcluder::VolumeExcluder(Acts::Volume* vol) : m_vol(vol)
{
}

// copy constructor
Acts::VolumeExcluder::VolumeExcluder(const VolumeExcluder& ex)
  : Acts::AreaExcluder(ex), m_vol(new Volume(*(ex.m_vol)))
{
}

// destructor
Acts::VolumeExcluder::~VolumeExcluder()
{
  delete m_vol;
}

/** Assignment operator */
Acts::VolumeExcluder&
Acts::VolumeExcluder::operator=(const VolumeExcluder& vol)
{
  if (&vol != this) {
    delete m_vol;
    AreaExcluder::operator=(vol);
    m_vol                 = new Volume(*(vol.m_vol));
  }
  return *this;
}

Acts::VolumeExcluder*
Acts::VolumeExcluder::clone() const
{
  return new Acts::VolumeExcluder(*this);
}
