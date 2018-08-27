// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Volume.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Volumes/Volume.hpp"
#include <iostream>
#include <utility>
#include "Acts/Volumes/VolumeBounds.hpp"

Acts::Volume::Volume()
  : GeometryObject()
  , m_transform(nullptr)
  , m_center(s_origin)
  , m_volumeBounds(nullptr)
{
}

Acts::Volume::Volume(const std::shared_ptr<const Transform3D>& htrans,
                     std::shared_ptr<const VolumeBounds>       volbounds)
  : GeometryObject()
  , m_transform(htrans)
  , m_center(s_origin)
  , m_volumeBounds(std::move(volbounds))
{
  if (htrans) {
    m_center = htrans->translation();
  }
}

Acts::Volume::Volume(const Volume& vol, const Transform3D* shift)
  : GeometryObject()
  , m_transform(vol.m_transform)
  , m_center(s_origin)
  , m_volumeBounds(vol.m_volumeBounds)
{
  // applyt he shift if it exists
  if (shift != nullptr) {
    m_transform = std::make_shared<const Transform3D>(transform() * (*shift));
  }
  // now set the center
  m_center = transform().translation();
}

Acts::Volume::~Volume() = default;

const Acts::Vector3D
Acts::Volume::binningPosition(Acts::BinningValue bValue) const
{
  // for most of the binning types it is actually the center,
  // just for R-binning types the
  if (bValue == Acts::binR || bValue == Acts::binRPhi) {
    // the binning Position for R-type may have an offset
    return (center() + m_volumeBounds->binningOffset(bValue));
  }
  // return the center
  return center();
}

// assignment operator
Acts::Volume&
Acts::Volume::operator=(const Acts::Volume& vol)
{
  if (this != &vol) {
    m_transform    = vol.m_transform;
    m_center       = vol.m_center;
    m_volumeBounds = vol.m_volumeBounds;
  }
  return *this;
}

Acts::Volume*
Acts::Volume::clone() const
{
  return new Acts::Volume(*this);
}

bool
Acts::Volume::inside(const Acts::Vector3D& gpos, double tol) const
{
  if (!m_transform) {
    return (volumeBounds()).inside(gpos, tol);
  }
  Acts::Vector3D posInVolFrame((transform().inverse()) * gpos);
  return (volumeBounds()).inside(posInVolFrame, tol);
}

std::ostream&
Acts::operator<<(std::ostream& sl, const Acts::Volume& vol)
{
  sl << "Voluem with " << vol.volumeBounds() << std::endl;
  return sl;
}
