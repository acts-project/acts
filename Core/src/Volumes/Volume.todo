// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Volume.cpp, ACTS project
///////////////////////////////////////////////////////////////////

// Geometry module
#include "ACTS/Volumes/Volume.hpp"
#include "ACTS/Volumes/VolumeBounds.hpp"
// STD/STL
#include <iostream>

// Default constructor
Acts::Volume::Volume()
  : GeometryObject()
  , m_transform(nullptr)
  , m_center(nullptr)
  , m_volumeBounds(nullptr)
{
}

// constructor with Transform3D
Acts::Volume::Volume(Acts::Transform3D*        htrans,
                     const Acts::VolumeBounds* volbounds)
  : GeometryObject()
  , m_transform(std::shared_ptr<Acts::Transform3D>(htrans))
  , m_center(nullptr)
  , m_volumeBounds(std::shared_ptr<const Acts::VolumeBounds>(volbounds))
{
}

// constructor with shared arguments Transform3D
Acts::Volume::Volume(std::shared_ptr<Acts::Transform3D>        htrans,
                     std::shared_ptr<const Acts::VolumeBounds> volbounds)
  : GeometryObject()
  , m_transform(htrans)
  , m_center(nullptr)
  , m_volumeBounds(volbounds)
{
}

// copy constructor - will up to now not copy the sub structure!
Acts::Volume::Volume(const Acts::Volume& vol)
  : GeometryObject()
  , m_transform(vol.m_transform)
  , m_center(nullptr)
  , m_volumeBounds(vol.m_volumeBounds)
{
}

// copy constructor with shift
Acts::Volume::Volume(const Acts::Volume& vol, const Acts::Transform3D& shift)
  : GeometryObject()
  , m_transform(std::shared_ptr<Acts::Transform3D>(
        new Acts::Transform3D(shift * vol.transform())))
  , m_center(nullptr)
  , m_volumeBounds(vol.m_volumeBounds)
{
}

// destructor
Acts::Volume::~Volume()
{
  delete m_center;
}

// The binning position method
Acts::Vector3D
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
    delete m_center;
    m_transform    = vol.m_transform;
    m_center       = nullptr;
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
Acts::Volume::inside(const Acts::Vector3D& gp, double tol) const
{
  if (!m_transform) return (volumeBounds()).inside(gp, tol);
  Acts::Vector3D posInVolFrame((transform().inverse()) * gp);
  return (volumeBounds()).inside(posInVolFrame, tol);
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream&
Acts::operator<<(std::ostream& sl, const Acts::Volume& vol)
{
  sl << "Voluem with " << vol.volumeBounds() << std::endl;
  return sl;
}
