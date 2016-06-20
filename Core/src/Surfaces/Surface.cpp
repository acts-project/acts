// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Surface.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/Surface.hpp"
#include <iomanip>
#include <iostream>

Acts::Surface::Surface(std::shared_ptr<Acts::Transform3D> tform)
  : m_transform(tform)
  , m_associatedDetElement(nullptr)
  , m_associatedDetElementId()
  , m_associatedLayer(nullptr)
  , m_associatedTrackingVolume(nullptr)
  , m_associatedMaterial(nullptr)
{}

Acts::Surface::Surface(const Acts::DetectorElementBase& detelement,
                       const Identifier&                id)
  : m_transform(nullptr)
  , m_associatedDetElement(&detelement)
  , m_associatedDetElementId(id)
  , m_associatedLayer(nullptr)
  , m_associatedTrackingVolume(nullptr)
  , m_associatedMaterial(nullptr)
{}

Acts::Surface::Surface(const Surface& sf)
  : m_transform(sf.m_transform)
  , m_associatedDetElement(nullptr)
  , m_associatedDetElementId()
  , m_associatedLayer(nullptr)
  , m_associatedTrackingVolume(nullptr)
  , m_associatedMaterial(sf.m_associatedMaterial)
{}

Acts::Surface::Surface(const Surface& sf, const Acts::Transform3D& shift)
  : m_transform(std::make_shared<Acts::Transform3D>(
        Acts::Transform3D(shift * sf.transform())))
  , m_associatedDetElement(nullptr)
  , m_associatedDetElementId()
  , m_associatedLayer(sf.m_associatedLayer)
  , m_associatedMaterial(sf.m_associatedMaterial)
{}

Acts::Surface::~Surface()
{}

Acts::Surface&
Acts::Surface::operator=(const Surface& sf)
{
  if (&sf != this){
    // detector element, identifier & layer association are unique
    m_transform               = m_transform;
    m_associatedDetElement    = nullptr;
    m_associatedDetElementId  = Identifier();
    m_associatedLayer         = nullptr;
    m_associatedMaterial      = sf.m_associatedMaterial;
  }
  return *this; 
}
  

bool
Acts::Surface::operator==(const Surface& sf) const
{
  // (a) fast exit for pointer comparison
  if (&sf == this) return true;
  // (b) fast exit for type
  if (sf.type() != type()) return false;
  // (c) fast exit for bounds
  if (sf.bounds() != bounds()) return false;
  // (d) comapre transform
  if (!sf.transform().isApprox(transform(), 10e-9)) return false;
  // we should be good
  return true; 
}
  
bool
Acts::Surface::isOnSurface(const Acts::Vector3D& gpos,
                           const BoundaryCheck&  bchk) const
{
  // create the local position
  Acts::Vector2D lpos;
  // global to local transformation
  bool g2L = globalToLocal(gpos, Acts::Vector3D::UnitX(), lpos);
  if (g2L) {
    // no boundary check, then return true
    if (!bchk) return true;
    // return what ever the bounds tell you
    return bounds().inside(lpos, bchk);
  }
  // did not succeed
  return false;
}

const Acts::RotationMatrix3D
Acts::Surface::measurementFrame(const Acts::Vector3D&,
                                const Acts::Vector3D&) const
{
  return transform().rotation();
}

// overload dump for stream operator
std::ostream&
Acts::Surface::dump(std::ostream& sl) const
{
  sl << std::setiosflags(std::ios::fixed);
  sl << std::setprecision(4);
  sl << name() << std::endl;
  sl << "     Center position  (x, y, z) = (" << center().x() << ", "
     << center().y() << ", " << center().z() << ")" << std::endl;
  Acts::RotationMatrix3D rot(transform().rotation());
  Acts::Vector3D         rotX(rot.col(0));
  Acts::Vector3D         rotY(rot.col(1));
  Acts::Vector3D         rotZ(rot.col(2));
  sl << std::setprecision(6);
  sl << "     Rotation:             colX = (" << rotX(0) << ", " << rotX(1)
     << ", " << rotX(2) << ")" << std::endl;
  sl << "                           colY = (" << rotY(0) << ", " << rotY(1)
     << ", " << rotY(2) << ")" << std::endl;
  sl << "                           colZ = (" << rotZ(0) << ", " << rotZ(1)
     << ", " << rotZ(2) << ")" << std::endl;
  sl << "     Bounds  : " << bounds();
  sl << std::setprecision(-1);
  return sl;
}

/**Overload of << operator for std::ostream for debug output*/
std::ostream&
Acts::operator<<(std::ostream& sl, const Acts::Surface& sf)
{
  return sf.dump(sl);
}
