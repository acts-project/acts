// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// Surface.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Surfaces/Surface.hpp"

#include <iomanip>
#include <iostream>
#include <utility>

Acts::Surface::Surface(std::shared_ptr<const Transform3D> tform)
  : GeometryObject(), m_transform(std::move(tform))
{
}

Acts::Surface::Surface(const DetectorElementBase& detelement)
  : GeometryObject(), m_transform(nullptr), m_associatedDetElement(&detelement)
{
}

Acts::Surface::Surface(const Surface& other)
  : GeometryObject(other)
  , std::enable_shared_from_this<Surface>()
  , m_transform(other.m_transform)
  , m_associatedMaterial(other.m_associatedMaterial)
{
}

Acts::Surface::Surface(const Surface& other, const Transform3D& shift)
  : GeometryObject()
  , m_transform(std::make_shared<const Transform3D>(
        Transform3D(shift * other.transform())))
  , m_associatedLayer(other.m_associatedLayer)
  , m_associatedMaterial(other.m_associatedMaterial)
{
}

Acts::Surface::~Surface() = default;

bool
Acts::Surface::isOnSurface(const Vector3D&      gpos,
                           const Vector3D&      gmom,
                           const BoundaryCheck& bcheck) const
{
  // create the local position
  Vector2D lpos;
  // global to local transformation
  bool gtlSuccess = globalToLocal(gpos, gmom, lpos);
  if (gtlSuccess) {
    return bcheck ? bounds().inside(lpos, bcheck) : true;
  }
  // did not succeed
  return false;
}

std::shared_ptr<Acts::Surface>
Acts::Surface::getSharedPtr()
{
  return shared_from_this();
}

std::shared_ptr<const Acts::Surface>
Acts::Surface::getSharedPtr() const
{
  return shared_from_this();
}

Acts::Surface&
Acts::Surface::operator=(const Surface& other)
{
  if (&other != this) {
    GeometryObject::operator=(other);
    // detector element, identifier & layer association are unique
    m_transform          = other.m_transform;
    m_associatedLayer    = other.m_associatedLayer;
    m_associatedMaterial = other.m_associatedMaterial;
    // assigning does invalidate the link to the detectore element
    // we want to have a unique association
    m_associatedDetElement = nullptr;
  }
  return *this;
}

bool
Acts::Surface::operator==(const Surface& other) const
{
  // (a) fast exit for pointer comparison
  if (&other == this) {
    return true;
  }
  // (b) fast exit for type
  if (other.type() != type()) {
    return false;
  }
  // (c) fast exit for bounds
  if (other.bounds() != bounds()) {
    return false;
  }
  // (d) comapre transform
  if (!other.transform().isApprox(transform(), 10e-9)) {
    return false;
  }
  // we should be good
  return true;
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
  Acts::RotationMatrix3D rot(transform().matrix().block<3, 3>(0, 0));
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

bool
Acts::Surface::operator!=(const Acts::Surface& sf) const
{
  return !(operator==(sf));
}
