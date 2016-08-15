// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceBounds.cpp, ACTS project
///////////////////////////////////////////////////////////////////

#include "ACTS/Surfaces/SurfaceBounds.hpp"

Acts::SurfaceBounds&
Acts::SurfaceBounds::operator=(const SurfaceBounds& sb)
{
  if (this != &sb) {
    m_valueStore = sb.m_valueStore;
  }
  return *this;
}

bool
Acts::SurfaceBounds::operator==(const SurfaceBounds& sb) const
{
  /// fast exit for pointer comparison
  if (&sb == this) return true;
  /// fast exit for type comparison
  if (sb.type() != type()) return false;
  // value comparison
  return (sb.m_valueStore == m_valueStore);
}

std::ostream&
Acts::operator<<(std::ostream& sl, const SurfaceBounds& sb)
{
  return sb.dump(sl);
}
