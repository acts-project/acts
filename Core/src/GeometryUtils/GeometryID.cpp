// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// GeometryIDCalculator.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Utilities/GeometryID.hpp"

bool
Acts::operator<(const Acts::GeometryID& one, const Acts::GeometryID& two)
{
  return (one.value() < two.value());
}

bool
Acts::operator<=(const Acts::GeometryID& one, const Acts::GeometryID& two)
{
  return (one.value() <= two.value());
}

bool
Acts::operator>(const Acts::GeometryID& one, const Acts::GeometryID& two)
{
  return (one.value() > two.value());
}

bool
Acts::operator>=(const Acts::GeometryID& one, const Acts::GeometryID& two)
{
  return (one.value() >= two.value());
}

std::ostream&
Acts::operator<<(std::ostream& sl, const Acts::GeometryID& tid)
{
  sl << tid.value();
  return sl;
}
