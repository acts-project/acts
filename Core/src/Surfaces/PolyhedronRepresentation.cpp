// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <sstream>

#include "Acts/Surfaces/PolyhedronRepresentation.hpp"

std::string
Acts::PolyhedronRepresentation::objString(size_t vtxOffset) const
{
  std::stringstream sstr;

  for (const auto& vtx : vertices) {
    sstr << "v " << vtx.x() << " " << vtx.y() << " " << vtx.z() << std::endl;
  }
  for (const auto& face : faces) {
    sstr << "f";
    for (const auto& idx : face) {
      sstr << " " << (vtxOffset + idx + 1);
    }
    sstr << std::endl;
  }
  return sstr.str();
}
