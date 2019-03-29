// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class IVisualization
{
public:
  using color_type = std::array<int, 3>;

  virtual void
  vertex(const Vector3D& vtx, color_type color = {120, 120, 120})
      = 0;

  virtual void
  face(const std::vector<Vector3D>& vtxs, color_type color = {120, 120, 120})
      = 0;

  virtual void
  line(const Vector3D& a, const Vector3D& b, color_type color = {120, 120, 120})
      = 0;

  virtual void
  write(std::ostream& os) const = 0;

  virtual void
  clear()
      = 0;

  void
  vertex(const Vector3F& vtx, color_type color = {120, 120, 120})
  {
    Vector3D vtxd = vtx.template cast<double>();
    vertex(vtxd, color);
  }

  void
  face(const std::vector<Vector3F>& vtxs, color_type color = {120, 120, 120})
  {
    std::vector<Vector3D> vtxsd;
    std::transform(vtxs.begin(),
                   vtxs.end(),
                   std::back_inserter(vtxsd),
                   [](auto& v) { return v.template cast<double>(); });
    face(vtxsd, color);
  }

  void
  line(const Vector3F& a, const Vector3F& b, color_type color = {120, 120, 120})
  {
    Vector3D ad = a.template cast<double>();
    Vector3D bd = b.template cast<double>();
    line(ad, bd, color);
  }
};

inline std::ostream&
operator<<(std::ostream& os, const IVisualization& hlp)
{
  hlp.write(os);
  return os;
}
}
