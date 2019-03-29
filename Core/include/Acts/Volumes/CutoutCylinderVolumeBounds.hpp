// This file is part of the Acts project.
//
// Copyright (C) 2019 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PolyhedronRepresentation.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Volumes/VolumeBounds.hpp"

namespace Acts {

class IVisualization;

class CutoutCylinderVolumeBounds : public VolumeBounds
{

public:
  CutoutCylinderVolumeBounds(double rmin,
                             double rmed,
                             double rmax,
                             double dz1,
                             double dz2)
    : m_rmin(rmin), m_rmed(rmed), m_rmax(rmax), m_dz1(dz1), m_dz2(dz2)
  {
  }

  virtual ~CutoutCylinderVolumeBounds() = default;

  VolumeBounds*
  clone() const override;

  bool
  inside(const Vector3D& gpos, double tol = 0) const override;

  std::vector<std::shared_ptr<const Surface>>
  decomposeToSurfaces(const Transform3D* transform) const override;
  std::ostream&
  toStream(std::ostream& sl) const override;

  void
  draw(IVisualization&    helper,
       const Transform3D& transform = Transform3D::Identity()) const;

  double
  rMin() const
  {
    return m_rmin;
  }

  double
  rMed() const
  {
    return m_rmed;
  }

  double
  rMax() const
  {
    return m_rmax;
  }

  double
  dZ1() const
  {
    return m_dz1;
  }

  double
  dZ2() const
  {
    return m_dz2;
  }

private:
  double m_rmin;
  double m_rmed;
  double m_rmax;
  double m_dz1;
  double m_dz2;
};

}  // namespace Acts
