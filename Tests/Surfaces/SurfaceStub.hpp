// This file is part of the ACTS project.
//
// Copyright (C) 2017-2018 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef ACTS_TEST_SURFACESTUB
#define ACTS_TEST_SURFACESTUB 1

#include "ACTS/Surfaces/InfiniteBounds.hpp"  //to get s_noBounds
#include "ACTS/Surfaces/PlanarBounds.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/VariantData.hpp"

namespace Acts {
/// Surface derived class stub
class SurfaceStub : public Surface
{
public:
  SurfaceStub(std::shared_ptr<const Transform3D> htrans = nullptr)
    : GeometryObject(), Surface(htrans)
  {
  }
  SurfaceStub(const SurfaceStub& sf, const Transform3D& transf)
    : GeometryObject(), Surface(sf, transf)
  {
  }
  SurfaceStub(const DetectorElementBase& detelement,
              const Identifier&          id = Identifier())
    : GeometryObject(), Surface(detelement, id)
  {
  }

  virtual ~SurfaceStub() { /*nop */}

  /// Implicit constructor
  Surface*
  clone(const Transform3D* /*shift = nullptr*/) const final
  {
    return nullptr;
  }

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType
  type() const final
  {
    return Surface::Other;
  }

  /// Return method for the normal vector of the surface
  const Vector3D
  normal(const Vector2D& /*lpos*/) const final
  {
    return Vector3D{0., 0., 0.};
  }

  /// Return method for SurfaceBounds
  const SurfaceBounds&
  bounds() const final
  {
    return s_noBounds;  // need to improve this for meaningful test
  }

  /// Local to global transformation
  void
  localToGlobal(const Vector2D& /*lpos*/,
                const Vector3D& /*gmom*/,
                Vector3D& /*gpos*/) const final
  {
    // nop
  }

  /// Global to local transformation
  bool
  globalToLocal(const Vector3D& /*gpos*/,
                const Vector3D& /*gmom*/,
                Vector2D& lpos) const final
  {
    lpos = Vector2D{20., 20.};
    return true;
  }

  /// Calculation of the path correction for incident
  double
  pathCorrection(const Vector3D& /*gpos*/, const Vector3D& /*gmom*/) const final
  {
    return 0.0;
  }

  /// Straight line intersection schema from parameters
  Intersection
  intersectionEstimate(const Vector3D& /*gpos*/,
                       const Vector3D& /*gdir*/,
                       bool /*forceDir = false*/,
                       const BoundaryCheck& /*bcheck   = false*/) const final
  {
    const Intersection is{Vector3D{1, 1, 1}, 20., true};
    return is;
  }

  /// Inherited from GeometryObject base
  const Vector3D binningPosition(BinningValue /*bValue*/) const final
  {
    const Vector3D v{0.0, 0.0, 0.0};
    return v;
  }

  /// Return properly formatted class name
  std::string
  name() const final
  {
    return std::string("SurfaceStub");
  }

  /// Simply return true to check a method can be called on a constructed object
  bool
  constructedOk() const
  {
    return true;
  }

  virtual variant_data
  toVariantData() const override
  {
    using namespace std::string_literals;
    variant_map data;
    data["type"] = "SurfaceStub"s;
    return data;
  }

private:
  /// the bounds of this surface
  std::shared_ptr<const PlanarBounds> m_bounds;
};
}
#endif
