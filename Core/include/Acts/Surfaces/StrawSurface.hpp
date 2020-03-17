// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {

class DetectorElementBase;
struct Polyhedron;

///  @class StrawSurface
///
///  Class for a StrawSurface in the TrackingGeometry
///  to describe dirft tube and straw like detectors.
///
/// @image html LineSurface.png
///
class StrawSurface : public LineSurface {
  friend Surface;

 protected:
  /// Constructor from Transform3D and bounds
  ///
  /// @param htrans is the transform that positions the surface in the global
  /// frame
  /// @param radius is the straw radius
  /// @param halez is the half length in z
  StrawSurface(std::shared_ptr<const Transform3D> htrans, double radius,
               double halez);

  /// Constructor from Transform3D and a shared bounds object
  ///
  /// @param htrans is the transform that positions the surface in the global
  /// frame
  /// @param lbounds are the bounds describing the straw dimensions, can be
  /// optionally nullptr
  StrawSurface(std::shared_ptr<const Transform3D> htrans,
               std::shared_ptr<const LineBounds> lbounds = nullptr);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param lbounds are the bounds describing the straw dimensions, they must
  /// not be nullptr
  /// @param detelement for which this surface is (at least) one representation
  StrawSurface(const std::shared_ptr<const LineBounds>& lbounds,
               const DetectorElementBase& detelement);

  /// Copy constructor
  ///
  /// @param slsf is the source surface for copying
  StrawSurface(const StrawSurface& other);

  /// Copy constructor - with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param other is the source cone surface
  /// @param transf is the additional transfrom applied after copying
  StrawSurface(const GeometryContext& gctx, const StrawSurface& other,
               const Transform3D& transf);

 public:
  /// Destructor - defaulted
  ~StrawSurface() override = default;

  /// Default Constructor - deleted
  StrawSurface() = delete;

  /// Assignment operator
  ///
  /// @param other is the source surface for copying
  StrawSurface& operator=(const StrawSurface& other);

  /// Clone method into a concrete type of StrawSurface with shift
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  std::shared_ptr<StrawSurface> clone(const GeometryContext& gctx,
                                      const Transform3D& shift) const;

  /// Return the surface type
  SurfaceType type() const final;

  /// Return properly formatted class name for screen output */
  std::string name() const final;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg Number of segments along curved lines, it represents
  /// the full 2*M_PI coverange, if lseg is set to 1 only the extrema
  /// are given @note if lseg is set to 1 then only the straw is created
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& gctx,
                                      size_t lseg) const final;

 private:
  /// Clone method implementation
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param shift applied to the surface
  StrawSurface* clone_impl(const GeometryContext& gctx,
                           const Transform3D& shift) const override;
};

inline Surface::SurfaceType StrawSurface::type() const {
  return Surface::Straw;
}

inline std::string Acts::StrawSurface::name() const {
  return "Acts::StrawSurface";
}

}  // namespace Acts
