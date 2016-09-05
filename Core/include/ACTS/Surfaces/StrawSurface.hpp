// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// StrawSurface.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_STRAWSURFACE_H
#define ACTS_SURFACES_STRAWSURFACE_H

#include "ACTS/Surfaces/LineBounds.hpp"
#include "ACTS/Surfaces/LineSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Identifier.hpp"

namespace Acts {

class DetectorElementBase;

///  @class StrawSurface
///
///  Class for a StrawSurface in the TrackingGeometry
///  to describe dirft tube and straw like detectors.

class StrawSurface : public LineSurface
{
public:
  /// Default Constructor - deleted
  StrawSurface() = delete;

  /// Constructor from Transform3D and bounds
  /// @param htrans is the transform that positions the surface in the global
  /// frame
  /// @param radius is the straw radius
  /// @param halex is the half length in z
  StrawSurface(std::shared_ptr<Transform3D> htrans, double radius, double halez)
    : LineSurface(htrans, radius, halez)
  {
  }

  /// Constructor from Transform3D and a shared bounds object
  /// @param htrans is the transform that positions the surface in the global
  /// frame
  /// @param lbounds are teh bounds describing the straw dimensions, can be
  /// optionally nullptr
  StrawSurface(std::shared_ptr<Transform3D>      htrans,
               std::shared_ptr<const LineBounds> lbounds = nullptr);

  /// Constructor from DetectorElementBase and Element identifier
  /// @param lbounds are teh bounds describing the straw dimensions, they must
  /// not be nullptr
  /// @param detelement for which this surface is (at least) one representation
  /// @param identifier
  StrawSurface(std::shared_ptr<const LineBounds> lbounds,
               const DetectorElementBase&        detelement,
               const Identifier&                 identifier = Identifier());

  /// Copy constructor
  /// @param slsf is teh source surface for copying
  StrawSurface(const StrawSurface& slsf) : LineSurface(slsf) {}
  /// Copy constructor with shift
  /// @param slsf is the source surface dor copying
  /// @param htrans is the additional transform applied after copying
  StrawSurface(const StrawSurface& slsf, const Transform3D& htrans)
    : LineSurface(slsf, htrans)
  {}

  /// Destructor
  virtual ~StrawSurface();

  /// Assignment operator
  StrawSurface&
  operator=(const StrawSurface& slsf);

  /// Implicit constructor - shift can be provided */
  virtual StrawSurface*
  clone(const Transform3D* shift = nullptr) const override;

  /// Return the surface type
  virtual SurfaceType
  type() const override
  {
    return Surface::Straw;
  }

  /// Return properly formatted class name for screen output */
  virtual std::string
  name() const override
  {
    return "Acts::StrawSurface";
  };
};

inline StrawSurface*
StrawSurface::clone(const Transform3D* shift) const
{
  if (shift) new StrawSurface(*this, *shift);
  return new StrawSurface(*this);
}

}  // end of namespace

#endif  // ACTS_SURFACESSTRAIGHTLINESURFACE_H
