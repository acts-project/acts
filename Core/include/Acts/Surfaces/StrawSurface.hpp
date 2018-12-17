// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

/////////////////////////////////////////////////////////////////
// StrawSurface.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Surfaces/LineBounds.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

class DetectorElementBase;
struct PolyhedronRepresentation;

///  @class StrawSurface
///
///  Class for a StrawSurface in the TrackingGeometry
///  to describe dirft tube and straw like detectors.
///
/// @image html LineSurface.png
///
class StrawSurface : public LineSurface
{
  friend Surface;

  /// Default Constructor - deleted
  StrawSurface() = delete;

protected:
  /// Constructor from Transform3D and bounds
  ///
  /// @param htrans is the transform that positions the surface in the global
  /// frame
  /// @param radius is the straw radius
  /// @param halez is the half length in z
  StrawSurface(std::shared_ptr<const Transform3D> htrans,
               double                             radius,
               double                             halez);

  /// Constructor from Transform3D and a shared bounds object
  ///
  /// @param htrans is the transform that positions the surface in the global
  /// frame
  /// @param lbounds are the bounds describing the straw dimensions, can be
  /// optionally nullptr
  StrawSurface(std::shared_ptr<const Transform3D> htrans,
               std::shared_ptr<const LineBounds>  lbounds = nullptr);

  /// Constructor from DetectorElementBase : Element proxy
  ///
  /// @param lbounds are the bounds describing the straw dimensions, they must
  /// not be nullptr
  /// @param detelement for which this surface is (at least) one representation
  StrawSurface(const std::shared_ptr<const LineBounds>& lbounds,
               const DetectorElementBase&               detelement);

  /// Copy constructor
  ///
  /// @param slsf is the source surface for copying
  StrawSurface(const StrawSurface& other);

  /// Copy constructor with shift
  ///
  /// @param other is the source surface dor copying
  /// @param htrans is the additional transform applied after copying
  StrawSurface(const StrawSurface& other, const Transform3D& htrans);

  /// Constructor which accepts @c variant_data
  ///
  /// @param vardata the @c variant_data to build from
  StrawSurface(const variant_data& vardata);

public:
  /// Destructor - defaulted
  ~StrawSurface() override = default;

  /// Assignment operator
  ///
  /// @param other is the source surface for copying
  StrawSurface&
  operator=(const StrawSurface& other);

  /// Clone method. Uses the copy constructor a new position can optionally be
  /// given a shift.
  ///
  /// @param shift additional, optional shift
  std::shared_ptr<StrawSurface>
  clone(const Transform3D* shift = nullptr) const;

  /// Return the surface type
  SurfaceType
  type() const final;

  /// Return properly formatted class name for screen output */
  std::string
  name() const final;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const override;

  /// Return a PolyhedronRepresentation for this object
  /// @param l0div Number of divisions along l0 (phi)
  /// @param l1div Number of divisions along l1 (z)
  virtual PolyhedronRepresentation
  polyhedronRepresentation(size_t l0div = 10, size_t l1div = 1) const;

private:
  /// Clone method. Uses the copy constructor a new position can optionally be
  /// given a shift.
  ///
  /// @param shift additional, optional shift
  StrawSurface*
  clone_impl(const Transform3D* shift = nullptr) const override;
};

inline Surface::SurfaceType
StrawSurface::type() const
{
  return Surface::Straw;
}

inline std::string
Acts::StrawSurface::name() const
{
  return "Acts::StrawSurface";
}

}  // end of namespace Acts
