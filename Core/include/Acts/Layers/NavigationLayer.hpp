// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// NavigationLayer.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include "Acts/Layers/Layer.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/GeometryStatics.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"

namespace Acts {

class Surface;
class BinUtility;

/// @class NavigationLayer
///
/// Class to be used for gaps in Volumes as a navigational link.
/// Navigation Layers have a surface representation, but should usually never be
/// propagated to.

class NavigationLayer : public Layer
{
public:
  ///  Factory Constructor - the surface representation is given by pointer
  /// (ownership passed)
  ///
  /// @param sRepresentation is the representation for extrapolation
  /// @param thickness is the thickness for the binning
  static LayerPtr
  create(std::shared_ptr<const Surface> sRepresentation, double thickness = 0.)
  {
    return LayerPtr(new NavigationLayer(std::move(sRepresentation), thickness));
  }

  /// Factory for shared Layer pointer, that accepts @c variant_data
  /// @param vardata The data to build from
  static LayerPtr
  create(const variant_data& vardata);

  /// Destructor
  ~NavigationLayer() override;

  /// The binning position method
  ///  - as default the center is given, but may be overloaded
  /// @return The return vector can be used for binning in a TrackingVolume
  const Vector3D
  binningPosition(BinningValue bValue) const final;

  /// Default Constructor - deleted
  NavigationLayer() = delete;

  /// Copy Constructor - deleted
  NavigationLayer(const NavigationLayer&) = delete;

  /// Assignment operator - deleted
  NavigationLayer&
  operator=(const NavigationLayer&)
      = delete;

  /// Transforms the layer into a Surface representation for extrapolation
  /// In general, extrapolation to a surface should be avoided
  const Surface&
  surfaceRepresentation() const final;

  // Non-const version
  Surface&
  surfaceRepresentation() final;

  /// Geometric isOnLayer() method
  /// using isOnSurface() with Layer specific tolerance
  ///
  /// @param gp is the global position for the check
  /// @param bcheck is the boundary check directive
  ///
  /// @return boolean that indicates if the position is on surface
  bool
  isOnLayer(const Vector3D& gp, const BoundaryCheck& bcheck = true) const final;

  /// Accept layer according to the following colelction directives
  ///
  /// @param resolveSensitive is the prescription to find the sensitive surfaces
  /// @param resolveMaterial is the precription to find material surfaces
  /// @param resolvePassive is the prescription to find all passive surfaces
  ///
  /// @note navigation layers are never accepted
  ///
  /// @return a boolean whether the layer is accepted for processing
  bool
  resolve(bool resolveSensitive,
          bool resolveMaterial,
          bool resolvePassive) const final;

  /// Produce a @c variant_data representation of this object
  /// @return The representation
  variant_data
  toVariantData() const;

protected:
  /// Private Constructor
  /// - this is called by the creat(args*) method
  /// passed spacer layer if needed
  ///
  /// @param surfaceRepresentation is the surface of the layer
  /// @param thickness ithe layer thickness
  NavigationLayer(std::shared_ptr<const Surface> surfaceRepresentation,
                  double                         thickness);

  /// for the navigation Volume the surface
  ///
  /// We will need to mutate this surface during the geometry building process,
  /// but the C++ type system has no const-correct way of expressing this.
  ///
  std::shared_ptr<const Surface> m_surfaceRepresentation;
};

inline const Surface&
NavigationLayer::surfaceRepresentation() const
{
  return (*m_surfaceRepresentation);
}

inline Surface&
NavigationLayer::surfaceRepresentation()
{
  return *(const_cast<Surface*>(m_surfaceRepresentation.get()));
}

inline const Vector3D
NavigationLayer::binningPosition(BinningValue bValue) const
{
  return m_surfaceRepresentation->binningPosition(bValue);
}

inline bool
NavigationLayer::isOnLayer(const Vector3D&      gp,
                           const BoundaryCheck& bcheck) const
{
  return m_surfaceRepresentation->isOnSurface(gp, s_origin, bcheck);
}

inline bool
NavigationLayer::resolve(bool /*resolveSensitive*/,
                         bool /*resolveMaterial*/,
                         bool /*reolvePassive*/) const
{
  return false;
}

}  // namespace Acts
