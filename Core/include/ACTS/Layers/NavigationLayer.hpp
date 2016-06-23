// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// NavigationLayer.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_LAYERS_NAVIGATIONLAYER_H
#define ACTS_LAYERS_NAVIGATIONLAYER_H

#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Definitions.hpp"

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
  /// @param sRepresenation is the represenation for extrapolation
  /// @param thickness is the thickness for the binning 
  static LayerPtr
  create(std::unique_ptr<const Surface> sRepresentation, double thickness = 0.)
  {
    return LayerPtr(new NavigationLayer(std::move(sRepresentation), thickness));
  }

  /// Clone with a shift 
  /// @param shift is the addition transform applied after cloning
  LayerPtr
  cloneWithShift(const Transform3D& shift) const override;

  /// Destructor
  virtual ~NavigationLayer();

  /// The binning position method 
  ///  - as default the center is given, but may be overloaded
  /// @return The return vector can be used for binning in a TrackingVolume
  virtual const Vector3D
  binningPosition(BinningValue bValue) const override;

  /// Default Constructor - deleted
  NavigationLayer() = delete;
  
  /// Copy Constructor - deleted 
  NavigationLayer(const NavigationLayer&) = delete;

  /// Assignment operator - deleted
  NavigationLayer&
  operator=(const NavigationLayer&) = delete;

  /// Transforms the layer into a Surface representation for extrapolation 
  /// In general, extrapolation to a surface should be avoided
  const Surface&
  surfaceRepresentation() const override;

  /// Geometric isOnLayer() method
  /// using isOnSurface() with Layer specific tolerance
  /// @param gpos is the global position for the check
  bool
  isOnLayer(const Vector3D&      gpos,
            const BoundaryCheck& bcheck = true) const override;

  ///  Boolean check method if layer has material:
  /// - checks if any of the layer surfaces has material:
  /// - can be approach surfaces or layer surface 
  bool
  hasMaterial() const override;

  ///  Boolean check method if layer has sensitive surfaces
  bool
  hasSensitive() const override;

protected:
  /// Private Constructor 
  /// - this is called by the creat(args*) method
  /// passed spacer layer if needed  */
  NavigationLayer(std::unique_ptr<const Surface> surfaceRepresentation, 
                  double thickness);
  
  //// for the navigation Volume the surface
  std::unique_ptr<const Surface>   m_surfaceRepresentation;
};

inline const Surface&
NavigationLayer::surfaceRepresentation() const
{
  return (*m_surfaceRepresentation);
}

inline const Vector3D
NavigationLayer::binningPosition(BinningValue bValue) const
{
  return m_surfaceRepresentation->binningPosition(bValue);
}

inline bool
NavigationLayer::hasMaterial() const
{
  return false;
}

inline bool
NavigationLayer::hasSensitive() const
{
  return false;
}

}  // end of namespace

#endif  // ACTS_LAYERS_NAVIGATIONLAYER_H
