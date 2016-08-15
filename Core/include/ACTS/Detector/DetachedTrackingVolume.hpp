// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetachedTrackingVolume.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DETECTOR_DETACHEDTRACKINGVOLUME_H
#define ACTS_DETECTOR_DETACHEDTRACKINGVOLUME_H 1

#include <memory>
#include "ACTS/Layers/Layer.hpp"
#include "ACTS/Layers/PlaneLayer.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/GeometrySignature.hpp"

namespace Acts {

class TrackingVolume;
class Surface;

typedef std::vector<LayerPtr> LayerVector;

// master typedefs
class DetachedTrackingVolume;
typedef std::shared_ptr<const DetachedTrackingVolume> DetachedTrackingVolumePtr;
typedef std::shared_ptr<const TrackingVolume>         TrackingVolumePtr;

/// @class DetachedTrackingVolume
///
/// Base Class for a navigation object (active/passive) in the Tracking
/// geometry.

class DetachedTrackingVolume
{
  friend class TrackingVolume;

public:
  /// Factory Constructor
  /// @param name is name identifier
  /// @param vol is the contained TrackingVolume
  /// @param layer is the contained layer
  /// @param multiLayer is the multi layer representation
  static DetachedTrackingVolumePtr
  create(const std::string& name,
         TrackingVolumePtr  vol,
         LayerPtr           layer      = nullptr,
         LayerVector        multiLayer = {})
  {
    return DetachedTrackingVolumePtr(
        new DetachedTrackingVolume(name, vol, layer, multiLayer));
  }

  /// Destructor
  ~DetachedTrackingVolume();

  /// returns the TrackingVolume
  const TrackingVolume*
  trackingVolume() const;

  /// returns the Name
  const std::string
  name() const;

  /// moving object around
  void
  move(Transform3D& shift) const;

  /// clone with shift
  DetachedTrackingVolumePtr
  clone(std::string name, Transform3D& shift) const;

  /// returns layer representation
  const Layer*
  layerRepresentation() const;

  /// returns (multi)layer representation
  const LayerVector
  multilayerRepresentation() const;

  /// sign the volume - the geometry builder has to do that
  void
  sign(GeometrySignature signat, GeometryType geotype) const;

  /// return the Signature
  GeometrySignature
  geometrySignature() const;

  /// return the Type
  GeometryType
  geometryType() const;

  /// set the simplified calculable components
  void
  saveConstituents(std::vector<std::pair<const Acts::Volume*, float>>*) const;

  /// get the simplified calculable components
  std::vector<std::pair<const Acts::Volume*, float>>*
  constituents() const;

  ///  alignment methods: set base transform
  /// default argument to current transform
  void
  setBaseTransform(Transform3D* transf = 0) const;

  /// alignment methods: realign  / default argument to base transform
  void
  realign(Transform3D* transf = 0) const;

protected:
  /// Default Constructor
  DetachedTrackingVolume();

  /// Constructor with name & layer representation
  DetachedTrackingVolume(const std::string&    name,
                         TrackingVolumePtr     vol,
                         LayerPtr              layer,
                         std::vector<LayerPtr> multilayer);

private:
  const std::string     m_name;
  TrackingVolumePtr     m_trkVolume;
  LayerPtr              m_layerRepresentation;
  std::vector<LayerPtr> m_multilayerRepresentation;
  mutable Transform3D*  m_baseTransform;
  mutable std::vector<std::pair<const Volume*, float>>* m_constituents;
};

inline const TrackingVolume*
DetachedTrackingVolume::trackingVolume() const
{
  return m_trkVolume.get();
}

inline const std::string
DetachedTrackingVolume::name() const
{
  return (m_name);
}

inline const Layer*
DetachedTrackingVolume::layerRepresentation() const
{
  return m_layerRepresentation.get();
}

inline const LayerVector
DetachedTrackingVolume::multilayerRepresentation() const
{
  return m_multilayerRepresentation;
}

inline void
DetachedTrackingVolume::saveConstituents(
    std::vector<std::pair<const Acts::Volume*, float>>* constituents) const
{
  m_constituents = constituents;
}

inline std::vector<std::pair<const Acts::Volume*, float>>*
DetachedTrackingVolume::constituents() const
{
  return m_constituents;
}

}  // end of namespace

#endif  // ACTS_DETECTOR_DETACHEDTRACKINGVOLUME_H
