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
  ///
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

  /// Returns the TrackingVolume
  /// @return a plain pointer to the representing TrackingVolume
  const TrackingVolume*
  trackingVolume() const;

  /// Returns the Name
  /// @return string identification of this DetachedTrackingVolume
  const std::string
  name() const;

  /// Moving object around post creation
  ///
  /// @note needs transform to be mutable
  /// @todo check if this needed
  /// @param shift is the applied transform to the detached volume
  void
  move(Transform3D& shift);

  /// Clone with shift
  ///
  /// @param name is the new name of the clone detached volume
  /// @param shift is the applid shift
  DetachedTrackingVolumePtr
  clone(std::string name, Transform3D& shift) const;

  /// Returns a layer representation
  ///
  /// @return pointer to a representation as a layer
  const Layer*
  layerRepresentation() const;

  /// Returns (multi)layer representation
  ///
  /// @return vector to representations as layers
  const LayerVector
  multilayerRepresentation() const;

  /// Sign the volume - the geometry builder has to do that
  ///
  /// @param signat is the volume signature
  /// @param geotype is the volume navigation type
  void
  sign(GeometrySignature signat, GeometryType geotype);

  /// Return the Signature
  /// @return geometry signature
  GeometrySignature
  geometrySignature() const;

  /// Return the Type
  /// @return geometry navigation type
  GeometryType
  geometryType() const;

  /// Set the simplified calculable components
  /// @todo check with Sharka
  ///
  /// @param consts are the constituents to be saved
  void
  saveConstituents(std::vector<std::pair<const Volume*, float>>* consts);

  /// Get the simplified calculable components
  ///
  /// @return the consituents
  std::vector<std::pair<const Acts::Volume*, float>>*
  constituents() const;

  /// Alignment methods: set base transform
  /// default argument to current transform
  ///
  /// @param transf is the relative transform for the alingment
  void
  setBaseTransform(Transform3D* transf = 0);

  /// Alignment methods: realign  / default argument to base transform
  ///
  /// @param transf is the relative transform for the alingment
  void
  realign(Transform3D* transf = 0);

protected:
  /// Default Constructor
  DetachedTrackingVolume();

  /// Constructor with name & layer representation
  ///
  /// @param name is name identifier
  /// @param vol is the contained TrackingVolume
  /// @param layer is the contained layer
  /// @param multiLayer is the multi layer representation
  DetachedTrackingVolume(const std::string&    name,
                         TrackingVolumePtr     vol,
                         LayerPtr              layer,
                         std::vector<LayerPtr> multiLayer);

private:
  const std::string     m_name;
  TrackingVolumePtr     m_trkVolume;
  LayerPtr              m_layerRepresentation;
  std::vector<LayerPtr> m_multilayerRepresentation;
  Transform3D*          m_baseTransform;
  std::vector<std::pair<const Volume*, float>>* m_constituents;
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
    std::vector<std::pair<const Acts::Volume*, float>>* constituents)
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
