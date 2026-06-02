// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/ApproachDescriptor.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GeometryObject.hpp"
#include "Acts/Geometry/Volume.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Propagator/NavigationTarget.hpp"
#include "Acts/Surfaces/BoundaryTolerance.hpp"
#include "Acts/Utilities/BinnedArray.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
#include <utility>

namespace Acts {

class Surface;
class SurfaceArray;
class ISurfaceMaterial;
class BinUtility;
class Volume;
class VolumeBounds;
class TrackingVolume;
class ApproachDescriptor;
class IMaterialDecorator;
template <typename object_t>
struct NavigationOptions;

class Layer;

/// @brief Type alias for a shared pointer to a layer
using LayerPtr = std::shared_ptr<const Layer>;
/// @brief Type alias for a mutable pointer to a layer
/// @details Used for non-const access to layer objects in the geometry
using MutableLayerPtr = std::shared_ptr<Layer>;
/// @brief Type alias for adjacent layer pointers
/// @details Stores pointers to the next inner and outer layers in the detector
using NextLayers = std::pair<const Layer*, const Layer*>;

/// @enum LayerType
///
/// For code readability, it distinguishes between different
/// type of layers, which steers the behaviour in the navigation
enum LayerType { navigation = -1, passive = 0, active = 1 };

/// @class Layer
///
/// Base Class for a Detector Layer in the Tracking Geometry
///
/// An actual implemented Detector Layer inheriting from this base
/// class has to inherit from a specific type of Surface as well.
/// In addition, a Layer can carry:
///
/// A SurfaceArray of Surfaces holding the actual detector elements or
/// subSurfaces.
/// A pointer to the TrackingVolume (can only be set by such)
/// An active/passive code :
/// 0      - active
/// 1      - passive
/// [....] - other
///
/// The search type for compatible surfaces on a layer is
///   [ the higher the number, the faster ]:
/// --------- Layer internal ------------------------------------------------
/// -1     - provide all intersection tested without boundary check
///  0     - provide all intersection tested with boundary check
///  1     - provide overlap descriptor based without boundary check
///  2     - provide overlap descriptor based with boundary check
///
/// A layer can have substructure regarding:
/// - m_ssRepresentingSurface -> always exists (1), can have material (2)
/// - m_ssSensitiveSurfaces   -> can or not exist (0,1), can have material (2)
/// - m_ssApproachSurfaces    -> can or not exist (0,1) cam have material (2)
///
class Layer : public virtual GeometryObject {
  /// Declare the TrackingVolume as a friend, to be able to register previous,
  /// next and set the enclosing TrackingVolume
  friend class TrackingVolume;
  friend class Gen1GeometryClosureVisitor;

 public:
  /// Destructor
  ~Layer() noexcept override;

  /// Return the entire SurfaceArray, returns a nullptr if no SurfaceArray
  /// @return Pointer to the surface array, or nullptr if not set
  const SurfaceArray* surfaceArray() const;

  /// Non-const version
  /// @return Mutable pointer to the surface array
  SurfaceArray* surfaceArray();

  /// Transforms the layer into a Surface representation for extrapolation
  /// @note the layer can be hosting many surfaces, but this is the global
  /// one to which one can extrapolate
  /// @return Reference to the layer's surface representation
  virtual const Surface& surfaceRepresentation() const = 0;

  /// Non-const version of surface representation access
  /// @return Mutable reference to the layer surface
  virtual Surface& surfaceRepresentation() = 0;

  /// Return the Thickness of the Layer
  /// this is by definition along the normal vector of the surfaceRepresentation
  /// @return The layer thickness value
  double layerThickness() const;

  /// geometrical isOnLayer() method
  ///
  /// @note using isOnSurface() with Layer specific tolerance
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the global position to be checked
  /// @param boundaryTolerance is the boundary check directive
  ///
  /// @return boolean that indicates success of the operation
  virtual bool isOnLayer(const GeometryContext& gctx, const Vector3& position,
                         const BoundaryTolerance& boundaryTolerance =
                             BoundaryTolerance::None()) const;

  /// Return method for the approach descriptor, can be nullptr
  /// @return Pointer to the approach descriptor, or nullptr if not set
  const ApproachDescriptor* approachDescriptor() const;

  /// Non-const version of the approach descriptor
  /// @return Mutable pointer to the approach descriptor
  ApproachDescriptor* approachDescriptor();

  /// Accept layer according to the following collection directives
  ///
  /// @tparam options_t Type of the options for navigation
  /// @param options Navigation options containing resolution settings
  ///
  /// @return a boolean whether the layer is accepted for processing
  template <typename options_t>
  bool resolve(const options_t& options) const {
    return resolve(options.resolveSensitive, options.resolveMaterial,
                   options.resolvePassive);
  }

  /// Accept layer according to the following collection directives
  ///
  /// @param resolveSensitive is the prescription to find the sensitive surfaces
  /// @param resolveMaterial is the precription to find material surfaces
  /// @param resolvePassive is the prescription to find all passive surfaces
  ///
  /// @return a boolean whether the layer is accepted for processing
  virtual bool resolve(bool resolveSensitive, bool resolveMaterial,
                       bool resolvePassive) const;

  /// @brief Decompose Layer into (compatible) surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Position parameter for searching
  /// @param direction Direction of the parameters for searching
  /// @param options The navigation options
  ///
  /// @return list of intersection of surfaces on the layer
  boost::container::small_vector<NavigationTarget, 10> compatibleSurfaces(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction,
      const NavigationOptions<Surface>& options) const;

  /// Surface seen on approach
  /// for layers without sub structure, this is the surfaceRepresentation
  /// for layers with sub structure, this is the approachSurface
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position Position for searching
  /// @param direction Direction for searching
  /// @param options The  navigation options
  ///
  /// @return the Surface intersection of the approach surface
  NavigationTarget surfaceOnApproach(
      const GeometryContext& gctx, const Vector3& position,
      const Vector3& direction, const NavigationOptions<Layer>& options) const;

  /// Fast navigation to next layer
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param position is the start position for the search
  /// @param direction is the direction for the search
  ///
  /// @return the pointer to the next layer
  const Layer* nextLayer(const GeometryContext& gctx, const Vector3& position,
                         const Vector3& direction) const;

  /// Get the confining TrackingVolume
  ///
  /// @return the pointer to the enclosing volume
  const TrackingVolume* trackingVolume() const;

  ///  Return the abstract volume that represents the layer
  ///
  /// @return the representing volume of the layer
  const Volume* representingVolume() const;

  /// return the LayerType
  /// @return The layer type (active, passive, or navigation)
  LayerType layerType() const;

 protected:
  /// Constructor with pointer to SurfaceArray (passing ownership)
  ///
  /// @param surfaceArray is the array of sensitive surfaces
  /// @param thickness is the normal thickness of the Layer
  /// @param ades oapproach descriptor
  /// @param laytyp is the layer type if active or passive
  explicit Layer(std::unique_ptr<SurfaceArray> surfaceArray,
                 double thickness = 0.,
                 std::unique_ptr<ApproachDescriptor> ades = nullptr,
                 LayerType laytyp = passive);

  ///  private method to set enclosing TrackingVolume, called by friend class
  /// only
  ///  optionally, the layer can be resized to the dimensions of the
  /// TrackingVolume
  ///  - Bounds of the Surface are resized
  ///  - MaterialSlab dimensions are resized
  ///  - SubSurface array boundaries are NOT resized
  ///
  /// @param tvol is the tracking volume the layer is confined
  void encloseTrackingVolume(const TrackingVolume& tvol);

  /// the previous Layer according to BinGenUtils
  NextLayers m_nextLayers;

  /// A binutility to find the next layer
  /// @todo check if this is needed
  const BinUtility* m_nextLayerUtility = nullptr;

  /// SurfaceArray on this layer Surface
  ///
  /// This array will be modified during signature and constant afterwards, but
  /// the C++ type system unfortunately cannot cleanly express this.
  ///
  std::unique_ptr<const SurfaceArray> m_surfaceArray;

  /// Thickness of the Layer
  double m_layerThickness = 0;

  /// descriptor for surface on approach
  ///
  /// The descriptor may need to be modified during geometry building, and will
  /// remain constant afterwards, but again C++ cannot currently express this.
  ///
  std::unique_ptr<const ApproachDescriptor> m_approachDescriptor;

  /// the enclosing TrackingVolume
  const TrackingVolume* m_trackingVolume = nullptr;

  /// Representing Volume
  /// can be used as approach surface sources
  std::unique_ptr<Volume> m_representingVolume;

  /// make a passive/active either way
  LayerType m_layerType;

  /// sub structure indication
  /// Substructure flag indicating representing surface configuration
  int m_ssRepresentingSurface = 0;
  /// Substructure flag indicating sensitive surface configuration
  int m_ssSensitiveSurfaces = 0;
  /// Substructure flag indicating approach surface configuration
  int m_ssApproachSurfaces = 0;

 private:
  /// Private helper method to close the geometry
  /// - it will assign material to the surfaces if needed
  /// - it will set the layer geometry ID for a unique identification
  /// - it will also register the internal sub structure
  ///
  /// @param materialDecorator is a decorator that assigns
  ///        optionally the surface material to where they belong
  /// @param layerID is the geometry id of the volume
  ///                as calculated by the TrackingGeometry
  /// @param hook Identifier hook to be applied to surfaces
  /// @param logger A @c Logger instance
  ///
  void closeGeometry(const IMaterialDecorator* materialDecorator,
                     const GeometryIdentifier& layerID,
                     const GeometryIdentifierHook& hook,
                     const Logger& logger = getDummyLogger());
};

/// Layers are constructed with shared_ptr factories, hence the layer array is
/// describes as:
using LayerArray = BinnedArray<LayerPtr>;

}  // namespace Acts

#include "Acts/Geometry/Layer.ipp"
