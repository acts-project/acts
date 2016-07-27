// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetExtension.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H 1

#include "ACTS/Plugins/DD4hepPlugins/IDetExtension.hpp"
// DD4hep
#include "DD4hep/Detector.h"
#include <vector>

namespace Acts {

/// @class DetExtension
///
/// Implementation of the IDetExtension class, which uses the extension
/// mechanism
/// of DD4hep, needed for the translation from the DD4hep geometry into the
/// tracking geometry of the ACTS package.
/// In this way, the segmentation of the sensitive detector elements can be
/// directly accessed from DD4hep to ensure consistency between the full and the
/// tracking geometry.
/// Since in DD4hep volumes used as a cylinder (detector layers are binned in r
/// and
/// z, e.g. central barrel volume) and discs (detector layers are binned in r
/// and
/// phi, e.g. end caps) are both described as a ROOT TGeoConeSeg one needs to
/// distinguish between these volume types by setting the shape.

class DetExtension : virtual public IDetExtension
{
public:
  /// Constructor
  DetExtension();
  /// Constructor for sensitive detector element
  /// @param segmentation DD4hep segmentation for the readout
  DetExtension(const DD4hep::Geometry::Segmentation segmentation);
  /// Constructor for volume with shape to distinguish between disc and cylinder
  /// volume
  /// @param shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  DetExtension(ShapeType shape);
  /// Constructor for layer with support structure
  /// @note the numer of bins determines the granularity of the material
  /// map of the layer
  /// @param bins1 The number of bins in first direction of the layer
  /// which is phi for both, cylinder and disc layers.
  /// @param bins2 The number of bins in second direction of the layer
  /// which is r in case of a disc layer and z in case of a cylinder layer
  /// @param layerMatPos states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  DetExtension(size_t bins1, size_t bins2, LayerMaterialPos layerMatPos);
  /// Constructor for layer with modules
  /// @param mod Possible sensitive modules contained by a layer
  DetExtension(std::vector<DD4hep::Geometry::DetElement> mod);
  /// Constructor for layer with support structure and modules
  /// @note the numer of bins determines the granularity of the material
  /// map of the layer
  /// @param bins1 The number of bins in first direction of the layer
  /// which is phi for both, cylinder and disc layers.
  /// @param bins2 The number of bins in second direction of the layer
  /// which is r in case of a disc layer and z in case of a cylinder layer
  /// @param layerMatPos states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  /// @param mod Possible sensitive modules contained by a layer
  DetExtension(size_t bins1, size_t bins2, LayerMaterialPos layerMatPos,
                 std::vector<DD4hep::Geometry::DetElement> mod);
  /// Copy constructor
  DetExtension(const DetExtension&, const DD4hep::Geometry::DetElement&);
  /// Virtual destructor
  virtual ~DetExtension() = default;
  /// Possibility to set shape of a volume to distinguish between disc and
  /// cylinder volume
  /// @param shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  void
  setShape(ShapeType shape) override;
  /// Access shape type of a volume to distinguish between disc and cylinder
  /// volume
  /// @return shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  ShapeType
  shape() const override;
  /// Method to set the DD4hep segmentation for the readout
  /// @param segmentation DD4hep segmentation for the readout
  void
  setSegmentation(const DD4hep::Geometry::Segmentation segmentation) override;
  /// Method to access the DD4hep segmentation for the readout
  /// @return segmentation DD4hep segmentation for the readout
  const DD4hep::Geometry::Segmentation
  segmentation() const override;
  /// possibility to mark layer to have support material
  /// @note the numer of bins determines the granularity of the material
  /// map of the layer
  /// @param bins1 The number of bins in first direction of the layer
  /// which is phi for both, cylinder and disc layers.
  /// @param bins2 The number of bins in second direction of the layer
  /// which is r in case of a disc layer and z in case of a cylinder layer
  /// @param layerMatPos states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  virtual void
  supportMaterial(size_t bins1, size_t bins2, LayerMaterialPos layerMatPos) override;
  /// Bool returning true if the layers should carry material
  bool
  hasSupportMaterial() const override;
  /// Access to the two bin numbers determining the granularity of the two dimensional grid
  /// on which the material of the layer should be mapped on
  std::pair<size_t,size_t> materialBins() const override;
  /// returns states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  virtual Acts::LayerMaterialPos layerMaterialPos() const override;
  /// Possibility to set contained detector modules of a layer
  /// @param mod Possible sensitive modules contained by a layer
  void
  setModules(std::vector<DD4hep::Geometry::DetElement> mod) override;
  /// Access modules detector module contained by a layer
  /// @return mod Possible sensitive modules contained by a layer
  virtual std::vector<DD4hep::Geometry::DetElement>
  modules() const final;

private:
  /// Segmentation of a sensitive detector module
  DD4hep::Geometry::Segmentation m_segmentation;
  /// Shape type of a volume defined in IDetExtension can be either disc or
  /// cylinder
  ShapeType m_shape;
  /// Stating if the layer will carry material
  bool m_supportMaterial;
  /// The number of bins in first direction of the layer
  /// which is phi for both, cylinder and disc layers.
  size_t m_bins1;
  /// The number of bins in second direction of the layer
  /// which is r in case of a disc layer and z in case of a cylinder layer
  size_t m_bins2;
  /// States if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  LayerMaterialPos m_layerMatPos;
  /// Possible contained modules of a layer
  std::vector<DD4hep::Geometry::DetElement> m_modules;
  
};
}

inline void
Acts::DetExtension::setShape(Acts::ShapeType type)
{
  m_shape = std::move(type);
}

inline Acts::ShapeType
Acts::DetExtension::shape() const
{
  return m_shape;
}

inline void
Acts::DetExtension::setSegmentation(const DD4hep::Geometry::Segmentation seg)
{
  m_segmentation = std::move(seg);
}

inline const DD4hep::Geometry::Segmentation
Acts::DetExtension::segmentation() const
{
  return m_segmentation;
}

inline void
Acts::DetExtension::supportMaterial(size_t bins1, size_t bins2, LayerMaterialPos layerMatPos)
{
  m_supportMaterial = true;
  m_bins1 = bins1;
  m_bins2 = bins2;
  m_layerMatPos = layerMatPos;
}

inline bool
Acts::DetExtension::hasSupportMaterial() const
{
  return m_supportMaterial;
}

inline Acts::LayerMaterialPos Acts::DetExtension::layerMaterialPos() const
{
    return m_layerMatPos;
}

inline std::pair<size_t,size_t> Acts::DetExtension::materialBins() const
{
    std::pair<size_t,size_t> bins(m_bins1,m_bins2);
    return (bins);
}

inline void
Acts::DetExtension::setModules(std::vector<DD4hep::Geometry::DetElement> mod)
{
  m_modules = std::move(mod);
}

inline std::vector<DD4hep::Geometry::DetElement>
Acts::DetExtension::modules() const
{
  return m_modules;
}



#endif  // ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
