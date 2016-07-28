// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DetExtension.hpp"

Acts::DetExtension::DetExtension()
  : Acts::IDetExtension()
  , m_segmentation(nullptr)
  , m_shape(Acts::ShapeType::None)
  , m_supportMaterial(false)
  , m_bins1(0)
  , m_bins2(0)
  , m_layerMatPos(LayerMaterialPos::inner)
  , m_modules()
{
}

Acts::DetExtension::DetExtension(
    const DD4hep::Geometry::Segmentation segmentation)
  : Acts::IDetExtension()
  , m_segmentation(segmentation)
  , m_shape(Acts::ShapeType::None)
  , m_supportMaterial(false)
  , m_bins1(0)
  , m_bins2(0)
  , m_layerMatPos(LayerMaterialPos::inner)
  , m_modules()
{
}

Acts::DetExtension::DetExtension(ShapeType shape)
  : Acts::IDetExtension()
  , m_segmentation(nullptr)
  , m_shape(shape)
  , m_supportMaterial(false)
  , m_bins1(0)
  , m_bins2(0)
  , m_layerMatPos(LayerMaterialPos::inner)
  , m_modules()
{
}

Acts::DetExtension::DetExtension(size_t                 bins1,
                                 size_t                 bins2,
                                 Acts::LayerMaterialPos layerMatPos)
  : Acts::IDetExtension()
  , m_segmentation(nullptr)
  , m_shape(Acts::ShapeType::None)
  , m_supportMaterial(true)
  , m_bins1(bins1)
  , m_bins2(bins2)
  , m_layerMatPos(layerMatPos)
  , m_modules()
{
}

Acts::DetExtension::DetExtension(std::vector<DD4hep::Geometry::DetElement> mod)
  : Acts::IDetExtension()
  , m_segmentation(nullptr)
  , m_shape(Acts::ShapeType::None)
  , m_supportMaterial(false)
  , m_bins1(0)
  , m_bins2(0)
  , m_layerMatPos(LayerMaterialPos::inner)
  , m_modules(mod)
{
}

Acts::DetExtension::DetExtension(size_t                 bins1,
                                 size_t                 bins2,
                                 Acts::LayerMaterialPos layerMatPos,
                                 std::vector<DD4hep::Geometry::DetElement> mod)
  : Acts::IDetExtension()
  , m_segmentation(nullptr)
  , m_shape(Acts::ShapeType::None)
  , m_supportMaterial(true)
  , m_bins1(bins1)
  , m_bins2(bins2)
  , m_layerMatPos(layerMatPos)
  , m_modules(mod)
{
}

Acts::DetExtension::DetExtension(const DetExtension& det,
                                 const DD4hep::Geometry::DetElement&)
  : Acts::IDetExtension()
  , m_segmentation(det.m_segmentation)
  , m_shape(det.m_shape)
  , m_bins1(det.m_bins1)
  , m_bins2(det.m_bins2)
  , m_layerMatPos(det.m_layerMatPos)
  , m_modules(det.m_modules)
{
}
