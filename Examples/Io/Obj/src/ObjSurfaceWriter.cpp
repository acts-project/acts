// This file is part of the Acts project.
//
// Copyright (C) 2017 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Plugins/Obj/ObjSurfaceWriter.hpp"

#include <Acts/Geometry/GeometryID.hpp>
#include <Acts/Geometry/Layer.hpp>
#include <Acts/Surfaces/CylinderBounds.hpp>
#include <Acts/Surfaces/PlanarBounds.hpp>
#include <Acts/Surfaces/RadialBounds.hpp>
#include <Acts/Surfaces/SurfaceBounds.hpp>
#include <ios>
#include <iostream>
#include <stdexcept>

FW::Obj::ObjSurfaceWriter::ObjSurfaceWriter(
    const FW::Obj::ObjSurfaceWriter::Config& cfg)
    : m_cfg(cfg) {
  // Validate the configuration
  if (!m_cfg.logger) {
    throw std::invalid_argument("Missing logger");
  } else if (m_cfg.name.empty()) {
    throw std::invalid_argument("Missing algorithm name");
  } else if (!m_cfg.outputStream) {
    throw std::invalid_argument("Missing output stream");
  }

  // Write down the file prefix
  (*(m_cfg.outputStream)) << m_cfg.filePrefix << '\n';
}

std::string FW::Obj::ObjSurfaceWriter::name() const {
  return m_cfg.name;
}

FW::ProcessCode FW::Obj::ObjSurfaceWriter::write(
    const AlgorithmContext& context, const Acts::Surface& surface) {
  std::lock_guard<std::mutex> lock(m_write_mutex);

  ACTS_DEBUG(">>Obj: Writer for Surface object called.");

  auto scalor = m_cfg.outputScalor;
  // let's get the bounds & the transform
  const Acts::SurfaceBounds& surfaceBounds = surface.bounds();
  auto sTransform = surface.transform(context.geoContext);

  // dynamic_cast to PlanarBounds
  const Acts::PlanarBounds* planarBounds =
      dynamic_cast<const Acts::PlanarBounds*>(&surfaceBounds);
  // only continue if the cast worked
  if (planarBounds && m_cfg.outputSensitive) {
    ACTS_VERBOSE(">>Obj: Writing out a PlaneSurface");
    // set the precision - just to be sure
    (*(m_cfg.outputStream)) << '\n';
    (*(m_cfg.outputStream)) << std::setprecision(m_cfg.outputPrecision);
    // get the vertices
    auto planarVertices = planarBounds->vertices();
    // loop over the vertices
    std::vector<Acts::Vector3D> vertices;
    vertices.reserve(planarVertices.size());
    for (auto pv : planarVertices) {
      // get the point in 3D
      Acts::Vector3D v3D(sTransform * Acts::Vector3D(pv.x(), pv.y(), 0.));
      vertices.push_back(v3D);
    }
    // get the thickness and vertical faces
    double thickness = 0.;
    std::vector<unsigned int> vfaces;
    if (surface.associatedDetectorElement() and m_cfg.outputThickness != 0.) {
      // get the thickness form the detector element
      thickness = surface.associatedDetectorElement()->thickness();
      vfaces = {1, 1, 1, 1};
    }
    // output to file
    Obj::writePlanarFace(*(m_cfg.outputStream), m_vtnCounter, scalor, vertices,
                         thickness, vfaces);
    (*(m_cfg.outputStream)) << '\n';
  }

  // check if you have layer and check what your have
  // dynamic cast to CylinderBounds work the same
  const Acts::CylinderBounds* cylinderBounds =
      dynamic_cast<const Acts::CylinderBounds*>(&surfaceBounds);
  if (cylinderBounds && m_cfg.outputLayerSurface) {
    ACTS_VERBOSE(">>Obj: Writing out a CylinderSurface with r = "
                 << cylinderBounds->get(Acts::CylinderBounds::eR));
    // name the object
    auto layerID = surface.geoID().layer();
    (*(m_cfg.outputStream))
        << " o Cylinder_" << std::to_string(layerID) << '\n';
    // output to the file
    Obj::writeTube(*(m_cfg.outputStream), m_vtnCounter, scalor,
                   m_cfg.outputPhiSegemnts, sTransform,
                   cylinderBounds->get(Acts::CylinderBounds::eR),
                   cylinderBounds->get(Acts::CylinderBounds::eHalfLengthZ),
                   m_cfg.outputThickness);
    (*(m_cfg.outputStream)) << '\n';
  }

  ////dynamic cast to RadialBounds or disc bounds work the same
  const Acts::RadialBounds* radialBounds =
      dynamic_cast<const Acts::RadialBounds*>(&surfaceBounds);
  if (radialBounds && m_cfg.outputLayerSurface) {
    ACTS_VERBOSE(">>Obj: Writing out a DiskSurface at z = "
                 << sTransform.translation().z());
    // name the object
    auto layerID = surface.geoID().layer();
    (*(m_cfg.outputStream)) << "o Disk_" << std::to_string(layerID) << '\n';
    // we use the tube writer in the other direction
    double rMin = radialBounds->rMin();
    double rMax = radialBounds->rMax();
    double thickness = rMax - rMin;
    // output to the file
    Obj::writeTube(*(m_cfg.outputStream), m_vtnCounter, scalor,
                   m_cfg.outputPhiSegemnts, sTransform, 0.5 * (rMin + rMax),
                   m_cfg.outputThickness, thickness);
    (*(m_cfg.outputStream)) << '\n';
  }

  // return success
  return FW::ProcessCode::SUCCESS;
}
