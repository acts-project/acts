// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/SurfaceGridSvgConverter.hpp"

#include "Acts/Experimental/DetectorVolume.hpp"
#include "Acts/Experimental/detail/SurfaceGridGenerator.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

std::tuple<std::vector<Acts::Svg::ProtoSurfaces>, Acts::Svg::ProtoGrid,
           std::vector<Acts::Svg::ProtoAssociations> >
Acts::Svg::SurfaceGridConverter::convert(
    const GeometryContext& gctx, const Experimental::DetectorVolume& volume,
    const SurfaceGridConverter::Options& cOptions) {
  // Prepare the return objects
  ProtoSurfaces pSurfaces;
  ProtoGrid pGrid;
  ProtoAssociations pAssociations;

  // The local logger
  ACTS_LOCAL_LOGGER(
      getDefaultLogger("SurfaceGridSvgConverter", cOptions.logLevel));

  const auto& surfaces = volume.surfaces();

  if (not surfaces.empty()) {
    enum ViewType { cylinder, polar, planar, none };
    ViewType vType = none;

    // Get access to the navigation state updator
    auto nStateUpdator = volume.navigationStateUpdator();
    auto nStateUpdatorImpl = nStateUpdator.implementation;
    if (nStateUpdatorImpl != nullptr) {
      // Get the grid parameters
      auto gParameters = Experimental::detail::extractGridParameters(
          *nStateUpdatorImpl, Experimental::detail::s_gridUpdatorTempaltes);

      // Associations
      pAssociations = gParameters.entries;

      if (gParameters.axes.size() == 2u) {
        // Helper method to convert from ACTS to Grid edges
        auto convertGridEdges =
            [](const std::vector<Acts::ActsScalar>& actsEdges)
            -> std::vector<actsvg::scalar> {
          std::vector<actsvg::scalar> svgEdges;
          svgEdges.reserve(actsEdges.size());
          for (const auto ae : actsEdges) {
            svgEdges.push_back(static_cast<actsvg::scalar>(ae));
          }
          return svgEdges;
        };

        // auto binning = surfaceArray.binningValues();
        auto fValue = gParameters.axes[0u].bValue;
        auto sValue = gParameters.axes[1u].bValue;

        if (fValue == binZ and sValue == binPhi) {
          vType = cylinder;
          pGrid._type = actsvg::proto::grid::e_z_phi;
        } else if (fValue == binZ and sValue == binPhi) {
          vType = polar;
          pGrid._type = actsvg::proto::grid::e_r_phi;
        }
        // Assign
        pGrid._edges_0 = convertGridEdges(gParameters.axes[0u].edges);
        pGrid._edges_1 = convertGridEdges(gParameters.axes[1u].edges);
      }

    }  // implementation of navigation state updator exists

    // Find the template surfaces & prepare tempalte objects to be assinged
    std::vector<actsvg::svg::object> templateObjects;
    std::vector<const SurfaceBounds*> templateBounds;

    for (const auto& sf : surfaces) {
      // Get bounds and check them
      const SurfaceBounds& sBounds = sf->bounds();
      // Helper to find bounds
      auto sameBounds = [&](const SurfaceBounds* test) {
        return ((*test) == sBounds);
      };
      // Check if you have this template object already
      auto tBounds = std::find_if(templateBounds.begin(), templateBounds.end(),
                                  sameBounds);
      // New reference bounds and new reference object
      if (tBounds == templateBounds.end()) {
        // Let's get the right style
        SurfaceConverter::Options sOptions;
        sOptions.templateSurface = true;
        // Find a corresponding file in the playbook
        auto sfStyle = cOptions.surfaceStyles.find(sf->geometryId());
        if (sfStyle != cOptions.surfaceStyles.end()) {
          sOptions.style = *sfStyle;
        }

        // Create a referese surface and reference object from it
        auto referenceSurface = SurfaceConverter::convert(gctx, *sf, sOptions);
        auto referenceObject =
            View::xy(referenceSurface,
                     "Template_" + std::to_string(templateObjects.size()));
        templateBounds.push_back(&sBounds);
        templateObjects.push_back(referenceObject);
      }
    }

    ACTS_VERBOSE("Found " << templateObjects.size()
                          << " templates for this layer");
    // Estimate a reference radius
    ActsScalar radius = 0.;

    // Now draw the surfaces from the correct template
    for (const auto& sf : surfaces) {
      radius += Acts::VectorHelpers::perp(sf->center(gctx));

      // Let's get the right style
      SurfaceConverter::Options sOptions;
      // Find a corresponding file in the playbook
      auto sfStyle = cOptions.surfaceStyles.find(sf->geometryId());
      if (sfStyle != cOptions.surfaceStyles.end()) {
        sOptions.style = *sfStyle;
      }

      // Convert the surface from ACTS to actsvg
      auto cSurface = Acts::Svg::SurfaceConverter::convert(gctx, *sf, sOptions);
      cSurface._name = "Module_n_" + std::to_string(pSurfaces.size());

      sOptions.templateSurface = vType != cylinder;

      cSurface._aux_info["grid_info"] = {
          "* module " + std::to_string(pSurfaces.size()) +
          ", surface = " + std::to_string(sf->geometryId().sensitive())};
      // Assign the template for cylinder layers
      if (vType == cylinder) {
        const SurfaceBounds& sBounds = sf->bounds();
        // Helper to find bounds
        auto sameBounds = [&](const SurfaceBounds* test) {
          return ((*test) == sBounds);
        };
        // Check if you have this template object already
        auto tBounds = std::find_if(templateBounds.begin(),
                                    templateBounds.end(), sameBounds);
        // New reference bounds and new reference object
        if (tBounds != templateBounds.end()) {
          size_t tObject = std::distance(templateBounds.begin(), tBounds);
          cSurface._template_object = templateObjects[tObject];
        }
      }
      // Correct view transfrom for disc/planar layers
      if (vType == planar or vType == polar) {
        // Get the transform and estimate the rotation of phi
        // Assumes x/y view
        const auto& sTransform = sf->transform(gctx);
        Vector3 localA = sTransform.rotation().col(0);
        Vector3 localZ = sTransform.rotation().col(2);
        // Find out orientation w.r.t. global transform
        ActsScalar projZ = localZ.dot(Vector3(0., 0., 1.));
        ActsScalar alpha = std::atan2(localA[1], localA[0]) / M_PI * 180.;
        if (projZ < 0.) {
          alpha += 180.;
        }
        auto surfaceCenter = sf->center(gctx);
        // Set the transform for an eventual placement
        cSurface._transform._tr = {
            static_cast<actsvg::scalar>(surfaceCenter[0]),
            static_cast<actsvg::scalar>(surfaceCenter[1])};
        cSurface._transform._rot = {static_cast<actsvg::scalar>(alpha), 0., 0.};
      }
      pSurfaces.push_back(cSurface);
    }
    radius /= surfaces.size();
  }
  // Return the surfaces and the grid
  std::vector<ProtoSurfaces> pSurfaceBatches = {pSurfaces};
  std::vector<ProtoAssociations> pAssociationBatchs = {pAssociations};
  return std::make_tuple(pSurfaceBatches, pGrid, pAssociationBatchs);
}
