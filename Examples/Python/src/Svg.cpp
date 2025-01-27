// This file is part of the Acts project.
//
// Copyright (C) 2022-2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/LayerSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/PortalSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Io/Svg/SvgPointWriter.hpp"
#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"

#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace {

// A cache object
using PortalCache = std::list<std::string>;

/// @brief helper tuple to define which views and which range to be used
///
/// @param view is the view type, e.g. 'xy', 'zr', ...
/// @param selection is the selection of the volume, e.g. 'all', 'sensitives', 'portals', 'materials'
using ViewAndRange =
    std::tuple<std::string, std::vector<std::string>, Acts::Extent>;

/// Helper function to be picked in different access patterns
///
/// @param pVolume is the proto volume to be drawn
/// @param identification is the identification of the volume
/// @param viewAndRange is the view, selection and range to be drawn
/// @param portalCache is a portal cache to avoid multiple drawings of the same portal
///
/// Returns an svg object in the right view
actsvg::svg::object viewDetectorVolume(const Svg::ProtoVolume& pVolume,
                                       const std::string& identification,
                                       const ViewAndRange& viewAndRange,
                                       PortalCache& portalCache) {
  actsvg::svg::object svgDet;
  svgDet._id = identification;
  svgDet._tag = "g";

  auto [view, selection, viewRange] = viewAndRange;

  // Translate selection into booleans
  bool all =
      std::find(selection.begin(), selection.end(), "all") != selection.end();
  bool sensitives = std::find(selection.begin(), selection.end(),
                              "sensitives") != selection.end();
  bool portals = std::find(selection.begin(), selection.end(), "portals") !=
                 selection.end();
  bool materials = std::find(selection.begin(), selection.end(), "materials") !=
                   selection.end();

  // Helper lambda for material selection
  auto materialSel = [&](const Svg::ProtoSurface& s) -> bool {
    return (materials &&
            s._decorations.find("material") != s._decorations.end());
  };

  // Helper lambda for view range selection
  auto viewRangeSel = [](const Svg::ProtoSurface& s,
                         const Extent& vRange) -> bool {
    bool accept = false;
    for (const auto& v : s._vertices) {
      if (vRange.contains(v)) {
        accept = true;
        break;
      }
    }
    return accept;
  };

  // -------------------- surface section
  // The surfaces to be drawn
  std::vector<Svg::ProtoVolume::surface_type> sSurfaces;
  sSurfaces.reserve(pVolume._v_surfaces.size());
  for (const auto& s : pVolume._v_surfaces) {
    if ((all || sensitives || materialSel(s)) && viewRangeSel(s, viewRange)) {
      sSurfaces.push_back(s);
    }
  }

  // Now draw all the surfaces
  for (const auto& vs : sSurfaces) {
    if (view == "xy") {
      svgDet.add_object(Svg::View::xy(vs, identification));
    } else if (view == "zr") {
      svgDet.add_object(Svg::View::zr(vs, identification));
    } else {
      throw std::invalid_argument("Unknown view type");
    }
  }

  // -------------------- portal section
  std::vector<Svg::ProtoPortal> sPortals;
  sPortals.reserve(pVolume._portals.size());
  for (const auto& vp : pVolume._portals) {
    if ((all || portals || materialSel(vp._surface)) &&
        viewRangeSel(vp._surface, viewRange)) {
      sPortals.push_back(vp);
    }
  }

  // Now draw all the portals - if not already in the cache
  for (const auto& vp : sPortals) {
    auto pgID = vp._surface._decorations.find("geo_id");
    std::string gpIDs = "";
    if (pgID != vp._surface._decorations.end()) {
      gpIDs = pgID->second._id;
    }

    if (std::find(portalCache.begin(), portalCache.end(), gpIDs) !=
        portalCache.end()) {
      continue;
    }

    // Register this portal to the cache
    portalCache.insert(portalCache.begin(), gpIDs);

    if (view == "xy") {
      svgDet.add_object(Svg::View::xy(vp, identification));
    } else if (view == "zr") {
      svgDet.add_object(Svg::View::zr(vp, identification));
    } else {
      throw std::invalid_argument("Unknown view type");
    }
  }
  return svgDet;
}

// Helper function to be picked in different access patterns
void viewDetector(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector& detector,
    const std::string& identification,
    const std::vector<std::tuple<int, Svg::DetectorVolumeConverter::Options>>&
        volumeIdxOpts,
    const std::vector<ViewAndRange>& viewAndRanges, const std::string& saveAs) {
  PortalCache portalCache;

  // The svg object to be returned
  std::vector<actsvg::svg::object> svgDetViews;
  svgDetViews.reserve(viewAndRanges.size());
  for (unsigned int i = 0; i < viewAndRanges.size(); ++i) {
    actsvg::svg::object svgDet;
    svgDet._id = identification;
    svgDet._tag = "g";
    svgDetViews.push_back(svgDet);
  }

  for (const auto& [vidx, vopts] : volumeIdxOpts) {
    // Get the volume and convert it
    const auto& v = detector.volumes()[vidx];
    auto [pVolume, pGrid] =
        Svg::DetectorVolumeConverter::convert(gctx, *v, vopts);

    for (auto [iv, var] : Acts::enumerate(viewAndRanges)) {
      auto [view, selection, range] = var;
      // Get the view and the range
      auto svgVolView = viewDetectorVolume(
          pVolume, identification + "_vol" + std::to_string(vidx) + "_" + view,
          var, portalCache);
      svgDetViews[iv].add_object(svgVolView);
    }
  }

  for (auto [iv, var] : Acts::enumerate(viewAndRanges)) {
    auto [view, selection, range] = var;
    Svg::toFile({svgDetViews[iv]}, saveAs + "_" + view + ".svg");
  }
}

}  // namespace

namespace Acts::Python {
void addSvg(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  auto svg = m.def_submodule("svg");

  // Some basics
  py::class_<actsvg::svg::object>(svg, "object");

  // Core components, added as an acts.svg submodule
  {
    auto c = py::class_<Svg::Style>(svg, "Style").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::Style);
    ACTS_PYTHON_MEMBER(fillColor);
    ACTS_PYTHON_MEMBER(fillOpacity);
    ACTS_PYTHON_MEMBER(highlightColor);
    ACTS_PYTHON_MEMBER(highlights);
    ACTS_PYTHON_MEMBER(strokeWidth);
    ACTS_PYTHON_MEMBER(strokeColor);
    ACTS_PYTHON_MEMBER(nSegments);
    ACTS_PYTHON_STRUCT_END();
  }

  // How surfaces should be drawn: Svg Surface options & drawning
  {
    auto c = py::class_<Svg::SurfaceConverter::Options>(svg, "SurfaceOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::SurfaceConverter::Options);
    ACTS_PYTHON_MEMBER(style);
    ACTS_PYTHON_MEMBER(templateSurface);
    ACTS_PYTHON_STRUCT_END();

    // Define the proto surface
    py::class_<Svg::ProtoSurface>(svg, "ProtoSurface");
    // Convert an Acts::Surface object into an acts::svg::proto::surface
    svg.def("convertSurface", &Svg::SurfaceConverter::convert);

    // Define the view functions
    svg.def("viewSurface", [](const Svg::ProtoSurface& pSurface,
                              const std::string& identification,
                              const std::string& view = "xy") {
      if (view == "xy") {
        return Svg::View::xy(pSurface, identification);
      } else if (view == "zr") {
        return Svg::View::zr(pSurface, identification);
      } else if (view == "zphi") {
        return Svg::View::zphi(pSurface, identification);
      } else if (view == "zrphi") {
        return Svg::View::zrphi(pSurface, identification);
      } else {
        throw std::invalid_argument("Unknown view type");
      }
    });
  }

  // How portals should be drawn: Svg Portal options & drawning
  {
    auto c = py::class_<Svg::PortalConverter::Options>(svg, "PortalOptions")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::PortalConverter::Options);
    ACTS_PYTHON_MEMBER(surfaceOptions);
    ACTS_PYTHON_MEMBER(linkLength);
    ACTS_PYTHON_MEMBER(volumeIndices);
    ACTS_PYTHON_STRUCT_END();

    // Define the proto portal
    py::class_<Svg::ProtoPortal>(svg, "ProtoPortal");
    // Convert an Acts::Experimental::Portal object into an
    // acts::svg::proto::portal
    svg.def("convertPortal", &Svg::PortalConverter::convert);

    // Define the view functions
    svg.def("viewPortal", [](const Svg::ProtoPortal& pPortal,
                             const std::string& identification,
                             const std::string& view = "xy") {
      if (view == "xy") {
        return Svg::View::xy(pPortal, identification);
      } else if (view == "zr") {
        return Svg::View::zr(pPortal, identification);
      } else {
        throw std::invalid_argument("Unknown view type");
      }
    });
  }

  // How detector volumes are drawn: Svg DetectorVolume options & drawning
  {
    auto c = py::class_<Svg::DetectorVolumeConverter::Options>(
                 svg, "DetectorVolumeOptions")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::DetectorVolumeConverter::Options);
    ACTS_PYTHON_MEMBER(portalIndices);
    ACTS_PYTHON_MEMBER(portalOptions);
    ACTS_PYTHON_MEMBER(surfaceOptions);
    ACTS_PYTHON_STRUCT_END();

    // Define the proto volume & indexed surface grid
    py::class_<Svg::ProtoVolume>(svg, "ProtoVolume");
    py::class_<Svg::ProtoIndexedSurfaceGrid>(svg, "ProtoIndexedSurfaceGrid");

    // Convert an Acts::Experimental::DetectorVolume object into an
    // acts::svg::proto::volume
    svg.def("convertDetectorVolume", &Svg::DetectorVolumeConverter::convert);

    // Define the view functions
    svg.def("viewDetectorVolume", &viewDetectorVolume);
  }

  // How a detector is drawn: Svg Detector options & drawning
  { svg.def("viewDetector", &viewDetector); }

  // Legacy geometry drawing
  {
    using DefinedStyle = std::pair<GeometryIdentifier, Svg::Style>;
    using DefinedStyleSet = std::vector<DefinedStyle>;

    auto sm = py::class_<GeometryHierarchyMap<Svg::Style>>(svg, "StyleMap")
                  .def(py::init<DefinedStyleSet>(), py::arg("elements"));

    auto c = py::class_<Svg::LayerConverter::Options>(svg, "LayerOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::LayerConverter::Options);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(surfaceStyles);
    ACTS_PYTHON_MEMBER(zRange);
    ACTS_PYTHON_MEMBER(phiRange);
    ACTS_PYTHON_MEMBER(gridInfo);
    ACTS_PYTHON_MEMBER(moduleInfo);
    ACTS_PYTHON_MEMBER(projectionInfo);
    ACTS_PYTHON_MEMBER(labelProjection);
    ACTS_PYTHON_MEMBER(labelGauge);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using DefinedLayerOptions =
        std::pair<GeometryIdentifier, Svg::LayerConverter::Options>;
    using DefinedLayerOptionsSet = std::vector<DefinedLayerOptions>;

    auto lom =
        py::class_<GeometryHierarchyMap<Svg::LayerConverter::Options>>(
            svg, "LayerOptionMap")
            .def(py::init<DefinedLayerOptionsSet>(), py::arg("elements"));

    auto c = py::class_<Svg::TrackingGeometryConverter::Options>(
                 svg, "TrackingGeometryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::TrackingGeometryConverter::Options);
    ACTS_PYTHON_MEMBER(prefix);
    ACTS_PYTHON_MEMBER(layerOptions);
    ACTS_PYTHON_STRUCT_END();
  }

  // Components from the ActsExamples - part of acts.examples

  {
    using Writer = ActsExamples::SvgTrackingGeometryWriter;
    auto w = py::class_<Writer, std::shared_ptr<Writer>>(
                 mex, "SvgTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const AlgorithmContext&,
                                                 const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_MEMBER(converterOptions);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Writer = ActsExamples::SvgPointWriter<ActsExamples::SimSpacePoint>;
    auto w =
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimSpacePointWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(writerName);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(inputCollection);
    ACTS_PYTHON_MEMBER(infoBoxTitle);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_STRUCT_END();
  }
  {
    using Writer =
        ActsExamples::SvgPointWriter<ActsExamples::SimHit,
                                     ActsExamples::AccessorPositionXYZ>;
    auto w =
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimHitWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Writer::Config);
    ACTS_PYTHON_MEMBER(writerName);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(inputCollection);
    ACTS_PYTHON_MEMBER(infoBoxTitle);
    ACTS_PYTHON_MEMBER(outputDir);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python
