// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Detector/Portal.hpp"
#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/ActSVG/DetectorSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/LayerSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/PortalSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Io/Svg/SvgPointWriter.hpp"
#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"
#include "ActsPython/Utilities/Context.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include <actsvg/core/draw.hpp>

#include <algorithm>
#include <memory>
#include <ranges>
#include <sstream>
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

// Helper lambda for view range selection
bool viewRangeSel(const Svg::ProtoSurface& s, const Extent& vRange) {
  for (const auto& v : s._vertices) {
    if (vRange.contains(v)) {
      return true;
    }
  }

  return false;
};

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
actsvg::svg::object drawDetectorVolume(const Svg::ProtoVolume& pVolume,
                                       const std::string& identification,
                                       const ViewAndRange& viewAndRange,
                                       PortalCache& portalCache) {
  actsvg::svg::object svgDet;
  svgDet._id = identification;
  svgDet._tag = "g";

  const auto& [view, selection, viewRange] = viewAndRange;

  // Translate selection into booleans
  const bool all = rangeContainsValue(selection, "all");
  const bool sensitives = rangeContainsValue(selection, "sensitives");
  const bool portals = rangeContainsValue(selection, "portals");
  const bool materials = rangeContainsValue(selection, "materials");

  // Helper lambda for material selection
  auto materialSel = [&](const Svg::ProtoSurface& s) -> bool {
    return (materials && s._decorations.contains("material"));
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

    if (rangeContainsValue(portalCache, gpIDs)) {
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
std::vector<actsvg::svg::object> drawDetector(
    const Acts::GeometryContext& gctx,
    const Acts::Experimental::Detector& detector,
    const std::string& identification,
    const std::vector<std::tuple<int, Svg::DetectorVolumeConverter::Options>>&
        volumeIdxOpts,
    const std::vector<ViewAndRange>& viewAndRanges) {
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
      auto svgVolView = drawDetectorVolume(
          pVolume, identification + "_vol" + std::to_string(vidx) + "_" + view,
          var, portalCache);
      svgDetViews[iv].add_object(svgVolView);
    }
  }
  return svgDetViews;
}

}  // namespace

namespace ActsPython {
void addSvgOutput(Context& ctx) {
  auto& mex = ctx.get("examples");

  // Components from the ActsExamples - part of acts.examples
  {
    using Writer = SvgTrackingGeometryWriter;
    auto w = py::class_<Writer, std::shared_ptr<Writer>>(
                 mex, "SvgTrackingGeometryWriter")
                 .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                      py::arg("config"), py::arg("level"))
                 .def("write", py::overload_cast<const AlgorithmContext&,
                                                 const Acts::TrackingGeometry&>(
                                   &Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, outputDir, converterOptions);
  }
  {
    using Writer = SvgPointWriter<SimSpacePoint>;
    auto w =
        py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimSpacePointWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, writerName, trackingGeometry, inputCollection,
                       infoBoxTitle, outputDir);
  }

  {
    using Writer =
        SvgPointWriter<SimHit,
                                     AccessorPositionXYZ>;
    auto w =
        py::class_<Writer, IWriter, std::shared_ptr<Writer>>(
            mex, "SvgSimHitWriter")
            .def(py::init<const Writer::Config&, Acts::Logging::Level>(),
                 py::arg("config"), py::arg("level"))
            .def("write",
                 py::overload_cast<const AlgorithmContext&>(&Writer::write));

    auto c = py::class_<Writer::Config>(w, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, writerName, trackingGeometry, inputCollection,
                       infoBoxTitle, outputDir);
  }
}
}  // namespace ActsPython
