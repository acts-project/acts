// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/DetectorVolume.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Plugins/ActSVG/IndexedSurfacesSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/LayerSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SvgUtils.hpp"
#include "Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/GeometryContainers.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/Io/Svg/SvgPointWriter.hpp"
#include "ActsExamples/Io/Svg/SvgTrackingGeometryWriter.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace Acts::Python {
void addSvg(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    auto c = py::class_<Svg::Style>(m, "SvgStyle").def(py::init<>());
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

  {
    auto c = py::class_<Svg::SurfaceConverter::Options>(m, "SvgSurfaceOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::SurfaceConverter::Options);
    ACTS_PYTHON_MEMBER(style);
    ACTS_PYTHON_MEMBER(templateSurface);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using DefinedStyle = std::pair<GeometryIdentifier, Svg::Style>;
    using DefinedStyleSet = std::vector<DefinedStyle>;

    auto sm = py::class_<GeometryHierarchyMap<Svg::Style>>(m, "SvgStyleMap")
                  .def(py::init<DefinedStyleSet>(), py::arg("elements"));

    auto c = py::class_<Svg::LayerConverter::Options>(m, "SvgLayerOptions")
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
            m, "SvgLayerOptionMap")
            .def(py::init<DefinedLayerOptionsSet>(), py::arg("elements"));

    auto c = py::class_<Svg::TrackingGeometryConverter::Options>(
                 m, "SvgTrackingGeometryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Svg::TrackingGeometryConverter::Options);
    ACTS_PYTHON_MEMBER(prefix);
    ACTS_PYTHON_MEMBER(layerOptions);
    ACTS_PYTHON_STRUCT_END();
  }

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

  {
    py::class_<actsvg::svg::object, std::shared_ptr<actsvg::svg::object>>(
        m, "SvgObject");
  }

  {
    m.def("svgIndexedSurfaceGrid",
          [](const GeometryContext& gctx,
             const Acts::Experimental::DetectorVolume& dvolume,
             const Acts::GeometryHierarchyMap<Acts::Svg::Style>& surfaceStyles)
              -> auto {
            if (!dvolume.surfaces().empty()) {
              // Create the draw options
              Acts::Svg::IndexedSurfacesConverter::Options iOptions;
              iOptions.surfaceStyles = surfaceStyles;

              auto protoIndexedSurfaceGrid =
                  Acts::Svg::IndexedSurfacesConverter::convert(
                      gctx, dvolume.surfaces(),
                      dvolume.surfaceCandidatesUpdater(), iOptions);

              return Acts::Svg::View::xy(protoIndexedSurfaceGrid,
                                         dvolume.name() + "_surface_grid");
            }
            return actsvg::svg::object{};
          });
  }

  {
    m.def("svgToFile", [](const std::vector<actsvg::svg::object>& objects,
                          const std::string& fileName) {
      Acts::Svg::toFile(objects, fileName);
    });
  }
}
}  // namespace Acts::Python
