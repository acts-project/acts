// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Arrow/ArrowUtil.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using ActsPlugins::ArrowUtil::ArrowSchemaHandle;
namespace ArrowUtil = ActsPlugins::ArrowUtil;

PYBIND11_MODULE(ActsPluginsPythonBindingsArrow, m) {
  // Opaque handle around arrow::Schema. The arrow library is built with
  // hidden visibility inside libActsPluginArrow, so arrow::Schema's
  // typeinfo isn't visible to this .so — pybind would fail to bind it
  // directly. ArrowSchemaHandle is an ACTS_ARROW_EXPORT wrapper whose
  // typeinfo IS visible, so pybind can register it normally. All
  // accessors funnel through methods compiled inside the arrow island.
  //
  // This handle is registered HERE in the plugin-level module so any
  // downstream binding (e.g. acts.examples.arrow.ParquetReader's
  // expectedSchemas) sees the same registered type. Examples consumers
  // that need the schemas import this module first via PLUGIN_IMPORT.
  py::class_<ArrowSchemaHandle>(m, "ArrowSchema")
      .def("__repr__",
           [](const ArrowSchemaHandle& s) {
             return "<ArrowSchema " + s.toString() + ">";
           })
      .def("__str__", &ArrowSchemaHandle::toString)
      .def("field_names", &ArrowSchemaHandle::fieldNames)
      .def("__len__", &ArrowSchemaHandle::numFields);

  auto wrap = [](std::shared_ptr<arrow::Schema> s) {
    return ArrowSchemaHandle{std::move(s)};
  };
  m.def(
      "particleSchema",
      [wrap]() { return wrap(ArrowUtil::particleSchema()); },
      "Schema produced by ArrowParticleOutputConverter.");
  m.def(
      "trackSchema", [wrap]() { return wrap(ArrowUtil::trackSchema()); },
      "Schema produced by ArrowTrackOutputConverter.");
  m.def(
      "simHitSchema", [wrap]() { return wrap(ArrowUtil::simHitSchema()); },
      "Schema produced by ArrowSimHitOutputConverter.");
  m.def(
      "caloHitSchema", [wrap]() { return wrap(ArrowUtil::caloHitSchema()); },
      "Schema produced by ArrowCaloHitOutputConverter.");
}
