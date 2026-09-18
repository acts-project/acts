// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Arrow/ArrowUtil.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <memory>
#include <stdexcept>

#include <arrow/c/abi.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using ActsPlugins::ArrowUtil::ArrowSchemaHandle;
using ActsPlugins::ArrowUtil::ArrowTable;
using ActsPython::WhiteBoardRegistry;
namespace ArrowUtil = ActsPlugins::ArrowUtil;

namespace {

// PyCapsule destructor for an ArrowSchema struct exported via the Arrow
// C Data Interface. The capsule owns the heap allocation; if the consumer
// (pyarrow et al.) hasn't already released the struct, we do it now to
// release the producer-side resources.
void releaseArrowSchemaCapsule(PyObject* capsule) {
  auto* c =
      static_cast<ArrowSchema*>(PyCapsule_GetPointer(capsule, "arrow_schema"));
  if (c == nullptr) {
    PyErr_Clear();
    return;
  }
  if (c->release != nullptr) {
    c->release(c);
  }
  delete c;
}

void releaseArrowArrayCapsule(PyObject* capsule) {
  auto* c =
      static_cast<ArrowArray*>(PyCapsule_GetPointer(capsule, "arrow_array"));
  if (c == nullptr) {
    PyErr_Clear();
    return;
  }
  if (c->release != nullptr) {
    c->release(c);
  }
  delete c;
}

}  // namespace

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
      .def("__len__", &ArrowSchemaHandle::numFields)
      // PyCapsule protocol — lets pyarrow / polars / duckdb consume our
      // schema directly, e.g. `pa.schema(acts.arrow.trackSchema())`.
      .def("__arrow_c_schema__", [](const ArrowSchemaHandle& self) {
        auto* c_schema = new ArrowSchema{};
        try {
          self.exportToC(c_schema);
        } catch (...) {
          delete c_schema;
          throw;
        }
        return py::reinterpret_steal<py::object>(
            PyCapsule_New(c_schema, "arrow_schema", releaseArrowSchemaCapsule));
      });

  // ArrowTable: opaque handle around shared_ptr<arrow::Table> that's also
  // the canonical WhiteBoard storage type for arrow tables. Bound with
  // smart_holder so WhiteBoardRegistry can register it for typed access
  // from Python algorithms; downstream `acts.examples.arrow.ParquetReader`
  // / `ParquetWriter` then talk in terms of this same type.
  //
  // Python users primarily get an ArrowTable off the WhiteBoard via a
  // typed `ReadDataHandle`. From there, `pa.record_batch(handle)` /
  // `handle.as_table()` cross into pyarrow via the Arrow C Data
  // Interface (PyCapsule protocol), zero-copy.
  auto arrowTableClass =
      py::class_<ArrowTable, py::smart_holder>(m, "ArrowTable")
          .def(py::init<>())
          .def("__repr__",
               [](const ArrowTable& t) {
                 return "<ArrowTable " + std::to_string(t.numRows()) +
                        " rows x " + std::to_string(t.numColumns()) + " cols>";
               })
          .def("__str__", &ArrowTable::toString)
          .def_property_readonly("num_rows", &ArrowTable::numRows)
          .def_property_readonly("num_columns", &ArrowTable::numColumns)
          .def_property_readonly("schema", &ArrowTable::schema)
          // PyCapsule protocol for the Arrow C Data Interface. The
          // signature matches what pyarrow / polars / duckdb expect:
          // returns (schema_capsule, array_capsule). The
          // requested_schema arg is ignored — we always export our
          // native schema.
          .def(
              "__arrow_c_array__",
              [](const ArrowTable& self,
                 [[maybe_unused]] const py::object& requested_schema) {
                auto* c_schema = new ArrowSchema{};
                auto* c_array = new ArrowArray{};
                try {
                  self.exportToC(c_schema, c_array);
                } catch (...) {
                  delete c_schema;
                  delete c_array;
                  throw;
                }
                py::object schemaCap =
                    py::reinterpret_steal<py::object>(PyCapsule_New(
                        c_schema, "arrow_schema", releaseArrowSchemaCapsule));
                py::object arrayCap =
                    py::reinterpret_steal<py::object>(PyCapsule_New(
                        c_array, "arrow_array", releaseArrowArrayCapsule));
                return py::make_tuple(std::move(schemaCap),
                                      std::move(arrayCap));
              },
              py::arg("requested_schema") = py::none())
          // Convenience: import pyarrow lazily and wrap the single
          // exported batch in a pa.Table. Fails with ImportError if
          // pyarrow isn't installed.
          .def("as_table",
               [](py::object self) {
                 auto pa = py::module_::import("pyarrow");
                 auto batch = pa.attr("record_batch")(self);
                 return pa.attr("Table").attr("from_batches")(
                     py::make_tuple(batch));
               })
          // Inverse of __arrow_c_array__: build an ArrowTable from any
          // object exposing the C Data Interface (pyarrow Table /
          // RecordBatch, polars DataFrame, duckdb relation, ...). For
          // multi-batch streams we read all batches via the stream
          // protocol and concatenate; for single-batch producers we
          // import directly. The producer's buffers are referenced
          // zero-copy via the C-Data release callbacks.
          .def_static(
              "from_arrow",
              [](const py::object& obj) {
                if (py::hasattr(obj, "__arrow_c_array__")) {
                  py::tuple capsules = obj.attr("__arrow_c_array__")();
                  if (capsules.size() != 2) {
                    throw py::type_error(
                        "__arrow_c_array__ returned a tuple of size " +
                        std::to_string(capsules.size()) + ", expected 2");
                  }
                  auto* sc = static_cast<ArrowSchema*>(
                      PyCapsule_GetPointer(capsules[0].ptr(), "arrow_schema"));
                  auto* ar = static_cast<ArrowArray*>(
                      PyCapsule_GetPointer(capsules[1].ptr(), "arrow_array"));
                  if (sc == nullptr || ar == nullptr) {
                    throw py::value_error(
                        "__arrow_c_array__ returned invalid PyCapsules");
                  }
                  // ImportRecordBatch consumes the structs; the
                  // capsule destructors then no-op when the capsules
                  // are GC'd.
                  return ArrowTable::importFromC(sc, ar);
                }
                if (py::hasattr(obj, "__arrow_c_stream__")) {
                  // Materialize the stream by routing through pyarrow:
                  // pa.table(obj) handles __arrow_c_stream__ uniformly,
                  // and the resulting Table exposes __arrow_c_array__
                  // (after combining chunks) which we then import.
                  auto pa = py::module_::import("pyarrow");
                  auto pa_table = pa.attr("table")(obj);
                  auto combined = pa_table.attr("combine_chunks")();
                  // pa.Table doesn't expose __arrow_c_array__ directly;
                  // it does expose __arrow_c_stream__. The cleanest
                  // single-batch path is to convert to a RecordBatch
                  // first via to_batches() (length 1 after combine).
                  py::list batches = combined.attr("to_batches")();
                  if (batches.size() != 1) {
                    throw py::value_error(
                        "expected 1 batch after combine_chunks, got " +
                        std::to_string(batches.size()));
                  }
                  py::tuple capsules = batches[0].attr("__arrow_c_array__")();
                  auto* sc = static_cast<ArrowSchema*>(
                      PyCapsule_GetPointer(capsules[0].ptr(), "arrow_schema"));
                  auto* ar = static_cast<ArrowArray*>(
                      PyCapsule_GetPointer(capsules[1].ptr(), "arrow_array"));
                  if (sc == nullptr || ar == nullptr) {
                    throw py::value_error(
                        "__arrow_c_array__ returned invalid PyCapsules");
                  }
                  return ArrowTable::importFromC(sc, ar);
                }
                throw py::type_error(
                    "ArrowTable.from_arrow: object does not implement the "
                    "Arrow C Data Interface (__arrow_c_array__ or "
                    "__arrow_c_stream__)");
              },
              py::arg("obj"),
              "Build an ArrowTable from any object implementing the Arrow "
              "C Data Interface (pyarrow Table/RecordBatch, polars "
              "DataFrame, duckdb relation, etc.). Zero-copy via the "
              "release-callback wiring; the producer's buffers stay alive "
              "until the resulting ArrowTable is destroyed.");

  WhiteBoardRegistry::registerClass(arrowTableClass);

  auto wrap = [](std::shared_ptr<arrow::Schema> s) {
    return ArrowSchemaHandle{std::move(s)};
  };
  m.def(
      "particleSchema", [wrap]() { return wrap(ArrowUtil::particleSchema()); },
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
