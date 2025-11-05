// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.


#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ActsPlugins/FpeMonitoring/FpeMonitor.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

namespace py = pybind11;
using namespace py::literals;

using namespace Acts;
using namespace ActsPlugins;
using namespace ActsPython;

namespace {

    void trigger_divbyzero() {
  volatile float j = 0.0;
  volatile float r = 123 / j;  // MARK: divbyzero
  (void)r;
}

void trigger_overflow() {
  volatile float j = std::numeric_limits<float>::max();
  volatile float r = j * j;  // MARK: overflow
  (void)r;
}

void trigger_invalid() {
  volatile float j = -1;
  volatile float r = std::sqrt(j);  // MARK: invalid
  (void)r;
}

}

PYBIND11_MODULE(ActsPluginsPythonBindingsFpe, fpe) {

  struct FpeMonitorContext {
    std::optional<FpeMonitor> mon;
  };

  auto fpeMon = py::class_<FpeMonitor>(fpe, "FpeMonitor")
                 .def_static("_trigger_divbyzero", &trigger_divbyzero)
                 .def_static("_trigger_overflow", &trigger_overflow)
                 .def_static("_trigger_invalid", &trigger_invalid)
                 .def_static("context", []() { return FpeMonitorContext(); });

  fpeMon.def_property_readonly("result", py::overload_cast<>(&FpeMonitor::result),
                            py::return_value_policy::reference_internal)
      .def("rearm", &FpeMonitor::rearm);

  py::class_<FpeMonitor::Result>(fpeMon, "Result")
      .def("merged", &FpeMonitor::Result::merged)
      .def("merge", &FpeMonitor::Result::merge)
      .def("count", &FpeMonitor::Result::count)
      .def("__str__", [](const FpeMonitor::Result& result) {
        std::stringstream os;
        result.summary(os);
        return os.str();
      });

  py::class_<FpeMonitorContext>(fpe, "_FpeMonitorContext")
      .def(py::init([]() { return std::make_unique<FpeMonitorContext>(); }))
      .def(
          "__enter__",
          [](FpeMonitorContext& fm) -> FpeMonitor& {
            fm.mon.emplace();
            return fm.mon.value();
          },
          py::return_value_policy::reference_internal)
      .def("__exit__", [](FpeMonitorContext& fm, py::object /*exc_type*/,
                          py::object /*exc_value*/,
                          py::object /*traceback*/) { fm.mon.reset(); });

  py::enum_<FpeType>(fpe, "FpeType")
      .value("INTDIV", FpeType::INTDIV)
      .value("INTOVF", FpeType::INTOVF)
      .value("FLTDIV", FpeType::FLTDIV)
      .value("FLTOVF", FpeType::FLTOVF)
      .value("FLTUND", FpeType::FLTUND)
      .value("FLTRES", FpeType::FLTRES)
      .value("FLTINV", FpeType::FLTINV)
      .value("FLTSUB", FpeType::FLTSUB)

      .def_property_readonly_static(
          "values", [](py::object /*self*/) -> const auto& {
            static const std::vector<FpeType> values = {
                FpeType::INTDIV, FpeType::INTOVF, FpeType::FLTDIV,
                FpeType::FLTOVF, FpeType::FLTUND, FpeType::FLTRES,
                FpeType::FLTINV, FpeType::FLTSUB};
            return values;
          });

  py::register_exception<FpeFailure>(fpe, "FpeFailure", PyExc_RuntimeError);

}