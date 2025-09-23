// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/AxisDefinitions.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <type_traits>

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace ActsPython {

class PythonLogger {
 public:
  PythonLogger(const std::string& name, Logging::Level level)
      : m_name{name}, m_logger{getDefaultLogger(m_name, level)} {}

  void log(Logging::Level level, const std::string& message) const {
    m_logger->log(level, message);
  }

  void setLevel(Logging::Level level) {
    m_logger = getDefaultLogger(m_name, level);
  }

 private:
  std::string m_name;
  std::unique_ptr<const Logger> m_logger;
};

/// @brief This adds the classes from Core/Utilities to the python module
/// @param m the pybind11 core module
void addUtilities(py::module_& m) {
  { py::class_<AnyBase<512>>(m, "AnyBase512").def(py::init<>()); }

  { py::class_<CalibrationContext>(m, "CalibrationContext").def(py::init<>()); }

  // Add l ogging infrastructure
  {
    auto logging = m.def_submodule("logging", "");

    auto levelEnum = py::enum_<Logging::Level>(logging, "Level")
                         .value("VERBOSE", Logging::VERBOSE)
                         .value("DEBUG", Logging::DEBUG)
                         .value("INFO", Logging::INFO)
                         .value("WARNING", Logging::WARNING)
                         .value("ERROR", Logging::ERROR)
                         .value("FATAL", Logging::FATAL)
                         .value("MAX", Logging::MAX)
                         .export_values();

    levelEnum
        .def("__lt__", [](Logging::Level self,
                          Logging::Level other) { return self < other; })
        .def("__gt__", [](Logging::Level self,
                          Logging::Level other) { return self > other; })
        .def("__le__", [](Logging::Level self,
                          Logging::Level other) { return self <= other; })
        .def("__ge__", [](Logging::Level self, Logging::Level other) {
          return self >= other;
        });

    auto makeLogFunction = [](Logging::Level level) {
      return
          [level](PythonLogger& logger, const std::string& fmt, py::args args) {
            auto locals = py::dict();
            locals["args"] = args;
            locals["fmt"] = fmt;
            py::exec(R"(
        message = fmt % args
    )",
                     py::globals(), locals);

            auto message = locals["message"].cast<std::string>();

            logger.log(level, message);
          };
    };

    py::class_<Logger>(m, "Logger");

    auto logger = py::class_<PythonLogger, std::shared_ptr<PythonLogger>>(
                      logging, "Logger")
                      .def("log", &PythonLogger::log)
                      .def("verbose", makeLogFunction(Logging::VERBOSE))
                      .def("debug", makeLogFunction(Logging::DEBUG))
                      .def("info", makeLogFunction(Logging::INFO))
                      .def("warning", makeLogFunction(Logging::WARNING))
                      .def("error", makeLogFunction(Logging::ERROR))
                      .def("fatal", makeLogFunction(Logging::FATAL))
                      .def("setLevel", &PythonLogger::setLevel);

    static std::unordered_map<std::string, std::shared_ptr<PythonLogger>>
        pythonLoggers = {
            {"root", std::make_shared<PythonLogger>("Python", Logging::INFO)}};

    logging.def(
        "getLogger",
        [](const std::string& name) {
          if (!pythonLoggers.contains(name)) {
            pythonLoggers[name] =
                std::make_shared<PythonLogger>(name, Logging::INFO);
          }
          return pythonLoggers[name];
        },
        py::arg("name") = "root");

    logging.def("setLevel", [](Logging::Level level) {
      pythonLoggers.at("root")->setLevel(level);
    });

    auto makeModuleLogFunction = [](Logging::Level level) {
      return [level](const std::string& fmt, py::args args) {
        auto locals = py::dict();
        locals["args"] = args;
        locals["fmt"] = fmt;
        py::exec(R"(
        message = fmt % args
    )",
                 py::globals(), locals);

        auto message = locals["message"].cast<std::string>();

        pythonLoggers.at("root")->log(level, message);
      };
    };

    logging.def("setFailureThreshold", &Logging::setFailureThreshold);
    logging.def("getFailureThreshold", &Logging::getFailureThreshold);

    struct ScopedFailureThresholdContextManager {
      std::optional<Logging::ScopedFailureThreshold> m_scopedFailureThreshold =
          std::nullopt;
      Logging::Level m_level;

      explicit ScopedFailureThresholdContextManager(Logging::Level level)
          : m_level(level) {}

      void enter() { m_scopedFailureThreshold.emplace(m_level); }

      void exit(const py::object& /*exc_type*/, const py::object& /*exc_value*/,
                const py::object& /*traceback*/) {
        m_scopedFailureThreshold.reset();
      }
    };

    py::class_<ScopedFailureThresholdContextManager>(logging,
                                                     "ScopedFailureThreshold")
        .def(py::init<Logging::Level>(), "level"_a)
        .def("__enter__", &ScopedFailureThresholdContextManager::enter)
        .def("__exit__", &ScopedFailureThresholdContextManager::exit);

    static py::exception<Logging::ThresholdFailure> exc(
        logging, "ThresholdFailure", PyExc_RuntimeError);
    // NOLINTNEXTLINE(performance-unnecessary-value-param)
    py::register_exception_translator([](std::exception_ptr p) {
      try {
        if (p) {
          std::rethrow_exception(p);
        }
      } catch (const std::exception& e) {
        std::string what = e.what();
        if (what.find("ACTS_LOG_FAILURE_THRESHOLD") != std::string::npos) {
          py::set_error(exc, e.what());
        } else {
          std::rethrow_exception(p);
        }
      }
    });

    logging.def("verbose", makeModuleLogFunction(Logging::VERBOSE));
    logging.def("debug", makeModuleLogFunction(Logging::DEBUG));
    logging.def("info", makeModuleLogFunction(Logging::INFO));
    logging.def("warning", makeModuleLogFunction(Logging::WARNING));
    logging.def("error", makeModuleLogFunction(Logging::ERROR));
    logging.def("fatal", makeModuleLogFunction(Logging::FATAL));
  }

  // Add axis related classes
  {
    auto binningValue = py::enum_<AxisDirection>(m, "AxisDirection")
                            .value("AxisX", AxisDirection::AxisX)
                            .value("AxisY", AxisDirection::AxisY)
                            .value("AxisZ", AxisDirection::AxisZ)
                            .value("AxisR", AxisDirection::AxisR)
                            .value("AxisPhi", AxisDirection::AxisPhi)
                            .value("AxisRPhi", AxisDirection::AxisRPhi)
                            .value("AxisTheta", AxisDirection::AxisTheta)
                            .value("AxisEta", AxisDirection::AxisEta)
                            .value("AxisMag", AxisDirection::AxisMag);

    auto boundaryType = py::enum_<AxisBoundaryType>(m, "AxisBoundaryType")
                            .value("Bound", AxisBoundaryType::Bound)
                            .value("Closed", AxisBoundaryType::Closed)
                            .value("Open", AxisBoundaryType::Open);

    auto axisType = py::enum_<AxisType>(m, "AxisType")
                        .value("equidistant", AxisType::Equidistant)
                        .value("variable", AxisType::Variable);
  }

  {
    // Be able to construct a proto binning
    py::class_<ProtoAxis>(m, "ProtoAxis")
        .def(py::init<AxisBoundaryType, const std::vector<double>&>(),
             "bType"_a, "e"_a)
        .def(py::init<AxisBoundaryType, double, double, std::size_t>(),
             "bType"_a, "minE"_a, "maxE"_a, "nbins"_a)
        .def(py::init<AxisBoundaryType, std::size_t>(), "bType"_a, "nbins"_a);

    py::class_<DirectedProtoAxis>(m, "DirectedProtoAxis")
        .def(py::init<AxisDirection, AxisBoundaryType,
                      const std::vector<double>&>(),
             "bValue"_a, "bType"_a, "e"_a)
        .def(py::init<AxisDirection, AxisBoundaryType, double, double,
                      std::size_t>(),
             "bValue"_a, "bType"_a, "minE"_a, "maxE"_a, "nbins"_a)
        .def(py::init<AxisDirection, AxisBoundaryType, std::size_t>(),
             "bValue"_a, "bType"_a, "nbins"_a);
  }
}

}  // namespace ActsPython
