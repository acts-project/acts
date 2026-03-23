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
#include "Acts/Utilities/Histogram.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/RangeXD.hpp"

#include <cmath>
#include <random>
#include <sstream>
#include <type_traits>
#include <unordered_map>

#include <boost/histogram.hpp>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace ActsPython {

namespace bh = boost::histogram;
using namespace Acts::Experimental;

/// Copy inner bin values (no flow) from a boost::histogram into a numpy array
/// in C (row-major) order. ValueFn extracts a double from the indexed element.
template <typename Hist, typename ValueFn>
static py::array_t<double> copyBins(const Hist& h, ValueFn fn) {
  const auto n = h.rank();
  std::vector<py::ssize_t> shape(n);
  for (std::size_t i = 0; i < n; ++i) {
    shape[i] = static_cast<py::ssize_t>(h.axis(i).size());
  }
  py::array_t<double> result(shape);
  double* ptr = result.mutable_data();
  for (auto&& x : bh::indexed(h, bh::coverage::inner)) {
    py::ssize_t flat = 0;
    py::ssize_t stride = 1;
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
      flat += static_cast<py::ssize_t>(x.index(i)) * stride;
      stride *= shape[i];
    }
    ptr[flat] = fn(x);
  }
  return result;
}

template <std::size_t Dim>
static void bindHistogram(py::module_& m) {
  using H = Histogram<Dim>;
  py::classh<H>(m, ("Histogram" + std::to_string(Dim)).c_str())
      .def_property_readonly("name", &H::name)
      .def_property_readonly("title", &H::title)
      .def_property_readonly("rank", [](const H&) { return Dim; })
      .def_property_readonly("histogram", &H::histogram,
                             py::return_value_policy::reference_internal);
}

template <std::size_t Dim>
static void bindEfficiency(py::module_& m) {
  using E = Efficiency<Dim>;
  py::classh<E>(m, ("Efficiency" + std::to_string(Dim)).c_str())
      .def_property_readonly("name", &E::name)
      .def_property_readonly("title", &E::title)
      .def_property_readonly("rank", [](const E&) { return Dim; })
      .def_property_readonly("accepted", &E::acceptedHistogram,
                             py::return_value_policy::reference_internal)
      .def_property_readonly("total", &E::totalHistogram,
                             py::return_value_policy::reference_internal);
}

/// @brief This adds the classes from Core/Utilities to the python module
/// @param m the pybind11 core module
void addUtilities(py::module_& m) {
  {
    py::class_<AnyBase<512>>(m, "AnyBase512").def(py::init<>());
  }

  {
    py::class_<CalibrationContext>(m, "CalibrationContext").def(py::init<>());
  }

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

    {
      auto makeLogFunction = [](Logging::Level level) {
        return
            [level](const Logger& logger, const py::str& fmt, py::args args) {
              if (!logger.doPrint(level)) {
                return;
              }

              logger.log(level, fmt.format(*args));
            };
      };

      auto logger =
          py::class_<Logger, py::smart_holder>(m, "Logger")
              .def("log", &Logger::log)
              .def_property_readonly("level", &Logger::level)
              .def_property_readonly("name", &Logger::name)
              .def("verbose", makeLogFunction(Logging::VERBOSE))
              .def("debug", makeLogFunction(Logging::DEBUG))
              .def("info", makeLogFunction(Logging::INFO))
              .def("warning", makeLogFunction(Logging::WARNING))
              .def("error", makeLogFunction(Logging::ERROR))
              .def("fatal", makeLogFunction(Logging::FATAL))
              .def("clone",
                   py::overload_cast<const std::optional<std::string>&,
                                     const std::optional<Logging::Level>&>(
                       &Logger::clone, py::const_),
                   py::arg("name") = py::none(), py::arg("level") = py::none())
              .def(
                  "clone",
                  py::overload_cast<Logging::Level>(&Logger::clone, py::const_),
                  py::arg("level"))
              .def("cloneWithSuffix", &Logger::cloneWithSuffix,
                   py::arg("suffix"), py::arg("level") = py::none());
    }

    static std::unordered_map<std::string, std::unique_ptr<const Logger>>
        loggerCache;

    m.def(
        "getDefaultLogger",
        [](const std::string& name, Logging::Level level) {
          return getDefaultLogger(name, level);
        },
        py::arg("name"), py::arg("level") = Logging::INFO,
        py::return_value_policy::take_ownership);

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

    logging.def("_consumeLoggerFunction",
                [](std::unique_ptr<const Logger> logger) {
                  logger->log(Logging::VERBOSE, "consumed logger logs");
                });

    struct ConfigWithLogger {
      std::unique_ptr<const Logger> logger;
    };

    py::class_<ConfigWithLogger>(logging, "_ConfigWithLogger")
        .def(py::init<std::unique_ptr<const Logger>>(), "logger"_a)
        .def_property_readonly(
            "logger", [](ConfigWithLogger& self) { return self.logger.get(); });
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

  {
    using RangeXDDim3 = RangeXD<3u, double>;

    py::class_<RangeXDDim3>(m, "RangeXDDim3")
        .def(py::init([](const std::array<double, 2u>& range0,
                         const std::array<double, 2u>& range1,
                         const std::array<double, 2u>& range2) {
          RangeXDDim3 range;
          range[0].shrink(range0[0], range0[1]);
          range[1].shrink(range1[0], range1[1]);
          range[2].shrink(range2[0], range2[1]);
          return range;
        }));
  }

  // Histogram bindings
  {
    // Single axis type covering regular, variable, and log axes
    py::class_<AxisVariant>(m, "Axis")
        .def_property_readonly("size", &AxisVariant::size)
        .def_property_readonly(
            "label",
            [](const AxisVariant& ax) -> std::string { return ax.metadata(); })
        .def_property_readonly(
            "edges", [](const AxisVariant& ax) { return extractBinEdges(ax); });

    // Core dense histogram (BoostHist — values as copied numpy array)
    py::class_<BoostHist>(m, "BoostHistogram")
        .def_property_readonly("rank", &BoostHist::rank)
        .def(
            "axis",
            [](const BoostHist& h, std::size_t i) -> AxisVariant {
              return h.axis(i);
            },
            "i"_a)
        .def("values", [](const BoostHist& h) {
          return copyBins(h, [](auto& x) { return static_cast<double>(*x); });
        });

    // Profile histogram (BoostProfileHist — means/variances as numpy arrays)
    py::class_<BoostProfileHist>(m, "BoostProfileHistogram")
        .def_property_readonly("rank", &BoostProfileHist::rank)
        .def(
            "axis",
            [](const BoostProfileHist& h, std::size_t i) -> AxisVariant {
              return h.axis(i);
            },
            "i"_a)
        .def("counts",
             [](const BoostProfileHist& h) {
               return copyBins(h, [](auto& x) {
                 return static_cast<double>((*x).count());
               });
             })
        .def("means",
             [](const BoostProfileHist& h) {
               return copyBins(h, [](auto& x) { return (*x).value(); });
             })
        .def("sum_of_deltas_squared", [](const BoostProfileHist& h) {
          return copyBins(h, [](auto& x) {
            const auto n = (*x).count();
            return n > 1.0 ? (*x).variance() * (n - 1.0) : 0.0;
          });
        });

    bindHistogram<1>(m);
    bindHistogram<2>(m);
    bindHistogram<3>(m);

    {
      using P = ProfileHistogram<1>;
      py::classh<P>(m, "ProfileHistogram1")
          .def_property_readonly("name", &P::name)
          .def_property_readonly("title", &P::title)
          .def_property_readonly("rank", [](const P&) { return 1u; })
          .def_property_readonly("sampleAxisTitle", &P::sampleAxisTitle)
          .def_property_readonly("histogram", &P::histogram,
                                 py::return_value_policy::reference_internal);
    }

    bindEfficiency<1>(m);
    bindEfficiency<2>(m);

    // Demo factory functions for testing and examples — not public API
    m.def("_demo_histogram1", []() -> Histogram1 {
      BoostRegularAxis ax(10, 0.0, 10.0, "x [a.u.]");
      Histogram1 h("demo", "Demo Histogram", {AxisVariant(ax)});
      std::mt19937 rng(42);
      std::uniform_real_distribution<double> x(0.0, 10.0);
      for (int i = 0; i < 1000; ++i) {
        h.fill({x(rng)});
      }
      return h;
    });

    m.def("_demo_profile1", []() -> ProfileHistogram1 {
      BoostRegularAxis ax(10, 0.0, 10.0, "x [a.u.]");
      ProfileHistogram1 h("demo", "Demo Profile", {AxisVariant(ax)},
                          "y [a.u.]");
      std::mt19937 rng(42);
      std::uniform_real_distribution<double> x(0.0, 10.0);
      std::normal_distribution<double> noise(0.0, 0.1);
      for (int i = 0; i < 1000; ++i) {
        double xi = x(rng);
        h.fill({xi}, std::sin(xi) + noise(rng));
      }
      return h;
    });

    m.def("_demo_efficiency1", []() -> Efficiency1 {
      BoostRegularAxis ax(10, 0.0, 10.0, "x [a.u.]");
      Efficiency1 h("demo", "Demo Efficiency", {AxisVariant(ax)});
      std::mt19937 rng(42);
      std::uniform_real_distribution<double> x(0.0, 10.0);
      for (int i = 0; i < 1000; ++i) {
        double xi = x(rng);
        std::bernoulli_distribution accepted(1.0 - std::exp(-xi));
        h.fill({xi}, accepted(rng));
      }
      return h;
    });
  }
}

}  // namespace ActsPython
