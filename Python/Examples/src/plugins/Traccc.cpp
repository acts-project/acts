// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Traccc/ActsMeasToTracccAlg.hpp"
#include "ActsExamples/Traccc/ActsSpToTracccAlg.hpp"
#include "ActsExamples/Traccc/TracccChain.hpp"
#include "ActsExamples/Traccc/TracccSeqAlg.hpp"
#include "ActsPython/Utilities/Helpers.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsTraccc, traccc) {
  // ── TracccSeqAlg ─────────────────────────────────────────────────────────
  {
    using Alg = Traccc::TracccSeqAlg;
    using Cfg = Alg::Config;
    using ChainCfg = Traccc::TracccChainConfig;

    auto cls =
        py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
            traccc, "TracccSeqAlgorithm")
            .def(py::init([](const Cfg& cfg, Acts::Logging::Level level) {
                   return std::make_shared<Alg>(cfg, level);
                 }),
                 "cfg"_a, "level"_a = Acts::Logging::INFO);

    py::enum_<ChainCfg::Backend>(traccc, "TracccBackend")
        .value("CPU", ChainCfg::Backend::CPU)
        .value("CUDA", ChainCfg::Backend::CUDA);

    cls.attr("Backend") = traccc.attr("TracccBackend");

    py::class_<Cfg>(cls, "Config")
        .def(py::init<>())
        // --- chain backend choice ---
        .def_property(
            "backend", [](const Cfg& c) { return c.chain.backend; },
            [](Cfg& c, ChainCfg::Backend v) { c.chain.backend = v; })
        // --- chain inputs ---
        .def_property(
            "detectorFile", [](const Cfg& c) { return c.chain.detectorFile; },
            [](Cfg& c, const std::string& v) { c.chain.detectorFile = v; })
        .def_property(
            "digitizationFile",
            [](const Cfg& c) { return c.chain.digitizationFile; },
            [](Cfg& c, const std::string& v) { c.chain.digitizationFile = v; })
        .def_property(
            "conditionsFile",
            [](const Cfg& c) { return c.chain.conditionsFile; },
            [](Cfg& c, const std::string& v) { c.chain.conditionsFile = v; })
        .def_property(
            "materialFile", [](const Cfg& c) { return c.chain.materialFile; },
            [](Cfg& c, const std::string& v) { c.chain.materialFile = v; })
        .def_property(
            "gridFile", [](const Cfg& c) { return c.chain.gridFile; },
            [](Cfg& c, const std::string& v) { c.chain.gridFile = v; })
        .def_property(
            "bfieldFile",
            [](const Cfg& c) { return c.chain.magneticFieldFile; },
            [](Cfg& c, const std::string& v) { c.chain.magneticFieldFile = v; })
        // --- direct input ---
        .def_readwrite("inputMeasurements", &Cfg::inputMeasurements)
        .def_readwrite("inputSpacepoints", &Cfg::inputSpacepoints);
  }

  // ── ActsMeasToTracccAlg ───────────────────────────────────────────────────
  {
    using Alg = ActsExamples::ActsMeasToTracccAlg;
    using Cfg = Alg::Config;

    auto cls =
        py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
            traccc, "ActsMeasToTracccAlg")
            .def(py::init([](const Cfg& cfg, Acts::Logging::Level level) {
              return std::make_shared<Alg>(
                  cfg, Acts::getDefaultLogger("ActsMeasToTracccAlg", level));
            }));

    py::class_<Cfg>(cls, "Config")
        .def(py::init<>())
        .def_readwrite("detectorFile", &Cfg::detectorFile)
        .def_readwrite("inputActsMeasurements", &Cfg::inputActsMeasurements)
        .def_readwrite("outputDetrayToActsMap", &Cfg::outputDetrayToActsMap)
        .def_readwrite("outputTracccMeasurements",
                       &Cfg::outputTracccMeasurements)
        .def_readwrite("trackingGeometry", &Cfg::trackingGeometry);
  }

  // ── ActsSpToTracccAlg ─────────────────────────────────────────────────────
  {
    using Alg = ActsExamples::ActsSpToTracccAlg;
    using Cfg = Alg::Config;

    auto cls =
        py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
            traccc, "ActsSpToTracccAlg")
            .def(py::init([](const Cfg& cfg, Acts::Logging::Level level) {
              return std::make_shared<Alg>(
                  cfg, Acts::getDefaultLogger("ActsSpToTracccAlg", level));
            }));

    py::class_<Cfg>(cls, "Config")
        .def(py::init<>())
        .def_readwrite("inputSpacePoints", &Cfg::inputSpacePoints)
        .def_readwrite("outputTracccSpacepoints",
                       &Cfg::outputTracccSpacepoints);
  }
}
