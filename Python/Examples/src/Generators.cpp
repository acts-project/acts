// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Utilities/ParametricParticleGenerator.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <cstddef>
#include <memory>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {

void addGenerators(py::module& mex) {
  {
    using Config = EventGenerator::Config;
    auto gen =
        py::class_<EventGenerator, IReader, std::shared_ptr<EventGenerator>>(
            mex, "EventGenerator")
            .def(py::init<const Config&, Logging::Level>(), py::arg("config"),
                 py::arg("level"))
            .def_property_readonly("config", &EventGenerator::config);

    py::class_<PrimaryVertexPositionGenerator,
               std::shared_ptr<PrimaryVertexPositionGenerator>>(
        gen, "VertexGenerator");
    py::class_<ParticlesGenerator, std::shared_ptr<ParticlesGenerator>>(
        gen, "ParticlesGenerator");
    py::class_<MultiplicityGenerator, std::shared_ptr<MultiplicityGenerator>>(
        gen, "MultiplicityGenerator");

    using EventGenerator = EventGenerator;
    using Generator = EventGenerator::Generator;
    py::class_<Generator>(gen, "Generator")
        .def(py::init<>())
        .def(py::init<std::shared_ptr<MultiplicityGenerator>,
                      std::shared_ptr<PrimaryVertexPositionGenerator>,
                      std::shared_ptr<ParticlesGenerator>>(),
             py::arg("multiplicity"), py::arg("vertex"), py::arg("particles"))
        .def_readwrite("multiplicity", &Generator::multiplicity)
        .def_readwrite("vertex", &Generator::vertex)
        .def_readwrite("particles", &Generator::particles);

    auto config = py::class_<Config>(gen, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(config, outputEvent, generators, randomNumbers,
                       printListing);
  }

  py::class_<GaussianPrimaryVertexPositionGenerator,
             PrimaryVertexPositionGenerator,
             std::shared_ptr<GaussianPrimaryVertexPositionGenerator>>(
      mex, "GaussianVertexGenerator")
      .def(py::init<>())
      .def(py::init([](const Vector4& stddev, const Vector4& mean) {
             GaussianPrimaryVertexPositionGenerator g;
             g.stddev = stddev;
             g.mean = mean;
             return g;
           }),
           py::arg("stddev"), py::arg("mean"))
      .def_readwrite("stddev", &GaussianPrimaryVertexPositionGenerator::stddev)
      .def_readwrite("mean", &GaussianPrimaryVertexPositionGenerator::mean);
  py::class_<GaussianDisplacedVertexPositionGenerator,
             PrimaryVertexPositionGenerator,
             std::shared_ptr<GaussianDisplacedVertexPositionGenerator>>(
      mex, "GaussianDisplacedVertexPositionGenerator")
      .def(py::init<>())
      .def(py::init([](double rMean, double rStdDev, double zMean,
                       double zStdDev, double tMean, double tStdDev) {
             GaussianDisplacedVertexPositionGenerator g;
             g.rMean = rMean;
             g.rStdDev = rStdDev;
             g.zMean = zMean;
             g.zStdDev = zStdDev;
             g.tMean = tMean;
             g.tStdDev = tStdDev;
             return g;
           }),
           py::arg("rMean"), py::arg("rStdDev"), py::arg("zMean"),
           py::arg("zStdDev"), py::arg("tMean"), py::arg("tStdDev"))
      .def_readwrite("rMean", &GaussianDisplacedVertexPositionGenerator::rMean)
      .def_readwrite("rStdDev",
                     &GaussianDisplacedVertexPositionGenerator::rStdDev)
      .def_readwrite("zMean", &GaussianDisplacedVertexPositionGenerator::zMean)
      .def_readwrite("zStdDev",
                     &GaussianDisplacedVertexPositionGenerator::zStdDev)
      .def_readwrite("tMean", &GaussianDisplacedVertexPositionGenerator::tMean)
      .def_readwrite("tStdDev",
                     &GaussianDisplacedVertexPositionGenerator::tStdDev);

  py::class_<UniformPrimaryVertexPositionGenerator,
             PrimaryVertexPositionGenerator,
             std::shared_ptr<UniformPrimaryVertexPositionGenerator>>(
      mex, "UniformVertexGenerator")
      .def(py::init<>())
      .def(py::init([](const Acts::Vector4& min, const Acts::Vector4& max) {
             UniformPrimaryVertexPositionGenerator g;
             g.min = min;
             g.max = max;
             return g;
           }),
           py::arg("min"), py::arg("max"))
      .def_readwrite("min", &UniformPrimaryVertexPositionGenerator::min)
      .def_readwrite("max", &UniformPrimaryVertexPositionGenerator::max);

  py::class_<FixedPrimaryVertexPositionGenerator,
             PrimaryVertexPositionGenerator,
             std::shared_ptr<FixedPrimaryVertexPositionGenerator>>(
      mex, "FixedVertexGenerator")
      .def(py::init<>())
      .def(py::init([](const Vector4& v) {
             FixedPrimaryVertexPositionGenerator g;
             g.fixed = v;
             return g;
           }),
           py::arg("fixed"))
      .def_readwrite("fixed", &FixedPrimaryVertexPositionGenerator::fixed);

  py::class_<SimParticle>(mex, "SimParticle");
  py::class_<SimParticleContainer>(mex, "SimParticleContainer");

  {
    using Config = ParametricParticleGenerator::Config;
    auto gen = py::class_<ParametricParticleGenerator, ParticlesGenerator,
                          std::shared_ptr<ParametricParticleGenerator>>(
                   mex, "ParametricParticleGenerator")
                   .def(py::init<const Config&>());

    py::class_<Config>(gen, "Config")
        .def(py::init<>())
        .def_readwrite("phiMin", &Config::phiMin)
        .def_readwrite("phiMax", &Config::phiMax)
        .def_readwrite("thetaMin", &Config::thetaMin)
        .def_readwrite("thetaMax", &Config::thetaMax)
        .def_readwrite("etaUniform", &Config::etaUniform)
        .def_readwrite("pMin", &Config::pMin)
        .def_readwrite("pMax", &Config::pMax)
        .def_readwrite("pTransverse", &Config::pTransverse)
        .def_readwrite("pLogUniform", &Config::pLogUniform)
        .def_readwrite("pdg", &Config::pdg)
        .def_readwrite("randomizeCharge", &Config::randomizeCharge)
        .def_readwrite("numParticles", &Config::numParticles)
        .def_readwrite("mass", &Config::mass)
        .def_readwrite("charge", &Config::charge)
        .def_property(
            "p", [](Config& cfg) { return std::pair{cfg.pMin, cfg.pMax}; },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.pMin = value.first;
              cfg.pMax = value.second;
            })
        .def_property(
            "phi",
            [](Config& cfg) { return std::pair{cfg.phiMin, cfg.phiMax}; },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.phiMin = value.first;
              cfg.phiMax = value.second;
            })
        .def_property(
            "theta",
            [](Config& cfg) { return std::pair{cfg.thetaMin, cfg.thetaMax}; },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.thetaMin = value.first;
              cfg.thetaMax = value.second;
            })
        .def_property(
            "eta",
            [](Config& cfg) {
              return std::pair{AngleHelpers::etaFromTheta(cfg.thetaMin),
                               AngleHelpers::etaFromTheta(cfg.thetaMax)};
            },
            [](Config& cfg, std::pair<double, double> value) {
              cfg.thetaMin = AngleHelpers::thetaFromEta(value.first);
              cfg.thetaMax = AngleHelpers::thetaFromEta(value.second);
            });
  }

  py::class_<FixedMultiplicityGenerator, MultiplicityGenerator,
             std::shared_ptr<FixedMultiplicityGenerator>>(
      mex, "FixedMultiplicityGenerator")
      .def(py::init<>())
      .def(py::init([](std::size_t n) {
             FixedMultiplicityGenerator g;
             g.n = n;
             return g;
           }),
           py::arg("n"))
      .def_readwrite("n", &FixedMultiplicityGenerator::n);

  py::class_<PoissonMultiplicityGenerator, MultiplicityGenerator,
             std::shared_ptr<PoissonMultiplicityGenerator>>(
      mex, "PoissonMultiplicityGenerator")
      .def(py::init<>())
      .def(py::init([](double mean) {
             PoissonMultiplicityGenerator g;
             g.mean = mean;
             return g;
           }),
           py::arg("mean"))
      .def_readwrite("mean", &PoissonMultiplicityGenerator::mean);
}

}  // namespace ActsPython
