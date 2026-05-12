// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Utilities/AngleHelpers.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/SimHit.hpp"
#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsExamples/Generators/EventGenerator.hpp"
#include "ActsExamples/Utilities/ParametricParticleGenerator.hpp"
#include "ActsExamples/Utilities/VertexGenerators.hpp"
#include "ActsFatras/EventData/Hit.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
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
        gen, "VertexGenerator")
        .def("__call__", &PrimaryVertexPositionGenerator::operator(),
             py::arg("rng"), py::arg("eventNumber"));
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

  py::class_<AdditiveVertexPositionGenerator, PrimaryVertexPositionGenerator,
             std::shared_ptr<AdditiveVertexPositionGenerator>>(
      mex, "AdditiveVertexGenerator")
      .def(py::init<>())
      .def(
          py::init([](const std::vector<
                       std::shared_ptr<PrimaryVertexPositionGenerator>>& gens) {
            AdditiveVertexPositionGenerator g;
            g.generators = gens;
            return g;
          }),
          py::arg("generators"))
      .def_readwrite("generators",
                     &AdditiveVertexPositionGenerator::generators);

  py::class_<LumiBlockVertexPositionGenerator, PrimaryVertexPositionGenerator,
             std::shared_ptr<LumiBlockVertexPositionGenerator>>(
      mex, "LumiBlockVertexGenerator")
      .def(py::init<>())
      .def(py::init([](std::size_t blockSize, const Vector4& stddev) {
             LumiBlockVertexPositionGenerator g;
             g.blockSize = blockSize;
             g.stddev = stddev;
             return g;
           }),
           py::arg("blockSize"), py::arg("stddev"))
      .def_readwrite("blockSize", &LumiBlockVertexPositionGenerator::blockSize)
      .def_readwrite("stddev", &LumiBlockVertexPositionGenerator::stddev);

  py::class_<LumiBlockRotationVertexPositionGenerator,
             PrimaryVertexPositionGenerator,
             std::shared_ptr<LumiBlockRotationVertexPositionGenerator>>(
      mex, "LumiBlockRotationVertexGenerator")
      .def(py::init<>())
      .def(py::init([](std::shared_ptr<PrimaryVertexPositionGenerator> base,
                       std::size_t blockSize, double xAngleStddev,
                       double yAngleStddev) {
             LumiBlockRotationVertexPositionGenerator g;
             g.base = std::move(base);
             g.blockSize = blockSize;
             g.xAngleStddev = xAngleStddev;
             g.yAngleStddev = yAngleStddev;
             return g;
           }),
           py::arg("base"), py::arg("blockSize"), py::arg("xAngleStddev"),
           py::arg("yAngleStddev"))
      .def_readwrite("base", &LumiBlockRotationVertexPositionGenerator::base)
      .def_readwrite("blockSize",
                     &LumiBlockRotationVertexPositionGenerator::blockSize)
      .def_readwrite("xAngleStddev",
                     &LumiBlockRotationVertexPositionGenerator::xAngleStddev)
      .def_readwrite("yAngleStddev",
                     &LumiBlockRotationVertexPositionGenerator::yAngleStddev);

  // Aliases for Fatras types mirroring C++
  auto fatras = py::module_::import("acts.fatras");
  mex.attr("SimBarcode") = fatras.attr("Barcode");
  mex.attr("GenerationProcess") = fatras.attr("GenerationProcess");
  mex.attr("SimulationOutcome") = fatras.attr("SimulationOutcome");
  mex.attr("SimParticleState") = fatras.attr("Particle");

  // SimParticle
  py::class_<SimParticle>(mex, "SimParticle")
      .def(py::init<>())
      .def(py::init<SimBarcode, Acts::PdgParticle, double, double>(),
           py::arg("particleId"), py::arg("pdg"), py::arg("charge"),
           py::arg("mass"))
      .def(py::init<SimBarcode, Acts::PdgParticle>(), py::arg("particleId"),
           py::arg("pdg"))
      .def(py::init<const SimParticleState&, const SimParticleState&>(),
           py::arg("initial"), py::arg("final"))
      .def_property("particleId", &SimParticle::particleId,
                    [](SimParticle& p, SimBarcode b) { p.setParticleId(b); })
      .def_property("pdg", &SimParticle::pdg,
                    [](SimParticle& p, Acts::PdgParticle v) { p.setPdg(v); })
      .def_property_readonly("absolutePdg", &SimParticle::absolutePdg)
      .def_property("charge", &SimParticle::charge,
                    [](SimParticle& p, double v) { p.setCharge(v); })
      .def_property_readonly("absoluteCharge", &SimParticle::absoluteCharge)
      .def_property("mass", &SimParticle::mass,
                    [](SimParticle& p, double v) { p.setMass(v); })
      .def_property("process", &SimParticle::process,
                    [](SimParticle& p, ActsFatras::GenerationProcess v) {
                      p.setProcess(v);
                    })
      .def_property_readonly("isSecondary", &SimParticle::isSecondary)
      .def_property("fourPosition", &SimParticle::fourPosition,
                    [](SimParticle& p, const Acts::Vector4& pos4) {
                      p.initialState().setPosition4(pos4);
                    })
      .def_property_readonly(
          "position",
          [](const SimParticle& p) { return Vector3(p.position()); })
      .def_property_readonly("time", &SimParticle::time)
      .def_property_readonly("fourMomentum", &SimParticle::fourMomentum)
      .def_property("direction", &SimParticle::direction,
                    [](SimParticle& p, const Acts::Vector3& dir) {
                      p.initialState().setDirection(dir);
                    })
      .def_property_readonly("theta", &SimParticle::theta)
      .def_property_readonly("phi", &SimParticle::phi)
      .def_property_readonly("transverseMomentum",
                             &SimParticle::transverseMomentum)
      .def_property("absoluteMomentum", &SimParticle::absoluteMomentum,
                    [](SimParticle& p, double v) {
                      p.initialState().setAbsoluteMomentum(v);
                    })
      .def_property_readonly("momentum", &SimParticle::momentum)
      .def_property_readonly("energy", &SimParticle::energy)
      .def_property_readonly("energyLoss", &SimParticle::energyLoss)
      .def_property_readonly("pathInX0", &SimParticle::pathInX0)
      .def_property_readonly("pathInL0", &SimParticle::pathInL0)
      .def(
          "setFinalMaterialPassed",
          [](SimParticle& p, double x0, double l0) -> SimParticle& {
            p.finalState().setMaterialPassed(x0, l0);
            return p;
          },
          py::arg("pathInX0"), py::arg("pathInL0"),
          py::return_value_policy::reference_internal)
      .def_property("numberOfHits", &SimParticle::numberOfHits,
                    [](SimParticle& p, std::uint32_t n) {
                      p.finalState().setNumberOfHits(n);
                    })
      .def_property("outcome", &SimParticle::outcome,
                    [](SimParticle& p, ActsFatras::SimulationOutcome o) {
                      p.finalState().setOutcome(o);
                    })
      .def_property_readonly(
          "initialState",
          py::overload_cast<>(&SimParticle::initialState, py::const_))
      .def_property_readonly(
          "finalState",
          py::overload_cast<>(&SimParticle::finalState, py::const_))
      .def("withParticleId", &SimParticle::withParticleId,
           py::arg("particleId"))
      .def("__repr__", [](const SimParticle& p) {
        std::ostringstream oss;
        oss << p;
        return oss.str();
      });

  auto simParticleContainer =
      py::classh<SimParticleContainer>(mex, "SimParticleContainer")
          .def(py::init<>())
          .def("__len__",
               [](const SimParticleContainer& c) { return c.size(); })
          .def(
              "__iter__",
              [](const SimParticleContainer& c) {
                return py::make_iterator(c.begin(), c.end());
              },
              py::keep_alive<0, 1>())
          .def("__contains__",
               [](const SimParticleContainer& c, SimBarcode barcode) {
                 return c.find(barcode) != c.end();
               })
          .def("__repr__",
               [](const SimParticleContainer& c) {
                 std::ostringstream oss;
                 oss << "SimParticleContainer(" << c.size() << " particles)";
                 return oss.str();
               })
          .def(
              "insert",
              [](SimParticleContainer& c, const SimParticle& p) {
                c.insert(p);
              },
              py::arg("particle"));

  WhiteBoardRegistry::registerClass(simParticleContainer);

  py::class_<SimHit>(mex, "SimHit")
      .def(py::init<Acts::GeometryIdentifier, SimBarcode, Acts::Vector4,
                    Acts::Vector4, Acts::Vector4, std::int32_t>(),
           py::arg("geometryId"), py::arg("particleId"), py::arg("pos4"),
           py::arg("before4"), py::arg("after4"), py::arg("index") = -1)
      .def_property_readonly("geometryId", &SimHit::geometryId)
      .def_property_readonly("particleId", &SimHit::particleId)
      .def_property_readonly("index", &SimHit::index)
      .def_property_readonly("fourPosition", &SimHit::fourPosition)
      .def_property_readonly("time", &SimHit::time)
      .def_property_readonly("momentum4Before", &SimHit::momentum4Before)
      .def_property_readonly("momentum4After", &SimHit::momentum4After)
      .def_property_readonly("depositedEnergy", &SimHit::depositedEnergy);

  auto simHitContainer =
      py::classh<SimHitContainer>(mex, "SimHitContainer")
          .def(py::init<>())
          .def("__len__", [](const SimHitContainer& c) { return c.size(); })
          .def(
              "__iter__",
              [](const SimHitContainer& c) {
                return py::make_iterator(c.begin(), c.end());
              },
              py::keep_alive<0, 1>())
          .def(
              "insert",
              [](SimHitContainer& c, const SimHit& h) { c.insert(h); },
              py::arg("hit"))
          .def("__repr__", [](const SimHitContainer& c) {
            std::ostringstream oss;
            oss << "SimHitContainer(" << c.size() << " hits)";
            return oss.str();
          });

  WhiteBoardRegistry::registerClass(simHitContainer);

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
