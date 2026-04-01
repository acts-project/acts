// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/PdgParticle.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/Particle.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <sstream>

#include <pybind11/pybind11.h>

namespace py = pybind11;

PYBIND11_MODULE(ActsFatrasPythonBindings, fatras) {
  using Barcode = ActsFatras::Barcode;
  using Particle = ActsFatras::Particle;

  py::class_<Barcode>(fatras, "Barcode")
      .def(py::init<>())
      .def_static("Invalid", &Barcode::Invalid)
      .def("isValid", [](const Barcode& b) { return b.isValid(); })
      .def_property(
          "vertexPrimary", [](const Barcode& b) { return b.vertexPrimary(); },
          [](Barcode& b, Barcode::PrimaryVertexId id) {
            b = b.withVertexPrimary(id);
          })
      .def_property(
          "vertexSecondary",
          [](const Barcode& b) { return b.vertexSecondary(); },
          [](Barcode& b, Barcode::SecondaryVertexId id) {
            b = b.withVertexSecondary(id);
          })
      .def_property(
          "particle", [](const Barcode& b) { return b.particle(); },
          [](Barcode& b, Barcode::ParticleId id) { b = b.withParticle(id); })
      .def_property(
          "generation", [](const Barcode& b) { return b.generation(); },
          [](Barcode& b, Barcode::GenerationId id) {
            b = b.withGeneration(id);
          })
      .def_property(
          "subParticle", [](const Barcode& b) { return b.subParticle(); },
          [](Barcode& b, Barcode::SubParticleId id) {
            b = b.withSubParticle(id);
          })
      .def("__repr__", [](const Barcode& b) {
        std::ostringstream oss;
        oss << b;
        return oss.str();
      });

  py::enum_<ActsFatras::ProcessType>(fatras, "ProcessType")
      .value("eUndefined", ActsFatras::ProcessType::eUndefined)
      .value("eDecay", ActsFatras::ProcessType::eDecay)
      .value("ePhotonConversion", ActsFatras::ProcessType::ePhotonConversion)
      .value("eBremsstrahlung", ActsFatras::ProcessType::eBremsstrahlung)
      .value("eNuclearInteraction",
             ActsFatras::ProcessType::eNuclearInteraction);

  py::enum_<ActsFatras::SimulationOutcome>(fatras, "SimulationOutcome")
      .value("Alive", ActsFatras::SimulationOutcome::Alive)
      .value("KilledInteraction",
             ActsFatras::SimulationOutcome::KilledInteraction)
      .value("KilledVolumeExit",
             ActsFatras::SimulationOutcome::KilledVolumeExit)
      .value("KilledTime", ActsFatras::SimulationOutcome::KilledTime)
      .value("KilledSecondaryParticle",
             ActsFatras::SimulationOutcome::KilledSecondaryParticle);

  py::class_<Particle>(fatras, "Particle")
      .def(py::init<>())
      .def(py::init<Barcode, Acts::PdgParticle, double, double>(),
           py::arg("particleId"), py::arg("pdg"), py::arg("charge"),
           py::arg("mass"))
      .def(py::init<Barcode, Acts::PdgParticle>(), py::arg("particleId"),
           py::arg("pdg"))
      .def_property_readonly("particleId", &Particle::particleId)
      .def_property_readonly("pdg", &Particle::pdg)
      .def_property_readonly("absolutePdg", &Particle::absolutePdg)
      .def_property_readonly("charge", &Particle::charge)
      .def_property_readonly("mass", &Particle::mass)
      .def_property_readonly("fourPosition", &Particle::fourPosition)
      .def_property_readonly("fourMomentum", &Particle::fourMomentum)
      .def_property_readonly("direction", &Particle::direction)
      .def_property_readonly("absoluteMomentum", &Particle::absoluteMomentum);
}
