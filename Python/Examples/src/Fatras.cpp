// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/SimParticle.hpp"
#include "ActsFatras/EventData/Barcode.hpp"
#include "ActsFatras/EventData/ParticleOutcome.hpp"
#include "ActsFatras/EventData/ProcessType.hpp"

#include <sstream>

#include <pybind11/pybind11.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addFatras(py::module& mex) {
  auto fatras = mex.def_submodule("fatras", "ACTS Fatras simulation types");

  // Barcode (ActsFatras::Barcode) - the canonical Fatras type
  py::class_<SimBarcode>(fatras, "Barcode")
      .def(py::init<>())
      .def_static("Invalid", &SimBarcode::Invalid)
      .def("isValid", [](const SimBarcode& b) { return b.isValid(); })
      .def_property(
          "vertexPrimary",
          [](const SimBarcode& b) { return b.vertexPrimary(); },
          [](SimBarcode& b, SimBarcode::PrimaryVertexId id) {
            b = b.withVertexPrimary(id);
          })
      .def_property(
          "vertexSecondary",
          [](const SimBarcode& b) { return b.vertexSecondary(); },
          [](SimBarcode& b, SimBarcode::SecondaryVertexId id) {
            b = b.withVertexSecondary(id);
          })
      .def_property(
          "particle", [](const SimBarcode& b) { return b.particle(); },
          [](SimBarcode& b, SimBarcode::ParticleId id) {
            b = b.withParticle(id);
          })
      .def_property(
          "generation", [](const SimBarcode& b) { return b.generation(); },
          [](SimBarcode& b, SimBarcode::GenerationId id) {
            b = b.withGeneration(id);
          })
      .def_property(
          "subParticle", [](const SimBarcode& b) { return b.subParticle(); },
          [](SimBarcode& b, SimBarcode::SubParticleId id) {
            b = b.withSubParticle(id);
          })
      .def("__repr__", [](const SimBarcode& b) {
        std::ostringstream oss;
        oss << b;
        return oss.str();
      });

  // ProcessType and ParticleOutcome enums
  py::enum_<ActsFatras::ProcessType>(fatras, "ProcessType")
      .value("eUndefined", ActsFatras::ProcessType::eUndefined)
      .value("eDecay", ActsFatras::ProcessType::eDecay)
      .value("ePhotonConversion", ActsFatras::ProcessType::ePhotonConversion)
      .value("eBremsstrahlung", ActsFatras::ProcessType::eBremsstrahlung)
      .value("eNuclearInteraction",
             ActsFatras::ProcessType::eNuclearInteraction);

  py::enum_<ActsFatras::ParticleOutcome>(fatras, "ParticleOutcome")
      .value("Alive", ActsFatras::ParticleOutcome::Alive)
      .value("KilledInteraction",
             ActsFatras::ParticleOutcome::KilledInteraction)
      .value("KilledVolumeExit", ActsFatras::ParticleOutcome::KilledVolumeExit)
      .value("KilledTime", ActsFatras::ParticleOutcome::KilledTime)
      .value("KilledSecondaryParticle",
             ActsFatras::ParticleOutcome::KilledSecondaryParticle);
}

}  // namespace ActsPython
