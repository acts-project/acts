// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/PdgParticle.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {

void addUnits(Context& ctx) {
  auto& m = ctx.get("main");
  auto u = m.def_submodule("UnitConstants");
  using namespace Acts::UnitConstants;

#define UNIT(x) u.attr(#x) = x;

  UNIT(fm)
  UNIT(pm)
  UNIT(um)
  UNIT(nm)
  UNIT(mm)
  UNIT(cm)
  UNIT(m)
  UNIT(km)
  UNIT(mm2)
  UNIT(cm2)
  UNIT(m2)
  UNIT(mm3)
  UNIT(cm3)
  UNIT(m3)
  UNIT(s)
  UNIT(fs)
  UNIT(ps)
  UNIT(ns)
  UNIT(us)
  UNIT(ms)
  UNIT(min)
  UNIT(h)
  UNIT(mrad)
  UNIT(rad)
  UNIT(degree)
  UNIT(eV)
  UNIT(keV)
  UNIT(MeV)
  UNIT(GeV)
  UNIT(TeV)
  UNIT(J)
  UNIT(u)
  UNIT(g)
  UNIT(kg)
  UNIT(e)
  UNIT(C)
  UNIT(T)
  UNIT(Gauss)
  UNIT(kGauss)
  UNIT(mol)

#undef UNIT
}

void addLogging(Acts::Python::Context& ctx) {
  auto& m = ctx.get("main");
  auto logging = m.def_submodule("logging", "");
  py::enum_<Acts::Logging::Level>(logging, "Level")
      .value("VERBOSE", Acts::Logging::VERBOSE)
      .value("DEBUG", Acts::Logging::DEBUG)
      .value("INFO", Acts::Logging::INFO)
      .value("WARNING", Acts::Logging::WARNING)
      .value("ERROR", Acts::Logging::ERROR)
      .value("FATAL", Acts::Logging::FATAL)
      .export_values();
}

void addPdgParticle(Acts::Python::Context& ctx) {
  auto& m = ctx.get("main");
  py::enum_<Acts::PdgParticle>(m, "PdgParticle")
      .value("eInvalid", Acts::PdgParticle::eInvalid)
      .value("eElectron", Acts::PdgParticle::eElectron)
      .value("eAntiElectron", Acts::PdgParticle::eAntiElectron)
      .value("ePositron", Acts::PdgParticle::ePositron)
      .value("eMuon", Acts::PdgParticle::eMuon)
      .value("eAntiMuon", Acts::PdgParticle::eAntiMuon)
      .value("eTau", Acts::PdgParticle::eTau)
      .value("eAntiTau", Acts::PdgParticle::eAntiTau)
      .value("eGamma", Acts::PdgParticle::eGamma)
      .value("ePionZero", Acts::PdgParticle::ePionZero)
      .value("ePionPlus", Acts::PdgParticle::ePionPlus)
      .value("ePionMinus", Acts::PdgParticle::ePionMinus)
      .value("eNeutron", Acts::PdgParticle::eNeutron)
      .value("eAntiNeutron", Acts::PdgParticle::eAntiNeutron)
      .value("eProton", Acts::PdgParticle::eProton)
      .value("eAntiProton", Acts::PdgParticle::eAntiProton);
}

void addAlgebra(Acts::Python::Context& ctx) {
  auto& m = ctx.get("main");
  py::class_<Acts::Vector3>(m, "Vector3")
      .def(py::init<double, double, double>())
      .def(py::init([](std::array<double, 3> a) {
        Acts::Vector3 v;
        v << a[0], a[1], a[2];
        return v;
      }));

  py::class_<Acts::Vector4>(m, "Vector4")
      .def(py::init<double, double, double, double>())
      .def(py::init([](std::array<double, 4> a) {
        Acts::Vector4 v;
        v << a[0], a[1], a[2], a[3];
        return v;
      }));
}

}  // namespace Acts::Python
