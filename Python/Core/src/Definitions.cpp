// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/PdgParticle.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/ParticleHypothesis.hpp"

#include <format>
#include <type_traits>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace ActsPython {

/// @brief This adds the classes from Core/Definitions to the python module
/// @param m the pybind11 core module
void addDefinitions(py::module_& m) {
  using namespace Acts;
  // Add algebra classes here
  py::class_<Vector2>(m, "Vector2")
      .def(py::init<double, double>())
      .def(py::init([](std::array<double, 2> a) {
        Vector2 v;
        v << a[0], a[1];
        return v;
      }))
      .def("__getitem__",
           [](const Vector2& self, Eigen::Index i) { return self[i]; })
      .def("__str__", [](const Vector2& self) {
        std::stringstream ss;
        ss << self.transpose();
        return ss.str();
      });

  py::class_<Vector3>(m, "Vector3")
      .def(py::init<double, double, double>())
      .def(py::init([](std::array<double, 3> a) {
        Vector3 v;
        v << a[0], a[1], a[2];
        return v;
      }))
      .def_static("UnitX", []() -> Vector3 { return Vector3::UnitX(); })
      .def_static("UnitY", []() -> Vector3 { return Vector3::UnitY(); })
      .def_static("UnitZ", []() -> Vector3 { return Vector3::UnitZ(); })
      .def("__add__",
           [](const Vector3& self, const Vector3& other) {
             return (self + other).eval();
           })
      .def("__sub__",
           [](const Vector3& self, const Vector3& other) {
             return (self - other).eval();
           })
      .def("__mul__",
           [](const Vector3& self, const Vector3& other) {
             return self.cwiseProduct(other).eval();
           })
      .def("cross",
           [](const Vector3& self, const Vector3& other) {
             return self.cross(other).eval();
           })
      .def("__getitem__",
           [](const Vector3& self, Eigen::Index i) { return self[i]; })
      .def("__str__", [](const Vector3& self) {
        return std::format("({}, {}, {})", self[0], self[1], self[2]);
      });

  py::class_<Vector4>(m, "Vector4")
      .def(py::init<double, double, double, double>())
      .def(py::init([](std::array<double, 4> a) {
        Vector4 v;
        v << a[0], a[1], a[2], a[3];
        return v;
      }))
      .def("__add__",
           [](const Vector4& self, const Vector4& other) {
             return (self + other).eval();
           })
      .def("__sub__",
           [](const Vector4& self, const Vector4& other) {
             return (self - other).eval();
           })
      .def("__mul__",
           [](const Vector4& self, const Vector4& other) {
             return self.cwiseProduct(other).eval();
           })
      .def("__getitem__",
           [](const Vector4& self, Eigen::Index i) { return self[i]; })
      .def("__str__", [](const Vector4& self) {
        return std::format("({}, {}, {}, {})", self[0], self[1], self[2],
                           self[3]);
      });

  py::class_<Translation3>(m, "Translation3")
      .def(py::init([](const Vector3& a) { return Translation3(a); }))
      .def(py::init([](std::array<double, 3> a) {
        return Translation3(Vector3(a[0], a[1], a[2]));
      }))
      .def("__getitem__", [](const Translation3& self,
                             Eigen::Index i) { return self.translation()[i]; })
      .def("__str__", [](const Translation3& self) {
        std::stringstream ss;
        ss << self.translation().transpose();
        return ss.str();
      });

  py::class_<RotationMatrix3>(m, "RotationMatrix3")
      .def(py::init([](const Vector3& u, const Vector3& v, const Vector3& w) {
        RotationMatrix3 r;
        r.col(0) = u;
        r.col(1) = v;
        r.col(2) = w;
        return r;
      }))
      .def("__getitem__",
           [](const RotationMatrix3& self, std::array<Eigen::Index, 2> ij) {
             return self.matrix()(ij[0], ij[1]);
           })
      .def("__str__", [](const RotationMatrix3& self) {
        std::stringstream ss;
        ss << self.matrix();
        return ss.str();
      });

  py::class_<AngleAxis3>(m, "AngleAxis3")
      .def(py::init([](double angle, const Vector3& axis) {
        return AngleAxis3(angle, axis);
      }))
      .def("__str__", [](const AngleAxis3& self) {
        std::stringstream ss;
        ss << "Angle: " << self.angle()
           << ", Axis: " << self.axis().transpose();
        return ss.str();
      });

  py::class_<Transform3>(m, "Transform3")
      .def(py::init([](const Vector3& translation) -> Transform3 {
        return Transform3{Translation3{translation}};
      }))
      .def(py::init([](const Vector3& translation,
                       const RotationMatrix3& rotation) -> Transform3 {
        Transform3 t;
        t.prerotate(rotation);
        t.pretranslate(translation);
        return t;
      }))
      .def_property_readonly(
          "translation",
          [](const Transform3& self) -> Vector3 { return self.translation(); })
      .def_property_readonly("rotation",
                             [](const Transform3& self) -> RotationMatrix3 {
                               return self.rotation();
                             })
      .def_static("Identity", &Transform3::Identity)
      .def("__mul__", [](const Transform3& self,
                         const Transform3& other) { return self * other; })
      .def("__mul__", [](const Transform3& self,
                         const Translation3& other) { return self * other; })
      .def("__mul__", [](const Transform3& self,
                         const AngleAxis3& other) { return self * other; })
      .def("__str__", [](const Transform3& self) {
        std::stringstream ss;
        ss << self.matrix();
        return ss.str();
      });

  // Add the unit constants
  auto u = m.def_submodule("UnitConstants");

#define UNIT(x) u.attr(#x) = UnitConstants::x;

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
  UNIT(T)
  UNIT(Gauss)
  UNIT(kGauss)
  UNIT(mol)

#undef UNIT

  // Add the PdgParticle enum
  py::enum_<PdgParticle>(m, "PdgParticle")
      .value("eInvalid", PdgParticle::eInvalid)
      .value("eElectron", PdgParticle::eElectron)
      .value("eAntiElectron", PdgParticle::eAntiElectron)
      .value("ePositron", PdgParticle::ePositron)
      .value("eMuon", PdgParticle::eMuon)
      .value("eAntiMuon", PdgParticle::eAntiMuon)
      .value("eTau", PdgParticle::eTau)
      .value("eAntiTau", PdgParticle::eAntiTau)
      .value("eGamma", PdgParticle::eGamma)
      .value("ePionZero", PdgParticle::ePionZero)
      .value("ePionPlus", PdgParticle::ePionPlus)
      .value("ePionMinus", PdgParticle::ePionMinus)
      .value("eKaonPlus", PdgParticle::eKaonPlus)
      .value("eKaonMinus", PdgParticle::eKaonMinus)
      .value("eNeutron", PdgParticle::eNeutron)
      .value("eAntiNeutron", PdgParticle::eAntiNeutron)
      .value("eProton", PdgParticle::eProton)
      .value("eAntiProton", PdgParticle::eAntiProton)
      .value("eLead", PdgParticle::eLead)
      .value("eJPsi", PdgParticle::eJPsi)
      .value("eB0", PdgParticle::eB0)
      .value("eBPlus", PdgParticle::eBPlus)
      .value("eD0", PdgParticle::eD0)
      .value("eDPlus", PdgParticle::eDPlus)
      .value("eAntiB0", PdgParticle::eAntiB0)
      .value("eAntiD0", PdgParticle::eAntiD0)
      .value("eNeutrinoE", PdgParticle::eNeutrinoE)
      .value("eNeutrinoMu", PdgParticle::eNeutrinoMu)
      .value("eNeutrinoTau", PdgParticle::eNeutrinoTau)
      .value("eAntiNeutrinoE", PdgParticle::eAntiNeutrinoE)
      .value("eAntiNeutrinoMu", PdgParticle::eAntiNeutrinoMu)
      .value("eAntiNeutrinoTau", PdgParticle::eAntiNeutrinoTau);

  // Add the parsePdgParticle function
  m.def("parsePdgParticle", &parsePdgParticle, py::arg("name"));

  py::class_<GenericParticleHypothesis<AnyCharge>>(
      m, "GenericParticleHypothesisAnyCharge");

  py::class_<ParticleHypothesis, GenericParticleHypothesis<AnyCharge>>(
      m, "ParticleHypothesis")
      .def(py::init([](PdgParticle absPdg, float mass, float absCharge) {
             return ParticleHypothesis(absPdg, mass, AnyCharge{absCharge});
           }),
           py::arg("pdg"), py::arg("mass"), py::arg("absCharge"))
      .def(py::init([](std::underlying_type_t<PdgParticle> absPdg, float mass,
                       float absCharge) {
             return ParticleHypothesis(static_cast<PdgParticle>(absPdg), mass,
                                       AnyCharge{absCharge});
           }),
           py::arg("absPdg"), py::arg("mass"), py::arg("absCharge"))
      .def("__str__",
           [](const ParticleHypothesis& particleHypothesis) {
             std::stringstream os;
             particleHypothesis.toStream(os);
             return os.str();
           })
      .def_property_readonly("absolutePdg", &ParticleHypothesis::absolutePdg)
      .def_property_readonly("mass", &ParticleHypothesis::mass)
      .def_property_readonly("absoluteCharge",
                             &ParticleHypothesis::absoluteCharge)
      .def_property_readonly_static(
          "muon",
          [](const py::object& /*self*/) { return ParticleHypothesis::muon(); })
      .def_property_readonly_static(
          "pion",
          [](const py::object& /*self*/) { return ParticleHypothesis::pion(); })
      .def_property_readonly_static("electron",
                                    [](const py::object& /*self*/) {
                                      return ParticleHypothesis::electron();
                                    })
      .def_property_readonly_static(
          "kaon",
          [](const py::object& /*self*/) { return ParticleHypothesis::kaon(); })
      .def_property_readonly_static("proton",
                                    [](const py::object& /*self*/) {
                                      return ParticleHypothesis::proton();
                                    })
      .def_property_readonly_static("geantino",
                                    [](const py::object& /*self*/) {
                                      return ParticleHypothesis::geantino();
                                    })
      .def_property_readonly_static(
          "chargedGeantino", [](const py::object& /*self*/) {
            return ParticleHypothesis::chargedGeantino();
          });
}

}  // namespace ActsPython
