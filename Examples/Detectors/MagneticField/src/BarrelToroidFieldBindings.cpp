#include <memory>
#include <pybind11/pybind11.h>

#include "ActsExamples/MagneticField/BarrelToroidField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

namespace py = pybind11;
using Acts::BarrelToroidField;

PYBIND11_MODULE(acts_barrel_field, m) {
  m.doc() = "Pybind for Acts::BarrelToroidField";

  // Bind Config first
  py::class_<BarrelToroidField::Config>(m, "Config")
      .def(py::init<>())
      .def_readwrite("R_in",   &BarrelToroidField::Config::R_in)
      .def_readwrite("R_out",  &BarrelToroidField::Config::R_out)
      .def_readwrite("c",      &BarrelToroidField::Config::c)
      .def_readwrite("b",      &BarrelToroidField::Config::b)
      .def_readwrite("I",      &BarrelToroidField::Config::I)
      .def_readwrite("Nturns", &BarrelToroidField::Config::Nturns)
      .def_readwrite("nArc",      &BarrelToroidField::Config::nArc)
      .def_readwrite("nStraight", &BarrelToroidField::Config::nStraight)
      .def_readwrite("closeLoop", &BarrelToroidField::Config::closeLoop)
      .def_readwrite("theta0_deg",    &BarrelToroidField::Config::theta0_deg)
      .def_readwrite("thetaStep_deg", &BarrelToroidField::Config::thetaStep_deg)
      .def_readwrite("nCoils",        &BarrelToroidField::Config::nCoils)
      .def_readwrite("eps",           &BarrelToroidField::Config::eps);

  // We expose it as a subclass of MagneticFieldProvider so it can be passed anywhere ACTS expects a field provider
  py::class_<BarrelToroidField, std::shared_ptr<BarrelToroidField>, Acts::MagneticFieldProvider>(m, "BarrelToroidField")
      .def(py::init<>())
      .def(py::init<const BarrelToroidField::Config&>())
      .def_property_readonly("config", &BarrelToroidField::config, py::return_value_policy::reference);

  // Convenience factory returning a base-class pointer (sometimes handy)
  m.def("make_barrel_toroid_field",
        [](const BarrelToroidField::Config& cfg) -> std::shared_ptr<Acts::MagneticFieldProvider> {
          return std::make_shared<BarrelToroidField>(cfg);
        },
        py::arg("cfg") = BarrelToroidField::Config{});
  
  // Also allow factory with default config
  m.def("make_barrel_toroid_field",
        []() -> std::shared_ptr<Acts::MagneticFieldProvider> {
          return std::make_shared<BarrelToroidField>();
        });
}

