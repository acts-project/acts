// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray core include(s)
#include "detray/core/detector.hpp"

// Detray propagation include(s)
#include "detray/propagator/propagation_config.hpp"

// Detray event generator include(s)
#include "detray/test/common/event_generator/uniform_track_generator_config.hpp"

// Detray test include(s)
#include "detray/test/cpu/material_scan.hpp"
#include "detray/test/framework/test_configuration.hpp"
#include "detray/test/validation/material_validation_config.hpp"

// Detray algebra plugin + detector metadata
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

namespace py = pybind11;

namespace {

using scalar_t = DETRAY_CUSTOM_SCALARTYPE;
using algebra_t = detray::array<scalar_t>;
using detector_t = detray::detector<detray::default_metadata<algebra_t>>;
using material_scan_config_t = detray::test::material_scan<detector_t>::config;
using propagation_config_t = detray::propagation::config;
using track_generator_config_t =
    detray::uniform_track_generator_config<scalar_t>;
using material_validation_config_t =
    detray::test::material_validation_config<algebra_t>;
using base_config_t = detray::test::configuration<scalar_t>;

template <typename T>
std::string to_string(const T &obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

}  // namespace

PYBIND11_MODULE(DetrayTestsPythonBindings, m) {
  m.doc() = "Detray tests bindings";

  py::class_<track_generator_config_t>(m, "TrackGeneratorConfig")
      .def(py::init<>())
      .def_property(
          "seed", [](const track_generator_config_t &c) { return c.seed(); },
          [](track_generator_config_t &c, std::uint64_t s) { c.seed(s); },
          "Monte-Carlo seed")
      .def(
          "nTracks",
          [](const track_generator_config_t &c) { return c.n_tracks(); },
          "Total number of tracks")
      .def_property(
          "phiRange",
          [](const track_generator_config_t &c) {
            const auto r = c.phi_range();
            return std::make_pair(r[0], r[1]);
          },
          [](track_generator_config_t &c,
             const std::pair<scalar_t, scalar_t> &r) {
            c.phi_range(r.first, r.second);
          },
          "Phi range (min, max) (native rad)")
      .def_property(
          "thetaRange",
          [](const track_generator_config_t &c) {
            const auto r = c.theta_range();
            return std::make_pair(r[0], r[1]);
          },
          [](track_generator_config_t &c,
             const std::pair<scalar_t, scalar_t> &r) {
            c.theta_range(r.first, r.second);
          },
          "Theta range (min, max) (native rad)")
      .def_property(
          "etaRange",
          [](const track_generator_config_t &c) {
            const auto r = c.eta_range();
            return std::make_pair(r[0], r[1]);
          },
          [](track_generator_config_t &c,
             const std::pair<scalar_t, scalar_t> &r) {
            c.eta_range(r.first, r.second);
          },
          "Eta range (min, max)")
      .def_property(
          "phiSteps",
          [](const track_generator_config_t &c) { return c.phi_steps(); },
          [](track_generator_config_t &c, std::size_t n) { c.phi_steps(n); },
          "Number of phi steps")
      .def_property(
          "thetaSteps",
          [](const track_generator_config_t &c) { return c.theta_steps(); },
          [](track_generator_config_t &c, std::size_t n) { c.theta_steps(n); },
          "Number of theta steps")
      .def_property(
          "etaSteps",
          [](const track_generator_config_t &c) { return c.eta_steps(); },
          [](track_generator_config_t &c, std::size_t n) { c.eta_steps(n); },
          "Number of eta steps")
      .def_property(
          "uniformEta",
          [](const track_generator_config_t &c) { return c.uniform_eta(); },
          [](track_generator_config_t &c, bool b) { c.uniform_eta(b); },
          "Whether to step uniformly in eta")
      .def_property(
          "origin",
          [](const track_generator_config_t &c) {
            const auto &o = c.origin();
            return std::make_tuple(o[0], o[1], o[2]);
          },
          [](track_generator_config_t &c,
             const std::tuple<scalar_t, scalar_t, scalar_t> &o) {
            c.origin(std::get<0>(o), std::get<1>(o), std::get<2>(o));
          },
          "Track origin")
      .def(
          "pT",
          [](track_generator_config_t &c,
             scalar_t p) -> track_generator_config_t & { return c.p_T(p); },
          py::arg("p"), py::return_value_policy::reference,
          "Set the transverse momentum magnitude")
      .def(
          "pTot",
          [](track_generator_config_t &c,
             scalar_t p) -> track_generator_config_t & { return c.p_tot(p); },
          py::arg("p"), py::return_value_policy::reference,
          "Set the total momentum magnitude")
      .def(
          "momRange",
          [](const track_generator_config_t &c) {
            const auto r = c.mom_range();
            return std::make_pair(r[0], r[1]);
          },
          "Momentum range")
      .def_property(
          "randomizeCharge",
          [](const track_generator_config_t &c) {
            return c.randomize_charge();
          },
          [](track_generator_config_t &c, bool b) { c.randomize_charge(b); },
          "Randomly flip the charge sign")
      .def_property(
          "time", [](const track_generator_config_t &c) { return c.time(); },
          [](track_generator_config_t &c, scalar_t t) { c.time(t); },
          "Track time")
      .def_property(
          "charge",
          [](const track_generator_config_t &c) { return c.charge(); },
          [](track_generator_config_t &c, scalar_t q) { c.charge(q); },
          "Track charge")
      .def(
          "isPT", [](const track_generator_config_t &c) { return c.is_pT(); },
          "Whether the momentum magnitude is interpreted as transverse")
      .def("__repr__", &to_string<track_generator_config_t>);

  py::class_<base_config_t>(m, "TestConfig")
      .def_property(
          "tol", [](const base_config_t &c) { return c.tol(); },
          [](base_config_t &c, scalar_t t) { c.tol(t); },
          "Tolerance to compare two floating point values")
      .def_property(
          "propagation",
          [](base_config_t &c) -> propagation_config_t & {
            return c.propagation();
          },
          [](base_config_t &c, const propagation_config_t &v) {
            c.propagation() = v;
          },
          py::return_value_policy::reference_internal,
          "Propagation configuration");

  py::class_<material_validation_config_t, base_config_t>(
      m, "MaterialValidationConfig")
      .def(py::init<>())
      .def_property(
          "name",
          [](const material_validation_config_t &c) { return c.name(); },
          [](material_validation_config_t &c, const std::string &n) {
            c.name(n);
          },
          "Name of the test")
      .def_property(
          "materialFile",
          [](const material_validation_config_t &c) {
            return c.material_file();
          },
          [](material_validation_config_t &c, const std::string &f) {
            c.material_file(f);
          },
          "Name of the output file with the navigation material traces")
      .def_property(
          "nTracks",
          [](const material_validation_config_t &c) { return c.n_tracks(); },
          [](material_validation_config_t &c, std::size_t n) { c.n_tracks(n); },
          "Maximal number of test tracks to run")
      .def_property(
          "relativeError",
          [](const material_validation_config_t &c) {
            return c.relative_error();
          },
          [](material_validation_config_t &c, scalar_t re) {
            c.relative_error(re);
          },
          "Allowed relative discrepancy between truth and navigation material")
      .def("__repr__", &to_string<material_validation_config_t>);

  py::class_<material_scan_config_t, base_config_t>(m, "MaterialScanConfig")
      .def(py::init<>())
      .def_property(
          "name", [](const material_scan_config_t &c) { return c.name(); },
          [](material_scan_config_t &c, const std::string &n) { c.name(n); },
          "Name of the test")
      .def_property(
          "materialFile",
          [](const material_scan_config_t &c) { return c.material_file(); },
          [](material_scan_config_t &c, const std::string &f) {
            c.material_file(f);
          },
          "Name of the output file with the material traces")
      .def_property(
          "overlapsRemoval",
          [](const material_scan_config_t &c) { return c.overlaps_removal(); },
          [](material_scan_config_t &c, bool o) { c.overlaps_removal(o); },
          "Perform overlaps removal")
      .def_property(
          "overlapsTol",
          [](const material_scan_config_t &c) { return c.overlaps_tol(); },
          [](material_scan_config_t &c, scalar_t t) { c.overlaps_tol(t); },
          "Tolerance for considering surfaces to be overlapping")
      .def_property(
          "trackGenerator",
          [](material_scan_config_t &c) -> track_generator_config_t & {
            return c.track_generator();
          },
          [](material_scan_config_t &c, const track_generator_config_t &v) {
            c.track_generator() = v;
          },
          py::return_value_policy::reference_internal,
          "Track generator configuration")
      .def("__repr__", &to_string<material_scan_config_t>);
}
