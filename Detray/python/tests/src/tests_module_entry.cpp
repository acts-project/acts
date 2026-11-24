// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray core include(s)
#include "detray/core/detector.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"
#include "detray/propagator/propagation_config.hpp"

// Detray test include(s)
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/common/event_generator/uniform_track_generator_config.hpp"
#include "detray/test/framework/test_configuration.hpp"
#include "detray/test/framework/types.hpp"
#include "detray/test/validation/material_validation.hpp"
#include "detray/test/validation/material_validation_config.hpp"

// Vecmem include(s)
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/memory/memory_resource.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>

namespace py = pybind11;

namespace {

using scalar_t = detray::test::scalar;
using algebra_t = detray::test::algebra;
using detector_t = detray::host::detector<detray::test::default_metadata>;
using propagation_config_t = detray::propagation::config;
using toy_det_config_t = detray::toy_det_config<scalar_t>;
using test_config_t = detray::test::configuration<scalar_t>;
using uniform_track_generator_config_t =
    detray::uniform_track_generator_config<scalar_t>;
using material_scan_config_t = detray::test::material_scan_config<algebra_t>;
using material_validation_result_t =
    detray::material_validator::material_validation_result;
using material_validation_config_t =
    detray::test::material_validation_config<algebra_t>;

template <typename T>
std::string to_string(const T &obj) {
  std::ostringstream os;
  os << obj;
  return os.str();
}

/// Run the CPU material validation for the given detector.
///
/// Same sequence as the material scan and material validation test fixtures.
material_validation_result_t run_material_validation(
    const detector_t &det, const detray::name_map &names,
    const material_scan_config_t &scan_cfg,
    const material_validation_config_t &val_cfg) {
  detector_t::geometry_context gctx{};
  vecmem::host_memory_resource host_mr{};

  const auto scan_data =
      detray::material_validator::material_scan(gctx, det, names, scan_cfg);

  return detray::material_validator::material_validation<
      detray::test::run_material_validation>(&host_mr, det, names, scan_data,
                                             val_cfg);
}

/// Build a toy detector using @p mr as configured by @p cfg .
///
/// The memory resource object is automatically kept alive at least as long as
/// the returned detector object.
std::pair<py::object, detray::name_map> build_toy_detector(
    std::shared_ptr<vecmem::memory_resource> mr, const toy_det_config_t &cfg) {
  auto [det, names] = detray::build_toy_detector<algebra_t>(*mr, cfg);

  py::object detector = py::cast(std::move(det));
  py::detail::keep_alive_impl(detector, py::cast(std::move(mr)));

  return {std::move(detector), std::move(names)};
}

}  // namespace

PYBIND11_MODULE(detray_tests_python_bindings, m) {
  m.doc() = "Detray tests bindings";

  py::class_<toy_det_config_t>(m, "ToyDetectorConfig")
      .def(py::init<>())
      .def_property(
          "n_brl_layers",
          [](const toy_det_config_t &c) { return c.n_brl_layers(); },
          [](toy_det_config_t &c, unsigned int n) { c.n_brl_layers(n); },
          "Number of barrel layers [0-4]")
      .def_property(
          "n_edc_layers",
          [](const toy_det_config_t &c) { return c.n_edc_layers(); },
          [](toy_det_config_t &c, unsigned int n) { c.n_edc_layers(n); },
          "Number of endcap layers on either side [0-7]")
      .def_property(
          "use_material_maps",
          [](const toy_det_config_t &c) { return c.use_material_maps(); },
          [](toy_det_config_t &c, bool b) { c.use_material_maps(b); },
          "Put material maps on the portals instead of homogeneous material "
          "on the modules")
      .def("__repr__", &to_string<toy_det_config_t>);

  m.def("build_toy_detector", &build_toy_detector, py::arg("memory_resource"),
        py::arg("config"),
        "Build a toy detector using the given memory resource, as configured "
        "by a ToyDetectorConfig. The returned detector keeps the memory "
        "resource alive for as long as it is used");

  py::class_<uniform_track_generator_config_t>(m, "UniformTrackGeneratorConfig")
      .def(py::init<>())
      .def_property(
          "seed",
          [](const uniform_track_generator_config_t &c) { return c.seed(); },
          [](uniform_track_generator_config_t &c, std::uint64_t s) {
            c.seed(s);
          },
          "Monte-Carlo seed")
      .def(
          "n_tracks",
          [](const uniform_track_generator_config_t &c) {
            return c.n_tracks();
          },
          "Total number of tracks")
      .def_property(
          "phi_range",
          [](const uniform_track_generator_config_t &c) {
            const auto r = c.phi_range();
            return std::make_pair(r[0], r[1]);
          },
          [](uniform_track_generator_config_t &c,
             const std::pair<scalar_t, scalar_t> &r) {
            c.phi_range(r.first, r.second);
          },
          "Phi range (min, max) (native rad)")
      .def_property(
          "theta_range",
          [](const uniform_track_generator_config_t &c) {
            const auto r = c.theta_range();
            return std::make_pair(r[0], r[1]);
          },
          [](uniform_track_generator_config_t &c,
             const std::pair<scalar_t, scalar_t> &r) {
            c.theta_range(r.first, r.second);
          },
          "Theta range (min, max) (native rad)")
      .def_property(
          "eta_range",
          [](const uniform_track_generator_config_t &c) {
            const auto r = c.eta_range();
            return std::make_pair(r[0], r[1]);
          },
          [](uniform_track_generator_config_t &c,
             const std::pair<scalar_t, scalar_t> &r) {
            c.eta_range(r.first, r.second);
          },
          "Eta range (min, max)")
      .def_property(
          "phi_steps",
          [](const uniform_track_generator_config_t &c) {
            return c.phi_steps();
          },
          [](uniform_track_generator_config_t &c, std::size_t n) {
            c.phi_steps(n);
          },
          "Number of phi steps")
      .def_property(
          "theta_steps",
          [](const uniform_track_generator_config_t &c) {
            return c.theta_steps();
          },
          [](uniform_track_generator_config_t &c, std::size_t n) {
            c.theta_steps(n);
          },
          "Number of theta steps")
      .def_property(
          "eta_steps",
          [](const uniform_track_generator_config_t &c) {
            return c.eta_steps();
          },
          [](uniform_track_generator_config_t &c, std::size_t n) {
            c.eta_steps(n);
          },
          "Number of eta steps")
      .def_property(
          "uniform_eta",
          [](const uniform_track_generator_config_t &c) {
            return c.uniform_eta();
          },
          [](uniform_track_generator_config_t &c, bool b) { c.uniform_eta(b); },
          "Whether to step uniformly in eta")
      .def_property(
          "origin",
          [](const uniform_track_generator_config_t &c) {
            const auto &o = c.origin();
            return std::make_tuple(o[0], o[1], o[2]);
          },
          [](uniform_track_generator_config_t &c,
             const std::tuple<scalar_t, scalar_t, scalar_t> &o) {
            c.origin(std::get<0>(o), std::get<1>(o), std::get<2>(o));
          },
          "Track origin")
      .def(
          "p_T",
          [](uniform_track_generator_config_t &c, scalar_t p)
              -> uniform_track_generator_config_t & { return c.p_T(p); },
          py::arg("p"), py::return_value_policy::reference,
          "Set the transverse momentum magnitude")
      .def(
          "p_tot",
          [](uniform_track_generator_config_t &c, scalar_t p)
              -> uniform_track_generator_config_t & { return c.p_tot(p); },
          py::arg("p"), py::return_value_policy::reference,
          "Set the total momentum magnitude")
      .def(
          "mom_range",
          [](const uniform_track_generator_config_t &c) {
            const auto r = c.mom_range();
            return std::make_pair(r[0], r[1]);
          },
          "Momentum range")
      .def_property(
          "randomize_charge",
          [](const uniform_track_generator_config_t &c) {
            return c.randomize_charge();
          },
          [](uniform_track_generator_config_t &c, bool b) {
            c.randomize_charge(b);
          },
          "Randomly flip the charge sign")
      .def_property(
          "time",
          [](const uniform_track_generator_config_t &c) { return c.time(); },
          [](uniform_track_generator_config_t &c, scalar_t t) { c.time(t); },
          "Track time")
      .def_property(
          "charge",
          [](const uniform_track_generator_config_t &c) { return c.charge(); },
          [](uniform_track_generator_config_t &c, scalar_t q) { c.charge(q); },
          "Track charge")
      .def(
          "is_pT",
          [](const uniform_track_generator_config_t &c) { return c.is_pT(); },
          "Whether the momentum magnitude is interpreted as transverse")
      .def("__repr__", &to_string<uniform_track_generator_config_t>);

  py::class_<test_config_t>(m, "TestConfig")
      .def_property(
          "tol", [](const test_config_t &c) { return c.tol(); },
          [](test_config_t &c, scalar_t t) { c.tol(t); },
          "Tolerance to compare two floating point values")
      .def_property(
          "propagation",
          [](test_config_t &c) -> propagation_config_t & {
            return c.propagation();
          },
          [](test_config_t &c, const propagation_config_t &v) {
            c.propagation() = v;
          },
          py::return_value_policy::reference_internal,
          "Propagation configuration");

  py::class_<material_validation_config_t, test_config_t>(
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
          "material_file",
          [](const material_validation_config_t &c) {
            return c.material_file();
          },
          [](material_validation_config_t &c, const std::string &f) {
            c.material_file(f);
          },
          "Name of the output file with the navigation material traces")
      .def_property(
          "n_tracks",
          [](const material_validation_config_t &c) { return c.n_tracks(); },
          [](material_validation_config_t &c, std::size_t n) { c.n_tracks(n); },
          "Maximal number of test tracks to run")
      .def_property(
          "relative_error",
          [](const material_validation_config_t &c) {
            return c.relative_error();
          },
          [](material_validation_config_t &c, scalar_t re) {
            c.relative_error(re);
          },
          "Allowed relative discrepancy between truth and navigation material")
      .def_property(
          "step_contingency",
          [](const material_validation_config_t &c) {
            return c.step_contingency();
          },
          [](material_validation_config_t &c, std::size_t n) {
            c.step_contingency(n);
          },
          "Additional material steps per track on top of the steps found by "
          "the material scan (device buffer capacity)")
      .def("__repr__", &to_string<material_validation_config_t>);

  py::class_<material_scan_config_t, test_config_t>(m, "MaterialScanConfig")
      .def(py::init<>())
      .def_property(
          "name", [](const material_scan_config_t &c) { return c.name(); },
          [](material_scan_config_t &c, const std::string &n) { c.name(n); },
          "Name of the test")
      .def_property(
          "material_file",
          [](const material_scan_config_t &c) { return c.material_file(); },
          [](material_scan_config_t &c, const std::string &f) {
            c.material_file(f);
          },
          "Name of the output file with the material traces")
      .def_property(
          "overlaps_removal",
          [](const material_scan_config_t &c) { return c.overlaps_removal(); },
          [](material_scan_config_t &c, bool o) { c.overlaps_removal(o); },
          "Perform overlaps removal")
      .def_property(
          "overlaps_tol",
          [](const material_scan_config_t &c) { return c.overlaps_tol(); },
          [](material_scan_config_t &c, scalar_t t) { c.overlaps_tol(t); },
          "Tolerance for considering surfaces to be overlapping")
      .def_property(
          "truth_from_file",
          [](const material_scan_config_t &c) { return c.truth_from_file(); },
          [](material_scan_config_t &c, bool b) { c.truth_from_file(b); },
          "Read the truth material from the material file instead of "
          "running the ray scan")
      .def_property(
          "track_generator",
          [](material_scan_config_t &c) -> uniform_track_generator_config_t & {
            return c.track_generator();
          },
          [](material_scan_config_t &c,
             const uniform_track_generator_config_t &v) {
            c.track_generator() = v;
          },
          py::return_value_policy::reference_internal,
          "Track generator configuration")
      .def("__repr__", &to_string<material_scan_config_t>);

  py::class_<material_validation_result_t>(m, "MaterialValidationResult")
      .def_readonly("n_tracks", &material_validation_result_t::n_tracks,
                    "Number of tracks that were compared")
      .def_readonly("n_failed", &material_validation_result_t::n_failed,
                    "Number of tracks with a material discrepancy")
      .def_property_readonly("passed", &material_validation_result_t::passed,
                             "No material discrepancy was found")
      .def("__bool__", &material_validation_result_t::passed)
      .def("__repr__", [](const material_validation_result_t &r) {
        return "MaterialValidationResult(n_tracks=" +
               std::to_string(r.n_tracks) +
               ", n_failed=" + std::to_string(r.n_failed) + ")";
      });

  m.def("run_material_validation", &run_material_validation,
        py::arg("detector"), py::arg("names"), py::arg("scan_config"),
        py::arg("validation_config"),
        "Run the CPU material validation for the given detector: Run the "
        "material scan (or read the truth material from file, if configured), "
        "then trace the material during navigation and compare it to the "
        "truth. The truth and the navigation material are written to csv "
        "files.");
}
