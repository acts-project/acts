// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

// Acts include(s)
#include "Acts/Plugins/Covfie/FieldConversion.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

// Pybind include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// Covfie include(s)
#include <covfie/core/backend/primitive/array.hpp>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {

namespace {

template <typename field_t>
typename field_t::view_t newView(const field_t& field) {
  typename field_t::view_t view(field);
  return view;
}

template <typename field_t>
void declareCovfieField(py::module& m, const std::string& fieldName) {
  using view_t = typename field_t::view_t;
  m.def("newView", static_cast<view_t (*)(const field_t&)>(&newView));
  py::class_<field_t, std::shared_ptr<field_t>>(m, fieldName.c_str());
  py::class_<view_t, std::shared_ptr<view_t>>(
      m, (fieldName + std::string("View")).c_str())
      .def("at", &view_t::template at<float, float, float>);
}

}  // namespace
void addCovfie(Context& ctx) {
  auto main = ctx.get("main");
  auto m = main.def_submodule("covfie", "Submodule for covfie conversion");

  using ConverterHost = Acts::CovfiePlugin::Converter<float, covfie::backend::array>;

  declareCovfieField<ConverterHost::ConstantField>(m,
                                                        "CovfieConstantField");
  declareCovfieField<ConverterHost::InterpolatedField>(
      m, "CovfieAffineLinearStridedField");


  // Only Python bindings for host.
  auto c = py::class_<ConverterHost, std::shared_ptr<ConverterHost>>(m, "CovfieConverterHost");

  c.def(py::init<>());

  c.def("covfieField",
        static_cast<ConverterHost::InterpolatedField (ConverterHost::*)(const Acts::InterpolatedMagneticField&) const>(
            &ConverterHost::covfieField));

  c.def("covfieField",
        static_cast<ConverterHost::ConstantField (ConverterHost::*)(const Acts::ConstantBField&) const>(
            &ConverterHost::covfieField));

  c.def("covfieField",
        static_cast<ConverterHost::InterpolatedField (ConverterHost::*)(const Acts::MagneticFieldProvider&,
                                                                        Acts::MagneticFieldProvider::Cache&,
                                                                        const std::vector<std::size_t>&,
                                                                        const std::vector<double>&,
                                                                        const std::vector<double>&) const>(
            &ConverterHost::covfieField));
}

}  // namespace Acts::Python
