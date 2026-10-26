// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsPython/Utilities/ProxyTether.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;

namespace ActsPython {

void addEventDataMeasurements(py::module& mex) {
  py::class_<IndexSourceLink>(mex, "IndexSourceLink")
      .def(py::init<Acts::GeometryIdentifier, Index>(), py::arg("geometryId"),
           py::arg("index"))
      .def_static(
          "fromSourceLink",
          [](Acts::SourceLink const& sl) { return sl.get<IndexSourceLink>(); })
      .def("toSourceLink",
           [](const IndexSourceLink& self) { return Acts::SourceLink(self); })
      .def("index", &IndexSourceLink::index)
      .def("geometryId", &IndexSourceLink::geometryId);

  // The measurement proxy is bound as a ProxyTether (see ProxyTether.hpp); the
  // type-erased alive-check lets both MeasurementContainer and
  // MeasurementSubset produce the same bound proxy type.
  using MeasTether = ProxyTether<ConstVariableBoundMeasurementProxy>;
  constexpr auto mcAlive = &ownerAlive<MeasurementContainer>;
  constexpr auto msAlive = &ownerAlive<MeasurementSubset>;

  // Register iterator types before the __iter__ bindings that return them.
  bindIndexIteratorTether<MeasurementContainer>(
      mex, "_MeasurementContainerIterator");
  bindIndexIteratorTether<MeasurementSubset>(mex, "_MeasurementSubsetIterator");

  using MeasProxy = ConstVariableBoundMeasurementProxy;
  py::class_<MeasTether>(mex, "ConstVariableBoundMeasurementProxy")
      .def_property_readonly(
          "geometryId", tetheredRead<MeasTether>(
                            [](const MeasProxy& s) { return s.geometryId(); }))
      .def_property_readonly(
          "size",
          tetheredRead<MeasTether>([](const MeasProxy& s) { return s.size(); }))
      .def_property_readonly("index",
                             tetheredRead<MeasTether>(
                                 [](const MeasProxy& s) { return s.index(); }))
      .def_property_readonly("fullParameters",
                             tetheredRead<MeasTether>([](const MeasProxy& s) {
                               return s.fullParameters();
                             }))
      .def_property_readonly("fullCovariance",
                             tetheredRead<MeasTether>([](const MeasProxy& s) {
                               return s.fullCovariance();
                             }))
      .def_property_readonly(
          "subspaceIndices",
          tetheredRead<MeasTether>([](const MeasProxy& self) {
            auto indices = self.subspaceHelper().indices();
            return std::vector<int>(indices.begin(), indices.end());
          }));

  auto measurementContainer =
      py::classh<MeasurementContainer>(mex, "MeasurementContainer")
          .def(py::init([]() { return MeasurementContainer(); }))
          .def("__len__", &MeasurementContainer::size)
          .def("reserve", &MeasurementContainer::reserve)
          .def(
              "emplaceMeasurement",
              [](py::object self, Acts::GeometryIdentifier geometryId,
                 const std::vector<int>& indices,
                 const std::vector<double>& par,
                 const std::vector<double>& cov) -> MeasTether {
                auto& container = self.cast<MeasurementContainer&>();
                if (indices.size() != par.size() ||
                    indices.size() != cov.size()) {
                  throw std::invalid_argument(
                      "Indices, parameters, and variances must have the same "
                      "size");
                }

                std::vector<Acts::BoundIndices> boundIndices;
                for (auto i : indices) {
                  if (i < 0 || i >= static_cast<int>(Acts::eBoundSize)) {
                    throw std::out_of_range("Subspace index out of range");
                  }
                  boundIndices.push_back(static_cast<Acts::BoundIndices>(i));
                }

                // Use existing helpers to convert the input to the measurement
                DigitizedParameters dParams;
                dParams.indices = boundIndices;
                dParams.values = par;
                dParams.variances = cov;
                return MeasTether{
                    self,
                    ConstVariableBoundMeasurementProxy{
                        createMeasurement(container, geometryId, dParams)},
                    mcAlive};
              },
              py::arg("geometryId"), py::arg("indices"), py::arg("parameters"),
              py::arg("covariance"))
          .def("__getitem__",
               [](py::object self, MeasurementContainer::Index idx) {
                 const auto& container =
                     self.cast<const MeasurementContainer&>();
                 return MeasTether{self, container.getMeasurement(idx),
                                   mcAlive};
               })
          .def("__iter__", [](py::object self) {
            return IndexIteratorTether<MeasurementContainer>{
                self, 0,
                [](const py::object& owner, MeasurementContainer& c,
                   std::size_t i) {
                  const MeasurementContainer& cc = c;
                  return py::cast(
                      MeasTether{owner, cc.getMeasurement(i),
                                 &ownerAlive<MeasurementContainer>});
                }};
          });

  WhiteBoardRegistry::registerClass(measurementContainer);

  // bind measurement subset
  auto measurementSubset =
      py::classh<MeasurementSubset>(mex, "MeasurementSubset")
          .def(py::init(
                   [](const MeasurementContainer& container,
                      const std::vector<MeasurementContainer::Index>& indices) {
                     return MeasurementSubset(container, indices);
                   }),
               py::keep_alive<0, 1>(), py::arg("container"), py::arg("indices"))
          .def("__len__",
               [](const MeasurementSubset& self) { return self.size(); })
          .def("__getitem__",
               [](py::object self, std::size_t i) {
                 const auto& subset = self.cast<const MeasurementSubset&>();
                 if (i >= subset.size()) {
                   throw py::index_error("index out of range");
                 }
                 return MeasTether{self, subset.at(i), msAlive};
               })
          .def("__iter__",
               [](py::object self) {
                 return IndexIteratorTether<MeasurementSubset>{
                     self, 0,
                     [](const py::object& owner, MeasurementSubset& c,
                        std::size_t i) {
                       const MeasurementSubset& cc = c;
                       return py::cast(MeasTether{
                           owner, cc.at(i), &ownerAlive<MeasurementSubset>});
                     }};
               })
          .def(
              "getMeasurement",
              [](py::object self, MeasurementContainer::Index idx) {
                const auto& subset = self.cast<const MeasurementSubset&>();
                return MeasTether{self, subset.getMeasurement(idx), msAlive};
              },
              py::arg("index"));

  WhiteBoardRegistry::registerClass(measurementSubset);
}

}  // namespace ActsPython
