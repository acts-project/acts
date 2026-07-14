// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/SpacePointContainer.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsPython/Utilities/ProxyTether.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// Prevent stl.h's list_caster-based type_caster<std::vector<T>> from matching
// ProtoTrackContainer, which would break py::cast<std::unique_ptr<T>> needed
// by WhiteBoardRegistry. The full specialization takes priority over stl.h's
// partial specialization regardless of include order.
PYBIND11_MAKE_OPAQUE(ActsExamples::ProtoTrackContainer)
// MeasurementSimHitsMap == SimHitMeasurementsMap at the C++ level (both are
// flat_multimap<std::uint32_t, std::uint32_t>), so only one MAKE_OPAQUE is
// needed.
PYBIND11_MAKE_OPAQUE(ActsExamples::MeasurementSimHitsMap)
PYBIND11_MAKE_OPAQUE(ActsExamples::MeasurementParticlesMap)
PYBIND11_MAKE_OPAQUE(ActsExamples::ParticleMeasurementsMap)

namespace py = pybind11;

using namespace ActsExamples;

namespace {

template <typename Map>
auto bindFlatMultimap(py::module& m, const char* name) {
  // Register the iterator type before the __iter__ binding that returns it.
  const std::string iterName = std::string(name) + "Iterator";
  ActsPython::bindIteratorTether<Map>(m, iterName.c_str());

  auto cls = py::classh<Map>(m, name)
                 .def(py::init<>())
                 .def("__len__", &Map::size)
                 .def("__iter__",
                      [](py::object self) {
                        const auto& map = self.cast<const Map&>();
                        return ActsPython::IteratorTether<Map>{
                            self, py::make_iterator(map.begin(), map.end())};
                      })
                 .def("__contains__",
                      [](const Map& self, const typename Map::key_type& key) {
                        return self.find(key) != self.end();
                      })
                 .def(
                     "values_for",
                     [](const Map& self, const typename Map::key_type& key) {
                       auto [first, last] = self.equal_range(key);
                       std::vector<typename Map::mapped_type> result;
                       for (auto it = first; it != last; ++it) {
                         result.push_back(it->second);
                       }
                       return result;
                     },
                     py::arg("key"))
                 .def(
                     "insert",
                     [](Map& self, const typename Map::key_type& key,
                        const typename Map::mapped_type& value) {
                       self.emplace(key, value);
                     },
                     py::arg("key"), py::arg("value"));
  ActsPython::WhiteBoardRegistry::registerClass(cls);
  return cls;
}

template <typename T>
void bindIndexMultimapPair(py::module& m, const char* forwardName,
                           const char* inverseName) {
  // Bind inverse first so its Python type is registered before being used as
  // the return type of .invert() on the forward map.
  bindFlatMultimap<ActsExamples::InverseMultimap<T>>(m, inverseName);
  auto fwd = bindFlatMultimap<ActsExamples::IndexMultimap<T>>(m, forwardName);
  fwd.def("invert", [](const ActsExamples::IndexMultimap<T>& self) {
    return ActsExamples::invertIndexMultimap(self);
  });
}

}  // namespace

namespace ActsPython {

void addEventData(py::module& mex) {
  py::class_<Acts::TrackStateType>(mex, "TrackStateTypeFlags")
      .def_property_readonly("hasMeasurement",
                             &Acts::TrackStateType::hasMeasurement)
      .def_property_readonly("isMeasurement",
                             &Acts::TrackStateType::isMeasurement)
      .def_property_readonly("isOutlier", &Acts::TrackStateType::isOutlier)
      .def_property_readonly("isHole", &Acts::TrackStateType::isHole)
      .def_property_readonly("hasMaterial", &Acts::TrackStateType::hasMaterial)
      .def_property_readonly("isSharedHit", &Acts::TrackStateType::isSharedHit)
      .def_property_readonly("hasParameters",
                             &Acts::TrackStateType::hasParameters);

  py::class_<Acts::ConstTrackStateTypeMap>(mex, "ConstTrackStateTypeFlags")
      .def_property_readonly("hasMeasurement",
                             &Acts::ConstTrackStateTypeMap::hasMeasurement)
      .def_property_readonly("isMeasurement",
                             &Acts::ConstTrackStateTypeMap::isMeasurement)
      .def_property_readonly("isOutlier",
                             &Acts::ConstTrackStateTypeMap::isOutlier)
      .def_property_readonly("isHole", &Acts::ConstTrackStateTypeMap::isHole)
      .def_property_readonly("hasMaterial",
                             &Acts::ConstTrackStateTypeMap::hasMaterial)
      .def_property_readonly("isSharedHit",
                             &Acts::ConstTrackStateTypeMap::isSharedHit)
      .def_property_readonly("hasParameters",
                             &Acts::ConstTrackStateTypeMap::hasParameters);

  py::class_<Acts::MutableTrackStateTypeMap>(mex, "MutableTrackStateTypeFlags")
      .def_property("hasMeasurement",
                    &Acts::MutableTrackStateTypeMap::hasMeasurement,
                    [](Acts::MutableTrackStateTypeMap& self, bool value) {
                      self.setHasMeasurement(value);
                    })
      .def_property("isMeasurement",
                    &Acts::MutableTrackStateTypeMap::isMeasurement,
                    [](Acts::MutableTrackStateTypeMap& self, bool value) {
                      if (value) {
                        self.setIsMeasurement();
                      } else if (self.isMeasurement()) {
                        self.setHasMeasurement(false);
                      }
                    })
      .def_property("isOutlier", &Acts::MutableTrackStateTypeMap::isOutlier,
                    [](Acts::MutableTrackStateTypeMap& self, bool value) {
                      self.setIsOutlier(value);
                    })
      .def_property("isHole", &Acts::MutableTrackStateTypeMap::isHole,
                    [](Acts::MutableTrackStateTypeMap& self, bool value) {
                      self.setIsHole(value);
                    })
      .def_property_readonly("hasMaterial",
                             &Acts::MutableTrackStateTypeMap::hasMaterial)
      .def_property("isSharedHit", &Acts::MutableTrackStateTypeMap::isSharedHit,
                    [](Acts::MutableTrackStateTypeMap& self, bool value) {
                      self.setIsSharedHit(value);
                    })
      .def_property_readonly("hasParameters",
                             &Acts::MutableTrackStateTypeMap::hasParameters);

  py::class_<ConstTrackStateProxy>(mex, "ConstTrackStateProxy")
      .def_property_readonly(
          "typeFlags",
          py::cpp_function(
              [](const ConstTrackStateProxy& self) { return self.typeFlags(); },
              py::keep_alive<0, 1>()))
      .def_property_readonly("hasPredicted",
                             &ConstTrackStateProxy::hasPredicted)
      .def_property_readonly("hasFiltered", &ConstTrackStateProxy::hasFiltered)
      .def_property_readonly("hasSmoothed", &ConstTrackStateProxy::hasSmoothed)
      .def_property_readonly("predicted",
                             [](const ConstTrackStateProxy& self) {
                               return Acts::BoundVector{self.predicted()};
                             })
      .def_property_readonly("filtered",
                             [](const ConstTrackStateProxy& self) {
                               return Acts::BoundVector{self.filtered()};
                             })
      .def_property_readonly("smoothed",
                             [](const ConstTrackStateProxy& self) {
                               return Acts::BoundVector{self.smoothed()};
                             })
      .def_property_readonly(
          "predictedCovariance",
          [](const ConstTrackStateProxy& self) {
            return Acts::BoundMatrix{self.predictedCovariance()};
          })
      .def_property_readonly(
          "filteredCovariance",
          [](const ConstTrackStateProxy& self) {
            return Acts::BoundMatrix{self.filteredCovariance()};
          })
      .def_property_readonly(
          "smoothedCovariance",
          [](const ConstTrackStateProxy& self) {
            return Acts::BoundMatrix{self.smoothedCovariance()};
          })
      .def_property_readonly("pathLength", &ConstTrackStateProxy::pathLength);

  py::class_<TrackStateProxy>(mex, "TrackStateProxy")
      .def_property_readonly("typeFlags", py::cpp_function(
                                              [](TrackStateProxy& self) {
                                                return self.typeFlags();
                                              },
                                              py::keep_alive<0, 1>()))
      .def_property_readonly("hasPredicted", &TrackStateProxy::hasPredicted)
      .def_property_readonly("hasFiltered", &TrackStateProxy::hasFiltered)
      .def_property_readonly("hasSmoothed", &TrackStateProxy::hasSmoothed)
      .def_property_readonly("hasReferenceSurface",
                             &TrackStateProxy::hasReferenceSurface)
      .def_property(
          "referenceSurface",
          [](const TrackStateProxy& self) -> const Acts::Surface& {
            return self.referenceSurface();
          },
          [](TrackStateProxy& self, std::shared_ptr<const Acts::Surface> srf) {
            self.setReferenceSurface(std::move(srf));
          })
      .def_property_readonly("predicted",
                             [](const TrackStateProxy& self) {
                               return Acts::BoundVector{self.predicted()};
                             })
      .def_property_readonly("filtered",
                             [](const TrackStateProxy& self) {
                               return Acts::BoundVector{self.filtered()};
                             })
      .def_property_readonly("smoothed",
                             [](const TrackStateProxy& self) {
                               return Acts::BoundVector{self.smoothed()};
                             })
      .def_property_readonly(
          "predictedCovariance",
          [](const TrackStateProxy& self) {
            return Acts::BoundMatrix{self.predictedCovariance()};
          })
      .def_property_readonly(
          "filteredCovariance",
          [](const TrackStateProxy& self) {
            return Acts::BoundMatrix{self.filteredCovariance()};
          })
      .def_property_readonly(
          "smoothedCovariance",
          [](const TrackStateProxy& self) {
            return Acts::BoundMatrix{self.smoothedCovariance()};
          })
      .def_property(
          "uncalibratedSourceLink", &TrackStateProxy::getUncalibratedSourceLink,
          [](TrackStateProxy& self, const Acts::SourceLink& sourceLink) {
            self.setUncalibratedSourceLink(Acts::SourceLink{sourceLink});
          })
      .def("allocateCalibrated",
           [](TrackStateProxy& self, std::size_t measdim) {
             self.allocateCalibrated(measdim);
           })
      .def_property_readonly("pathLength", [](const TrackStateProxy& self) {
        return self.pathLength();
      });

  auto constTrackProxy =
      py::class_<ConstTrackProxy>(mex, "ConstTrackProxy")
          .def_property_readonly("index", &ConstTrackProxy::index)
          .def_property_readonly(
              "tipIndex",
              [](const ConstTrackProxy& self) { return self.tipIndex(); })
          .def_property_readonly(
              "stemIndex",
              [](const ConstTrackProxy& self) { return self.stemIndex(); })
          .def_property_readonly("referenceSurface",
                                 &ConstTrackProxy::referenceSurface)
          .def_property_readonly("hasReferenceSurface",
                                 &ConstTrackProxy::hasReferenceSurface)
          // Convert the parameters to a BoundVector so they're accessible by
          // value
          .def_property_readonly("parameters",
                                 [](const ConstTrackProxy& self) {
                                   return Acts::BoundVector{self.parameters()};
                                 })
          // Convert the covariance to a BoundMatrix so they're accessible by
          // value
          .def_property_readonly("covariance",
                                 [](const ConstTrackProxy& self) {
                                   return Acts::BoundMatrix{self.covariance()};
                                 })
          .def_property_readonly("particleHypothesis",
                                 &ConstTrackProxy::particleHypothesis)
          .def_property_readonly(
              "nMeasurements",
              [](const ConstTrackProxy& self) { return self.nMeasurements(); })
          .def_property_readonly(
              "nHoles",
              [](const ConstTrackProxy& self) { return self.nHoles(); })
          .def_property_readonly("isForwardLinked",
                                 [](const ConstTrackProxy& self) {
                                   return self.isForwardLinked();
                                 })
          .def_property_readonly("trackStatesReversed",
                                 py::cpp_function(
                                     [](const ConstTrackProxy& self) {
                                       auto range = self.trackStatesReversed();
                                       return py::make_iterator(range.begin(),
                                                                range.end());
                                     },
                                     py::keep_alive<0, 1>()))
          .def_property_readonly("trackStates",
                                 py::cpp_function(
                                     [](const ConstTrackProxy& self) {
                                       auto range = self.trackStates();
                                       return py::make_iterator(range.begin(),
                                                                range.end());
                                     },
                                     py::keep_alive<0, 1>()));

  // Mark a numpy array as non-writeable and return it.
  const auto readOnly = [](auto arr) {
    arr.attr("flags").attr("writeable") = py::bool_(false);
    return arr;
  };

  // Factory for zero-copy 1-D views over a plain SoA column.
  // `accessor` is called with the backend and must return a const-ref to the
  // desired std::vector member.  dtype and stride are derived automatically.
  const auto col1D = [&readOnly](auto accessor) {
    return [accessor, &readOnly](const py::object& self_py) {
      const auto& backend =
          self_py.cast<const ConstTrackContainer&>().container();
      const auto N = static_cast<py::ssize_t>(backend.size_impl());
      const auto& vec = accessor(backend);
      using T = typename std::decay_t<decltype(vec)>::value_type;
      if (N == 0) {
        return readOnly(py::array_t<T>(py::ssize_t{0}));
      }
      return readOnly(py::array_t<T>({N}, {static_cast<py::ssize_t>(sizeof(T))},
                                     vec.data(), self_py));
    };
  };

  auto constTrackContainer =
      py::classh<ConstTrackContainer>(mex, "ConstTrackContainer")
          .def("__len__", &ConstTrackContainer::size)
          .def(
              "__iter__",
              [](const ConstTrackContainer& self) {
                return py::make_iterator(self.begin(), self.end());
              },
              py::keep_alive<0, 1>())
          .def("__getitem__",
               py::overload_cast<ConstTrackContainer::IndexType>(
                   &ConstTrackContainer::getTrack, py::const_),
               py::keep_alive<0, 1>())

          // Zero-copy numpy array views of the underlying SoA columns.
          // The returned arrays are read-only and keep the container alive via
          // the numpy base-object mechanism as long as any array is alive.

          // shape (N, 6), float64 — bound track parameters [loc0, loc1, phi,
          // theta, q/p, t]. Strides account for potential Eigen alignment
          // padding between entries.
          .def_property_readonly(
              "parameters",
              [&readOnly](const py::object& self_py) -> py::array_t<double> {
                using CoeffsType = Acts::detail_tsp::FixedSizeTypes<
                    Acts::eBoundSize>::Coefficients;
                const auto& backend =
                    self_py.cast<const ConstTrackContainer&>().container();
                const auto N = static_cast<py::ssize_t>(backend.size_impl());
                if (N == 0) {
                  return readOnly(py::array_t<double>({N, py::ssize_t{6}}));
                }
                return readOnly(py::array_t<double>(
                    {N, py::ssize_t{6}},
                    {static_cast<py::ssize_t>(sizeof(CoeffsType)),
                     static_cast<py::ssize_t>(sizeof(double))},
                    backend.m_params[0].data(), self_py));
              })

          // shape (N, 6, 6), float64, column-major sub-matrices (Eigen
          // default). arr[k, i, j] is row i, column j of track k's covariance.
          .def_property_readonly(
              "covariance",
              [&readOnly](const py::object& self_py) -> py::array_t<double> {
                using CovType = Acts::detail_tsp::FixedSizeTypes<
                    Acts::eBoundSize>::Covariance;
                const auto& backend =
                    self_py.cast<const ConstTrackContainer&>().container();
                const auto N = static_cast<py::ssize_t>(backend.size_impl());
                if (N == 0) {
                  return readOnly(
                      py::array_t<double>({N, py::ssize_t{6}, py::ssize_t{6}}));
                }
                // Eigen column-major: row stride = 1 double, col stride = 6
                // doubles.
                constexpr py::ssize_t dbl = sizeof(double);
                return readOnly(py::array_t<double>(
                    {N, py::ssize_t{6}, py::ssize_t{6}},
                    {static_cast<py::ssize_t>(sizeof(CovType)), dbl, 6 * dbl},
                    backend.m_cov[0].data(), self_py));
              })

          .def_property_readonly(
              "tipIndex",
              col1D([](const auto& b) -> const auto& { return b.m_tipIndex; }))
          .def_property_readonly(
              "stemIndex",
              col1D([](const auto& b) -> const auto& { return b.m_stemIndex; }))
          .def_property_readonly("nMeasurements",
                                 col1D([](const auto& b) -> const auto& {
                                   return b.m_nMeasurements;
                                 }))
          .def_property_readonly(
              "nHoles",
              col1D([](const auto& b) -> const auto& { return b.m_nHoles; }))
          .def_property_readonly(
              "chi2",
              col1D([](const auto& b) -> const auto& { return b.m_chi2; }))
          .def_property_readonly("ndf", col1D([](const auto& b) -> const auto& {
                                   return b.m_ndf;
                                 }));

  WhiteBoardRegistry::registerClass(constTrackContainer);

  py::bind_vector<ProtoTrack>(mex, "ProtoTrack");

  auto protoTrackContainer =
      py::bind_vector<ProtoTrackContainer, py::smart_holder>(
          mex, "ProtoTrackContainer");
  WhiteBoardRegistry::registerClass(protoTrackContainer);

  mex.attr("kTrackIndexInvalid") = Acts::kTrackIndexInvalid;

  py::class_<IndexSourceLink>(mex, "IndexSourceLink")
      .def("FromSourceLink",
           [](Acts::SourceLink const& sl) { return sl.get<IndexSourceLink>(); })
      .def("index", &IndexSourceLink::index)
      .def("geometryId", &IndexSourceLink::geometryId);

  py::class_<TrackProxy>(mex, "TrackProxy")
      .def_property(
          "referenceSurface",
          [](const TrackProxy& self) -> const Acts::Surface& {
            return self.referenceSurface();
          },
          [](TrackProxy& self, std::shared_ptr<const Acts::Surface> srf) {
            self.setReferenceSurface(std::move(srf));
          })
      .def_property(
          "parameters",
          [](const TrackProxy& self) {
            return Acts::BoundVector{self.parameters()};
          },
          [](TrackProxy& self, const Acts::BoundVector& v) {
            self.parameters() = v;
          })
      .def_property(
          "covariance",
          [](const TrackProxy& self) {
            return Acts::BoundMatrix{self.covariance()};
          },
          [](TrackProxy& self, const Acts::BoundMatrix& m) {
            self.covariance() = m;
          })
      .def_property(
          "particleHypothesis",
          [](const TrackProxy& self) { return self.particleHypothesis(); },
          [](TrackProxy& self, const Acts::ParticleHypothesis& hyp) {
            self.setParticleHypothesis(hyp);
          })
      .def_property(
          "nMeasurements",
          [](const TrackProxy& self) -> std::uint32_t {
            return self.nMeasurements();
          },
          [](TrackProxy& self, std::uint32_t n) { self.nMeasurements() = n; })
      .def_property(
          "nHoles",
          [](const TrackProxy& self) -> std::uint32_t { return self.nHoles(); },
          [](TrackProxy& self, std::uint32_t n) { self.nHoles() = n; })
      .def_property(
          "chi2", [](const TrackProxy& self) -> float { return self.chi2(); },
          [](TrackProxy& self, float v) { self.chi2() = v; })
      .def(
          "appendTrackState",
          [](TrackProxy& self) { return self.appendTrackState(); },
          py::keep_alive<0, 1>());

  py::class_<TrackContainer>(mex, "TrackContainer")
      .def(py::init([]() {
        return TrackContainer{std::make_shared<Acts::VectorTrackContainer>(),
                              std::make_shared<Acts::VectorMultiTrajectory>()};
      }))
      .def("__len__", &TrackContainer::size)
      .def("makeTrack", &TrackContainer::makeTrack, py::keep_alive<0, 1>())
      .def("makeConst", [](TrackContainer& self) {
        return ConstTrackContainer{
            std::make_shared<Acts::ConstVectorTrackContainer>(
                std::move(self.container())),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(
                std::move(self.trackStateContainer()))};
      });

  // bind measurements
  // The measurement proxy is bound as a ProxyTether (see ProxyTether.hpp). The
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
  // MeasurementSimHitsMap and SimHitMeasurementsMap are the same C++ type
  // (flat_multimap<std::uint32_t, std::uint32_t>) because SimHitIndex == Index
  // == std::uint32_t. Bind once and alias the second name.
  auto simHitsMap =
      bindFlatMultimap<MeasurementSimHitsMap>(mex, "MeasurementSimHitsMap");
  simHitsMap.def("invert", [](const MeasurementSimHitsMap& self) {
    return invertIndexMultimap(self);
  });
  mex.attr("SimHitMeasurementsMap") = mex.attr("MeasurementSimHitsMap");

  bindIndexMultimapPair<SimBarcode>(mex, "MeasurementParticlesMap",
                                    "ParticleMeasurementsMap");
}

}  // namespace ActsPython
