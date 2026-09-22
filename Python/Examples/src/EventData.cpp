// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/EventData/AnyTrackProxy.hpp"
#include "Acts/EventData/AnyTrackStateProxy.hpp"
#include "Acts/EventData/SpacePointContainer.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "ActsExamples/Digitization/MeasurementCreation.hpp"
#include "ActsExamples/EventData/Cluster.hpp"
#include "ActsExamples/EventData/IndexSourceLink.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/EventData/ProtoTrack.hpp"
#include "ActsExamples/EventData/SeedSpacePointSelection.hpp"
#include "ActsExamples/EventData/Track.hpp"
#include "ActsExamples/EventData/TruthMatching.hpp"
#include "ActsPython/Utilities/ProxyTether.hpp"
#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <string>

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

// Prevent stl.h's list_caster-based type_caster<std::vector<T>> from matching,
// which would break py::cast<std::unique_ptr<T>> needed
// by WhiteBoardRegistry. The full specialization takes priority over stl.h's
// partial specialization regardless of include order.
PYBIND11_MAKE_OPAQUE(ActsExamples::ClusterContainer)
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
                     "valuesFor",
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

/// Bind the bit-flag accessors shared by TrackStateType/ConstTrackStateTypeMap/
/// MutableTrackStateTypeMap. Getters are always available; setters only exist
/// when `!Map::IsReadOnly`. isMeasurement mirrors the C++ API: it's a
/// computed flag with only a no-argument setIsMeasurement() (no inverse), so
/// it stays a read-only property with setIsMeasurement() as a plain method,
/// rather than inventing bool semantics for it.
template <typename Map>
void addTrackStateTypeFlags(py::class_<Map>& cls) {
  cls.def_property_readonly("hasMaterial", &Map::hasMaterial)
      .def_property_readonly("hasParameters", &Map::hasParameters)
      .def_property_readonly("hasMeasurement", &Map::hasMeasurement)
      .def_property_readonly("isMeasurement", &Map::isMeasurement)
      .def_property_readonly("isOutlier", &Map::isOutlier)
      .def_property_readonly("isHole", &Map::isHole)
      .def_property_readonly("isSharedHit", &Map::isSharedHit);

  if constexpr (!Map::IsReadOnly) {
    cls.def_property("hasMeasurement", &Map::hasMeasurement,
                     &Map::setHasMeasurement)
        .def("setIsMeasurement", &Map::setIsMeasurement)
        .def_property("isOutlier", &Map::isOutlier, &Map::setIsOutlier)
        .def_property("isHole", &Map::isHole, &Map::setIsHole)
        .def_property("isSharedHit", &Map::isSharedHit, &Map::setIsSharedHit);
  }
}

/// The plain Eigen type that `access(self)` evaluates to, e.g.
/// `Eigen::Map<const Acts::BoundVector>::PlainObject` is `Acts::BoundVector`.
template <typename Proxy, typename Access>
using EigenValueT = typename std::decay_t<decltype(std::declval<Access>()(
    std::declval<const Proxy&>()))>::PlainObject;

/// Bind a property whose accessor returns an Eigen::Map, always read-only
/// (for accessors with no mutable overload, e.g. TrackStateProxy's "best
/// available" parameters()/covariance()).
template <typename Cls, typename Access>
void bindEigenReadonly(Cls& cls, const char* name, Access access) {
  using Proxy = typename Cls::type;
  using Value = EigenValueT<Proxy, Access>;
  cls.def_property_readonly(
      name, [access](const Proxy& self) { return Value{access(self)}; });
}

/// Bind a property whose accessor returns an Eigen::Map, read-write when
/// `!Proxy::ReadOnly`. `access` must be a generic lambda `(auto& self) ->
/// decltype(auto)`, since these accessors are overloaded on const vs. mutable
/// with different Eigen::Map return types; the constness of `self` picks it.
template <typename Cls, typename Access>
void bindEigen(Cls& cls, const char* name, Access access) {
  using Proxy = typename Cls::type;
  using Value = EigenValueT<Proxy, Access>;
  auto get = [access](const Proxy& self) { return Value{access(self)}; };
  if constexpr (Proxy::ReadOnly) {
    cls.def_property_readonly(name, get);
  } else {
    cls.def_property(
        name, get, [access](Proxy& self, const Value& v) { access(self) = v; });
  }
}

/// Bind a property for an accessor overloaded on const (by value) vs. mutable
/// (by reference) with the same name, e.g. tipIndex, nMeasurements, chi2.
/// Read-write when `!Proxy::ReadOnly`; see @ref bindEigen for `access`.
template <typename Cls, typename Access>
void bindScalar(Cls& cls, const char* name, Access access) {
  using Proxy = typename Cls::type;
  using T = std::decay_t<decltype(access(std::declval<const Proxy&>()))>;
  auto get = [access](const Proxy& self) -> T { return access(self); };
  if constexpr (Proxy::ReadOnly) {
    cls.def_property_readonly(name, get);
  } else {
    cls.def_property(name, get,
                     [access](Proxy& self, T v) { access(self) = v; });
  }
}

/// Copy a dynamic-size Eigen vector/matrix map into std::vector(s).
template <typename Map>
std::vector<double> toVector(const Map& map) {
  return std::vector<double>(map.data(), map.data() + map.size());
}
template <typename Map>
std::vector<std::vector<double>> toNestedVector(const Map& map) {
  std::vector<std::vector<double>> out(static_cast<std::size_t>(map.rows()));
  for (Eigen::Index i = 0; i < map.rows(); ++i) {
    out[static_cast<std::size_t>(i)].assign(map.row(i).begin(),
                                            map.row(i).end());
  }
  return out;
}

/// Assign a std::vector/nested std::vector into a dynamic-size Eigen map,
/// raising ValueError on a size mismatch.
template <typename Map>
void assignVector(Map map, const std::vector<double>& values) {
  if (static_cast<Eigen::Index>(values.size()) != map.size()) {
    throw std::invalid_argument("size mismatch");
  }
  std::copy(values.begin(), values.end(), map.begin());
}
template <typename Map>
void assignNestedVector(Map map,
                        const std::vector<std::vector<double>>& values) {
  if (static_cast<Eigen::Index>(values.size()) != map.rows()) {
    throw std::invalid_argument("size mismatch");
  }
  for (Eigen::Index i = 0; i < map.rows(); ++i) {
    const auto& row = values[static_cast<std::size_t>(i)];
    if (static_cast<Eigen::Index>(row.size()) != map.cols()) {
      throw std::invalid_argument("size mismatch");
    }
    std::copy(row.begin(), row.end(), map.row(i).begin());
  }
}

/// Bind the surface shared by the concrete and Any track proxies: all four
/// (`ConstTrackProxy`/`TrackProxy`/`AnyConstTrackProxy`/`AnyMutableTrackProxy`)
/// publicly inherit the same `TrackProxyCommon` base.
template <typename Cls>
void bindTrackProxyCommon(Cls& cls) {
  using Proxy = typename Cls::type;
  cls.def_property_readonly("index", &Proxy::index)
      .def_property_readonly("referenceSurface", &Proxy::referenceSurface)
      .def_property_readonly("hasReferenceSurface", &Proxy::hasReferenceSurface)
      .def_property_readonly("particleHypothesis", &Proxy::particleHypothesis)
      .def_property_readonly("nTrackStates", &Proxy::nTrackStates)
      .def_property_readonly("isForwardLinked", &Proxy::isForwardLinked)
      .def_property_readonly("loc0", &Proxy::loc0)
      .def_property_readonly("loc1", &Proxy::loc1)
      .def_property_readonly("phi", &Proxy::phi)
      .def_property_readonly("theta", &Proxy::theta)
      .def_property_readonly("time", &Proxy::time)
      .def_property_readonly("qOverP", &Proxy::qOverP)
      .def_property_readonly("charge", &Proxy::charge)
      .def_property_readonly("absoluteMomentum", &Proxy::absoluteMomentum)
      .def_property_readonly("transverseMomentum", &Proxy::transverseMomentum)
      .def_property_readonly("direction", &Proxy::direction)
      .def_property_readonly("momentum", &Proxy::momentum)
      .def_property_readonly("fourMomentum", &Proxy::fourMomentum)
      .def(
          "hasColumn",
          [](const Proxy& self, const std::string& key) {
            return self.hasColumn(Acts::hashStringDynamic(key));
          },
          py::arg("key"))
      .def_property_readonly("trackStates",
                             py::cpp_function(
                                 [](Proxy& self) {
                                   auto range = self.trackStates();
                                   return py::make_iterator(range.begin(),
                                                            range.end());
                                 },
                                 py::keep_alive<0, 1>()))
      .def_property_readonly("trackStatesReversed",
                             py::cpp_function(
                                 [](Proxy& self) {
                                   auto range = self.trackStatesReversed();
                                   return py::make_iterator(range.begin(),
                                                            range.end());
                                 },
                                 py::keep_alive<0, 1>()))
      .def_property_readonly(
          "innermostTrackState",
          [](Proxy& self) { return self.innermostTrackState(); })
      .def_property_readonly("outermostTrackState", [](Proxy& self) {
        return self.outermostTrackState();
      });

  // tipIndex/stemIndex/nMeasurements/nHoles/nOutliers/nSharedHits/chi2/nDoF
  // are overloaded on const (by value) vs. mutable (by reference); the
  // accessor lambdas below select the overload via self's constness.
  bindScalar(cls, "tipIndex",
             [](auto& self) -> decltype(auto) { return self.tipIndex(); });
  bindScalar(cls, "stemIndex",
             [](auto& self) -> decltype(auto) { return self.stemIndex(); });
  bindScalar(cls, "nMeasurements",
             [](auto& self) -> decltype(auto) { return self.nMeasurements(); });
  bindScalar(cls, "nHoles",
             [](auto& self) -> decltype(auto) { return self.nHoles(); });
  bindScalar(cls, "nOutliers",
             [](auto& self) -> decltype(auto) { return self.nOutliers(); });
  bindScalar(cls, "nSharedHits",
             [](auto& self) -> decltype(auto) { return self.nSharedHits(); });
  bindScalar(cls, "chi2",
             [](auto& self) -> decltype(auto) { return self.chi2(); });
  bindScalar(cls, "nDoF",
             [](auto& self) -> decltype(auto) { return self.nDoF(); });
  bindEigen(cls, "parameters",
            [](auto& self) -> decltype(auto) { return self.parameters(); });
  bindEigen(cls, "covariance",
            [](auto& self) -> decltype(auto) { return self.covariance(); });
}

/// Bind the extra surface only the concrete `ConstTrackProxy`/`TrackProxy`
/// have: `AnyTrackProxy` is a value-mutate handle, not a track-building API,
/// so it lacks even the `referenceSurface`/`particleHypothesis` setters.
template <typename Cls>
void bindTrackProxyVector(Cls& cls) {
  using Proxy = typename Cls::type;
  cls.def_property_readonly("createParametersAtReference",
                            &Proxy::createParametersAtReference);

  if constexpr (!Proxy::ReadOnly) {
    cls.def_property("referenceSurface", &Proxy::referenceSurface,
                     &Proxy::setReferenceSurface)
        .def_property("particleHypothesis", &Proxy::particleHypothesis,
                      &Proxy::setParticleHypothesis)
        .def("appendTrackState", &Proxy::appendTrackState,
             py::arg("mask") = Acts::TrackStatePropMask::All,
             py::keep_alive<0, 1>())
        .def("linkForward", &Proxy::linkForward)
        .def("reverseTrackStates", &Proxy::reverseTrackStates,
             py::arg("invertJacobians") = false)
        // copyFrom*() are member function templates, so explicit
        // instantiation per source type replaces what would otherwise need a
        // lambda (a template has no single address to bind).
        .def("copyFrom", &Proxy::template copyFrom<ConstTrackProxy>)
        .def("copyFrom", &Proxy::template copyFrom<TrackProxy>)
        .def("copyFromWithoutStates",
             &Proxy::template copyFromWithoutStates<ConstTrackProxy>)
        .def("copyFromWithoutStates",
             &Proxy::template copyFromWithoutStates<TrackProxy>)
        .def("copyFromShallow",
             &Proxy::template copyFromShallow<ConstTrackProxy>)
        .def("copyFromShallow", &Proxy::template copyFromShallow<TrackProxy>);
  }
}

/// Bind the surface shared by the concrete and Any track state proxies (same
/// rationale as @ref bindTrackProxyCommon). `effectiveCalibrated`/
/// `effectiveCalibratedCovariance` are always read-only here; see
/// @ref bindTrackStateProxyVector for why.
template <typename Cls>
void bindTrackStateProxyCommon(Cls& cls) {
  using Proxy = typename Cls::type;
  cls.def_property_readonly("index", &Proxy::index)
      .def_property_readonly("hasPrevious", &Proxy::hasPrevious)
      .def_property_readonly("referenceSurface", &Proxy::referenceSurface)
      .def_property_readonly("hasReferenceSurface", &Proxy::hasReferenceSurface)
      .def_property_readonly("mask", &Proxy::getMask)
      .def_property_readonly("hasPredicted", &Proxy::hasPredicted)
      .def_property_readonly("hasFiltered", &Proxy::hasFiltered)
      .def_property_readonly("hasSmoothed", &Proxy::hasSmoothed)
      .def_property_readonly("hasJacobian", &Proxy::hasJacobian)
      .def_property_readonly("hasProjector", &Proxy::hasProjector)
      .def_property_readonly("hasUncalibratedSourceLink",
                             &Proxy::hasUncalibratedSourceLink)
      .def_property_readonly("hasCalibrated", &Proxy::hasCalibrated)
      .def_property_readonly("uncalibratedSourceLink",
                             &Proxy::getUncalibratedSourceLink)
      .def_property_readonly("calibratedSize", &Proxy::calibratedSize)
      .def_property_readonly("projectorSubspaceIndices",
                             [](const Proxy& self) {
                               auto indices = self.projectorSubspaceIndices();
                               return std::vector<int>(
                                   indices.begin(),
                                   indices.begin() + self.calibratedSize());
                             })
      .def(
          "hasColumn",
          [](const Proxy& self, const std::string& key) {
            return self.hasColumn(Acts::hashStringDynamic(key));
          },
          py::arg("key"))
      // effectiveCalibrated(+Covariance) write access is Vector-only; see
      // bindTrackStateProxyVector.
      .def_property_readonly("effectiveCalibrated",
                             [](const Proxy& self) {
                               return toVector(self.effectiveCalibrated());
                             })
      .def_property_readonly(
          "effectiveCalibratedCovariance", [](const Proxy& self) {
            return toNestedVector(self.effectiveCalibratedCovariance());
          });

  // typeFlags returns genuinely different types depending on constness
  // (ConstTrackStateTypeMap vs. MutableTrackStateTypeMap), so it can't share
  // a getter+setter pair with the generic accessors below.
  cls.def_property_readonly(
      "typeFlags",
      py::cpp_function([](const Proxy& self) { return self.typeFlags(); },
                       py::keep_alive<0, 1>()));

  bindScalar(cls, "previous",
             [](auto& self) -> decltype(auto) { return self.previous(); });
  bindScalar(cls, "chi2",
             [](auto& self) -> decltype(auto) { return self.chi2(); });
  bindScalar(cls, "pathLength",
             [](auto& self) -> decltype(auto) { return self.pathLength(); });
  bindEigen(cls, "predicted",
            [](auto& self) -> decltype(auto) { return self.predicted(); });
  bindEigen(cls, "filtered",
            [](auto& self) -> decltype(auto) { return self.filtered(); });
  bindEigen(cls, "smoothed",
            [](auto& self) -> decltype(auto) { return self.smoothed(); });
  bindEigen(cls, "predictedCovariance", [](auto& self) -> decltype(auto) {
    return self.predictedCovariance();
  });
  bindEigen(cls, "filteredCovariance", [](auto& self) -> decltype(auto) {
    return self.filteredCovariance();
  });
  bindEigen(cls, "smoothedCovariance", [](auto& self) -> decltype(auto) {
    return self.smoothedCovariance();
  });
  bindEigen(cls, "jacobian",
            [](auto& self) -> decltype(auto) { return self.jacobian(); });
  // "best available" of predicted/filtered/smoothed; no mutable overload
  // exists at all, so always read-only regardless of Proxy.
  bindEigenReadonly(cls, "parameters",
                    [](const Proxy& self) { return self.parameters(); });
  bindEigenReadonly(cls, "covariance",
                    [](const Proxy& self) { return self.covariance(); });

  if constexpr (!Proxy::ReadOnly) {
    cls.def_property("referenceSurface", &Proxy::referenceSurface,
                     &Proxy::setReferenceSurface)
        // overrides the read-only common binding above: the mutable overload
        // returns MutableTrackStateTypeMap, not ConstTrackStateTypeMap.
        .def_property_readonly(
            "typeFlags",
            py::cpp_function([](Proxy& self) { return self.typeFlags(); },
                             py::keep_alive<0, 1>()))
        .def_property(
            "uncalibratedSourceLink", &Proxy::getUncalibratedSourceLink,
            [](Proxy& self, const Acts::SourceLink& sourceLink) {
              // Some concrete setUncalibratedSourceLink overloads take
              // SourceLink&&; pybind's argument caster only ever produces an
              // lvalue, so a direct member-pointer binding doesn't compile.
              self.setUncalibratedSourceLink(Acts::SourceLink{sourceLink});
            })
        .def(
            "setProjectorSubspaceIndices",
            [](Proxy& self, const std::vector<int>& indices) {
              self.setProjectorSubspaceIndices(indices);
            },
            py::arg("indices"))
        .def(
            "allocateCalibrated",
            [](Proxy& self, std::size_t measdim) {
              self.allocateCalibrated(measdim);
            },
            py::arg("measdim"))
        .def("unset", &Proxy::unset, py::arg("mask"))
        .def("addComponents", &Proxy::addComponents, py::arg("mask"));
  }
}

/// Bind the extra surface only the concrete Vector-backed proxies have:
/// `copyFrom` (a member function template, absent from `AnyTrackStateProxy`)
/// and write access to `effectiveCalibrated`/`effectiveCalibratedCovariance`.
///
/// The latter is kept off `Any*` deliberately, not because it's unsafe today:
/// these accessors always return a plain (non-strided) `Eigen::Map`, which
/// Eigen only allows to be contiguous, so `toVector`'s raw `.data()` walk is
/// safe regardless of backend. But that's a property of today's concrete
/// backends, not something `AnyTrackStateProxy`'s interface enforces — keep
/// the write path Vector-only rather than relying on that going forward.
template <typename Cls>
void bindTrackStateProxyVector(Cls& cls) {
  using Proxy = typename Cls::type;
  if constexpr (!Proxy::ReadOnly) {
    cls.def("copyFrom", &Proxy::template copyFrom<ConstTrackStateProxy>,
            py::arg("other"), py::arg("mask") = Acts::TrackStatePropMask::All,
            py::arg("onlyAllocated") = true)
        .def("copyFrom", &Proxy::template copyFrom<TrackStateProxy>,
             py::arg("other"), py::arg("mask") = Acts::TrackStatePropMask::All,
             py::arg("onlyAllocated") = true)
        .def_property(
            "effectiveCalibrated",
            [](const Proxy& self) {
              return toVector(self.effectiveCalibrated());
            },
            [](Proxy& self, const std::vector<double>& v) {
              assignVector(self.effectiveCalibrated(), v);
            })
        .def_property(
            "effectiveCalibratedCovariance",
            [](const Proxy& self) {
              return toNestedVector(self.effectiveCalibratedCovariance());
            },
            [](Proxy& self, const std::vector<std::vector<double>>& v) {
              assignNestedVector(self.effectiveCalibratedCovariance(), v);
            });
  }
}

}  // namespace

namespace ActsPython {

void addEventData(py::module& mex) {
  // TrackStatePropMask is bound first: several TrackProxy/TrackStateProxy
  // methods below use it as a default argument, which requires the type
  // caster to already be registered when those .def() calls execute.
  py::enum_<Acts::TrackStatePropMask>(mex, "TrackStatePropMask")
      .value("None_", Acts::TrackStatePropMask::None)
      .value("Predicted", Acts::TrackStatePropMask::Predicted)
      .value("Filtered", Acts::TrackStatePropMask::Filtered)
      .value("Smoothed", Acts::TrackStatePropMask::Smoothed)
      .value("Jacobian", Acts::TrackStatePropMask::Jacobian)
      .value("Calibrated", Acts::TrackStatePropMask::Calibrated)
      .value("All", Acts::TrackStatePropMask::All)
      .def("__or__", [](Acts::TrackStatePropMask a,
                        Acts::TrackStatePropMask b) { return a | b; })
      .def("__and__", [](Acts::TrackStatePropMask a,
                         Acts::TrackStatePropMask b) { return a & b; });

  py::class_<Acts::TrackStateType> trackStateTypeFlags(mex,
                                                       "TrackStateTypeFlags");
  addTrackStateTypeFlags(trackStateTypeFlags);

  py::class_<Acts::ConstTrackStateTypeMap> constTrackStateTypeFlags(
      mex, "ConstTrackStateTypeFlags");
  addTrackStateTypeFlags(constTrackStateTypeFlags);

  py::class_<Acts::MutableTrackStateTypeMap> mutableTrackStateTypeFlags(
      mex, "MutableTrackStateTypeFlags");
  addTrackStateTypeFlags(mutableTrackStateTypeFlags);

  py::class_<ConstTrackStateProxy> constTrackState(mex, "ConstTrackStateProxy");
  bindTrackStateProxyCommon(constTrackState);
  bindTrackStateProxyVector(constTrackState);

  py::class_<TrackStateProxy> trackState(mex, "TrackStateProxy");
  bindTrackStateProxyCommon(trackState);
  bindTrackStateProxyVector(trackState);

  py::class_<ConstTrackProxy> constTrackProxy(mex, "ConstTrackProxy");
  bindTrackProxyCommon(constTrackProxy);
  bindTrackProxyVector(constTrackProxy);

  py::class_<TrackProxy> trackProxy(mex, "TrackProxy");
  bindTrackProxyCommon(trackProxy);
  bindTrackProxyVector(trackProxy);

  py::classh<Acts::AnyConstTrackStateProxy> anyConstTrackState(
      mex, "AnyConstTrackStateProxy");
  bindTrackStateProxyCommon(anyConstTrackState);
  anyConstTrackState.def(py::init<ConstTrackStateProxy&>());
  py::implicitly_convertible<ConstTrackStateProxy,
                             Acts::AnyConstTrackStateProxy>();
  // Kept for backwards compatibility with the previous (misspelled) name.
  mex.attr("AnyConstTrackState") = mex.attr("AnyConstTrackStateProxy");

  py::classh<Acts::AnyMutableTrackStateProxy> anyMutableTrackState(
      mex, "AnyMutableTrackStateProxy");
  bindTrackStateProxyCommon(anyMutableTrackState);
  anyMutableTrackState.def(py::init<TrackStateProxy&>());
  py::implicitly_convertible<TrackStateProxy,
                             Acts::AnyMutableTrackStateProxy>();

  py::classh<Acts::AnyConstTrackProxy> anyConstTrack(mex, "AnyConstTrackProxy");
  bindTrackProxyCommon(anyConstTrack);
  anyConstTrack.def(py::init<ConstTrackProxy>());
  py::implicitly_convertible<ConstTrackProxy, Acts::AnyConstTrackProxy>();

  py::classh<Acts::AnyMutableTrackProxy> anyMutableTrack(
      mex, "AnyMutableTrackProxy");
  bindTrackProxyCommon(anyMutableTrack);
  anyMutableTrack.def(py::init<TrackProxy>());
  py::implicitly_convertible<TrackProxy, Acts::AnyMutableTrackProxy>();

  // Mark a numpy array as non-writeable and return it.
  const auto readOnly = [](auto arr) {
    arr.attr("flags").attr("writeable") = py::bool_(false);
    return arr;
  };

  // Factory for zero-copy 1-D views over a plain SoA column.
  // `accessor` is called with the backend and must return a const-ref to the
  // desired std::vector member.  dtype and stride are derived automatically.
  const auto col1D = [readOnly](auto accessor) {
    return [accessor, readOnly](const py::object& self_py) {
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
          .def("getTrack",
               py::overload_cast<ConstTrackContainer::IndexType>(
                   &ConstTrackContainer::getTrack, py::const_),
               py::arg("index"), py::keep_alive<0, 1>())
          .def(
              "hasColumn",
              [](const ConstTrackContainer& self, const std::string& key) {
                return self.hasColumn(key);
              },
              py::arg("key"))
          // Build an independent, fully mutable copy of this container. The
          // source container is left untouched (unlike @c makeConst on
          // TrackContainer, which moves out of its backends), since a
          // ConstTrackContainer may be shared, e.g. via the whiteboard.
          .def("makeMutable",
               [](const ConstTrackContainer& self) {
                 return TrackContainer{
                     std::make_shared<Acts::VectorTrackContainer>(
                         self.container()),
                     std::make_shared<Acts::VectorMultiTrajectory>(
                         self.trackStateContainer())};
               })

          // Zero-copy numpy array views of the underlying SoA columns.
          // The returned arrays are read-only and keep the container alive via
          // the numpy base-object mechanism as long as any array is alive.

          // shape (N, 6), float64 — bound track parameters [loc0, loc1, phi,
          // theta, q/p, t]. Strides account for potential Eigen alignment
          // padding between entries.
          .def_property_readonly(
              "parameters",
              [readOnly](const py::object& self_py) -> py::array_t<double> {
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
              [readOnly](const py::object& self_py) -> py::array_t<double> {
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
                                 }))
          .def_property_readonly(
              "nOutliers",
              col1D([](const auto& b) -> const auto& { return b.m_nOutliers; }))
          .def_property_readonly("nSharedHits",
                                 col1D([](const auto& b) -> const auto& {
                                   return b.m_nSharedHits;
                                 }));

  WhiteBoardRegistry::registerClass(constTrackContainer);

  py::class_<TrackContainer>(mex, "TrackContainer")
      .def(py::init([]() {
        return TrackContainer{std::make_shared<Acts::VectorTrackContainer>(),
                              std::make_shared<Acts::VectorMultiTrajectory>()};
      }))
      .def(
          py::init([](const ConstTrackContainer& other) {
            return TrackContainer{
                std::make_shared<Acts::VectorTrackContainer>(other.container()),
                std::make_shared<Acts::VectorMultiTrajectory>(
                    other.trackStateContainer())};
          }),
          py::arg("other"))
      .def("__len__", &TrackContainer::size)
      .def(
          "__iter__",
          [](TrackContainer& self) {
            return py::make_iterator(self.begin(), self.end());
          },
          py::keep_alive<0, 1>())
      .def("__getitem__",
           py::overload_cast<TrackContainer::IndexType>(
               &TrackContainer::getTrack),
           py::keep_alive<0, 1>())
      .def("getTrack",
           py::overload_cast<TrackContainer::IndexType>(
               &TrackContainer::getTrack),
           py::arg("index"), py::keep_alive<0, 1>())
      .def("makeTrack", &TrackContainer::makeTrack, py::keep_alive<0, 1>())
      .def("removeTrack", &TrackContainer::removeTrack, py::arg("index"))
      .def("clear", &TrackContainer::clear)
      .def(
          "hasColumn",
          [](const TrackContainer& self, const std::string& key) {
            return self.hasColumn(key);
          },
          py::arg("key"))
      .def(
          "ensureDynamicColumns",
          [](TrackContainer& self, const ConstTrackContainer& other) {
            self.ensureDynamicColumns(other);
          },
          py::arg("other"))
      .def(
          "ensureDynamicColumns",
          [](TrackContainer& self, const TrackContainer& other) {
            self.ensureDynamicColumns(other);
          },
          py::arg("other"))
      .def("makeConst", [](TrackContainer& self) {
        return ConstTrackContainer{
            std::make_shared<Acts::ConstVectorTrackContainer>(
                std::move(self.container())),
            std::make_shared<Acts::ConstVectorMultiTrajectory>(
                std::move(self.trackStateContainer()))};
      });

  mex.attr("kTrackIndexInvalid") = Acts::kTrackIndexInvalid;

  py::bind_vector<ProtoTrack>(mex, "ProtoTrack");

  auto protoTrackContainer =
      py::bind_vector<ProtoTrackContainer, py::smart_holder>(
          mex, "ProtoTrackContainer");
  WhiteBoardRegistry::registerClass(protoTrackContainer);

  py::class_<Cluster>(mex, "Cluster")
      .def(py::init<>())
      .def_readwrite("sizeLoc0", &Cluster::sizeLoc0)
      .def_readwrite("sizeLoc1", &Cluster::sizeLoc1)
      .def_readwrite("globalPosition", &Cluster::globalPosition)
      .def_readwrite("localDirection", &Cluster::localDirection)
      .def_readwrite("lengthDirection", &Cluster::lengthDirection)
      .def_readwrite("localEta", &Cluster::localEta)
      .def_readwrite("localPhi", &Cluster::localPhi)
      .def_readwrite("globalEta", &Cluster::globalEta)
      .def_readwrite("globalPhi", &Cluster::globalPhi)
      .def_readwrite("etaAngle", &Cluster::etaAngle)
      .def_readwrite("phiAngle", &Cluster::phiAngle)
      .def("sumActivations", &Cluster::sumActivations);

  auto clusterContainer = py::bind_vector<ClusterContainer, py::smart_holder>(
      mex, "ClusterContainer");
  WhiteBoardRegistry::registerClass(clusterContainer);

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

  py::enum_<SeedSpacePointSelection>(mex, "SeedSpacePointSelection")
      .value("FirstThree", SeedSpacePointSelection::FirstThree)
      .value("InnermostTriplet", SeedSpacePointSelection::InnermostTriplet)
      .value("SpreadTriplet", SeedSpacePointSelection::SpreadTriplet);
}

}  // namespace ActsPython
