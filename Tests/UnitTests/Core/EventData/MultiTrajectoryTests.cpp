// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/data/test_case.hpp>
#include <boost/test/tools/context.hpp>
#include <boost/test/unit_test.hpp>

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/ProxyAccessor.hpp"
#include "Acts/EventData/TrackContainer.hpp"
#include "Acts/EventData/TrackStatePropMask.hpp"
#include "Acts/EventData/TrackStateType.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/MultiTrajectoryTestsCommon.hpp"
#include "Acts/EventData/detail/TestTrackState.hpp"
#include "Acts/Geometry/GeometryContext.hpp"

#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <type_traits>
#include <typeinfo>
#include <utility>

namespace {

using namespace Acts;
using namespace Acts::UnitLiterals;
using namespace Acts::detail::Test;
using namespace Acts::HashedStringLiteral;
namespace bd = boost::unit_test::data;

using ParametersVector = BoundVector;
using CovarianceMatrix = BoundMatrix;
using Jacobian = BoundMatrix;

const auto gctx = GeometryContext::dangerouslyDefaultConstruct();
// fixed seed for reproducible tests
std::default_random_engine rng(31415);

struct Factory {
  using trajectory_t = VectorMultiTrajectory;
  using const_trajectory_t = ConstVectorMultiTrajectory;

  VectorMultiTrajectory create() { return {}; }
  ConstVectorMultiTrajectory createConst() { return {}; }
};

using CommonTests = MultiTrajectoryTestsCommon<Factory>;

}  // namespace

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(EventDataSuite)

BOOST_AUTO_TEST_CASE(Build) {
  CommonTests ct;
  ct.testBuild();
}

BOOST_AUTO_TEST_CASE(TrackStateTypeIsReadOnly) {
  static_assert(!TrackStateType::IsReadOnly);
  static_assert(ConstTrackStateTypeMap::IsReadOnly);
  static_assert(!MutableTrackStateTypeMap::IsReadOnly);
}

BOOST_AUTO_TEST_CASE(TrackStateTypeSharedSplitExclusive) {
  TrackStateType type;

  type.setIsSharedHit();
  BOOST_CHECK(type.hasMeasurement());
  BOOST_CHECK(type.isSharedHit());
  BOOST_CHECK(!type.isSplitHit());

  type.setIsSplitHit();
  BOOST_CHECK(type.hasMeasurement());
  BOOST_CHECK(type.isSplitHit());
  BOOST_CHECK(!type.isSharedHit());

  type.setIsSharedHit();
  BOOST_CHECK(type.isSharedHit());
  BOOST_CHECK(!type.isSplitHit());

  type.setIsSplitHit();
  type.setIsHole();
  BOOST_CHECK(type.isHole());
  BOOST_CHECK(!type.isSplitHit());
  BOOST_CHECK(!type.isSharedHit());
}

BOOST_AUTO_TEST_CASE(ConstCorrectness) {
  // make mutable
  VectorMultiTrajectory t;
  auto i0 = t.addTrackState();

  BOOST_CHECK(!IsReadOnlyMultiTrajectory<decltype(t)>::value);

  {
    VectorMultiTrajectory::TrackStateProxy tsp = t.getTrackState(i0);
    static_cast<void>(tsp);
    VectorMultiTrajectory::ConstTrackStateProxy ctsp{t.getTrackState(i0)};
    static_cast<void>(ctsp);

    tsp.predicted().setRandom();
    // const auto& tsp_const = tsp;
    // tsp_const.predicted().setRandom();
    // ctsp.predicted().setRandom();
  }

  // is this something we actually want?
  ConstVectorMultiTrajectory ct{t};
  BOOST_CHECK_EQUAL(ct.size(), t.size());

  ConstVectorMultiTrajectory ctm{std::move(t)};
  BOOST_CHECK_EQUAL(ctm.size(), ct.size());

  {
    static_assert(
        std::is_same_v<ConstVectorMultiTrajectory::ConstTrackStateProxy,
                       decltype(ctm.getTrackState(i0))>,
        "Got mutable track state proxy");
    ConstVectorMultiTrajectory::ConstTrackStateProxy ctsp =
        ctm.getTrackState(i0);
    static_cast<void>(ctsp);

    // doesn't compile:
    // ctsp.predictedCovariance().setIdentity();
  }

  // doesn't compile:
  // ct.clear();
  // ct.addTrackState();
}

BOOST_AUTO_TEST_CASE(ConstMutableRoundTrip) {
  VectorMultiTrajectory mtj;
  mtj.template addColumn<int>("extra_column");

  TestTrackState pc(rng, 2u);
  auto ts = mtj.getTrackState(
      mtj.addTrackState(TrackStatePropMask::All, kTrackIndexInvalid));
  fillTrackState<VectorMultiTrajectory>(pc, TrackStatePropMask::All, ts);
  ts.template component<int>("extra_column") = 42;

  ConstVectorMultiTrajectory cmtj{mtj};

  // convert back to mutable
  VectorMultiTrajectory roundTripped{cmtj};
  BOOST_CHECK_EQUAL(roundTripped.size(), mtj.size());

  auto orig = mtj.getTrackState(ts.index());
  auto copy = roundTripped.getTrackState(ts.index());

  BOOST_CHECK_EQUAL(copy.predicted(), orig.predicted());
  BOOST_CHECK_EQUAL(copy.predictedCovariance(), orig.predictedCovariance());
  BOOST_CHECK_EQUAL(copy.filtered(), orig.filtered());
  BOOST_CHECK_EQUAL(copy.smoothed(), orig.smoothed());
  BOOST_CHECK_EQUAL(copy.jacobian(), orig.jacobian());
  BOOST_CHECK_EQUAL(copy.chi2(), orig.chi2());
  BOOST_CHECK_EQUAL(copy.pathLength(), orig.pathLength());
  BOOST_CHECK_EQUAL(&copy.referenceSurface(), &orig.referenceSurface());
  BOOST_CHECK(copy.hasUncalibratedSourceLink());
  BOOST_CHECK_EQUAL((copy.template component<int>("extra_column")),
                    (orig.template component<int>("extra_column")));

  // the round-tripped trajectory is independently mutable
  copy.chi2() = orig.chi2() + 1.0f;
  BOOST_CHECK_NE(copy.chi2(), orig.chi2());
}

BOOST_AUTO_TEST_CASE(Clear) {
  CommonTests ct;
  ct.testClear();
}

BOOST_AUTO_TEST_CASE(ApplyWithAbort) {
  CommonTests ct;
  ct.testApplyWithAbort();
}

BOOST_AUTO_TEST_CASE(AddTrackStateWithBitMask) {
  CommonTests ct;
  ct.testAddTrackStateWithBitMask();
}

BOOST_AUTO_TEST_CASE(AddTrackStateComponents) {
  CommonTests ct;
  ct.testAddTrackStateComponents();
}

BOOST_AUTO_TEST_CASE(AddTrackStateComponentsAfterShareAndUnset) {
  CommonTests ct;
  ct.testAddTrackStateComponentsAfterShareAndUnset();
}

// assert expected "cross-talk" between trackstate proxies
BOOST_AUTO_TEST_CASE(TrackStateProxyCrossTalk) {
  CommonTests ct;
  ct.testTrackStateProxyCrossTalk(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateReassignment) {
  CommonTests ct;
  ct.testTrackStateReassignment(rng);
}

BOOST_DATA_TEST_CASE(TrackStateProxyStorage, bd::make({1u, 2u}),
                     nMeasurements) {
  CommonTests ct;
  ct.testTrackStateProxyStorage(rng, nMeasurements);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyAllocations) {
  CommonTests ct;
  ct.testTrackStateProxyAllocations(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateProxyGetMask) {
  CommonTests ct;
  ct.testTrackStateProxyGetMask();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopy) {
  CommonTests ct;
  ct.testTrackStateProxyCopy(rng);
}

BOOST_AUTO_TEST_CASE(TrackStateCopyDynamicColumns) {
  CommonTests ct;
  ct.testTrackStateCopyDynamicColumns();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyCopyDiffMTJ) {
  CommonTests ct;
  ct.testTrackStateProxyCopyDiffMTJ();
}

BOOST_AUTO_TEST_CASE(ProxyAssignment) {
  CommonTests ct;
  ct.testProxyAssignment();
}

BOOST_AUTO_TEST_CASE(CopyFromConst) {
  CommonTests ct;
  ct.testCopyFromConst();
}

BOOST_AUTO_TEST_CASE(TrackStateProxyShare) {
  CommonTests ct;
  ct.testTrackStateProxyShare(rng);
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumns) {
  CommonTests ct;
  ct.testMultiTrajectoryExtraColumns();
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryExtraColumnsRuntime) {
  CommonTests ct;
  ct.testMultiTrajectoryExtraColumnsRuntime();
}

BOOST_AUTO_TEST_CASE(MultiTrajectoryAllocateCalibratedInit) {
  CommonTests ct;
  ct.testMultiTrajectoryAllocateCalibratedInit(rng);
}

BOOST_AUTO_TEST_CASE(MemoryStats) {
  using namespace boost::histogram;
  using cat = axis::category<std::string>;

  VectorMultiTrajectory mt;

  auto stats = mt.statistics();

  std::stringstream ss;
  stats.toStream(ss);
  std::string out = ss.str();
  BOOST_CHECK(!out.empty());
  BOOST_CHECK_NE(out.find("total"), std::string::npos);

  const auto& h = stats.hist;

  auto column_axis = axis::get<cat>(h.axis(0));
  auto type_axis = axis::get<axis::category<>>(h.axis(1));

  for (int t = 0; t < type_axis.size(); t++) {
    for (int c = 0; c < column_axis.size(); c++) {
      double v = h.at(c, t);
      BOOST_CHECK_EQUAL(v, 0.0);
    }
  }

  TestTrackState pc(rng, 2u);
  auto ts = mt.makeTrackState();
  fillTrackState<VectorMultiTrajectory>(pc, TrackStatePropMask::All, ts);

  stats = mt.statistics();

  for (int t = 0; t < type_axis.size(); t++) {
    BOOST_TEST_CONTEXT((type_axis.bin(t) == 1 ? "meas" : "other"))
    for (int c = 0; c < column_axis.size(); c++) {
      std::string key = column_axis.bin(c);
      BOOST_TEST_CONTEXT("column: " << key) {
        double v = h.at(c, t);
        if (t == 0) {
          BOOST_CHECK_NE(v, 0.0);
        } else {
          BOOST_CHECK_EQUAL(v, 0.0);
        }
      }
    }
  }
}

BOOST_AUTO_TEST_CASE(Accessors) {
  VectorMultiTrajectory mtj;
  mtj.addColumn<unsigned int>("ndof");
  mtj.addColumn<double>("super_chi2");

  auto ts = mtj.makeTrackState();

  ProxyAccessor<unsigned int> ndof("ndof");
  ConstProxyAccessor<unsigned int> ndofConst("ndof");
  ProxyAccessor<double> superChi2("super_chi2");
  ConstProxyAccessor<double> superChi2Const("super_chi2");

  ndof(ts) = 65;
  BOOST_CHECK_EQUAL((ts.component<unsigned int, "ndof"_hash>()), 65);
  BOOST_CHECK_EQUAL(ndofConst(ts), 65);

  // should not compile
  // ndofConst(ts) = 66;

  superChi2(ts) = 123.45;
  BOOST_CHECK_EQUAL((ts.component<double, "super_chi2"_hash>()), 123.45);
  BOOST_CHECK_EQUAL(superChi2Const(ts), 123.45);

  // should not compile
  // superChi2Const(ts) = 66.66;
}

BOOST_AUTO_TEST_CASE(ChangeSourceLinkType) {
  VectorMultiTrajectory mtj;
  auto ts = mtj.makeTrackState();

  int value = 5;
  ts.setUncalibratedSourceLink(SourceLink{value});

  BOOST_CHECK_EQUAL(ts.getUncalibratedSourceLink().get<int>(), value);
  BOOST_CHECK_THROW(ts.getUncalibratedSourceLink().get<double>(),
                    std::bad_any_cast);

  double otherValue = 42.42;

  // this changes the stored type
  ts.setUncalibratedSourceLink(SourceLink{otherValue});
  BOOST_CHECK_EQUAL(ts.getUncalibratedSourceLink().get<double>(), otherValue);
  BOOST_CHECK_THROW(ts.getUncalibratedSourceLink().get<int>(),
                    std::bad_any_cast);
}

BOOST_AUTO_TEST_CASE(EnsureDynamicColumnsCopyOnTrackState) {
  TrackContainer source{VectorTrackContainer{}, VectorMultiTrajectory{}};
  source.trackStateContainer().template addColumn<std::uint64_t>("counter");
  source.trackStateContainer().template addColumn<std::uint8_t>("odd");

  auto sourceTrack = source.makeTrack();
  auto sourceState = sourceTrack.appendTrackState();
  sourceState.template component<std::uint64_t>("counter") = 42;
  sourceState.template component<std::uint8_t>("odd") = 1;

  TrackContainer target{VectorTrackContainer{}, VectorMultiTrajectory{}};
  target.ensureDynamicColumns(source);

  BOOST_CHECK(target.trackStateContainer().hasColumn(hashString("counter")));
  BOOST_CHECK(target.trackStateContainer().hasColumn(hashString("odd")));

  auto targetTrack = target.makeTrack();
  auto targetState = targetTrack.appendTrackState();
  targetState.copyFrom(sourceState);

  BOOST_CHECK_EQUAL(targetState.template component<std::uint64_t>("counter"),
                    sourceState.template component<std::uint64_t>("counter"));
  BOOST_CHECK_EQUAL(targetState.template component<std::uint8_t>("odd"),
                    sourceState.template component<std::uint8_t>("odd"));
}

/// Check that the backend knows the column @p key, that it maps the key back
/// to @p expName and that it reports @p expType as data type of the column.
/// The same information must be available via the MultiTrajectory frontend.
template <typename backend_t>
void checkTrackStateColumnInfo(const backend_t& backend, HashedString key,
                               std::string_view expName,
                               const std::type_info& expType) {
  BOOST_TEST_CONTEXT("Column '" << expName << "'") {
    BOOST_CHECK(backend.hasColumn(key));
    BOOST_CHECK_EQUAL(backend.columnName_impl(key), expName);
    BOOST_CHECK_EQUAL(backend.columnName(key), expName);
    const std::type_info& type = backend.columnDataType_impl(key);
    BOOST_CHECK_MESSAGE(type == expType, "Expected data type "
                                             << expType.name() << " but got "
                                             << type.name());
    BOOST_CHECK(backend.columnDataType(key) == type);
  }
}

template <typename backend_t>
void checkStaticTrackStateColumns(const backend_t& backend) {
  using ConstParameters = TrackStateTraits<eBoundSize, true>::Parameters;
  using ConstCovariance = TrackStateTraits<eBoundSize, true>::Covariance;
  using IndexType = TrackIndexType;

  checkTrackStateColumnInfo(backend, "predicted"_hash, "predicted",
                            typeid(ConstParameters));
  checkTrackStateColumnInfo(backend, "filtered"_hash, "filtered",
                            typeid(ConstParameters));
  checkTrackStateColumnInfo(backend, "smoothed"_hash, "smoothed",
                            typeid(ConstParameters));
  checkTrackStateColumnInfo(backend, "jacobian"_hash, "jacobian",
                            typeid(ConstCovariance));
  checkTrackStateColumnInfo(backend, "projector"_hash, "projector",
                            typeid(SerializedSubspaceIndices));
  checkTrackStateColumnInfo(backend, "previous"_hash, "previous",
                            typeid(IndexType));
  checkTrackStateColumnInfo(backend, "next"_hash, "next", typeid(IndexType));
  checkTrackStateColumnInfo(backend, "uncalibratedSourceLink"_hash,
                            "uncalibratedSourceLink",
                            typeid(std::optional<SourceLink>));
  checkTrackStateColumnInfo(backend, "referenceSurface"_hash,
                            "referenceSurface", typeid(const Surface*));
  checkTrackStateColumnInfo(backend, "measdim"_hash, "measdim",
                            typeid(IndexType));
  checkTrackStateColumnInfo(backend, "chi2"_hash, "chi2", typeid(float));
  checkTrackStateColumnInfo(backend, "pathLength"_hash, "pathLength",
                            typeid(double));
  checkTrackStateColumnInfo(backend, "typeFlags"_hash, "typeFlags",
                            typeid(TrackStateType::raw_type));

  // The calibrated columns have a measurement-dependent size. Only the name
  // mapping is checked here.
  BOOST_CHECK_EQUAL(backend.columnName("calibrated"_hash), "calibrated");
  BOOST_CHECK_EQUAL(backend.columnName("calibratedCov"_hash), "calibratedCov");
}

template <typename backend_t>
void checkDynamicTrackStateColumns(const backend_t& backend) {
  checkTrackStateColumnInfo(backend, hashString("col_float"), "col_float",
                            typeid(float));
  checkTrackStateColumnInfo(backend, hashString("col_bool"), "col_bool",
                            typeid(bool));
  checkTrackStateColumnInfo(backend, hashString("col_u64"), "col_u64",
                            typeid(std::uint64_t));
}

void addDynamicTrackStateColumns(VectorMultiTrajectory& mtj) {
  mtj.addColumn<float>("col_float");
  mtj.addColumn<bool>("col_bool");
  mtj.addColumn<std::uint64_t>("col_u64");
}

BOOST_AUTO_TEST_CASE(VectorMultiTrajectoryStaticColumnInfo) {
  VectorMultiTrajectory mtj{};
  checkStaticTrackStateColumns(mtj);

  // The column information must not depend on the number of track states
  auto i0 = mtj.addTrackState();
  mtj.addTrackState(TrackStatePropMask::All, i0);
  checkStaticTrackStateColumns(mtj);

  ConstVectorMultiTrajectory cmtj{std::move(mtj)};
  checkStaticTrackStateColumns(cmtj);
}

BOOST_AUTO_TEST_CASE(VectorMultiTrajectoryDynamicColumnInfo) {
  VectorMultiTrajectory mtj{};
  addDynamicTrackStateColumns(mtj);
  checkDynamicTrackStateColumns(mtj);

  // Columns added after track states have been created
  mtj.addTrackState();
  mtj.addColumn<double>("late");
  checkTrackStateColumnInfo(mtj, hashString("late"), "late", typeid(double));
  checkDynamicTrackStateColumns(mtj);

  // Name and type survive copies and conversions to the const backend
  VectorMultiTrajectory copy{mtj};
  checkDynamicTrackStateColumns(copy);
  checkTrackStateColumnInfo(copy, hashString("late"), "late", typeid(double));

  ConstVectorMultiTrajectory cmtj{mtj};
  checkDynamicTrackStateColumns(cmtj);
  checkTrackStateColumnInfo(cmtj, hashString("late"), "late", typeid(double));

  ConstVectorMultiTrajectory cmtjMoved{std::move(copy)};
  checkDynamicTrackStateColumns(cmtjMoved);

  // The dynamic columns must not shadow the static ones
  checkStaticTrackStateColumns(cmtj);
}

BOOST_AUTO_TEST_CASE(VectorMultiTrajectoryEnsureDynamicColumnInfo) {
  VectorMultiTrajectory source{};
  addDynamicTrackStateColumns(source);

  VectorMultiTrajectory target{};
  target.ensureDynamicColumns_impl(source);
  checkDynamicTrackStateColumns(target);

  // Also via the track container frontend
  TrackContainer tcSource{VectorTrackContainer{}, VectorMultiTrajectory{}};
  addDynamicTrackStateColumns(tcSource.trackStateContainer());
  TrackContainer tcTarget{VectorTrackContainer{}, VectorMultiTrajectory{}};
  tcTarget.ensureDynamicColumns(tcSource);
  checkDynamicTrackStateColumns(tcTarget.trackStateContainer());
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
