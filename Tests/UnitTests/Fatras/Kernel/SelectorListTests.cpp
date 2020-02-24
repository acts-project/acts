// This file is part of the Acts project.
//
// Copyright (C) 2018-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include <string>

#include "ActsFatras/Kernel/SelectorList.hpp"

namespace {
// mock-up objects for particle and detector.

struct Object {
  int feature = 0;
  std::string name = "";
};

struct Environment {
  int feature = 1;
};

/// Selector that selects on the feature. Only acts on the object.
struct FeatureSelector {
  int select_on = 0;

  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &object) const {
    return object.feature == select_on;
  }
};

/// Selector that selects on the name. Only acts on the object.
struct NameSelector {
  std::string select_on = "";

  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &, const particle_t &object) const {
    return object.name == select_on;
  }
};

/// Selector that selects on the feature given an environment.
struct EnvironmentSelector {
  template <typename detector_t, typename particle_t>
  bool operator()(const detector_t &environment,
                  const particle_t &object) const {
    return object.feature == environment.feature;
  }
};

struct Fixture {
  FeatureSelector selectObjectFeature = {1};
  NameSelector selectObjectName = {"same_name"};
  EnvironmentSelector selectObjectEnvironmentMatch;

  Object obj = {1, "same_name"};
  Environment env = {1};
  Object objWrongFeature = {2, "same_name"};
  Environment envWrongFeature = {2};
};
}  // namespace

BOOST_AUTO_TEST_SUITE(FatrasSelectorList)

BOOST_AUTO_TEST_CASE(And1) {
  Fixture f;

  ActsFatras::SelectorListAND<FeatureSelector> select;
  select.get<FeatureSelector>() = f.selectObjectFeature;
  BOOST_TEST(select(f.env, f.obj));
  BOOST_TEST(not select(f.env, f.objWrongFeature));
  BOOST_TEST(select(f.envWrongFeature, f.obj));
  BOOST_TEST(not select(f.envWrongFeature, f.objWrongFeature));
}

BOOST_AUTO_TEST_CASE(And2) {
  Fixture f;

  ActsFatras::SelectorListAND<FeatureSelector, NameSelector> select;
  select.get<FeatureSelector>() = f.selectObjectFeature;
  select.get<NameSelector>() = f.selectObjectName;
  BOOST_TEST(select(f.env, f.obj));
  BOOST_TEST(not select(f.env, f.objWrongFeature));
  BOOST_TEST(select(f.envWrongFeature, f.obj));
  BOOST_TEST(not select(f.envWrongFeature, f.objWrongFeature));
}

BOOST_AUTO_TEST_CASE(And3) {
  Fixture f;

  ActsFatras::SelectorListAND<FeatureSelector, NameSelector,
                              EnvironmentSelector>
      select;
  select.get<FeatureSelector>() = f.selectObjectFeature;
  select.get<NameSelector>() = f.selectObjectName;
  select.get<EnvironmentSelector>() = f.selectObjectEnvironmentMatch;
  BOOST_TEST(select(f.env, f.obj));
  BOOST_TEST(not select(f.env, f.objWrongFeature));
  BOOST_TEST(not select(f.envWrongFeature, f.obj));
  BOOST_TEST(not select(f.envWrongFeature, f.objWrongFeature));
}

BOOST_AUTO_TEST_CASE(Or1) {
  Fixture f;

  ActsFatras::SelectorListOR<FeatureSelector> select;
  select.get<FeatureSelector>() = f.selectObjectFeature;
  BOOST_TEST(select(f.env, f.obj));
  BOOST_TEST(not select(f.env, f.objWrongFeature));
  BOOST_TEST(select(f.envWrongFeature, f.obj));
  BOOST_TEST(not select(f.envWrongFeature, f.objWrongFeature));
}

BOOST_AUTO_TEST_CASE(Or2) {
  Fixture f;

  ActsFatras::SelectorListOR<FeatureSelector, EnvironmentSelector> select;
  select.get<FeatureSelector>() = f.selectObjectFeature;
  select.get<EnvironmentSelector>() = f.selectObjectEnvironmentMatch;
  BOOST_TEST(select(f.env, f.obj));
  BOOST_TEST(not select(f.env, f.objWrongFeature));
  BOOST_TEST(select(f.envWrongFeature, f.obj));
  BOOST_TEST(select(f.envWrongFeature, f.objWrongFeature));
}

BOOST_AUTO_TEST_CASE(Or3) {
  Fixture f;

  ActsFatras::SelectorListOR<FeatureSelector, NameSelector, EnvironmentSelector>
      select;
  select.get<FeatureSelector>() = f.selectObjectFeature;
  select.get<NameSelector>() = f.selectObjectName;
  select.get<EnvironmentSelector>() = f.selectObjectEnvironmentMatch;
  BOOST_TEST(select(f.env, f.obj));
  BOOST_TEST(select(f.env, f.objWrongFeature));
  BOOST_TEST(select(f.envWrongFeature, f.obj));
  BOOST_TEST(select(f.envWrongFeature, f.objWrongFeature));
}

BOOST_AUTO_TEST_SUITE_END()
