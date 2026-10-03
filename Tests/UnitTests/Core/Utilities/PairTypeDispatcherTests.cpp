// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/PairTypeDispatcher.hpp"

#include <memory>
#include <string>

using namespace Acts;

namespace ActsTests {

// Test hierarchy
class BaseClass {
 public:
  virtual ~BaseClass() = default;
  virtual std::string getType() const = 0;
};

class DerivedA : public BaseClass {
 public:
  std::string getType() const override { return "DerivedA"; }
};

class DerivedB : public BaseClass {
 public:
  std::string getType() const override { return "DerivedB"; }
};

class DerivedC : public BaseClass {
 public:
  std::string getType() const override { return "DerivedC"; }
};

// Intermediate class for inheritance tests
class DerivedAA : public DerivedA {
 public:
  std::string getType() const override { return "DerivedAA"; }
};

// Test functions
std::string processAB(const DerivedA& a, const DerivedB& b,
                      const std::string& prefix) {
  return prefix + a.getType() + "+" + b.getType();
}

std::string processAA(const DerivedA& a1, const DerivedA& a2) {
  return a1.getType() + "&" + a2.getType();
}

BOOST_AUTO_TEST_SUITE(PairTypeDispatcherTests)

BOOST_AUTO_TEST_CASE(SymmetricDispatch) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);  // (DerivedA, DerivedB) auto-detected

  DerivedA objA;
  DerivedB objB;

  // Both argument orders resolve to the same handler, and the handler always
  // receives its arguments in the registered (A, B) order.
  BOOST_CHECK_EQUAL(dispatcher(objA, objB, "p_"), "p_DerivedA+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(objB, objA, "p_"), "p_DerivedA+DerivedB");
}

BOOST_AUTO_TEST_CASE(ExplicitTemplateParameterRegistration) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction<DerivedA, DerivedB>(processAB);

  DerivedA objA;
  DerivedB objB;

  BOOST_CHECK_EQUAL(dispatcher(objA, objB, "e_"), "e_DerivedA+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(objB, objA, "e_"), "e_DerivedA+DerivedB");
}

BOOST_AUTO_TEST_CASE(SameTypePair) {
  PairTypeDispatcher<BaseClass, std::string()> dispatcher;
  dispatcher.registerFunction(processAA);

  DerivedA objA;
  BOOST_CHECK_EQUAL(dispatcher(objA, objA), "DerivedA&DerivedA");
}

BOOST_AUTO_TEST_CASE(HasFunctionCheck) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);

  DerivedA objA;
  DerivedB objB;
  DerivedC objC;

  // hasFunction on objects, both orientations
  BOOST_CHECK(dispatcher.hasFunction(objA, objB));
  BOOST_CHECK(dispatcher.hasFunction(objB, objA));
  BOOST_CHECK(!dispatcher.hasFunction(objA, objC));
  BOOST_CHECK(!dispatcher.hasFunction(objA, objA));

  // hasFunction on types is unordered
  BOOST_CHECK((dispatcher.hasFunction<DerivedA, DerivedB>()));
  BOOST_CHECK((dispatcher.hasFunction<DerivedB, DerivedA>()));
  BOOST_CHECK(!(dispatcher.hasFunction<DerivedA, DerivedC>()));
}

BOOST_AUTO_TEST_CASE(UnregisteredPairThrows) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);

  DerivedA objA;
  DerivedC objC;

  BOOST_CHECK_THROW(dispatcher(objA, objC, "p_"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(UnregisteredPairErrorMessage) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);

  DerivedA objA;
  DerivedC objC;
  try {
    dispatcher(objA, objC, "p_");
    BOOST_FAIL("Expected std::runtime_error to be thrown");
  } catch (const std::runtime_error& e) {
    std::string message = e.what();
    BOOST_CHECK(message.find("No function registered for type pair") !=
                std::string::npos);
    BOOST_CHECK(message.find("DerivedC") != std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(DuplicateRegistrationPrevention) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);

  // Registering the same unordered pair again must throw, regardless of the
  // order the derived types are given in.
  auto reversed = [](const DerivedB& b, const DerivedA& a,
                     const std::string& prefix) {
    return prefix + b.getType() + "+" + a.getType();
  };
  std::string (*funcReversed)(const DerivedB&, const DerivedA&,
                              const std::string&) = reversed;

  BOOST_CHECK_THROW(dispatcher.registerFunction(funcReversed),
                    std::runtime_error);
}

BOOST_AUTO_TEST_CASE(RegistrationTimeConflictDetection) {
  auto handleAB = [](const DerivedA& a, const DerivedB& b) {
    return a.getType() + b.getType();
  };
  auto handleAAB = [](const DerivedAA& a, const DerivedB& b) {
    return a.getType() + b.getType();
  };
  std::string (*funcAB)(const DerivedA&, const DerivedB&) = handleAB;
  std::string (*funcAAB)(const DerivedAA&, const DerivedB&) = handleAAB;

  PairTypeDispatcher<BaseClass, std::string()> dispatcher;
  dispatcher.registerFunction(funcAB);

  // (DerivedAA, DerivedB) is default constructible and would be claimed by the
  // (DerivedA, DerivedB) handler -> conflict detected at registration time.
  BOOST_CHECK_THROW(dispatcher.registerFunction(funcAAB), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(AmbiguousDispatchThrows) {
  // Two handlers whose subtype overlap cannot be detected at registration
  // time (non-default-constructible), but a concrete call is ambiguous.
  class Nd : public BaseClass {
   public:
    explicit Nd(int) {}
    std::string getType() const override { return "Nd"; }
  };
  class NdSub : public Nd {
   public:
    explicit NdSub(int v) : Nd(v) {}
    std::string getType() const override { return "NdSub"; }
  };

  auto handleNdB = [](const Nd&, const DerivedB&) { return std::string("NdB"); };
  auto handleSubB = [](const NdSub&, const DerivedB&) {
    return std::string("SubB");
  };
  std::string (*funcNdB)(const Nd&, const DerivedB&) = handleNdB;
  std::string (*funcSubB)(const NdSub&, const DerivedB&) = handleSubB;

  PairTypeDispatcher<BaseClass, std::string()> dispatcher;
  dispatcher.registerFunction(funcNdB);
  BOOST_CHECK_NO_THROW(dispatcher.registerFunction(funcSubB));

  NdSub obj(1);
  DerivedB objB;
  // Both (Nd, B) and (NdSub, B) can handle (NdSub, B).
  BOOST_CHECK_THROW(dispatcher(obj, objB), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(PolymorphicDispatch) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);

  std::unique_ptr<BaseClass> objA = std::make_unique<DerivedA>();
  std::unique_ptr<BaseClass> objB = std::make_unique<DerivedB>();

  BOOST_CHECK_EQUAL(dispatcher(*objA, *objB, "poly_"), "poly_DerivedA+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(*objB, *objA, "poly_"), "poly_DerivedA+DerivedB");
}

BOOST_AUTO_TEST_CASE(InheritanceMatching) {
  // A handler registered for a base subtype should claim derived objects.
  auto handleAB = [](const DerivedA& a, const DerivedB& b) {
    return a.getType() + "+" + b.getType();
  };
  std::string (*funcAB)(const DerivedA&, const DerivedB&) = handleAB;

  PairTypeDispatcher<BaseClass, std::string()> dispatcher;
  dispatcher.registerFunction(funcAB);

  DerivedAA objAA;  // is-a DerivedA
  DerivedB objB;

  // getType() reports the concrete type, dispatch matches via the base subtype.
  BOOST_CHECK_EQUAL(dispatcher(objAA, objB), "DerivedAA+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(objB, objAA), "DerivedAA+DerivedB");
}

BOOST_AUTO_TEST_CASE(ClearAndSize) {
  PairTypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processAB);
  BOOST_CHECK_EQUAL(dispatcher.size(), 1u);

  dispatcher.clear();
  BOOST_CHECK_EQUAL(dispatcher.size(), 0u);

  DerivedA objA;
  DerivedB objB;
  BOOST_CHECK(!dispatcher.hasFunction(objA, objB));
}

BOOST_AUTO_TEST_CASE(ConstructorWithMultipleFunctions) {
  auto handleAC = [](const DerivedA& a, const DerivedC& c) {
    return a.getType() + "+" + c.getType();
  };
  std::string (*funcAC)(const DerivedA&, const DerivedC&) = handleAC;
  std::string (*funcAA)(const DerivedA&, const DerivedA&) = processAA;

  PairTypeDispatcher<BaseClass, std::string()> dispatcher(funcAC, funcAA);

  DerivedA objA;
  DerivedC objC;

  BOOST_CHECK_EQUAL(dispatcher(objA, objC), "DerivedA+DerivedC");
  BOOST_CHECK_EQUAL(dispatcher(objC, objA), "DerivedA+DerivedC");
  BOOST_CHECK_EQUAL(dispatcher(objA, objA), "DerivedA&DerivedA");
  BOOST_CHECK_EQUAL(dispatcher.size(), 2u);
}

// A base-typed slot matches any object (dynamic_cast to a base always
// succeeds), so it acts as a wildcard: (BaseClass, DerivedB) behaves as
// "a DerivedB paired with anything".
BOOST_AUTO_TEST_CASE(BaseTypedSlotIsWildcard) {
  auto handle = [](const BaseClass& any, const DerivedB& b) {
    return any.getType() + "+" + b.getType();
  };
  std::string (*func)(const BaseClass&, const DerivedB&) = handle;

  PairTypeDispatcher<BaseClass, std::string()> dispatcher;
  dispatcher.registerFunction(func);

  DerivedA objA;
  DerivedB objB;
  DerivedC objC;

  // The DerivedB object always lands in the B slot and the other object fills
  // the wildcard, regardless of call order.
  BOOST_CHECK_EQUAL(dispatcher(objA, objB), "DerivedA+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(objB, objA), "DerivedA+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(objC, objB), "DerivedC+DerivedB");
  BOOST_CHECK_EQUAL(dispatcher(objB, objB), "DerivedB+DerivedB");

  // No DerivedB present in either slot -> no match.
  BOOST_CHECK_THROW(dispatcher(objA, objC), std::runtime_error);
}

// An overlap between a wildcard handler and a more specific handler cannot be
// resolved by specificity (there is no "most derived wins"). Whether the clash
// is caught at registration or only at call time depends on registration order,
// because the abstract root base cannot be sampled for conflict detection.
BOOST_AUTO_TEST_CASE(WildcardOverlapAmbiguityIsOrderDependent) {
  auto wide = [](const BaseClass& x, const DerivedB& b) {
    return "wide:" + x.getType() + "+" + b.getType();
  };
  auto spec = [](const DerivedA& a, const DerivedB& b) {
    return "spec:" + a.getType() + "+" + b.getType();
  };
  std::string (*funcWide)(const BaseClass&, const DerivedB&) = wide;
  std::string (*funcSpec)(const DerivedA&, const DerivedB&) = spec;

  DerivedA objA;
  DerivedB objB;
  DerivedC objC;

  // Wide first, then spec: registering spec materialises samples of
  // (DerivedA, DerivedB) -- both default constructible -- which the wide
  // checker claims, so the conflict is detected at registration time.
  {
    PairTypeDispatcher<BaseClass, std::string()> dispatcher;
    dispatcher.registerFunction(funcWide);
    BOOST_CHECK_THROW(dispatcher.registerFunction(funcSpec),
                      std::runtime_error);
  }

  // Spec first, then wide: registering wide would need to sample the abstract
  // BaseClass, which is impossible, so detection is skipped and registration
  // succeeds. The overlap then surfaces only at call time, and only for the
  // pairs that actually trigger it.
  {
    PairTypeDispatcher<BaseClass, std::string()> dispatcher;
    dispatcher.registerFunction(funcSpec);
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction(funcWide));

    // (objA, objB) is claimed by both handlers -> ambiguous, not resolved in
    // favour of the more specific (DerivedA, DerivedB) handler.
    BOOST_CHECK_THROW(dispatcher(objA, objB), std::runtime_error);

    // (objC, objB) is claimed only by the wildcard handler -> still works.
    BOOST_CHECK_EQUAL(dispatcher(objC, objB), "wide:DerivedC+DerivedB");
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
