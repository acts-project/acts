// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/TypeDispatcher.hpp"

#include <functional>
#include <memory>
#include <string>

namespace Acts::Test {

// Test hierarchy
class BaseClass {
 public:
  virtual ~BaseClass() = default;
  virtual std::string getType() const = 0;
};

class DerivedA : public BaseClass {
 public:
  std::string getType() const override { return "DerivedA"; }
  int getValue() const { return 42; }
};

class DerivedB : public BaseClass {
 public:
  std::string getType() const override { return "DerivedB"; }
  double getValue() const { return 3.14; }
};

class DerivedC : public BaseClass {
 public:
  std::string getType() const override { return "DerivedC"; }
};

// Extended hierarchy for inheritance testing
class DerivedAA : public DerivedA {
 public:
  std::string getType() const override { return "DerivedAA"; }
  std::string getExtraInfo() const { return "extra"; }
};

class DerivedAAA : public DerivedAA {
 public:
  std::string getType() const override { return "DerivedAAA"; }
};

// Test functions
std::string processA(const DerivedA& obj, const std::string& prefix) {
  return prefix + "A:" + std::to_string(obj.getValue());
}

std::string processB(const DerivedB& obj, const std::string& prefix) {
  return prefix + "B:" + std::to_string(obj.getValue());
}

BOOST_AUTO_TEST_SUITE(TypeDispatcherTests)

BOOST_AUTO_TEST_CASE(RegisterAndCallFunctions) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction<DerivedA>(processA);
  dispatcher.registerFunction<DerivedB>(processB);

  DerivedA objA;
  DerivedB objB;

  BOOST_CHECK_EQUAL(dispatcher(objA, "test_"), "test_A:42");
  BOOST_CHECK_EQUAL(dispatcher(objB, "test_"), "test_B:3.140000");
}

BOOST_AUTO_TEST_CASE(AutoDetectedFunctionRegistration) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;

  // Auto-detection: no need for explicit template parameters!
  dispatcher.registerFunction(processA);  // Type auto-detected as DerivedA
  dispatcher.registerFunction(processB);  // Type auto-detected as DerivedB

  DerivedA objA;
  DerivedB objB;

  BOOST_CHECK_EQUAL(dispatcher(objA, "auto_"), "auto_A:42");
  BOOST_CHECK_EQUAL(dispatcher(objB, "auto_"), "auto_B:3.140000");
}

BOOST_AUTO_TEST_CASE(HasFunctionCheck) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction<DerivedA>(processA);
  dispatcher.registerFunction<DerivedB>(processB);

  DerivedA objA;
  DerivedB objB;
  DerivedC objC;

  BOOST_CHECK(dispatcher.hasFunction(objA));
  BOOST_CHECK(dispatcher.hasFunction(objB));
  BOOST_CHECK(!dispatcher.hasFunction(objC));

  BOOST_CHECK(dispatcher.hasFunction<DerivedA>());
  BOOST_CHECK(dispatcher.hasFunction<DerivedB>());
  BOOST_CHECK(!dispatcher.hasFunction<DerivedC>());
}

BOOST_AUTO_TEST_CASE(UnregisteredTypeThrows) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction<DerivedA>(processA);
  dispatcher.registerFunction<DerivedB>(processB);

  DerivedC objC;

  BOOST_CHECK_THROW(dispatcher(objC, "test_"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(RegisterLambda) {
  TypeDispatcher<BaseClass, int()> lambdaDispatcher;

  lambdaDispatcher.registerFunction<DerivedA>(
      [](const DerivedA& obj) { return obj.getValue() * 2; });

  DerivedA objA;
  BOOST_CHECK_EQUAL(lambdaDispatcher(objA), 84);
}

BOOST_AUTO_TEST_CASE(RegisterFunctionObject) {
  struct Multiplier {
    int factor;
    explicit Multiplier(int f) : factor(f) {}

    int operator()(const DerivedA& obj) const {
      return obj.getValue() * factor;
    }
  };

  TypeDispatcher<BaseClass, int()> funcObjDispatcher;
  Multiplier mult(3);
  funcObjDispatcher.registerFunction<DerivedA>(mult);

  DerivedA objA;
  BOOST_CHECK_EQUAL(funcObjDispatcher(objA), 126);
}

BOOST_AUTO_TEST_CASE(RegisterStdFunctionByRvalue) {
  TypeDispatcher<BaseClass, int()> stdFuncDispatcher;

  std::function<int(const DerivedA&)> func = [](const DerivedA& obj) {
    return obj.getValue() + 10;
  };

  stdFuncDispatcher.registerFunction<DerivedA>(std::move(func));

  DerivedA objA;
  BOOST_CHECK_EQUAL(stdFuncDispatcher(objA), 52);
}

BOOST_AUTO_TEST_CASE(ClearAndSize) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction<DerivedA>(processA);
  dispatcher.registerFunction<DerivedB>(processB);

  BOOST_CHECK_EQUAL(dispatcher.size(), 2u);

  dispatcher.clear();
  BOOST_CHECK_EQUAL(dispatcher.size(), 0u);

  DerivedA objA;
  BOOST_CHECK(!dispatcher.hasFunction(objA));
}

BOOST_AUTO_TEST_CASE(PolymorphicDispatch) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction<DerivedA>(processA);
  dispatcher.registerFunction<DerivedB>(processB);

  std::unique_ptr<BaseClass> objA = std::make_unique<DerivedA>();
  std::unique_ptr<BaseClass> objB = std::make_unique<DerivedB>();

  BOOST_CHECK_EQUAL(dispatcher(*objA, "poly_"), "poly_A:42");
  BOOST_CHECK_EQUAL(dispatcher(*objB, "poly_"), "poly_B:3.140000");
}

BOOST_AUTO_TEST_CASE(MultipleArguments) {
  TypeDispatcher<BaseClass, std::string(int, double, const std::string&)>
      multiArgDispatcher;

  multiArgDispatcher.registerFunction<DerivedA>(
      [](const DerivedA& obj, int i, double d, const std::string& s) {
        return std::format("A:{}:{}:{:.6f}:{}", obj.getValue(), i, d, s);
      });

  DerivedA objA;
  BOOST_CHECK_EQUAL(multiArgDispatcher(objA, 10, 2.5, "hello"),
                    "A:42:10:2.500000:hello");
}

BOOST_AUTO_TEST_CASE(VoidReturnType) {
  TypeDispatcher<BaseClass, void(std::string&)> voidDispatcher;

  voidDispatcher.registerFunction<DerivedA>(
      [](const DerivedA& obj, std::string& result) {
        result = "Modified by A: " + std::to_string(obj.getValue());
      });

  DerivedA objA;
  std::string result;
  voidDispatcher(objA, result);
  BOOST_CHECK_EQUAL(result, "Modified by A: 42");
}

BOOST_AUTO_TEST_CASE(ImprovedErrorMessages) {
  TypeDispatcher<BaseClass, std::string()> dispatcher;

  // Register a function for DerivedA but try to call with DerivedC
  dispatcher.registerFunction<DerivedA>(
      [](const DerivedA&) { return std::string("A"); });

  DerivedC objC;
  try {
    dispatcher(objC);
    BOOST_FAIL("Expected std::runtime_error to be thrown");
  } catch (const std::runtime_error& e) {
    std::string message = e.what();
    // Check that the error message contains demangled type name
    BOOST_CHECK(message.find("No function registered for type:") !=
                std::string::npos);
    // The actual demangled name will depend on the compiler, but it should be
    // readable (not mangled like "N4Acts4Test8DerivedCE")
    BOOST_CHECK(message.find("DerivedC") != std::string::npos);
  }
}

BOOST_AUTO_TEST_CASE(InheritanceTreeHandling) {
  // Test objects at different levels of the hierarchy
  DerivedA objA;      // Level 1: inherits from BaseClass
  DerivedAA objAA;    // Level 2: inherits from DerivedA
  DerivedAAA objAAA;  // Level 3: inherits from DerivedAA -> DerivedA
  DerivedC objC;      // Different branch, inherits from BaseClass

  // Scenario 1: Register only for intermediate class DerivedA
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction<DerivedA>([](const DerivedA&) {
      return std::string("handled by DerivedA function");
    });

    // All objects that inherit from DerivedA should work
    BOOST_CHECK_EQUAL(dispatcher(objA), "handled by DerivedA function");
    BOOST_CHECK_EQUAL(dispatcher(objAA), "handled by DerivedA function");
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "handled by DerivedA function");

    // Object from different branch should fail
    BOOST_CHECK_THROW(dispatcher(objC), std::runtime_error);
  }

  // Scenario 2: Register for base class - handles everything
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction<BaseClass>([](const BaseClass&) {
      return std::string("handled by BaseClass function");
    });

    // All objects should work since they all inherit from BaseClass
    BOOST_CHECK_EQUAL(dispatcher(objA), "handled by BaseClass function");
    BOOST_CHECK_EQUAL(dispatcher(objAA), "handled by BaseClass function");
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "handled by BaseClass function");
    BOOST_CHECK_EQUAL(dispatcher(objC), "handled by BaseClass function");
  }

  // Scenario 3: Register for leaf class - only handles exact type
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction<DerivedAAA>([](const DerivedAAA&) {
      return std::string("handled by DerivedAAA function");
    });

    // Only DerivedAAA should work
    BOOST_CHECK_THROW(dispatcher(objA), std::runtime_error);
    BOOST_CHECK_THROW(dispatcher(objAA), std::runtime_error);
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "handled by DerivedAAA function");
    BOOST_CHECK_THROW(dispatcher(objC), std::runtime_error);
  }

  // Scenario 4: Register for multiple branches - each handles its subtree
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction<DerivedA>(
        [](const DerivedA&) { return std::string("A branch"); });

    dispatcher.registerFunction<DerivedC>(
        [](const DerivedC&) { return std::string("C branch"); });

    // DerivedA subtree uses A function
    BOOST_CHECK_EQUAL(dispatcher(objA), "A branch");
    BOOST_CHECK_EQUAL(dispatcher(objAA), "A branch");
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "A branch");

    // DerivedC uses C function
    BOOST_CHECK_EQUAL(dispatcher(objC), "C branch");
  }
}

BOOST_AUTO_TEST_CASE(RegistrationTimeConflictDetection) {
  // Helper functions for testing auto-detection
  auto funcA = [](const DerivedA&) { return std::string("A function"); };
  auto funcAA = [](const DerivedAA&) { return std::string("AA function"); };
  auto funcB = [](const DerivedB&) { return std::string("B function"); };

  // Test with explicit template parameters
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register for DerivedA
    dispatcher.registerFunction<DerivedA>(funcA);

    // This should throw because DerivedAA is default constructible and would be
    // handled by DerivedA
    BOOST_CHECK_THROW(dispatcher.registerFunction<DerivedAA>(funcAA),
                      std::runtime_error);

    // But registering for a different branch should work fine
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction<DerivedB>(funcB));
  }

  // Test with auto-detected types using function pointers
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register for DerivedA using auto-detection
    std::string (*ptrA)(const DerivedA&) = [](const DerivedA&) {
      return std::string("A function");
    };
    dispatcher.registerFunction(ptrA);  // DerivedA auto-detected

    // This should throw for the same reason - auto-detection doesn't change
    // conflict logic
    std::string (*ptrAA)(const DerivedAA&) = [](const DerivedAA&) {
      return std::string("AA function");
    };
    BOOST_CHECK_THROW(dispatcher.registerFunction(ptrAA), std::runtime_error);

    // Different branch should still work fine with auto-detection
    std::string (*ptrB)(const DerivedB&) = [](const DerivedB&) {
      return std::string("B function");
    };
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction(ptrB));
  }
}

BOOST_AUTO_TEST_CASE(NonDefaultConstructibleTypes) {
  // Test with a non-default-constructible type
  class NonDefaultConstructible : public BaseClass {
    int value;

   public:
    explicit NonDefaultConstructible(int v) : value(v) {
      static_cast<void>(value);
    }
    std::string getType() const override { return "NonDefaultConstructible"; }
  };

  class DerivedFromNonDefault : public NonDefaultConstructible {
   public:
    explicit DerivedFromNonDefault(int v) : NonDefaultConstructible(v) {}
    std::string getType() const override { return "DerivedFromNonDefault"; }
  };

  // Test with explicit template parameters
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register for non-default-constructible type
    dispatcher.registerFunction<NonDefaultConstructible>(
        [](const NonDefaultConstructible&) {
          return std::string("Non-default function");
        });

    // Since DerivedFromNonDefault is not default constructible, we can't detect
    // the conflict at registration time
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction<DerivedFromNonDefault>(
        [](const DerivedFromNonDefault&) {
          return std::string("Derived function");
        }));

    // But calling with a DerivedFromNonDefault object should detect the
    // ambiguity at runtime
    DerivedFromNonDefault obj(42);
    BOOST_CHECK_THROW(dispatcher(obj), std::runtime_error);
  }

  // Test with auto-detected types (should behave the same)
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register using auto-detection for non-default-constructible type
    std::string (*ptrBase)(const NonDefaultConstructible&) =
        [](const NonDefaultConstructible&) {
          return std::string("Non-default function");
        };
    dispatcher.registerFunction(
        ptrBase);  // NonDefaultConstructible auto-detected

    // Since DerivedFromNonDefault is not default constructible, registration
    // succeeds
    std::string (*ptrDerived)(const DerivedFromNonDefault&) =
        [](const DerivedFromNonDefault&) {
          return std::string("Derived function");
        };
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction(
        ptrDerived));  // DerivedFromNonDefault auto-detected

    // But calling should still detect the runtime ambiguity
    DerivedFromNonDefault obj(42);
    BOOST_CHECK_THROW(dispatcher(obj), std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(DuplicateRegistrationPrevention) {
  // Test with explicit template parameters
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // First registration should work
    dispatcher.registerFunction<DerivedA>(
        [](const DerivedA&) { return std::string("first"); });

    // Second registration for the same type should throw
    BOOST_CHECK_THROW(
        dispatcher.registerFunction<DerivedA>(
            [](const DerivedA&) { return std::string("second"); }),
        std::runtime_error);

    // Original registration should still work
    DerivedA objA;
    BOOST_CHECK_EQUAL(dispatcher(objA), "first");
  }

  // Test with auto-detected types
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // First registration using auto-detection should work
    std::string (*ptr1)(const DerivedA&) = [](const DerivedA&) {
      return std::string("first");
    };
    dispatcher.registerFunction(ptr1);  // DerivedA auto-detected

    // Second registration for same auto-detected type should throw
    std::string (*ptr2)(const DerivedA&) = [](const DerivedA&) {
      return std::string("second");
    };
    BOOST_CHECK_THROW(dispatcher.registerFunction(ptr2), std::runtime_error);

    // Original registration should still work
    DerivedA objA;
    BOOST_CHECK_EQUAL(dispatcher(objA), "first");
  }

  // Test mixing explicit and auto-detected (should also conflict)
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register with explicit template
    dispatcher.registerFunction<DerivedA>(
        [](const DerivedA&) { return std::string("explicit"); });

    // Try to register same type with auto-detection - should throw
    std::string (*ptr)(const DerivedA&) = [](const DerivedA&) {
      return std::string("auto");
    };
    BOOST_CHECK_THROW(dispatcher.registerFunction(ptr), std::runtime_error);

    DerivedA objA;
    BOOST_CHECK_EQUAL(dispatcher(objA), "explicit");
  }
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace Acts::Test
