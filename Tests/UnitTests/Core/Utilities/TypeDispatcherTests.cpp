// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/TypeDispatcher.hpp"

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
  dispatcher.registerFunction(processA);  // Auto-detected as DerivedA
  dispatcher.registerFunction(processB);  // Auto-detected as DerivedB

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

BOOST_AUTO_TEST_CASE(ExplicitTemplateParameterRegistration) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;

  // Explicit template parameters can still be used if needed
  dispatcher.registerFunction<DerivedA>(processA);
  dispatcher.registerFunction<DerivedB>(processB);

  DerivedA objA;
  DerivedB objB;

  BOOST_CHECK_EQUAL(dispatcher(objA, "explicit_"), "explicit_A:42");
  BOOST_CHECK_EQUAL(dispatcher(objB, "explicit_"), "explicit_B:3.140000");
}

BOOST_AUTO_TEST_CASE(HasFunctionCheck) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processA);  // Auto-detected as DerivedA
  dispatcher.registerFunction(processB);  // Auto-detected as DerivedB

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
  dispatcher.registerFunction(processA);  // Auto-detected as DerivedA
  dispatcher.registerFunction(processB);  // Auto-detected as DerivedB

  DerivedC objC;

  BOOST_CHECK_THROW(dispatcher(objC, "test_"), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ClearAndSize) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processA);  // Auto-detected as DerivedA
  dispatcher.registerFunction(processB);  // Auto-detected as DerivedB

  BOOST_CHECK_EQUAL(dispatcher.size(), 2u);

  dispatcher.clear();
  BOOST_CHECK_EQUAL(dispatcher.size(), 0u);

  DerivedA objA;
  BOOST_CHECK(!dispatcher.hasFunction(objA));
}

BOOST_AUTO_TEST_CASE(PolymorphicDispatch) {
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher;
  dispatcher.registerFunction(processA);  // Auto-detected as DerivedA
  dispatcher.registerFunction(processB);  // Auto-detected as DerivedB

  std::unique_ptr<BaseClass> objA = std::make_unique<DerivedA>();
  std::unique_ptr<BaseClass> objB = std::make_unique<DerivedB>();

  BOOST_CHECK_EQUAL(dispatcher(*objA, "poly_"), "poly_A:42");
  BOOST_CHECK_EQUAL(dispatcher(*objB, "poly_"), "poly_B:3.140000");
}

// Test function for multiple arguments
std::string processMultiArg(const DerivedA& obj, int i, double d,
                            const std::string& s) {
  return std::format("A:{}:{}:{:.6f}:{}", obj.getValue(), i, d, s);
}

BOOST_AUTO_TEST_CASE(MultipleArguments) {
  TypeDispatcher<BaseClass, std::string(int, double, const std::string&)>
      multiArgDispatcher;

  multiArgDispatcher.registerFunction(processMultiArg);

  DerivedA objA;
  BOOST_CHECK_EQUAL(multiArgDispatcher(objA, 10, 2.5, "hello"),
                    "A:42:10:2.500000:hello");
}

// Test function for void return type
void processVoid(const DerivedA& obj, std::string& result) {
  result = "Modified by A: " + std::to_string(obj.getValue());
}

BOOST_AUTO_TEST_CASE(VoidReturnType) {
  TypeDispatcher<BaseClass, void(std::string&)> voidDispatcher;

  voidDispatcher.registerFunction(processVoid);

  DerivedA objA;
  std::string result;
  voidDispatcher(objA, result);
  BOOST_CHECK_EQUAL(result, "Modified by A: 42");
}

// Test function for error message testing
std::string processErrorTest(const DerivedA& /*arg*/) {
  return std::string("A");
}

BOOST_AUTO_TEST_CASE(ImprovedErrorMessages) {
  TypeDispatcher<BaseClass, std::string()> dispatcher;

  // Register a function for DerivedA but try to call with DerivedC
  dispatcher.registerFunction(processErrorTest);

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

  // Test functions for inheritance scenarios
  auto handleDerivedA = [](const DerivedA&) -> std::string {
    return "handled by DerivedA function";
  };
  std::string (*funcA)(const DerivedA&) = handleDerivedA;

  // Scenario 1: Register only for intermediate class DerivedA
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction(funcA);

    // All objects that inherit from DerivedA should work
    BOOST_CHECK_EQUAL(dispatcher(objA), "handled by DerivedA function");
    BOOST_CHECK_EQUAL(dispatcher(objAA), "handled by DerivedA function");
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "handled by DerivedA function");

    // Object from different branch should fail
    BOOST_CHECK_THROW(dispatcher(objC), std::runtime_error);
  }

  auto handleBaseClass = [](const BaseClass&) -> std::string {
    return "handled by BaseClass function";
  };
  std::string (*funcBase)(const BaseClass&) = handleBaseClass;

  // Scenario 2: Register for base class - handles everything
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction(funcBase);

    // All objects should work since they all inherit from BaseClass
    BOOST_CHECK_EQUAL(dispatcher(objA), "handled by BaseClass function");
    BOOST_CHECK_EQUAL(dispatcher(objAA), "handled by BaseClass function");
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "handled by BaseClass function");
    BOOST_CHECK_EQUAL(dispatcher(objC), "handled by BaseClass function");
  }

  auto handleDerivedAAA = [](const DerivedAAA&) -> std::string {
    return "handled by DerivedAAA function";
  };
  std::string (*funcAAA)(const DerivedAAA&) = handleDerivedAAA;

  // Scenario 3: Register for leaf class - only handles exact type
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction(funcAAA);

    // Only DerivedAAA should work
    BOOST_CHECK_THROW(dispatcher(objA), std::runtime_error);
    BOOST_CHECK_THROW(dispatcher(objAA), std::runtime_error);
    BOOST_CHECK_EQUAL(dispatcher(objAAA), "handled by DerivedAAA function");
    BOOST_CHECK_THROW(dispatcher(objC), std::runtime_error);
  }

  auto handleABranch = [](const DerivedA&) -> std::string {
    return "A branch";
  };
  auto handleCBranch = [](const DerivedC&) -> std::string {
    return "C branch";
  };
  std::string (*funcABranch)(const DerivedA&) = handleABranch;
  std::string (*funcCBranch)(const DerivedC&) = handleCBranch;

  // Scenario 4: Register for multiple branches - each handles its subtree
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    dispatcher.registerFunction(funcABranch);
    dispatcher.registerFunction(funcCBranch);

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
  auto lambdaA = [](const DerivedA&) { return std::string("A function"); };
  auto lambdaAA = [](const DerivedAA&) { return std::string("AA function"); };
  auto lambdaB = [](const DerivedB&) { return std::string("B function"); };

  std::string (*funcA)(const DerivedA&) = lambdaA;
  std::string (*funcAA)(const DerivedAA&) = lambdaAA;
  std::string (*funcB)(const DerivedB&) = lambdaB;

  // Test with auto-detected types using function pointers
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register for DerivedA using auto-detection
    dispatcher.registerFunction(funcA);  // DerivedA auto-detected

    // This should throw because DerivedAA is default constructible and would be
    // handled by DerivedA
    BOOST_CHECK_THROW(dispatcher.registerFunction(funcAA), std::runtime_error);

    // But registering for a different branch should work fine
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction(funcB));
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

  // Helper functions for non-default constructible test
  auto lambdaNonDefault = [](const NonDefaultConstructible&) {
    return std::string("Non-default function");
  };
  auto lambdaDerived = [](const DerivedFromNonDefault&) {
    return std::string("Derived function");
  };
  std::string (*funcNonDefault)(const NonDefaultConstructible&) =
      lambdaNonDefault;
  std::string (*funcDerived)(const DerivedFromNonDefault&) = lambdaDerived;

  // Test with auto-detected types using function pointers
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // Register for non-default-constructible type
    dispatcher.registerFunction(funcNonDefault);

    // Since DerivedFromNonDefault is not default constructible, we can't detect
    // the conflict at registration time
    BOOST_CHECK_NO_THROW(dispatcher.registerFunction(funcDerived));

    // But calling with a DerivedFromNonDefault object should detect the
    // ambiguity at runtime
    DerivedFromNonDefault obj(42);
    BOOST_CHECK_THROW(dispatcher(obj), std::runtime_error);
  }
}

BOOST_AUTO_TEST_CASE(DuplicateRegistrationPrevention) {
  // Helper functions for duplicate registration test
  auto lambda1 = [](const DerivedA&) { return std::string("first"); };
  auto lambda2 = [](const DerivedA&) { return std::string("second"); };
  std::string (*func1)(const DerivedA&) = lambda1;
  std::string (*func2)(const DerivedA&) = lambda2;

  // Test with auto-detected types using function pointers
  {
    TypeDispatcher<BaseClass, std::string()> dispatcher;

    // First registration should work
    dispatcher.registerFunction(func1);  // DerivedA auto-detected

    // Second registration for the same type should throw
    BOOST_CHECK_THROW(dispatcher.registerFunction(func2), std::runtime_error);

    // Original registration should still work
    DerivedA objA;
    BOOST_CHECK_EQUAL(dispatcher(objA), "first");
  }
}

BOOST_AUTO_TEST_CASE(ConstructorWithMultipleFunctionPointers) {
  // Test constructor with multiple function pointers - types auto-detected!
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher(
      processA, processB);

  DerivedA objA;
  DerivedB objB;
  DerivedC objC;

  // Both functions should be registered and working
  BOOST_CHECK_EQUAL(dispatcher(objA, "ctor_"), "ctor_A:42");
  BOOST_CHECK_EQUAL(dispatcher(objB, "ctor_"), "ctor_B:3.140000");

  // Unregistered type should throw
  BOOST_CHECK_THROW(dispatcher(objC, "ctor_"), std::runtime_error);

  // Verify both functions are registered
  BOOST_CHECK(dispatcher.hasFunction<DerivedA>());
  BOOST_CHECK(dispatcher.hasFunction<DerivedB>());
  BOOST_CHECK(!dispatcher.hasFunction<DerivedC>());
  BOOST_CHECK_EQUAL(dispatcher.size(), 2u);
}

BOOST_AUTO_TEST_CASE(ConstructorWithSingleFunction) {
  // Test constructor with just one function pointer
  TypeDispatcher<BaseClass, std::string(const std::string&)> dispatcher(
      processA);

  DerivedA objA;
  DerivedB objB;

  BOOST_CHECK_EQUAL(dispatcher(objA, "single_"), "single_A:42");
  BOOST_CHECK_THROW(dispatcher(objB, "single_"), std::runtime_error);
  BOOST_CHECK_EQUAL(dispatcher.size(), 1u);
}

BOOST_AUTO_TEST_CASE(ConstructorWithDifferentSignatures) {
  // Test constructor with functions that have different signatures
  int (*funcInt)(const DerivedA&) = [](const DerivedA& obj) {
    return obj.getValue();
  };
  int (*funcIntB)(const DerivedB&) = [](const DerivedB& obj) {
    return static_cast<int>(obj.getValue());
  };

  TypeDispatcher<BaseClass, int()> intDispatcher(funcInt, funcIntB);

  DerivedA objA;
  DerivedB objB;

  BOOST_CHECK_EQUAL(intDispatcher(objA), 42);
  BOOST_CHECK_EQUAL(intDispatcher(objB), 3);  // 3.14 truncated to int
  BOOST_CHECK_EQUAL(intDispatcher.size(), 2u);
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
