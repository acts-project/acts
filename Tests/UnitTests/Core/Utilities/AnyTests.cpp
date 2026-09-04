// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <boost/test/unit_test.hpp>

#include "Acts/Utilities/Any.hpp"

#include <any>
#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <utility>

using namespace Acts;

#if defined(_ACTS_ANY_ENABLE_TRACK_ALLOCATIONS)
#define CHECK_ANY_ALLOCATIONS()                         \
  do {                                                  \
    detail::_AnyAllocationReporter::checkAllocations(); \
  } while (0)
#else
#define CHECK_ANY_ALLOCATIONS() \
  do {                          \
  } while (0)
#endif

namespace ActsTests {

BOOST_AUTO_TEST_SUITE(UtilitiesSuite)

BOOST_AUTO_TEST_CASE(AnyConstructPrimitive) {
  {
    // small type
    Any a;
    BOOST_CHECK(!a);

    int v = 5;
    a = Any{v};
    BOOST_CHECK(!!a);

    BOOST_CHECK_EQUAL(a.as<int>(), v);
    BOOST_CHECK_NE(a.as<int>(), v + 1);

    BOOST_CHECK_THROW(a.as<float>(), std::bad_any_cast);
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // type that is large
    Any a;
    BOOST_CHECK(!a);

    std::array<int, 2> v{1, 2};
    a = Any{v};
    BOOST_CHECK(!!a);

    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<decltype(v)>().begin(),
                                  a.as<decltype(v)>().end(), v.begin(),
                                  v.end());
    BOOST_CHECK_THROW(a.as<float>(), std::bad_any_cast);
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // type that is large
    Any a;
    BOOST_CHECK(!a);

    std::array<unsigned long, 5> v{1, 2, 3, 4, 5};
    a = Any{v};
    BOOST_CHECK(!!a);

    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<decltype(v)>().begin(),
                                  a.as<decltype(v)>().end(), v.begin(),
                                  v.end());
    BOOST_CHECK_THROW(a.as<float>(), std::bad_any_cast);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyAsPtr) {
  {
    // small type: correct type returns non-null pointer
    Any a{42};
    int* p = a.asPtr<int>();
    BOOST_REQUIRE_NE(p, static_cast<int*>(nullptr));
    BOOST_CHECK_EQUAL(*p, 42);

    // wrong type returns nullptr
    BOOST_CHECK_EQUAL(a.asPtr<float>(), static_cast<float*>(nullptr));
    BOOST_CHECK_EQUAL(a.asPtr<double>(), static_cast<double*>(nullptr));

    // mutation through pointer
    *p = 99;
    BOOST_CHECK_EQUAL(a.as<int>(), 99);
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // large (heap-allocated) type: correct type returns non-null pointer
    std::array<unsigned long, 5> v{10, 20, 30, 40, 50};
    Any a{v};
    auto* p = a.asPtr<std::array<unsigned long, 5>>();
    BOOST_REQUIRE_NE(p, static_cast<decltype(p)>(nullptr));
    BOOST_CHECK_EQUAL_COLLECTIONS(p->begin(), p->end(), v.begin(), v.end());

    // wrong type returns nullptr
    BOOST_CHECK_EQUAL(a.asPtr<int>(), static_cast<int*>(nullptr));
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // empty Any returns nullptr
    Any a;
    BOOST_CHECK_EQUAL(a.asPtr<int>(), static_cast<int*>(nullptr));
    BOOST_CHECK_EQUAL(a.asPtr<float>(), static_cast<float*>(nullptr));
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // const overload
    const Any a{3.14f};
    const float* p = a.asPtr<float>();
    BOOST_REQUIRE_NE(p, static_cast<const float*>(nullptr));
    BOOST_CHECK_EQUAL(*p, 3.14f);

    // wrong type on const Any
    BOOST_CHECK_EQUAL(a.asPtr<int>(), static_cast<const int*>(nullptr));
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // const overload with empty Any
    const Any a;
    BOOST_CHECK_EQUAL(a.asPtr<int>(), static_cast<const int*>(nullptr));
  }
  CHECK_ANY_ALLOCATIONS();

  {
    // const overload with heap-allocated type
    std::array<int, 64> v{};
    v.fill(7);
    const Any a{v};
    const auto* p = a.asPtr<std::array<int, 64>>();
    BOOST_REQUIRE_NE(p, static_cast<decltype(p)>(nullptr));
    for (const auto& elem : *p) {
      BOOST_CHECK_EQUAL(elem, 7);
    }
    BOOST_CHECK_EQUAL(a.asPtr<float>(), static_cast<const float*>(nullptr));
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyAssignConstructEmpty) {
  Any a;
  Any b;
  a = b;
  Any c{a};
  a = std::move(b);
  Any d{std::move(a)};

  BOOST_CHECK(!a);
  BOOST_CHECK(!b);
  BOOST_CHECK(!c);
  BOOST_CHECK(!d);

  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyConstructCustom) {
  struct A {
    int value;
    A() { value = 76; }
  };

  Any a;
  BOOST_CHECK(!a);
  a = Any{A{}};

  BOOST_CHECK(!!a);

  BOOST_CHECK_EQUAL(a.as<A>().value, 76);

  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyConstructCustomInPlace) {
  struct A {
    int value;
    explicit A(int v) { value = v; }
  };

  Any a{std::in_place_type<A>, 42};
  BOOST_CHECK(!!a);
  BOOST_CHECK_EQUAL(a.as<A>().value, 42);

  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyMove) {
  {
    // small type
    Any a;
    BOOST_CHECK(!a);

    int v = 5;
    a = Any{v};
    BOOST_CHECK(!!a);

    Any b = std::move(a);
    BOOST_CHECK(!!b);
    BOOST_CHECK_EQUAL(b.as<int>(), 5);

    Any c;
    c = std::move(b);
    BOOST_CHECK(!!c);
    BOOST_CHECK_EQUAL(c.as<int>(), 5);
  }

  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyCopy) {
  {
    // small type
    Any a;
    BOOST_CHECK(!a);

    int v = 5;
    a = Any{v};
    BOOST_CHECK(!!a);

    Any b = a;
    BOOST_CHECK(!!b);
    BOOST_CHECK_EQUAL(b.as<int>(), 5);

    Any c;
    c = a;
    BOOST_CHECK(!!c);
    BOOST_CHECK_EQUAL(c.as<int>(), 5);
  }
  CHECK_ANY_ALLOCATIONS();
}

struct D {
  bool* destroyed;
  explicit D(bool* d) : destroyed{d} {}
  ~D() { *destroyed = true; }
};

struct D2 {
  bool* destroyed{nullptr};
  std::array<char, 512> blob{};

  explicit D2(bool* d) : destroyed{d} {}

  ~D2() { *destroyed = true; }
};

BOOST_AUTO_TEST_CASE(AnyEmplace) {
  {
    Any a;
    auto& value = a.emplace<int>(42);
    BOOST_CHECK_EQUAL(value, 42);
    BOOST_CHECK_EQUAL(a.as<int>(), 42);
    value = 84;
    BOOST_CHECK_EQUAL(a.as<int>(), 84);
  }
  CHECK_ANY_ALLOCATIONS();

  {
    bool destroyed = false;
    Any a{std::in_place_type<D>, &destroyed};
    BOOST_CHECK(!destroyed);
    a.emplace<int>(7);
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(a.as<int>(), 7);
  }
  CHECK_ANY_ALLOCATIONS();

  {
    bool destroyed = false;
    Any a{std::in_place_type<D2>, &destroyed};
    BOOST_CHECK(!destroyed);
    bool destroyed2 = false;
    auto& ref = a.emplace<D2>(&destroyed2);
    BOOST_CHECK(destroyed);
    BOOST_CHECK(!destroyed2);
    BOOST_CHECK_EQUAL(ref.destroyed, &destroyed2);
    BOOST_CHECK_EQUAL(a.as<D2>().destroyed, &destroyed2);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyMoveTypeChange) {
  BOOST_TEST_CONTEXT("Small type") {
    bool destroyed = false;
    D d{&destroyed};
    Any a{std::move(d)};
    BOOST_CHECK(!destroyed);

    int value = 5;
    Any b{value};
    a = std::move(b);
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(a.as<int>(), value);
  }

  bool destroyed = false;
  BOOST_TEST_CONTEXT("Large type") {
    D2 d{&destroyed};
    Any a{std::move(d)};
    BOOST_CHECK(!destroyed);

    int value = 5;
    Any b{value};
    a = std::move(b);
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(a.as<int>(), value);
  }
}

BOOST_AUTO_TEST_CASE(AnyCopyTypeChange) {
  BOOST_TEST_CONTEXT("Small type") {
    bool destroyed = false;
    D d{&destroyed};
    Any a{std::move(d)};
    BOOST_CHECK(!destroyed);

    int value = 5;
    Any b{value};
    a = b;
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(a.as<int>(), value);
  }

  bool destroyed = false;
  BOOST_TEST_CONTEXT("Large type") {
    D2 d{&destroyed};
    Any a{std::move(d)};
    BOOST_CHECK(!destroyed);

    int value = 5;
    Any b{value};
    a = b;
    BOOST_CHECK(destroyed);
    BOOST_CHECK_EQUAL(a.as<int>(), value);
  }
}

BOOST_AUTO_TEST_CASE(AnyCopyTypeChangeToHeap) {
  // Copy-assigning a heap-allocated value over an Any that already holds a
  // value of a different type. Before the fix this read the stale contents of
  // the internal buffer as the copy destination (a freed pointer for the
  // heap->heap case, or the previous local value's bytes for local->heap),
  // resulting in a wild write / use-after-free.
  using Large = std::array<unsigned long, 5>;

  BOOST_TEST_CONTEXT("local -> heap") {
    Any a{42};
    Large v{1, 2, 3, 4, 5};
    Any b{v};
    a = b;
    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<Large>().begin(), a.as<Large>().end(),
                                  v.begin(), v.end());
    // source is left intact
    BOOST_CHECK_EQUAL_COLLECTIONS(b.as<Large>().begin(), b.as<Large>().end(),
                                  v.begin(), v.end());
  }
  CHECK_ANY_ALLOCATIONS();

  BOOST_TEST_CONTEXT("heap -> heap (different type)") {
    std::array<int, 4> v1{1, 2, 3, 4};
    Large v2{6, 7, 8, 9, 10};
    Any a{v1};
    Any b{v2};
    a = b;
    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<Large>().begin(), a.as<Large>().end(),
                                  v2.begin(), v2.end());
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyDestroy) {
  {  // small type
    bool destroyed = false;
    D d{&destroyed};
    BOOST_CHECK(!destroyed);
    {
      Any a{std::move(d)};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }
  CHECK_ANY_ALLOCATIONS();

  {  // large type
    bool destroyed = false;
    D2 d{&destroyed};
    BOOST_CHECK(!destroyed);
    {
      Any a{std::move(d)};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyDestroyCopy) {
  {  // small type
    bool destroyed = false;

    {
      Any b;
      {
        Any a{std::in_place_type<D>, &destroyed};
        BOOST_CHECK(!destroyed);
        b = a;
        BOOST_CHECK(!destroyed);
      }
      BOOST_CHECK(destroyed);  // a destroyed, should be true
      destroyed = false;
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);  // b destroyed, should be true again
  }
  CHECK_ANY_ALLOCATIONS();

  {  // large type
    bool destroyed = false;

    {
      Any b;
      {
        Any a{std::in_place_type<D2>, &destroyed};
        BOOST_CHECK(!destroyed);
        b = a;
        BOOST_CHECK(!destroyed);
      }
      BOOST_CHECK(destroyed);  // a destroyed, should be true
      destroyed = false;
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);  // b destroyed, should be true again
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyDestroyInPlace) {
  {  // small type
    bool destroyed = false;
    BOOST_CHECK(!destroyed);
    {
      Any a{std::in_place_type<D>, &destroyed};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }
  CHECK_ANY_ALLOCATIONS();

  {  // large type
    bool destroyed = false;
    BOOST_CHECK(!destroyed);
    {
      Any a{std::in_place_type<D2>, &destroyed};
      BOOST_CHECK(!destroyed);
    }
    BOOST_CHECK(destroyed);
  }
  CHECK_ANY_ALLOCATIONS();
}

struct D3 {
  std::size_t* destroyed{nullptr};
  std::array<char, 512> blob{};

  explicit D3(std::size_t* d) : destroyed{d} {}

  ~D3() { (*destroyed)++; }
};

BOOST_AUTO_TEST_CASE(LeakCheck) {
  std::size_t destroyed = 0;
  for (std::size_t i = 0; i < 10000; i++) {
    {
      BOOST_CHECK_EQUAL(destroyed, i);
      Any a;
      BOOST_CHECK_EQUAL(destroyed, i);
      a = Any{std::in_place_type<D3>, &destroyed};
      BOOST_CHECK_EQUAL(destroyed, i);
    }
    BOOST_CHECK_EQUAL(destroyed, i + 1);
  }
  CHECK_ANY_ALLOCATIONS();
}

struct LifecycleCounters {
  std::size_t nDestroy = 0;
  std::size_t nCopyConstruct = 0;
  std::size_t nCopy = 0;
  std::size_t nMoveConstruct = 0;
  std::size_t nMove = 0;
};

template <std::size_t PADDING>
struct Lifecycle;

template <>
struct Lifecycle<0> {
  LifecycleCounters* counters;

  explicit Lifecycle(LifecycleCounters* _counters) : counters{_counters} {}

  Lifecycle(Lifecycle&& o) {
    counters = o.counters;
    counters->nMoveConstruct++;
  }

  Lifecycle& operator=(Lifecycle&& o) {
    counters = o.counters;
    counters->nMove++;
    return *this;
  }

  Lifecycle(const Lifecycle& o) {
    counters = o.counters;
    counters->nCopyConstruct++;
  }

  Lifecycle& operator=(const Lifecycle& o) {
    counters = o.counters;
    counters->nCopy++;
    return *this;
  }

  ~Lifecycle() { counters->nDestroy++; }
};

template <std::size_t PADDING>
struct Lifecycle : public Lifecycle<0> {
  std::array<char, PADDING> m_padding{};

  explicit Lifecycle(LifecycleCounters* _counters) : Lifecycle<0>(_counters) {}
};

template <std::size_t PADDING>
struct LifecycleHandle {
  LifecycleCounters counters;
  Lifecycle<PADDING> inner;

  LifecycleHandle() : counters{}, inner{&counters} {}
};

#define checkCounters()                                                      \
  do {                                                                       \
    BOOST_REQUIRE_EQUAL(l.counters.nCopy, counters.nCopy);                   \
    BOOST_REQUIRE_EQUAL(l.counters.nCopyConstruct, counters.nCopyConstruct); \
    BOOST_REQUIRE_EQUAL(l.counters.nMove, counters.nMove);                   \
    BOOST_REQUIRE_EQUAL(l.counters.nMoveConstruct, counters.nMoveConstruct); \
    BOOST_REQUIRE_EQUAL(l.counters.nDestroy, counters.nDestroy);             \
  } while (0)

#define makeCounter(counter, n) \
  do {                          \
    counter += n;               \
    checkCounters();            \
  } while (0)

#define incCopyConstruct(n) makeCounter(counters.nCopyConstruct, n)
#define incCopy(n) makeCounter(counters.nCopy, n)
#define incMoveConstruct(n) makeCounter(counters.nMoveConstruct, n)
#define incMove(n) makeCounter(counters.nMove, n)
#define incDestroy(n) makeCounter(counters.nDestroy, n)

BOOST_AUTO_TEST_CASE(LifeCycleSmall) {
  LifecycleCounters counters;
  LifecycleHandle<0> l;

  checkCounters();

  {
    const auto& o = l.inner;  // force copy
    Any a{o};
    incCopyConstruct(1);

    const auto& _a = a;  // force copy later
    {
      Any b{_a};
      incCopyConstruct(1);
    }
    incDestroy(1);

    {
      Any b;
      b = _a;
      incCopyConstruct(1);
      b = _a;
      incCopy(1);
    }
    incDestroy(1);

    {
      Any b{a};
      incCopyConstruct(1);
      b = a;
      incCopy(1);
    }
    incDestroy(1);

    {
      auto _a2 = a;
      incCopyConstruct(1);
      Any b;
      b = std::move(_a2);
      incMoveConstruct(1);
      auto _a3 = a;
      incCopyConstruct(1);
      b = std::move(_a3);
      incMove(1);
    }
    incDestroy(3);
  }
  incDestroy(1);

  checkCounters();

  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(LifeCycleHeap) {
  LifecycleCounters counters;
  LifecycleHandle<512> l;

  checkCounters();

  {
    const auto& o = l.inner;  // force copy
    Any a{o};
    incCopyConstruct(1);

    const auto& _a = a;  // force copy later
    {
      Any b{_a};
      incCopyConstruct(1);
    }
    incDestroy(1);

    {
      Any b;
      b = _a;
      incCopyConstruct(1);
      b = _a;
      incCopy(1);
    }
    incDestroy(1);

    {
      Any b{a};
      incCopyConstruct(1);
      b = a;
      incCopy(1);
    }
    incDestroy(1);

    {
      Any _a2 = a;
      incCopyConstruct(1);
      Any b;
      b = std::move(_a2);
      // no actual move

      Any _a3 = a;
      incCopyConstruct(1);
      b = std::move(_a3);
      // no actual move
      incDestroy(1);
    }
    incDestroy(1);
  }
  incDestroy(1);

  checkCounters();

  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyMoveOnlyMoveOnlyTypes) {
  using MoveOnlyAny = Acts::AnyMoveOnly;

  using Ptr = std::unique_ptr<int>;

  // AnyMoveOnly can store move-only types
  {
    auto ptr = std::make_unique<int>(42);
    MoveOnlyAny a{std::move(ptr)};
    BOOST_CHECK(!!a);
    Ptr const* storedPtr = a.asPtr<Ptr>();
    BOOST_REQUIRE_NE(storedPtr, nullptr);
    BOOST_CHECK_NE(storedPtr->get(), nullptr);
    BOOST_CHECK_EQUAL(**storedPtr, 42);
  }

  // AnyMoveOnly is moveable
  {
    auto ptr = std::make_unique<int>(7);
    MoveOnlyAny a{std::move(ptr)};
    MoveOnlyAny b = std::move(a);
    // Note: moved-from Any may still report non-empty for local storage
    BOOST_CHECK(!!b);
    Ptr const* bPtr = b.asPtr<Ptr>();
    BOOST_REQUIRE_NE(bPtr, nullptr);
    int val = **bPtr;
    BOOST_CHECK_EQUAL(val, 7);
  }

  // AnyMoveOnly is not copyable
  static_assert(!std::is_copy_constructible_v<MoveOnlyAny>);
  static_assert(!std::is_copy_assignable_v<MoveOnlyAny>);
}

BOOST_AUTO_TEST_CASE(AnyTake) {
  // take() moves value out and leaves Any empty
  {
    Any a{42};
    BOOST_CHECK(!!a);
    int val = a.take<int>();
    BOOST_CHECK_EQUAL(val, 42);
    BOOST_CHECK(!a);
  }

  // take() with move-only type
  {
    auto ptr = std::make_unique<int>(99);
    Acts::AnyMoveOnly a{std::move(ptr)};
    BOOST_CHECK(!!a);
    auto taken = a.take<std::unique_ptr<int>>();
    BOOST_REQUIRE_NE(taken.get(), nullptr);
    BOOST_CHECK_EQUAL(*taken, 99);
    BOOST_CHECK(!a);
  }

  // take() throws on wrong type
  {
    Any a{42};
    BOOST_CHECK_THROW(a.take<float>(), std::bad_any_cast);
    BOOST_CHECK(!!a);
    BOOST_CHECK_EQUAL(a.as<int>(), 42);
  }

  // take() throws on empty
  {
    Any a;
    BOOST_CHECK_THROW(a.take<int>(), std::bad_any_cast);
  }
}

struct ThrowOnDemand {
  // larger than the small-buffer size, so this type is heap-allocated
  std::array<char, 64> blob{};
  explicit ThrowOnDemand(bool doThrow) {
    if (doThrow) {
      throw std::runtime_error{"boom"};
    }
  }
};

BOOST_AUTO_TEST_CASE(AnyEmplaceThrowingConstructor) {
  // If the constructor of the emplaced type throws, the previously held value
  // has already been destroyed, so the Any must be left empty rather than
  // claiming to hold a value whose storage was never initialised (which would
  // make the subsequent destruction operate on stale/freed memory).
  BOOST_TEST_CONTEXT("local") {
    struct SmallThrow {
      int x{0};
      explicit SmallThrow(bool doThrow) {
        if (doThrow) {
          throw std::runtime_error{"boom"};
        }
      }
    };
    Any a{std::in_place_type<SmallThrow>, false};
    BOOST_CHECK(!!a);
    BOOST_CHECK_THROW(a.emplace<SmallThrow>(true), std::runtime_error);
    BOOST_CHECK(!a);
  }
  CHECK_ANY_ALLOCATIONS();

  BOOST_TEST_CONTEXT("heap") {
    Any a{std::in_place_type<ThrowOnDemand>, false};
    BOOST_CHECK(!!a);
    BOOST_CHECK_THROW(a.emplace<ThrowOnDemand>(true), std::runtime_error);
    BOOST_CHECK(!a);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnySelfAssign) {
  using Large = std::array<unsigned long, 5>;

  BOOST_TEST_CONTEXT("self copy-assign local") {
    Any a{7};
    Any& r = a;
    a = r;
    BOOST_CHECK_EQUAL(a.as<int>(), 7);
  }
  BOOST_TEST_CONTEXT("self copy-assign heap") {
    Large v{1, 2, 3, 4, 5};
    Any a{v};
    Any& r = a;
    a = r;
    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<Large>().begin(), a.as<Large>().end(),
                                  v.begin(), v.end());
  }
  BOOST_TEST_CONTEXT("self move-assign local") {
    Any a{7};
    Any& r = a;
    a = std::move(r);
    BOOST_CHECK(!!a);
    BOOST_CHECK_EQUAL(a.as<int>(), 7);
  }
  BOOST_TEST_CONTEXT("self move-assign heap") {
    Large v{1, 2, 3, 4, 5};
    Any a{v};
    Any& r = a;
    a = std::move(r);
    BOOST_CHECK(!!a);
    BOOST_CHECK_EQUAL_COLLECTIONS(a.as<Large>().begin(), a.as<Large>().end(),
                                  v.begin(), v.end());
  }
  CHECK_ANY_ALLOCATIONS();
}

// Copying a heap-allocated value allocates and may run a throwing copy
// constructor, so the copy operations must not be marked noexcept (otherwise
// such a throw would call std::terminate). Moves only steal a pointer for
// heap values, so they are noexcept (which lets std::vector<Any> move rather
// than copy on reallocation) -- except in the allocation-tracking debug build
// (_ACTS_ANY_ENABLE_TRACK_ALLOCATIONS), where the tracking hooks may throw and
// kAnyNoexcept is false.
static_assert(!std::is_nothrow_copy_constructible_v<Any>);
static_assert(!std::is_nothrow_copy_assignable_v<Any>);
static_assert(std::is_nothrow_move_constructible_v<Any> ==
              detail::kAnyNoexcept);
static_assert(std::is_nothrow_move_assignable_v<Any> == detail::kAnyNoexcept);

struct HeapCopyThrows {
  // larger than the small-buffer size, so this type is heap-allocated
  std::array<int, 8> blob{};
  HeapCopyThrows() = default;
  HeapCopyThrows(const HeapCopyThrows& /*unused*/) {
    throw std::runtime_error{"copy boom"};
  }
  HeapCopyThrows& operator=(const HeapCopyThrows&) = default;
  HeapCopyThrows(HeapCopyThrows&&) = default;
  HeapCopyThrows& operator=(HeapCopyThrows&&) = default;
};

BOOST_AUTO_TEST_CASE(AnyCopyExceptionPropagates) {
  {
    Any a{std::in_place_type<HeapCopyThrows>};

    // copy construction propagates the stored type's throwing copy constructor
    BOOST_CHECK_THROW(Any b{a}, std::runtime_error);

    // copy assignment likewise propagates and leaves the target empty
    Any c;
    BOOST_CHECK_THROW(c = a, std::runtime_error);
    BOOST_CHECK(!c);
  }
  CHECK_ANY_ALLOCATIONS();
}

// ---------------------------------------------------------------------------
// AnyOf<Base>: typed type-erasure over a common polymorphic base
// ---------------------------------------------------------------------------

namespace {

struct AnyOfBase {
  virtual ~AnyOfBase() = default;
  virtual int value() const = 0;
};

// Small enough to live in the inline buffer (vtable ptr + int).
struct DerivedSmall : public AnyOfBase {
  int v;
  explicit DerivedSmall(int x) : v{x} {}
  int value() const override { return v; }
};

// Large enough to force heap storage for the default sb_size of AnyOf (48).
struct DerivedLarge : public AnyOfBase {
  std::array<std::int64_t, 16> pad{};
  int v;
  explicit DerivedLarge(int x) : v{x} { pad.fill(x); }
  int value() const override { return v; }
};

// AnyOfBase is the SECOND base, so its subobject sits at a non-zero offset
// within the complete object. Exercises the upcast offset adjustment.
struct OtherBase {
  virtual ~OtherBase() = default;
  std::int64_t tag = 0;
  virtual std::int64_t getTag() const { return tag; }
};

struct DerivedMI : public OtherBase, public AnyOfBase {
  int v;
  explicit DerivedMI(int x) : v{x} { tag = x; }
  int value() const override { return v; }
};

struct Unrelated {
  int v;
};

template <typename T>
concept HasAsBase = requires(T t) { t.asBase(); };

}  // namespace

// Plain Any (Base == void) exposes no base accessors; AnyOf<Base> does.
static_assert(!HasAsBase<Any>);
static_assert(HasAsBase<AnyOf<AnyOfBase>>);

// Storable iff convertible to AnyOfBase*; unrelated types are rejected.
static_assert(std::is_constructible_v<AnyOf<AnyOfBase>, DerivedSmall>);
static_assert(std::is_constructible_v<AnyOf<AnyOfBase>, DerivedLarge>);
static_assert(std::is_constructible_v<AnyOf<AnyOfBase>, DerivedMI>);
static_assert(!std::is_constructible_v<AnyOf<AnyOfBase>, Unrelated>);

BOOST_AUTO_TEST_CASE(AnyOfBasicAccess) {
  {
    AnyOf<AnyOfBase> a{DerivedSmall{42}};
    BOOST_CHECK(static_cast<bool>(a));

    // base access via asBase / operator* / operator->
    BOOST_REQUIRE(a.asBase() != nullptr);
    BOOST_CHECK_EQUAL(a.asBase()->value(), 42);
    BOOST_CHECK_EQUAL((*a).value(), 42);
    BOOST_CHECK_EQUAL(a->value(), 42);

    // concrete type is still recoverable
    BOOST_REQUIRE(a.asPtr<DerivedSmall>() != nullptr);
    BOOST_CHECK_EQUAL(a.asPtr<DerivedSmall>()->v, 42);
    BOOST_CHECK(a.asPtr<DerivedLarge>() == nullptr);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyOfConstAccess) {
  {
    const AnyOf<AnyOfBase> a{DerivedSmall{7}};
    const AnyOfBase* base = a.asBase();
    BOOST_REQUIRE(base != nullptr);
    BOOST_CHECK_EQUAL(base->value(), 7);
    BOOST_CHECK_EQUAL((*a).value(), 7);
    BOOST_CHECK_EQUAL(a->value(), 7);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyOfEmpty) {
  {
    AnyOf<AnyOfBase> a;
    BOOST_CHECK(!a);
    BOOST_CHECK(a.asBase() == nullptr);

    const AnyOf<AnyOfBase> b;
    BOOST_CHECK(b.asBase() == nullptr);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyOfCopyMoveLocal) {
  {
    AnyOf<AnyOfBase> a{DerivedSmall{1}};

    // copy construct
    AnyOf<AnyOfBase> b{a};
    BOOST_CHECK_EQUAL(a->value(), 1);
    BOOST_CHECK_EQUAL(b->value(), 1);
    // independent copies
    BOOST_CHECK(a.asBase() != b.asBase());

    // copy assign
    AnyOf<AnyOfBase> c;
    c = a;
    BOOST_CHECK_EQUAL(c->value(), 1);

    // move construct
    AnyOf<AnyOfBase> d{std::move(b)};
    BOOST_CHECK_EQUAL(d->value(), 1);

    // move assign
    AnyOf<AnyOfBase> e;
    e = std::move(c);
    BOOST_CHECK_EQUAL(e->value(), 1);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyOfCopyMoveHeap) {
  {
    AnyOf<AnyOfBase> a{DerivedLarge{5}};
    BOOST_CHECK_EQUAL(a->value(), 5);

    // copy construct allocates a fresh, independent heap object
    AnyOf<AnyOfBase> b{a};
    BOOST_CHECK_EQUAL(b->value(), 5);
    BOOST_CHECK(a.asBase() != b.asBase());

    // copy assign
    AnyOf<AnyOfBase> c;
    c = a;
    BOOST_CHECK_EQUAL(c->value(), 5);

    // move construct steals the pointer
    AnyOf<AnyOfBase> d{std::move(b)};
    BOOST_CHECK_EQUAL(d->value(), 5);
    BOOST_CHECK(!b);

    // move assign
    AnyOf<AnyOfBase> e;
    e = std::move(c);
    BOOST_CHECK_EQUAL(e->value(), 5);
    BOOST_CHECK(!c);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_CASE(AnyOfMultipleInheritanceOffset) {
  {
    AnyOf<AnyOfBase> a{DerivedMI{99}};

    DerivedMI* concrete = a.asPtr<DerivedMI>();
    BOOST_REQUIRE(concrete != nullptr);

    // Sanity: the AnyOfBase subobject is genuinely at a non-zero offset, so
    // this test actually exercises the offset adjustment.
    BOOST_REQUIRE(static_cast<void*>(static_cast<AnyOfBase*>(concrete)) !=
                  static_cast<void*>(concrete));

    // The function-pointer upcast must apply that same offset.
    BOOST_CHECK_EQUAL(a.asBase(), static_cast<AnyOfBase*>(concrete));
    BOOST_CHECK_EQUAL(a->value(), 99);
  }
  CHECK_ANY_ALLOCATIONS();
}

BOOST_AUTO_TEST_SUITE_END()

}  // namespace ActsTests
