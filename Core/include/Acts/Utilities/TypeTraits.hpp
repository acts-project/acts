// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <functional>
#include <string>
#include <type_traits>

namespace Acts {

namespace detail {

/// This file contains an implementation of the detection idiom in C++.
/// It's not currently in the standard, but can be implemented using
/// standard features. This implementation is largely taken from the C++
/// [technical specifications, library fundamentals
/// v2](https://en.cppreference.com/w/cpp/experimental/is_detected)
///
/// The detector pattern works like this: there is a default type, that
/// accepts an "operation" that can be basically anything. It also accepts
/// variadic arguments for that operation. The default type
/// has a member type that is std::false_type to indicate success or
/// failure. It also has a member type "type" which captures a type result.
/// Then there is a specialization which attempts to instantiate the operation
/// with the given parameters, and tries to assign it into std::void_t. If the
/// operation fails to instantiate (say, the checked for type does not exist),
/// the specialization will not be instantiated, and the compiler falls back to
/// the default type which contains std::false_type. Since it happens inside
/// while the compiler tries to find a better matching template specialization
/// than the default one (so basically overload resolution), a compile error
/// inside the operation is handled as a substitution failure, and is not an
/// error. If the instantiation succeeds, the specialization contains a
/// std::true_type, and an alias to the result of the operation.
///
/// Essentially, it provides a convenient way to "lift" operations into this
/// overload resolution, allowing testing expressions and evaluating them into
/// compile time booleans (instead of compilation failures).

/// Helper struct which cannot be constructed (or destroyed) at all.
struct nonesuch {
  ~nonesuch() = delete;
  nonesuch(nonesuch const&) = delete;
  void operator=(nonesuch const&) = delete;
};

/// This is the default specialization.
/// It does not attempt to instantiate ``Op<Args...>`` at all.
/// @tparam Default The default type to set
/// @tparam AlwaysVoid Helper type that accepts the void instantiation
/// @tparam Op The operation to test
/// @tparam Args Arguments to the operation
template <class Default, class AlwaysVoid, template <class...> class Op,
          class... Args>
struct detector {
  using value_t = std::false_type;
  using type = Default;
};

/// This is the specialization which attempts to instantiate ``Op<Args...``.
/// @tparam Default Default type to set if substitution fails
/// @tparam Op The operation to test
/// @tparam Args Arguments to the operation
template <class Default, template <class...> class Op, class... Args>
struct detector<Default, std::void_t<Op<Args...>>, Op, Args...> {
  // Note that std::void_t is a C++17 feature
  using value_t = std::true_type;
  using type = Op<Args...>;
};

}  // namespace detail

namespace Concepts {

/// This type ties together the detection idiom. It instantiates the
/// ``detector`` template with the ``Op`` and ``Args`` and resolves to the exact
/// value type. In essence, if ``Op<Args...>`` succeeds, this will evaluate to
/// ``std::true_type``, and if not, it will evaluate to ``std::false_type``.
/// @tparam Op The operation to test
/// @tparam Args The arguments to the operation
template <template <class...> class Op, class... Args>
using is_detected =
    typename detail::detector<detail::nonesuch, void, Op, Args...>::value_t;

/// This type calls into the detector (same as ``is_detected``) but it extracts
/// the return type of ``Op<Args...>``.
/// @tparam Op The operation
/// @tparam Args The arguments to the operation
template <template <class...> class Op, class... Args>
using detected_t =
    typename detail::detector<detail::nonesuch, void, Op, Args...>::type;

/// This invokes ``detected_t``, and checks whether its result matches
/// ``Expected``.
/// @tparam Expected The expected result of the operation.
/// @tparam Op The operation
/// @tparam Args The arguments to the operation
template <class Expected, template <class...> class Op, class... Args>
using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;

/// This evaluates ``Op`` inside the detector, and checks whether the resolved
/// type
/// is convertible to ``To``.
/// @tparam To The type to check convertibility to.
/// @tparam Op The operation
/// @tparam Args The arguments to the operation
template <class To, template <class...> class Op, class... Args>
using is_detected_convertible =
    std::is_convertible<detected_t<Op, Args...>, To>;

/// Helper which invokes the detector with a default type, and resolves to the
/// type.
/// @tparam Default The type to resolve to if ``Op<Args...>`` does not resolve.
/// @tparam Op The operation
/// @tparam Args The argument to the operation
template <class Default, template <class...> class Op, class... Args>
using detected_or = detail::detector<Default, void, Op, Args...>;

/// Define some sort of "Domain Specific Language" to declare concepts a little
/// more naturally. These are taken from
/// https://izzys.casa/2016/09/implementing-concepts-in-cxx/

/// Helper which combines a set of predicates (constexpr bools) with a logical
/// AND. Converts to ``std::bool_constant``.
/// @tparam Bs The booleans to combine
template <bool... Bs>
constexpr bool require = std::conjunction<std::bool_constant<Bs>...>::value;

/// Helper which forms the logical OR of its arguments.
/// Converts to ``std::bool_constant``.
/// @tparam Bs The booleans to combine.
template <bool... Bs>
constexpr bool either = std::disjunction<std::bool_constant<Bs>...>::value;

/// Alias for the negation of a ``require``. This is essentially a NOT ANY test.
/// @tparam Bs The booleans.
template <bool... Bs>
constexpr bool disallow = !require<Bs...>;

/// Alias to ``is_detected`` which unpacks the constexpr boolean value.
/// @tparam Op The operation
/// @tparam Args The arguments to the operation.
template <template <class...> class Op, class... Args>
constexpr bool exists = is_detected<Op, Args...>::value;

/// Alias to conversion check, which also extracts the constexpr boolean value.
/// @tparam To The type to check convertibility to.
/// @tparam Op The operation
/// @tparam Args The arguments to the operation.
template <class To, template <class...> class Op, class... Args>
constexpr bool converts_to = is_detected_convertible<To, Op, Args...>::value;

/// Unpacks the constexpr boolean value from ``is_detected_exact``
/// @tparam Exact The type to check identity against
/// @tparam Op The operation
/// @tparam Args The arguments to the operation.
template <class Exact, template <class...> class Op, class... Args>
constexpr bool identical_to = is_detected_exact<Exact, Op, Args...>::value;

/// Helper which evaluates whether the type ``T`` has a method with a given
/// signature.
/// @tparam T The type to check on. This can contain a const qualifier if you
/// want to check on that.
/// @tparam R The return type
/// @tparam M The method trait, as generated by METHOD_TRAIT
/// @tparam Arguments The argument types that make up the signature.
template <typename T, typename R, template <class...> class M,
          typename... Arguments>
constexpr bool has_method = M<T, R, Arguments...>::template tv<T>::value;

/// Helper to assert if a member of a given type exists. Basically only calls
/// into ``identical_to`` but is nicer to read.
/// @tparam T The type to check existence of member on.
/// @tparam M The member type trait
/// @tparam V The type that the member is supposed to have.
template <typename T, template <class...> class M, typename V>
constexpr bool has_member = identical_to<V, M, T>;

/// Have a look at ``TypeTraitsTest.cpp`` to see most of this in action.
}  // namespace Concepts
}  // namespace Acts

/// These helpers allow writing checks. The missing piece is something that you
/// can put into these. This is the point where this turns into type checking.
/// ``Op`` from above can be anything. For instance, we can write an expression,
/// and have the compiler try to calculate that expression's resulting type,
/// like so:
///
/// .. code-block:: cpp
///
///    decltype(std::declval<T>().member_a)
///
/// ``std::declval<T>()`` constructs a **pseudo-value** of type ``T`` which
/// works even if ``T`` is not actually constructible like this. Then we access
/// a member called ``member_a`` inside and instruct the compiler to calculate
/// the resulting type using ``decltype``. This will only work if the expression
/// is valid. If not, this would normally cause a compilation error. If we wrap
/// it into the
/// **detection idiom**, however, we can turn that compilation error into a
/// simple constexpr false!
///
/// To do that, we need to put that expression above into a **metafunction**. If
/// we simply wrote
///
/// .. code-block:: cpp
///
///    is_detected<decltype(std::declval<T>().member_a)>
///
/// where ``decltype(std::declval<T>().member_a)`` is ``Op``, and ``Args`` is
/// empty, we still get a compilation error. This is because the compiler
/// evaluates the first argument before even passing it into ``is_detected``,
/// and that means
/// **outside** the overload resolution context. This is why we need to pass a
/// metafunction around the expression as ``Op``, and the concrete type ``T`` as
/// ``Args`` so that the actual evaluation of the expression happens **inside
/// the overload resolution context**. Luckily, writing a metafunction around
/// the expression is as easy as
///
/// .. code-block:: cpp
///
///    template <typename T>
///    using member_a_t = decltype(std::declval<T>().member_a>);
///
/// and we can then use it like ``is_detected<member_a_t, T>``.
///
/// Basically, what you have to do to write type assertions using this pattern
/// is metafunctions like the one described above. Those can be member checks,
/// like shown before, nested type checks like
///
/// .. code-block:: cpp
///
///    template <typename T>
///    using nested_a_t = typename T::NestedA;
///
/// and checks for contained templates like
///
/// .. code-block:: cpp
///
///    template <typename T>
///    using meta_t = typename T::template meta<void, void>;
///
/// but also arbitrary expressions, like operators, or even operators between
/// types. For instance, say we have two types ``U`` and ``V``. We can then
/// write a metafunction that attempts to instantiate a binary operation between
/// the two:
///
/// .. code-block:: cpp
///
///    template <typename U, typename V>
///    using binary_op_t = decltype(std::declval<U>() + std::declval<V>());
///
/// and simply takes both ``U`` and ``V`` as meta-arguments.
/// Invoked like ``is_detected<binary_op_t, TypeA, TypeB>`` on some types
/// ``TypeA`` and ``TypeB``, this will tell you whether you can call that binary
/// operation on values of those types.
///
/// Implementing method checks is a little more involved. A simple method check
/// can be written like
///
/// .. code-block:: cpp
///
///    template <typename T>
///    using foo_method_t = decltype(std::declval<T>().foo(double, int)>);
///
/// This only checks if the given expression is valid. That means this
/// will evaluate to true even if the actual arguments of that function are
/// references, as the compiler will figure that out behind the scenes and still
/// allow you to call it. That can be fine, if you're really only interested if
/// that specific call is allowed. Remember
/// that ``decltype`` calculates the type of the expression. In this context, it
/// will evaluate to **the return type** of the called function. Using
/// ``identical_to`` you can therefore assert the return type to be a given
/// value.
///
/// You might, however, want to constrain the *exact* types of the arguments,
/// and the method's const qualifier. This is less straightforward. To achieve
/// this, a macro is provided below, which will generate the required code for a
/// trait predicate. You *cannot* use the result of this macro directly with the
/// helpers defined above. What the macro provides is a metafunction, which
/// still takes the type, the return type and the arguments as metaarguments. So
/// the result of the macro only encodes the name of the method, not the rest of
/// its signature. That means you can use the same method trait to check for any
/// signature involving the same method name. Say you generated
///
/// .. code-block:: cpp
///
///    // foo_method_t is the name of the struct generated by this macro,
///    // foo is the name of the method you want to check for.
///    METHOD_TRAIT(foo_method_t, foo);
///
/// An additional helper ``has_method`` is provided,
/// which wraps a little boilerplate to check a specific signature.
/// You can then write
///
/// .. code-block:: cpp
///
///    has_method<T, R, foo_method_t, bool, const int&>
///
/// to check for a signature of the form
///
/// .. code-block:: cpp
///
///    R T::foo(bool, const int&)
///
/// If you want to check for a const method you can modify this to
///
/// .. code-block:: cpp
///
///    has_method<const T, R, foo_method_t, bool, const int&>
///
/// Note that both will only evaluate to true if the const qualifier matches
/// what you gave exactly. If you want to check for a method of a given
/// specifier and you don't care if the method is const or not, you have to
/// write out both variants explicitly, and combine them with ``either``.

/// @cond

/// This macro generates some boilerplate code that is necessary to correctly
/// implement a method type trait that evaluates to constexpr bools correctly
/// @param trait_name The resulting name of the trait.
/// @param method_name The name of the method the trait shall check.
#define METHOD_TRAIT(trait_name, method_name)                                  \
  template <class T, typename R, typename... Arguments>                        \
  struct trait_name {                                                          \
    /* Meta function to check if a type has a const qualifier*/                \
    /* (by stripping it and seeing if something changed */                     \
    template <typename T_>                                                     \
    static constexpr bool is_const =                                           \
        !std::is_same_v<std::remove_const_t<T_>, T_>;                          \
                                                                               \
    /*These following meta-functions basically to this: they check whether or  \
     * not the actual function pointer extracted through ``&T::method_name``   \
     * can                                                                     \
     * be assigned to a prepared function pointer type with the given          \
     * signature. This checks the exact signature, and not just callability    \
     * and validity of the expression. */                                      \
                                                                               \
    /* Meta function which constructs the right type to check a function       \
     * pointer, non-const version*/                                            \
    template <typename T_, typename = int>                                     \
    struct fptr_meta {                                                         \
      template <typename... Arguments_>                                        \
      using type = typename std::integral_constant<                            \
          decltype(std::declval<T_>().method_name(                             \
              std::declval<Arguments_>()...)) (T_::*)(Arguments_...),          \
          &T_::method_name>::value_type;                                       \
    };                                                                         \
                                                                               \
    /* Meta function which constructs the right type to check a function       \
     * pointer, const version*/                                                \
    /* The ``const`` needs to be put in there in one specific spot, that's why \
     * the metafunction is needed*/                                            \
    template <typename T_>                                                     \
    struct fptr_meta<T_, std::enable_if_t<is_const<T_>, int>> {                \
      template <typename... Arguments_>                                        \
      using type = typename std::integral_constant<                            \
          decltype(std::declval<T_>().method_name(                             \
              std::declval<Arguments_>()...)) (T_::*)(Arguments_...) const,    \
          &T_::method_name>::value_type;                                       \
    };                                                                         \
                                                                               \
    /* Helper on top of the function pointer metafunction */                   \
    template <typename T_, typename... Arguments_>                             \
    using fptr_meta_t = typename fptr_meta<T_>::template type<Arguments_...>;  \
                                                                               \
    /* Trait check for the qualifier and the return type of the function */    \
    /* This does not check the const qualifier at all */                       \
    template <typename T_, typename... Arguments_>                             \
    using qual_ret = decltype(std::declval<T_>().method_name(                  \
        std::declval<Arguments_>()...));                                       \
                                                                               \
    /* The problem is this: while the above is fine with and without const,    \
     * and with and without exact argument type match, the assignment to the   \
     * function pointer fails hard if there is no method at all with the given \
     * name. That is undesirable. The following first uses the expression      \
     * validity check to assert that there is in fact a method of the given    \
     * name, and only if that is the case, try to compile the function pointer \
     * based signature check. That way, there is no hard failures, only        \
     * substitution failures and we're happy. */                               \
    template <typename T_, typename = int>                                     \
    struct tv {                                                                \
      static constexpr bool value = false;                                     \
    };                                                                         \
    template <typename T_>                                                     \
    struct tv<T_, std::enable_if_t<Acts::Concepts::is_detected_exact<          \
                                       R, qual_ret, T_, Arguments...>::value,  \
                                   int>> {                                     \
      /* This is only ever evaluate if the method exists!*/                    \
      static constexpr bool value =                                            \
          Acts::Concepts::is_detected<fptr_meta_t, T, Arguments...>::value;    \
    };                                                                         \
  }

/// @endcond
