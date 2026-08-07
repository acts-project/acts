// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/concepts.hpp"
#include "detray/algebra/common/ksm/detail/math_helpers.hpp"
#include "detray/algebra/common/ksm/substructure.hpp"
#include "detray/algebra/common/ksm/value.hpp"
#include "detray/algebra/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/assert.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

// System include(s)
#include <array>
#include <concepts>
#include <cstddef>
#include <type_traits>
#include <utility>

// This file defines the KSM matrix, which materializes a substructure into a
// type that can actually be used.
namespace detray::ksm {

/// @brief The primary KSM matrix type.
///
/// This type is defined as a substructure and a scalar type, and manifests as
/// a real, usable C++ type that can be used to perform matrix arithmetic.
template <typename substructure_t, typename scalar_t>
struct matrix {
  /// @brief The substructure types, both as-is and the canonical variant
  ///
  /// @{
  using substructure_type = substructure_t;
  using canonical_substructure_type = typename substructure_t::canonical_type;
  /// @}

  // Protect against absolutely disastrous failures in the canonicalization
  // algorithm
  static_assert(substructure_type::columns ==
                canonical_substructure_type::columns);
  static_assert(substructure_type::rows == canonical_substructure_type::rows);
  static_assert(std::same_as<canonical_substructure_type,
                             typename detail::canonicalize_substructure<
                                 canonical_substructure_type>::type>);
  static_assert(concepts::is_canonical<canonical_substructure_type>);

  /// @brief The scalar type stored inside this matrix.
  using scalar_type = scalar_t;

  /// @brief The shape of the matrix
  ///
  /// @{
  static constexpr std::size_t rows = substructure_t::rows;
  static constexpr std::size_t columns = substructure_t::columns;
  /// @}

  /// @brief Constant accessor for a value at a given index.
  ///
  /// Because substructures are compile-time, we must accept the indices I and
  /// J as template arguments (unless we want to emit an enormous jumble of
  /// if-else clauses).
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_t at() const {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    static_assert(concepts::is_canonical_value<value_t>);

    // Decide whether the value is structural or not; if it is, just return
    // the value as a scalar. If it is a runtime value, load it from memory
    // using the element_index accessor.
    if constexpr (concepts::is_index_variable<value_t>) {
      return m_elems[canonical_substructure_type::template element_index<I, J>];
    } else {
      return static_cast<scalar_t>(
          canonical_substructure_type::template value_at<I, J>::value);
    }
  }

  /// @brief Mutable accessor for a value at a given index.
  ///
  /// Because a mutable accessor to a structural value makes no sense, this
  /// function carries a requires clause which makes it a compile-time error
  /// to try to assign a value to a structural cell. This may make the API of
  /// this type less ergonomic at first, but it allows us to enforce structure
  /// at compile time which is an enormously powerful tool.
  template <std::size_t I, std::size_t J>
    requires(concepts::is_index_variable<
             typename canonical_substructure_type::template value_at<I, J>>)
  DETRAY_HOST_DEVICE scalar_t &at() {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    static_assert(concepts::is_canonical_value<value_t>);
    return m_elems[canonical_substructure_type::template element_index<I, J>];
  }

  /// @brief Constant accessor for a value at a given index, with detray
  /// terminology.
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_t element() const {
    return this->template at<I, J>();
  }

  /// @brief Mutable accessor for a value at a given index, with detray
  /// terminology.
  template <std::size_t I, std::size_t J>
    requires(concepts::is_index_variable<
             typename canonical_substructure_type::template value_at<I, J>>)
  DETRAY_HOST_DEVICE scalar_t &element() {
    return this->template at<I, J>();
  }

  /// @returns The identity matrix in this substructure.
  ///
  /// For this function to be usable, the substructure must have only
  /// variables or static ones on the diagonal, and static zeroes or variables
  /// on the other cells.
  DETRAY_HOST_DEVICE static matrix identity()
    requires(concepts::can_represent_identity<canonical_substructure_type>)
  {
    matrix rv{};
    // fill_with already skips a cell without storage, and writes a shared one
    // only through the cell that owns it. The requires clause above is what
    // makes that the right thing to do: it rules out a constant cell that
    // disagrees with the identity, and cells sharing a variable that want
    // different values.
    fill_with(rv, []<std::size_t I, std::size_t J>() {
      return (I == J) ? scalar_t{1} : scalar_t{0};
    });
    return rv;
  }

  /// @brief Turn this matrix into a dense detray matrix with the same values.
  template <::detray::concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE dmatrix<algebra_t, rows, columns> to_dense() const {
    dmatrix<algebra_t, rows, columns> rv{};
    to_dense_cells<algebra_t>(rv, std::make_index_sequence<rows * columns>{});
    return rv;
  }

  /// @brief Create a structured matrix from a dense detray matrix.
  template <::detray::concepts::algebra algebra_t>
  DETRAY_HOST_DEVICE static matrix from_dense(
      const dmatrix<algebra_t, rows, columns> &m) {
    matrix rv{};
    // Load first, then check. fill_with writes each stored value once, through
    // the cell that owns it, and skips the cells the substructure fixes.
    fill_with(rv, [&m]<std::size_t I, std::size_t J>() {
      return static_cast<scalar_t>(getter::element<I, J>(m));
    });
    rv.check_dense_cells(m, std::make_index_sequence<rows * columns>{});
    return rv;
  }

  /// @brief Copy every element of another matrix into this matrix, at a given
  /// offset.
  template <std::size_t RowOffset, std::size_t ColOffset, typename src_t>
  DETRAY_HOST_DEVICE void copy_from(const src_t &src) {
    constexpr auto src_rows =
        static_cast<std::size_t>(::detray::traits::rows<src_t>);
    constexpr auto src_columns =
        static_cast<std::size_t>(::detray::traits::columns<src_t>);

    static_assert(RowOffset + src_rows <= rows,
                  "the copied elements do not fit the destination");
    static_assert(ColOffset + src_columns <= columns,
                  "the copied elements do not fit the destination");

    const auto element_of =
        [&src]<std::size_t I, std::size_t J>() -> decltype(auto) {
      // Read cell (I, J) of src, which may be a dense matrix or a vector. A
      // vector stands for a single column, so only its row index is meaningful.
      if constexpr (::detray::concepts::matrix<src_t>) {
        return ::detray::getter::element<I, J>(src);
      } else {
        static_assert(J == 0u, "a vector stands for a single column");
        return ::detray::getter::element<I>(src);
      }
    };

    [this, &element_of]<std::size_t... Ks>(std::index_sequence<Ks...>) {
      ((this->template at<RowOffset + Ks / src_columns,
                          ColOffset + Ks % src_columns>() =
            static_cast<scalar_t>(
                element_of.template
                operator()<Ks / src_columns, Ks % src_columns>())),
       ...);
    }(std::make_index_sequence<src_rows * src_columns>{});
  }

  /// @brief Compute the product of this structured matrix with another.
  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto operator*(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type = matrix<
        typename canonical_substructure_type::template multiplication_type<
            typename other_substructure_t::canonical_type>,
        scalar_t>;

    result_type rv{};
    // For multiplication we fully rely on the `product_dot` helper function
    // to compute the value at each cell.
    fill_with(rv, [this, &o]<std::size_t I, std::size_t J>() {
      return this->template product_dot<I, J>(
          o, std::make_index_sequence<canonical_substructure_type::columns>{});
    });
    return rv;
  }

  /// @brief Compute the addition of this structured matrix with another.
  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto operator+(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type =
        matrix<typename canonical_substructure_type::template addition_type<
                   typename other_substructure_t::canonical_type>,
               scalar_t>;

    result_type rv{};
    fill_with(rv, [this, &o]<std::size_t I, std::size_t J>() {
      // Here we essentially re-encode the algebra system defined in the
      // substructure: A + 0 = A, 0 + A = A, and anything else we just add up
      // as usual. This is important so we do not emit useless additions with
      // zero.
      using A_val =
          typename canonical_substructure_type::template value_at<I, J>;
      using B_val =
          typename other_substructure_t::canonical_type::template value_at<I,
                                                                           J>;
      if constexpr (std::same_as<A_val, zero>) {
        return o.template at<I, J>();
      } else if constexpr (std::same_as<B_val, zero>) {
        return this->template at<I, J>();
      } else {
        return this->template at<I, J>() + o.template at<I, J>();
      }
    });
    return rv;
  }

  /// @brief Compute the subtraction of this structured matrix with another.
  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto operator-(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type =
        matrix<typename canonical_substructure_type::template subtraction_type<
                   typename other_substructure_t::canonical_type>,
               scalar_t>;

    result_type rv{};
    fill_with(rv, [this, &o]<std::size_t I, std::size_t J>() {
      // Same logic as with addition. Note that we don't need to account for
      // the case where we subtract something from itself, as the substructure
      // will already have reduced that to a structural zero!
      using A_val =
          typename canonical_substructure_type::template value_at<I, J>;
      using B_val =
          typename other_substructure_t::canonical_type::template value_at<I,
                                                                           J>;
      if constexpr (std::same_as<A_val, zero>) {
        return -o.template at<I, J>();
      } else if constexpr (std::same_as<B_val, zero>) {
        return this->template at<I, J>();
      } else {
        return this->template at<I, J>() - o.template at<I, J>();
      }
    });
    return rv;
  }

  /// @brief Compute the congruence of this structured matrix with another,
  /// i.e. compute this * other * this^T.
  template <typename other_substructure_t>
  DETRAY_HOST_DEVICE auto congruence(
      const matrix<other_substructure_t, scalar_t> &o) const {
    using result_type =
        matrix<typename canonical_substructure_type::template congruence_type<
                   typename other_substructure_t::canonical_type>,
               scalar_t>;

    // First we compute the product this * other, and then this^T separately.
    const auto ab = *this * o;
    const auto at = this->transpose();

    result_type rv{};
    // Note that, in this case, fill_with is smart enough to recognise that
    // the result type (as decided by the substructure) is symmetric if the
    // other matrix is symmetric, and it will compute and write the shared
    // cells only once.
    fill_with(rv, [this, &o, &ab, &at]<std::size_t I, std::size_t J>() {
      return ab.template product_dot<I, J>(
          at, std::make_index_sequence<canonical_substructure_type::columns>{});
    });
    return rv;
  }

  /// @brief Compute the inverse of this structured matrix.
  ///
  /// @note Currently only works for dimensions 1 and 2.
  DETRAY_HOST_DEVICE auto inverse() const {
    using result_type = matrix<typename detail::inverted_substructure<
                                   canonical_substructure_type>::type,
                               scalar_t>;

    result_type rv{};

    // Check some constraints that are relevant for inversion.
    static_assert(rows == 1 || rows == 2);
    static_assert(concepts::is_invertible<canonical_substructure_type>);

    using determinant_type =
        detail::get_determinant_helper<canonical_substructure_type>::type;

    static_assert(!std::same_as<determinant_type, zero>);

    // HACK: This code is a big mess; hopefully we can refactor it in the
    // future to be more readable.

    scalar_t det;

    // Compute the determinant first.
    if constexpr (std::same_as<determinant_type, one>) {
      // This should never be used, as we will use the structural knowledge
      // to avoid division entirely.
      det = scalar_t{0};
    } else {
      // If the determinant is not structurally one, compute it for the 1D and
      // 2D cases.
      if constexpr (rows == 1) {
        det = this->template at<0, 0>();
      } else {
        det = this->template at<0, 0>() * this->template at<1, 1>() -
              this->template at<1, 0>() * this->template at<0, 1>();
      }
    }

    // Then compute the inverse.
    if constexpr (rows == 1) {
      // The 1D case is easy, we just divide by the determinant. Special case
      // if the value of the sole element in the matrix is one.
      fill_with(rv, [this, det]<std::size_t I, std::size_t J>() {
        if constexpr (std::same_as<determinant_type, one>) {
          return this->template at<I, J>();
        } else {
          return scalar_t{1} / det;
        }
      });
    } else {
      // My apologies for anyone reading this. This is the logic for turning
      // the matrix [[a, b], [c, d]] into the matrix [[d, -b], [-c, a]].
      fill_with(rv, [this, det]<std::size_t I, std::size_t J>() {
        if constexpr (std::same_as<determinant_type, one>) {
          if constexpr (I == 0 && J == 0) {
            return this->template at<1, 1>();
          } else if constexpr (I == 1 && J == 1) {
            return this->template at<0, 0>();
          } else {
            return -this->template at<I, J>();
          }
        } else {
          if constexpr (I == 0 && J == 0) {
            return this->template at<1, 1>() / det;
          } else if constexpr (I == 1 && J == 1) {
            return this->template at<0, 0>() / det;
          } else {
            return -this->template at<I, J>() / det;
          }
        }
      });
    }

    return rv;
  }

  /// @brief Compute the transpose of this structured matrix.
  DETRAY_HOST_DEVICE auto transpose() const {
    using result_type = matrix<typename detail::transposed_substructure<
                                   canonical_substructure_type>::type,
                               scalar_t>;

    result_type rv{};
    // This one is delightfully simple, we just get the position at the
    // swapped indices!
    fill_with(rv, [this]<std::size_t I, std::size_t J>() {
      return this->template at<J, I>();
    });
    return rv;
  }

  /// @brief Compute the negative of this structured matrix.
  DETRAY_HOST_DEVICE auto operator-() const {
    using result_type = matrix<typename detail::negated_substructure<
                                   canonical_substructure_type>::type,
                               scalar_t>;

    result_type rv{};
    // This one is also easy.
    fill_with(rv, [this]<std::size_t I, std::size_t J>() {
      return -this->template at<I, J>();
    });
    return rv;
  }

 private:
  /// @brief Write every cell, structural ones included, into a dense matrix.
  template <typename algebra_t, typename dense_t, std::size_t... Ks>
  DETRAY_HOST_DEVICE void to_dense_cells(dense_t &rv,
                                         std::index_sequence<Ks...>) const {
    ((getter::element<Ks / columns, Ks % columns>(rv) =
          static_cast<dscalar<algebra_t>>(
              this->template at<Ks / columns, Ks % columns>())),
     ...);
  }

  /// @brief Check that cell (I, J) of the dense matrix holds what this
  /// substructure needs.
  ///
  /// This functions as a debug assertion for `from_dense`.
  ///
  /// @{
  template <std::size_t I, std::size_t J, typename dense_t>
  DETRAY_HOST_DEVICE void check_dense_cell(
      [[maybe_unused]] const dense_t &m) const {
    using value_t =
        typename canonical_substructure_type::template value_at<I, J>;

    if constexpr (concepts::is_index_variable<value_t>) {
      assert((static_cast<scalar_t>(getter::element<I, J>(m)) ==
              this->template at<I, J>()));
    } else {
      assert((static_cast<scalar_t>(getter::element<I, J>(m)) ==
              static_cast<scalar_t>(value_t::value)));
    }
  }

  template <typename dense_t, std::size_t... Ks>
  DETRAY_HOST_DEVICE void check_dense_cells(const dense_t &m,
                                            std::index_sequence<Ks...>) const {
    (check_dense_cell<Ks / columns, Ks % columns>(m), ...);
  }
  /// @}

  /// @brief Helper function for the dot product, where both terms are
  /// constants.
  ///
  /// Computes the product between A_IK and B_KJ, for use in the product of
  /// row I of A and column J of B.
  ///
  /// @warning This returns 0, the additive identity, when either value is not
  /// an integral value.
  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t K>
  DETRAY_HOST_DEVICE static constexpr long dot_term_constant() {
    using A_val = typename canonical_substructure_type::template value_at<I, K>;
    using B_val = typename OtherSub::canonical_type::template value_at<K, J>;
    if constexpr (concepts::is_integral_value<A_val> &&
                  concepts::is_integral_value<B_val>) {
      return static_cast<long>(A_val::value) * static_cast<long>(B_val::value);
    } else {
      return 0;
    }
  }

  /// @brief Compute the sum of constants in the dot product of two rows.
  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t... Ks>
  DETRAY_HOST_DEVICE static constexpr long dot_constant_offset(
      std::index_sequence<Ks...>) {
    // Unlike everything we do at runtime on floating point numbers, we can
    // add a zero at the front of this fold expression without any cost, as
    // the compiler can elide that without risk.
    //
    // Recall that dot_term_constant computes the product of row I of A and
    // column J of B, and Ks are the indices along those rows.
    return (0 + ... + dot_term_constant<I, J, OtherSub, Ks>());
  }

  /// @brief Compute the value of the dot product for variables.
  ///
  /// This function uses the index K to determine the position along row I
  /// of matrix A and column J of matrix B, and recurses with the tail Ks.
  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t K,
            std::size_t... Ks>
  DETRAY_HOST_DEVICE auto product_dot_vars(
      const matrix<OtherSub, scalar_t> &o,
      std::index_sequence<K, Ks...>) const {
    using A_val = typename canonical_substructure_type::template value_at<I, K>;
    using B_val = typename OtherSub::canonical_type::template value_at<K, J>;

    // First, we compute if this multiplication will result in a zero value
    // because either operand is zero.
    //
    // NOTE: The substructure has already done this computation, but we do it
    // again for the runtime logic here.
    constexpr bool is_zero_term =
        std::same_as<A_val, zero> || std::same_as<B_val, zero>;

    // Second, check whether the result of this computation will be constant.
    constexpr bool both_const = concepts::is_integral_value<A_val> &&
                                concepts::is_integral_value<B_val>;

    // Finally, check whether this will actually produce a variable.
    constexpr bool is_var_term = !is_zero_term && !both_const;

    // Check whether we have any additional tail to process. This logic is so
    // complex simply to minimize the number of addition instructions we emit.
    // At compile time we can emit as much nonsense as we want and it will not
    // affect runtimes, but here we don't have that luxury so we have to
    // carefully break down all cases.
    if constexpr (sizeof...(Ks) > 0) {
      // First, the case where we have a tail. First, recurse to find the
      // value of the remainder of this dot product. Note that this might be
      // a scalar or a compile-time constant!
      const auto rem = product_dot_vars<I, J>(o, std::index_sequence<Ks...>{});
      using rem_type = std::decay_t<decltype(rem)>;

      if constexpr (!is_var_term) {
        // If this case returns a non-variable, we simply return the
        // remainder of the dot product.
        //
        // NOTE: The attentive reader will ask: "what about the constant value
        // product of the values at K?!". The answer is that those structural
        // constants will be pooled together by dot_term_constant in order to,
        // again, minimize runtime flops.
        return rem;
      } else if constexpr (std::same_as<A_val, one>) {
        // If we do return a variable but A is one, the value of this
        // multiplication pair is determined entirely by B.
        if constexpr (std::same_as<rem_type, zero>) {
          // If the remainder is structurally zero, don't emit B + 0 but just
          // B.
          return o.template at<K, J>();
        } else {
          // If the remainder is non-zero, emit an addition flop.
          return o.template at<K, J>() + rem;
        }
      } else if constexpr (std::same_as<B_val, one>) {
        // Similar case to above, but this time B is one, so the product is
        // determined by A.
        if constexpr (std::same_as<rem_type, zero>) {
          // If the remainder is zero, don't emit A + 0, but just A.
          return this->template at<I, K>();
        } else {
          // Otherwise, emit the full addition.
          return this->template at<I, K>() + rem;
        }
      } else {
        // In all other cases, we have to compute the product of A and B,
        // because we know them only at runtime.
        if constexpr (std::same_as<rem_type, zero>) {
          // Same logic as before. If the remainder is zero, emit just (A * B)
          return this->template at<I, K>() * o.template at<K, J>();
        } else {
          // Otherwise, emit (A * B) + R.
          return this->template at<I, K>() * o.template at<K, J>() + rem;
        }
      }
    } else {
      // If there is no tail to the dot product, the logic becomes simpler.
      // This logic is mostly a mirror of the logic above without the branch
      // on the value of the remainder.
      if constexpr (!is_var_term) {
        // If we produce a constant value, return a _structural_ zero. As
        // before, the dot_term_constant function will collect the value of
        // the constant products.
        return zero{};
      } else if constexpr (std::same_as<A_val, one>) {
        // If A is one, just emit B.
        return o.template at<K, J>();
      } else if constexpr (std::same_as<B_val, one>) {
        // If B is one, just emit A.
        return this->template at<I, K>();
      } else {
        // And if neither is one, just emit the product (A * B).
        return this->template at<I, K>() * o.template at<K, J>();
      }
    }
  }

  /// @brief Compute the value of the dot product between row I of this and
  /// column J of the other matrix.
  template <std::size_t I, std::size_t J, typename OtherSub, std::size_t... Ks>
  DETRAY_HOST_DEVICE auto product_dot(const matrix<OtherSub, scalar_t> &o,
                                      std::index_sequence<Ks...> seq) const {
    // First, compute the sums of the constants as well as the sum of the
    // variables.
    constexpr long offset = dot_constant_offset<I, J, OtherSub>(seq);
    const auto var_sum = product_dot_vars<I, J>(o, seq);
    using vs_type = std::decay_t<decltype(var_sum)>;

    // Check what the type of the variable sum is.
    if constexpr (std::same_as<vs_type, zero>) {
      // If the variable sum is somehow zero, we simply convert the offset
      // to a scalar. If the offset is zero, we still need to do this because
      // this function _must_ return a runtime scalar.
      return scalar_t(offset);
    } else if constexpr (offset == 0) {
      // If the offset is zero but the variable sum is not, just emit the
      // variable sum.
      return var_sum;
    } else {
      // If neither is zero, emit an addition flop.
      return var_sum + scalar_t(offset);
    }
  }

  /// @brief Helper logic for all cases in which we have to fill cells in KSM
  /// matrices.
  ///
  /// @{
  ///
  /// @brief Fill cell (I, J) of rv with the value returned by f.
  template <std::size_t I, std::size_t J, typename Result, typename F>
  DETRAY_HOST_DEVICE static void fill_cell(Result &rv, F &f) {
    // Remember that we cannot write to structural cells, and anonymous
    // variables don't exist in canonical substructures, so we only write to
    // non-anonymous cells.
    if constexpr (concepts::is_index_variable<
                      typename Result::canonical_substructure_type::
                          template value_at<I, J>>) {
      // We also don't want to do unnecessary work writing to non-anonymous
      // shared variables twice, so we check that this is the first time the
      // value appears.
      if constexpr (detail::first_equal_flat_index_v<
                        typename Result::canonical_substructure_type,
                        typename Result::canonical_substructure_type::
                            template value_at<I, J>> ==
                    I * Result::columns + J) {
        // If this is a non-anonymous cell, and it is the first time it
        // appears, call into f and set the value of the cell.
        rv.template at<I, J>() = f.template operator()<I, J>();
      }
    }
  }

  /// @brief Helper function that takes a flat pack of indices Ks and turns
  /// them into (I, J) pairs.
  template <typename Result, typename F, std::size_t... Ks>
  DETRAY_HOST_DEVICE static void fill_all(Result &rv, F f,
                                          std::index_sequence<Ks...>) {
    constexpr std::size_t C = Result::canonical_substructure_type::columns;
    (fill_cell<Ks / C, Ks % C>(rv, f), ...);
  }

  /// @brief Fill matrix cells using a callable object f.
  template <typename Result, typename F>
  DETRAY_HOST_DEVICE static void fill_with(Result &rv, F f) {
    using result_sub = typename Result::canonical_substructure_type;
    fill_all(
        rv, f,
        std::make_index_sequence<result_sub::rows * result_sub::columns>{});
  }
  /// @}

  /// @brief The storage of the matrix elements.
  std::array<scalar_t, canonical_substructure_type::num_variables> m_elems{};

  // All KSM matrices are friends. <3
  template <class, class>
  friend class matrix;
};
}  // namespace detray::ksm
