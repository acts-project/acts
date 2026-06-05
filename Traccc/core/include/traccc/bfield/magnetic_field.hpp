/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Covfie include(s).
#include <covfie/core/concepts.hpp>
#include <covfie/core/field.hpp>

// System include(s).
#include <any>
#include <sstream>

namespace traccc {

/// Typeless, owning, host-only magnetic field object
class magnetic_field {

    public:
    /// Default constructor
    magnetic_field() = default;

    /// Constructor from a specific b-field object
    ///
    /// @tparam bfield_backend_t The backend type of the b-field object
    /// @param obj The b-field object to construct from
    ///
    template <covfie::concepts::field_backend bfield_backend_t>
    explicit magnetic_field(covfie::field<bfield_backend_t>&& obj);

    /// Set a specific b-field object
    ///
    /// @tparam bfield_backend_t The backend type of the b-field object
    /// @param obj The b-field object to set
    ///
    template <covfie::concepts::field_backend bfield_backend_t>
    void set(covfie::field<bfield_backend_t>&& obj);

    /// Check if the b-field is of a certain type
    ///
    /// @tparam bfield_backend_t The covfie backend type to check
    /// @return @c true if the b-field is of the specified type,
    ///         @c false otherwise
    ///
    template <covfie::concepts::field_backend bfield_backend_t>
    bool is() const;

    /// Get a b-field view object as a specific type
    ///
    /// @tparam bfield_backend_t The covfie backend type to use
    /// @return The b-field view object of the specified type
    ///
    template <covfie::concepts::field_backend bfield_backend_t>
    typename covfie::field<bfield_backend_t>::view_t as_view() const;

    /// Get the b-field object as a specific type
    ///
    /// @tparam bfield_backend_t The covfie backend type to use
    /// @return The b-field object cast to the specified type
    ///
    template <covfie::concepts::field_backend bfield_backend_t>
    const covfie::field<bfield_backend_t>& as_field() const;

    /// @brief Return type information about the contained magnetic field.
    const std::type_info& type() const { return m_field.type(); }

    private:
    /// The actualy covfie b-field object
    std::any m_field;

};  // class magnetic_field

/// @brief Helper function for `bfield_visitor`
template <typename callable_t, typename bfield_t, typename... bfield_ts>
auto magnetic_field_visitor_helper(const magnetic_field& bfield,
                                   callable_t&& callable,
                                   std::tuple<bfield_t, bfield_ts...>*)
    requires(covfie::concepts::field_backend<bfield_t> &&
             (covfie::concepts::field_backend<bfield_ts> && ...))
{
    if (bfield.is<bfield_t>()) {
        return callable(bfield.as_view<bfield_t>());
    } else {
        if constexpr (sizeof...(bfield_ts) > 0) {
            return magnetic_field_visitor_helper(
                bfield, std::forward<callable_t>(callable),
                static_cast<std::tuple<bfield_ts...>*>(nullptr));
        } else {
            std::stringstream exception_message;

            exception_message
                << "Invalid B-field type (" << bfield.type().name()
                << ") received, but this type is not supported" << std::endl;

            throw std::invalid_argument(exception_message.str());
        }
    }
}

/// @brief Visitor for polymorphic magnetic field types
///
/// This function takes a list of supported magnetic field types and checks
/// if the provided field is one of them. If it is, it will call the provided
/// callable on a view of it and otherwise it will throw an exception.
template <typename bfield_list_t, typename callable_t>
auto magnetic_field_visitor(const magnetic_field& bfield,
                            callable_t&& callable) {
    return magnetic_field_visitor_helper(bfield,
                                         std::forward<callable_t>(callable),
                                         static_cast<bfield_list_t*>(nullptr));
}

}  // namespace traccc

// Include the implementation.
#include "traccc/bfield/impl/magnetic_field.ipp"
