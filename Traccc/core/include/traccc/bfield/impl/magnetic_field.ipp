/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

namespace traccc {

template <covfie::concepts::field_backend bfield_backend_t>
magnetic_field::magnetic_field(covfie::field<bfield_backend_t>&& obj)
    : m_field(std::move(obj)) {}

template <covfie::concepts::field_backend bfield_backend_t>
void magnetic_field::set(covfie::field<bfield_backend_t>&& obj) {
    m_field = std::move(obj);
}

template <covfie::concepts::field_backend bfield_backend_t>
bool magnetic_field::is() const {
    return (m_field.type() == typeid(covfie::field<bfield_backend_t>));
}

template <covfie::concepts::field_backend bfield_backend_t>
typename covfie::field<bfield_backend_t>::view_t magnetic_field::as_view()
    const {
    return typename covfie::field<bfield_backend_t>::view_t{
        std::any_cast<const covfie::field<bfield_backend_t>&>(m_field)};
}

template <covfie::concepts::field_backend bfield_backend_t>
const covfie::field<bfield_backend_t>& magnetic_field::as_field() const {
    return std::any_cast<const covfie::field<bfield_backend_t>&>(m_field);
}

}  // namespace traccc
