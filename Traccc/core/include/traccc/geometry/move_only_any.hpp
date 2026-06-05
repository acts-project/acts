/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <any>

#include "traccc/definitions/primitives.hpp"

namespace traccc {

class move_only_any {
    public:
    move_only_any() = default;
    move_only_any(const move_only_any &) = delete;
    move_only_any &operator=(const move_only_any &) = delete;

    template <typename obj_t>
    explicit move_only_any(obj_t &&obj)
        requires(!std::same_as<std::decay_t<obj_t>, move_only_any>)
        : m_obj(std::malloc(sizeof(obj_t))),
          m_type(&typeid(obj_t)),
          m_destructor(get_destructor<obj_t>()) {
        new (m_obj) obj_t(std::forward<obj_t>(obj));
    }

    move_only_any(move_only_any &&other) noexcept
        : m_obj(other.m_obj),
          m_type(other.m_type),
          m_destructor(other.m_destructor) {
        other.m_obj = nullptr;
        other.m_type = nullptr;
        other.m_destructor = nullptr;
    }

    move_only_any &operator=(move_only_any &&other) noexcept {
        reset();

        m_obj = other.m_obj;
        other.m_obj = nullptr;
        m_type = other.m_type;
        other.m_type = nullptr;
        m_destructor = other.m_destructor;
        other.m_destructor = nullptr;

        return *this;
    }

    ~move_only_any() { reset(); }

    void reset() {
        if (m_obj != nullptr) {
            assert(m_destructor != nullptr);
            m_destructor(m_obj);
            std::free(m_obj);
        }
    }

    template <typename obj_t>
    void set(obj_t &&obj)
        requires(!std::same_as<std::decay_t<obj_t>, move_only_any>)
    {
        if (m_obj != nullptr) {
            assert(m_destructor != nullptr);
            m_destructor(m_obj);
            std::free(m_obj);
        }

        m_obj = std::malloc(sizeof(obj_t));
        new (m_obj) obj_t(std::forward<obj_t>(obj));

        m_type = &typeid(obj_t);
        m_destructor = get_destructor<obj_t>();
    }

    template <typename obj_t>
    bool is() const {
        if (m_type == nullptr) {
            return false;
        } else {
            return (*m_type == typeid(obj_t));
        }
    }

    bool has_value() const {
        if (m_type != nullptr) {
            assert(m_obj != nullptr && m_destructor != nullptr);
            return true;
        } else {
            assert(m_obj == nullptr && m_destructor == nullptr);
            return false;
        }
    }

    const std::type_info &type() const {
        if (!has_value()) {
            throw std::logic_error(
                "Type ID for `traccc::move_only_any` requested, but no value "
                "exists.");
        }

        return *m_type;
    }

    template <typename obj_t>
    obj_t &as() const {
        if (!has_value()) {
            throw std::logic_error(
                "Value for `traccc::move_only_any` requested, but no value "
                "exists.");
        }

        return *static_cast<obj_t *>(m_obj);
    }

    private:
    template <typename obj_t>
    void (*get_destructor() const)(void *) {
        return static_cast<void (*)(void *)>(
            [](void *ptr) { static_cast<obj_t *>(ptr)->~obj_t(); });
    }

    void *m_obj = nullptr;
    const std::type_info *m_type = nullptr;
    void (*m_destructor)(void *) = nullptr;
};

}  // namespace traccc
