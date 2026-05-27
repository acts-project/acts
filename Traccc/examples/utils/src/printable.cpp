/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/examples/utils/printable.hpp"

// System include(s).
#include <sstream>

namespace traccc {
std::string configuration_printable::print() const {
    std::size_t mkw = get_max_key_width_impl();
    std::ostringstream oss;
    print_impl("", "", 0, mkw, oss);
    return oss.str();
};

configuration_printable::~configuration_printable() = default;

configuration_category::configuration_category(std::string_view n) : name(n) {}

void configuration_category::add_child(
    std::unique_ptr<configuration_printable> elem) {
    elements.push_back(std::move(elem));
}

void configuration_category::print_impl(const std::string& self_prefix,
                                        const std::string& child_prefix,
                                        std::size_t prefix_len,
                                        std::size_t max_key_width,
                                        std::ostream& out) const {
    out << self_prefix << name << ":\n";

    for (std::size_t i = 0; i < elements.size(); i++) {
        if (i == elements.size() - 1) {
            elements[i]->print_impl(child_prefix + "└ ", child_prefix + "  ",
                                    prefix_len + 2, max_key_width, out);
        } else {
            elements[i]->print_impl(child_prefix + "├ ", child_prefix + "│ ",
                                    prefix_len + 2, max_key_width, out);
        }
    }
}

std::size_t configuration_category::get_max_key_width_impl() const {
    std::size_t res = 0;

    for (auto& i : elements) {
        res = std::max(res, 2 + i->get_max_key_width_impl());
    }

    return res;
}

configuration_category::~configuration_category() = default;

configuration_list::configuration_list() = default;

void configuration_list::add_child(
    std::unique_ptr<configuration_printable> elem) {
    elements.push_back(std::move(elem));
}

void configuration_list::print_impl(const std::string& self_prefix,
                                    const std::string& child_prefix,
                                    std::size_t prefix_len,
                                    std::size_t max_key_width,
                                    std::ostream& out) const {
    for (auto& i : elements) {
        i->print_impl(self_prefix, child_prefix, prefix_len, max_key_width,
                      out);
    }
}

std::size_t configuration_list::get_max_key_width_impl() const {
    std::size_t res = 0;

    for (auto& i : elements) {
        res = std::max(res, i->get_max_key_width_impl());
    }

    return res;
}

configuration_list::~configuration_list() = default;

configuration_kv_pair::configuration_kv_pair(std::string_view k,
                                             std::string_view v)
    : key(k), value(v) {}

void configuration_kv_pair::print_impl(const std::string& self_prefix,
                                       const std::string&,
                                       std::size_t prefix_len,
                                       std::size_t max_key_width,
                                       std::ostream& out) const {
    out << self_prefix << key << ": "
        << std::string((max_key_width - (prefix_len + key.length())), ' ')
        << value << "\n";
}

std::size_t configuration_kv_pair::get_max_key_width_impl() const {
    return key.length();
}

configuration_kv_pair::~configuration_kv_pair() = default;
}  // namespace traccc
