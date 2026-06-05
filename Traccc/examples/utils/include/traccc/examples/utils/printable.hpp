/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// System include(s).
#include <iosfwd>
#include <memory>
#include <string>
#include <string_view>
#include <vector>

namespace traccc {
class configuration_printable {
    public:
    virtual std::string print() const;

    virtual void print_impl(const std::string& self_prefix,
                            const std::string& child_prefix,
                            std::size_t prefix_len, std::size_t max_key_width,
                            std::ostream& out) const = 0;

    virtual std::size_t get_max_key_width_impl() const = 0;

    virtual ~configuration_printable();
};

class configuration_category final : public configuration_printable {
    public:
    explicit configuration_category(std::string_view n);

    void add_child(std::unique_ptr<configuration_printable> elem);

    void print_impl(const std::string& self_prefix,
                    const std::string& child_prefix, std::size_t prefix_len,
                    std::size_t max_key_width,
                    std::ostream& out) const override;

    std::size_t get_max_key_width_impl() const override;

    ~configuration_category() override;

    private:
    std::string name;
    std::vector<std::unique_ptr<configuration_printable>> elements;
};

class configuration_list final : public configuration_printable {
    public:
    configuration_list();

    void add_child(std::unique_ptr<configuration_printable> elem);

    void print_impl(const std::string& self_prefix,
                    const std::string& child_prefix, std::size_t prefix_len,
                    std::size_t max_key_width,
                    std::ostream& out) const override;

    std::size_t get_max_key_width_impl() const override;

    ~configuration_list() override;

    private:
    std::vector<std::unique_ptr<configuration_printable>> elements;
};

class configuration_kv_pair final : public configuration_printable {
    public:
    configuration_kv_pair(std::string_view k, std::string_view v);

    void print_impl(const std::string& self_prefix, const std::string&,
                    std::size_t prefix_len, std::size_t max_key_width,
                    std::ostream& out) const override;

    std::size_t get_max_key_width_impl() const override;

    ~configuration_kv_pair() override;

    private:
    std::string key;
    std::string value;
};
}  // namespace traccc
