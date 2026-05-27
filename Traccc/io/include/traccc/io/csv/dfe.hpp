// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2015,2018-2020 Moritz Kiehn

#pragma once

#include <Acts/Utilities/Concepts.hpp>
#include <Acts/Utilities/Helpers.hpp>
#include <algorithm>
#include <array>
#include <cassert>
#include <fstream>
#include <ostream>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#define DFE_NAMEDTUPLE(name, ...)                                           \
    using Tuple = decltype(::std::make_tuple(__VA_ARGS__));                 \
    static ::std::array<::std::string, ::std::tuple_size<Tuple>::value>     \
    names() {                                                               \
        return ::traccc::io::dfe::unstringify<                              \
            ::std::tuple_size<Tuple>::value>((#__VA_ARGS__));               \
    }                                                                       \
    template <typename... U>                                                \
    name& operator=(const ::std::tuple<U...>& other) {                      \
        ::std::tie(__VA_ARGS__) = other;                                    \
        return *this;                                                       \
    }                                                                       \
    template <typename... U>                                                \
    name& operator=(::std::tuple<U...>&& other) {                           \
        ::std::tie(__VA_ARGS__) = ::std::forward<std::tuple<U...>>(other);  \
        return *this;                                                       \
    }                                                                       \
    operator Tuple() const {                                                \
        return ::std::make_tuple(__VA_ARGS__);                              \
    }                                                                       \
    Tuple tuple() const {                                                   \
        return ::std::make_tuple(__VA_ARGS__);                              \
    }                                                                       \
    template <std::size_t I>                                                \
    constexpr ::std::tuple_element_t<I, Tuple>& get() {                     \
        return ::std::get<I>(std::tie(__VA_ARGS__));                        \
    }                                                                       \
    template <::std::size_t I>                                              \
    constexpr const ::std::tuple_element_t<I, Tuple>& get() const {         \
        return ::std::get<I>(std::tie(__VA_ARGS__));                        \
    }                                                                       \
    template <::std::size_t I>                                              \
    friend constexpr ::std::tuple_element_t<I, Tuple>& get(name& nt) {      \
        return nt.template get<I>();                                        \
    }                                                                       \
    template <::std::size_t I>                                              \
    friend constexpr const ::std::tuple_element_t<I, Tuple>& get(           \
        const name& nt) {                                                   \
        return nt.template get<I>();                                        \
    }                                                                       \
    friend inline ::std::ostream& operator<<(::std::ostream& os,            \
                                             const name& nt) {              \
        return ::traccc::io::dfe::print_tuple(                              \
            os, nt.names(), nt.tuple(),                                     \
            ::std::make_index_sequence<::std::tuple_size<Tuple>::value>{}); \
    }

namespace traccc::io {
namespace dfe {

template <std::size_t N>
constexpr std::array<std::string, N> unstringify(const char* str) {
    assert(str && "Input string must be non-null");

    std::array<std::string, N> out;

    for (std::size_t idx = 0; idx < N; ++idx) {
        while ((*str != '\0') && (*str == ' ')) {
            ++str;
        }
        const char* sep = str;
        while ((*sep != '\0') && (*sep != ',')) {
            ++sep;
        }
        out[idx].assign(str, static_cast<std::size_t>(sep - str));
        if (*sep == '\0') {
            break;
        }
        str = ++sep;
    }
    return out;
}

template <typename Names, typename Values, std::size_t... I>
inline std::ostream& print_tuple(std::ostream& os, const Names& n,
                                 const Values& v,
                                 std::index_sequence<I...> /*seq*/) {
    using std::get;
    using Vacuum = int[];
    (void)Vacuum{
        (os << ((0 < I) ? " " : "") << get<I>(n) << "=" << get<I>(v), 0)...};
    return os;
}

template <char Delimiter>
class DsvWriter {
    public:
    DsvWriter() = delete;
    DsvWriter(const DsvWriter&) = delete;
    DsvWriter(DsvWriter&&) = default;
    ~DsvWriter() = default;
    DsvWriter& operator=(const DsvWriter&) = delete;
    DsvWriter& operator=(DsvWriter&&) = default;

    DsvWriter(const std::vector<std::string>& columns, const std::string& path,
              int precision = std::numeric_limits<double>::max_digits10);

    template <typename Arg0, typename... Args>
    void append(Arg0&& arg0, Args&&... args);

    private:
    std::ofstream m_file;
    std::size_t m_num_columns;

    template <typename T>
    unsigned write(T&& x, std::ostream& os);

    template <typename T, typename Allocator>
    static unsigned write(const std::vector<T, Allocator>& xs,
                          std::ostream& os);
};

template <char Delimiter>
class DsvReader {
    public:
    DsvReader() = delete;
    DsvReader(const DsvReader&) = delete;
    DsvReader(DsvReader&&) = default;
    ~DsvReader() = default;
    DsvReader& operator=(const DsvReader&) = delete;
    DsvReader& operator=(DsvReader&&) = default;

    explicit DsvReader(const std::string& path);

    bool read(std::vector<std::string>& columns);

    std::size_t num_lines() const { return m_num_lines; }

    private:
    std::ifstream m_file;
    std::string m_line;
    std::size_t m_num_lines = 0;
};

template <char Delimiter, typename NamedTuple>
class NamedTupleDsvWriter {
    public:
    NamedTupleDsvWriter() = delete;
    NamedTupleDsvWriter(const NamedTupleDsvWriter&) = delete;
    NamedTupleDsvWriter(NamedTupleDsvWriter&&) = default;
    ~NamedTupleDsvWriter() = default;
    NamedTupleDsvWriter& operator=(const NamedTupleDsvWriter&) = delete;
    NamedTupleDsvWriter& operator=(NamedTupleDsvWriter&&) = default;

    NamedTupleDsvWriter(
        const std::string& path,
        int precision = std::numeric_limits<double>::max_digits10)
        : m_writer(colum_names(), path, precision) {}

    void append(const NamedTuple& record) {
        append_impl(record,
                    std::make_index_sequence<
                        std::tuple_size_v<typename NamedTuple::Tuple>>{});
    }

    private:
    DsvWriter<Delimiter> m_writer;

    static std::vector<std::string> colum_names() {
        const auto& from_record = NamedTuple::names();
        return {from_record.begin(), from_record.end()};
    }
    template <std::size_t... I>
    void append_impl(const NamedTuple& values,
                     std::index_sequence<I...> /*seq*/) {
        using std::get;
        m_writer.append(get<I>(values)...);
    }
};

template <typename T>
static void parse(const std::string& str, T& value) {
    std::istringstream is(str);
    is >> value;
}

template <char Delimiter, typename NamedTuple>
class NamedTupleDsvReader {
    public:
    NamedTupleDsvReader() = delete;
    NamedTupleDsvReader(const NamedTupleDsvReader&) = delete;
    NamedTupleDsvReader(NamedTupleDsvReader&&) = default;
    ~NamedTupleDsvReader() = default;
    NamedTupleDsvReader& operator=(const NamedTupleDsvReader&) = delete;
    NamedTupleDsvReader& operator=(NamedTupleDsvReader&&) = default;

    NamedTupleDsvReader(const std::string& path,
                        const std::vector<std::string>& optional_columns = {},
                        bool verify_header = true);
    bool read(NamedTuple& record);

    template <typename T>
    bool read(NamedTuple& record, std::vector<T>& extra);

    std::size_t num_extra_columns() const { return m_extra_columns.size(); }
    std::size_t num_records() const { return m_reader.num_lines() - 1u; }

    private:
    using Tuple = typename NamedTuple::Tuple;

    DsvReader<Delimiter> m_reader;
    std::vector<std::string> m_columns;
    std::size_t m_num_columns = SIZE_MAX;
    std::array<std::size_t, std::tuple_size<Tuple>::value> m_tuple_column_map;
    std::vector<std::size_t> m_extra_columns;

    void use_default_columns();
    void parse_header(const std::vector<std::string>& optional_columns);
    template <std::size_t... I>
    void parse_record(NamedTuple& record,
                      std::index_sequence<I...> /*seq*/) const {
        using Vacuum = int[];
        (void)Vacuum{(parse_element<I>(record), 0)...};
    }
    template <std::size_t I>
    void parse_element(NamedTuple& record) const {
        using std::get;
        if (m_tuple_column_map[I] != SIZE_MAX) {
            parse(m_columns[m_tuple_column_map[I]], get<I>(record));
        }
    }
};

template <char Delimiter>
inline DsvWriter<Delimiter>::DsvWriter(const std::vector<std::string>& columns,
                                       const std::string& path, int precision)
    : m_file(path,
             std::ios_base::binary | std::ios_base::out | std::ios_base::trunc),
      m_num_columns(columns.size()) {
    if (!m_file.is_open() || m_file.fail()) {
        throw std::runtime_error("Could not open file '" + path + "'");
    }
    m_file.precision(precision);
    if (m_num_columns == 0) {
        throw std::invalid_argument("No columns were specified");
    }
    append(columns);
}

template <char Delimiter>
template <typename Arg0, typename... Args>
inline void DsvWriter<Delimiter>::append(Arg0&& arg0, Args&&... args) {
    std::stringstream line;
    unsigned written_columns[] = {
        write(std::forward<Arg0>(arg0), line),
        (line << Delimiter, write(std::forward<Args>(args), line))...,
    };
    line << '\n';
    unsigned total_columns = 0;
    for (auto nc : written_columns) {
        total_columns += nc;
    }
    if (total_columns < m_num_columns) {
        throw std::invalid_argument("Not enough columns");
    }
    if (m_num_columns < total_columns) {
        throw std::invalid_argument("Too many columns");
    }
    m_file << line.rdbuf();
    if (!m_file.good()) {
        throw std::runtime_error("Could not write data to file");
    }
}

template <char Delimiter>
template <typename T>
inline unsigned DsvWriter<Delimiter>::write(T&& x, std::ostream& os) {
    os << x;
    return 1u;
}

template <char Delimiter>
template <typename T, typename Allocator>
inline unsigned DsvWriter<Delimiter>::write(const std::vector<T, Allocator>& xs,
                                            std::ostream& os) {
    unsigned n = 0;
    for (const auto& x : xs) {
        if (0 < n) {
            os << Delimiter;
        }
        os << x;
        n += 1;
    }
    return n;
}

template <char Delimiter>
inline DsvReader<Delimiter>::DsvReader(const std::string& path)
    : m_file(path, std::ios_base::binary | std::ios_base::in) {
    if (!m_file.is_open() || m_file.fail()) {
        throw std::runtime_error("Could not open file '" + path + "'");
    }
}

template <char Delimiter>
inline bool DsvReader<Delimiter>::read(std::vector<std::string>& columns) {
    std::getline(m_file, m_line);
    if (m_file.eof()) {
        return false;
    }
    if (m_file.fail()) {
        throw std::runtime_error(std::string("Could not read line ") +
                                 std::to_string(m_num_lines));
    }
    m_num_lines += 1;

    columns.clear();
    for (std::string::size_type pos = 0; pos < m_line.size();) {
        auto del = m_line.find_first_of(Delimiter, pos);
        if (del == std::string::npos) {
            columns.emplace_back(m_line, pos);
            break;
        } else {
            columns.emplace_back(m_line, pos, del - pos);
            pos = del + 1;
        }
    }
    return true;
}

template <char Delimiter, typename NamedTuple>
inline NamedTupleDsvReader<Delimiter, NamedTuple>::NamedTupleDsvReader(
    const std::string& path, const std::vector<std::string>& optional_columns,
    bool verify_header)
    : m_reader(path) {
    if ((!optional_columns.empty()) && (!verify_header)) {
        throw std::runtime_error(
            "Optional columns can not be used without header verification");
    }
    if (!m_reader.read(m_columns)) {
        throw std::runtime_error("Could not read header from '" + path + "'");
    }
    if (verify_header) {
        parse_header(optional_columns);
    } else {
        use_default_columns();
    }
}

template <char Delimiter, typename NamedTuple>
inline bool NamedTupleDsvReader<Delimiter, NamedTuple>::read(
    NamedTuple& record) {
    if (!m_reader.read(m_columns)) {
        return false;
    }
    if (m_columns.size() < m_num_columns) {
        throw std::runtime_error("Too few columns in line " +
                                 std::to_string(m_reader.num_lines()));
    }
    if (m_num_columns < m_columns.size()) {
        throw std::runtime_error("Too many columns in line " +
                                 std::to_string(m_reader.num_lines()));
    }
    parse_record(record, std::make_index_sequence<std::tuple_size_v<Tuple>>{});
    return true;
}

template <char Delimiter, typename NamedTuple>
template <typename T>
inline bool NamedTupleDsvReader<Delimiter, NamedTuple>::read(
    NamedTuple& record, std::vector<T>& extra) {
    if (!read(record)) {
        return false;
    }
    extra.resize(m_extra_columns.size());
    for (std::size_t i = 0; i < m_extra_columns.size(); ++i) {
        parse(m_columns[m_extra_columns[i]], extra[i]);
    }
    return true;
}

template <char Delimiter, typename NamedTuple>
inline void NamedTupleDsvReader<Delimiter, NamedTuple>::use_default_columns() {
    m_num_columns = std::tuple_size_v<Tuple>;
    for (std::size_t i = 0; i < m_tuple_column_map.size(); ++i) {
        m_tuple_column_map[i] = i;
    }
    m_extra_columns.clear();
}

template <char Delimiter, typename NamedTuple>
inline void NamedTupleDsvReader<Delimiter, NamedTuple>::parse_header(
    const std::vector<std::string>& optional_columns) {
    const auto& names = NamedTuple::names();

    m_num_columns = m_columns.size();

    for (const auto& name : names) {
        if (::Acts::rangeContainsValue(optional_columns, name)) {
            continue;
        }
        if (!::Acts::rangeContainsValue(m_columns, name)) {
            throw std::runtime_error("Missing header column '" + name + "'");
        }
    }

    m_tuple_column_map.fill(SIZE_MAX);

    m_extra_columns.clear();
    for (std::size_t i = 0; i < m_columns.size(); ++i) {
        auto it = std::ranges::find(names, m_columns[i]);
        if (it != names.end()) {
            m_tuple_column_map[static_cast<std::size_t>(
                std::distance(names.begin(), it))] = i;
        } else {
            m_extra_columns.push_back(i);
        }
    }
}

using CsvWriter = dfe::DsvWriter<','>;

using TsvWriter = dfe::DsvWriter<'\t'>;

template <typename T>
using NamedTupleCsvWriter = dfe::NamedTupleDsvWriter<',', T>;

template <typename T>
using NamedTupleCsvReader = dfe::NamedTupleDsvReader<',', T>;

template <typename T>
using NamedTupleTsvWriter = dfe::NamedTupleDsvWriter<'\t', T>;

template <typename T>
using NamedTupleTsvReader = dfe::NamedTupleDsvReader<'\t', T>;
}  // namespace dfe
}  // namespace traccc::io
