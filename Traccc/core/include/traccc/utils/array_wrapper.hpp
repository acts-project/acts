/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <traccc/definitions/qualifiers.hpp>
#include <traccc/utils/functor.hpp>
#include <type_traits>
#include <utility>
#include <vecmem/memory/memory_resource.hpp>
#include <vecmem/memory/unique_ptr.hpp>

namespace traccc {
/**
 * @brief Tuple of plain-old-data types.
 *
 * This is just a super-simple re-implementation of `std::tuple` that ensures
 * that the type is trivially constructible.
 */
template <typename... Ts>
struct pod {};

template <>
struct pod<> {};

template <typename T, typename... Ts>
struct pod<T, Ts...> {
    static_assert(std::is_standard_layout_v<T> && std::is_trivial_v<T>,
                  "Types in POD-tuple must each be POD.");

    T v;
    pod<Ts...> r;
};

/**
 * @brief Get a given value from a POD tuple.
 *
 * @tparam I The index in the POD type to retrieve.
 * @tparam Ts The types of the POD tuple.
 *
 * @param p A valid of the POD tuple type described by `Ts...`.
 *
 * @return A value of type `Ts[I]` as contained in `p`.
 */
template <std::size_t I, typename... Ts>
constexpr auto& pod_get(pod<Ts...>& p) {
    if constexpr (I == 0) {
        return p.v;
    } else {
        return pod_get<I - 1>(p.r);
    }
}

/**
 * @brief An array wrapper that can swap between SoA and AoS layouts almost
 * transparently.
 *
 * @tparam F The higher-order layout type, either `soa` or `aos`.
 * @tparam T The struct type to embed in memory.
 */
template <template <typename...> typename F,
          template <template <typename> typename> typename T>
struct array_wrapper {
    struct owner {
        owner(vecmem::memory_resource& mr, std::size_t n) : data(mr, n) {}

        typename details::functor::reapply<
            F, typename T<details::functor::identity>::tuple_t>::type::owner
            data;
    };

    struct handle {
        using handle_t = typename details::functor::reapply<
            F, typename T<details::functor::identity>::tuple_t>::type::handle;

        handle(const owner& o) : data(o.data) {}

        TRACCC_HOST_DEVICE std::size_t size() const { return data.size(); }

        template <std::size_t I>
        constexpr TRACCC_HOST_DEVICE auto& get(std::size_t i) {
            return data.get<I>(i);
        }

        template <std::size_t I>
        constexpr TRACCC_HOST_DEVICE auto get(std::size_t i) const {
            return data.get<I>(i);
        }

        template <std::size_t... Ns>
        constexpr TRACCC_HOST_DEVICE T<details::functor::identity>
        _construct_helper_identity(std::index_sequence<Ns...>,
                                   std::size_t i) const {
            return T<details::functor::identity>{get<Ns>(i)...};
        }

        template <std::size_t... Ns>
        constexpr TRACCC_HOST_DEVICE T<details::functor::reference>
        _construct_helper_reference(std::index_sequence<Ns...>, std::size_t i) {
            return T<details::functor::reference>{get<Ns>(i)...};
        }

        constexpr TRACCC_HOST_DEVICE T<details::functor::identity> operator[](
            std::size_t i) const {
            return _construct_helper_identity(
                std::make_index_sequence<std::tuple_size_v<
                    typename T<details::functor::identity>::tuple_t>>(),
                i);
        }

        constexpr TRACCC_HOST_DEVICE T<details::functor::reference> operator[](
            std::size_t i) {
            return _construct_helper_reference(
                std::make_index_sequence<std::tuple_size_v<
                    typename T<details::functor::identity>::tuple_t>>(),
                i);
        }

        handle_t data;
    };
};

/**
 * @brief Helper function to convert a list of unique pointers to raw pointers.
 */
template <std::size_t... Ns, typename... Ts>
std::tuple<Ts*...> _get_ptrs(
    std::index_sequence<Ns...>,
    const std::tuple<vecmem::unique_alloc_ptr<Ts[]>...>& o) {
    return {std::get<Ns>(o).get()...};
}

/**
 * @brief A struct-of-arrays (SoA) layout scheme.
 *
 * Struct-of-arrays is a common approach in vectorized programming both on CPUs
 * as well as CPUs. If only parts of the struct are accessed at the same time,
 * an SoA layout may increase locality and allow for coalescing of accesses.
 *
 * @tparam Ts The types contained in the POD struct.
 */
template <typename... Ts>
struct soa {
    struct owner {
        owner(vecmem::memory_resource& mr, std::size_t n)
            : _size(n), _ptrs{vecmem::make_unique_alloc<Ts[]>(mr, n)...} {}

        constexpr std::size_t size() const { return _size; }

        const std::tuple<vecmem::unique_alloc_ptr<Ts[]>...>& pointers() const {
            return _ptrs;
        }

        private:
        std::size_t _size;
        std::tuple<vecmem::unique_alloc_ptr<Ts[]>...> _ptrs;
    };

    struct handle {
        private:
        using tuple_t = std::tuple<Ts*...>;

        public:
        handle(const owner& o)
            : _size(o.size()),
              _ptrs(_get_ptrs(std::make_index_sequence<sizeof...(Ts)>(),
                              o.pointers())) {}

        constexpr TRACCC_HOST_DEVICE std::size_t size() const { return _size; }

        template <std::size_t I>
        constexpr TRACCC_HOST_DEVICE auto& get(std::size_t i) {
            return std::get<I>(_ptrs)[i];
        }

        template <std::size_t I>
        constexpr TRACCC_HOST_DEVICE const auto& get(std::size_t i) const {
            return std::get<I>(_ptrs)[i];
        }

        private:
        std::size_t _size;
        std::tuple<Ts*...> _ptrs;
    };
};

/**
 * @brief An array-of-structs (Aos) layout scheme.
 *
 * Array-of-structs is the default layout scheme in C++. It is useful in cases
 * where all elements of a struct need to be accessed simultaneously.
 *
 * @tparam Ts The types contained in the POD struct.
 */
template <typename... Ts>
struct aos {
    struct owner {
        constexpr owner(vecmem::memory_resource& mr, std::size_t n)
            : _size(n), _ptr{vecmem::make_unique_alloc<pod<Ts...>[]>(mr, n)} {}

        constexpr std::size_t size() const { return _size; }

        constexpr const vecmem::unique_alloc_ptr<pod<Ts...>[]>& pointer()
            const {
            return _ptr;
        }

        private:
        std::size_t _size;
        vecmem::unique_alloc_ptr<pod<Ts...>[]> _ptr;
    };

    struct handle {
        private:
        using tuple_t = pod<Ts...>;

        public:
        constexpr handle(const owner& o)
            : _size(o.size()), _ptr(o.pointer().get()) {}

        constexpr TRACCC_HOST_DEVICE std::size_t size() const { return _size; }

        template <std::size_t I>
        constexpr TRACCC_HOST_DEVICE auto& get(std::size_t i) {
            return pod_get<I>(_ptr[i]);
        }

        template <std::size_t I>
        constexpr TRACCC_HOST_DEVICE const auto& get(std::size_t i) const {
            return pod_get<I>(_ptr[i]);
        }

        private:
        std::size_t _size;
        tuple_t* _ptr;
    };
};
}  // namespace traccc
