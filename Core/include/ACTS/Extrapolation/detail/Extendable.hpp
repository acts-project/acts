#ifndef ACTS_EXTENDABLE_HPP
#define ACTS_EXTENDABLE_HPP 1

#include <tuple>

namespace Acts {

namespace detail {

  template <typename... Extensions>
  struct Extendable
  {
    template <typename R>
    const R&
    get() const
    {
      return std::get<R>(m_tExtensions);
    }

    template <typename R>
    R&
    get()
    {
      return std::get<R>(m_tExtensions);
    }

    const std::tuple<Extensions...>&
    tuple() const
    {
      return m_tExtensions;
    }

    std::tuple<Extensions...>&
    tuple()
    {
      return m_tExtensions;
    }

  private:
    std::tuple<Extensions...> m_tExtensions;
  };

}  // namespace detail

}  // namespace Acts
#endif  // ACTS_EXTENDABLE_HPP
