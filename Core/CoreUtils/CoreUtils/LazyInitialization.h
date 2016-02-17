#ifndef ATSUTILS_LAZYINITIALISATION
#define ATSUTILS_LAZYINITIALISATION 1

// STL include(s)
#include <mutex>
#include <atomic>
#include <memory>

namespace Ats
{
  /**
   * @struct LazyInit
   *
   * @brief Helper struct for lazy and thread-safe initialization of pointer variables
   *
   * This helper struct should simplify the thread-safe lazy initialization of pointer
   * variables. It employs the double checked locking pattern (<a href="https://en.wikipedia.org/wiki/Double-checked_locking">DCLP</a>).
   *   * @attention Lazy initialization of @c std::shared_ptr may not be thread-safe due to implementation deficits in g++ 4.9.3.
   *
   * @todo The lazy initialization of @c std::shared_ptr is currently not thread-safe on
   *       CPUs with weakly-ordered memory instructions due to atomic operations for
   *       @c std::shared_ptr are not fully implemented in g++ 4.9.3 (<a href="https://gcc.gnu.org/onlinedocs/gcc-4.9.3/libstdc++/manual/manual/status.html#status.iso.2011">link</a>).
   *       This is fixed in g++ 5.3 and the code should be updated once we switch to this
   *       compiler version.
   *
   * @tparam T       type of the object which gets lazily-initialised
   * @tparam Mutex   type of mutex used for locking
   *
   * @example LazyInitExample.cxx Example for using the lazy initialization helper
   */
  template<typename T,typename Mutex = std::mutex>
  struct LazyInit
  {
    template<typename... Args>
    static T* init(std::atomic<T*>& atomic,Mutex& m,Args&&... args);

    template<typename... Args,typename FUNC>
    static T* init(std::atomic<T*>& atomic,Mutex& m,const FUNC& creator,Args&&... args);
    
    template<typename... Args>
    static std::shared_ptr<T>& init(std::shared_ptr<T>& rShared,Mutex& m,Args&&... args);

    template<typename... Args,typename FUNC>
    static std::shared_ptr<T>& init(std::shared_ptr<T>& rShared,Mutex& m,const FUNC& creator,Args&&... args);
  }; 
} // end of namespace Ats

#include "CoreUtils/LazyInitialization.icc"

#endif // ATSUTILS_LAZYINITIALISATION
