#ifndef ATSUTILS_LAZYINITIALISATION
#define ATSUTILS_LAZYINITIALISATION 1

// STL include(s)
#include <utility>
#include <mutex>
#include <atomic>

namespace Ats
{
  template<typename T,typename Mutex = std::mutex>
  struct LazyInit
  {
    template<typename... Args>
    static T* init(std::atomic<T*>& atomic,Mutex& m,Args&&... args)
    {
      T* tmp = atomic.load(std::memory_order_acquire);
      if(!tmp)
      {
	std::lock_guard<Mutex> lock(m);
	tmp = atomic.load(std::memory_order_relaxed);
	if(!tmp)
	{
	  tmp = new T(std::forward<Args>(args)...);
	  atomic.store(tmp,std::memory_order_release);
	} // if not initialized
      } // end of DCLP

      return tmp;
    }
  };
  
} // end of namespace ATS

#endif // ATSUTILS_LAZYINITIALISATION
