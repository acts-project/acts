// STL include(s)
#include <atomic>
#include <mutex>
#include <functional>

// package include(s)
#include "AtsUtils/LazyInitialization.h"

// some (arbitrary) function returning an pointer to int
int* someFunction(int i);

class A
{
public:
  // Option 1: create new integer object and set it to 3
  int* getInt1() const
  {
    return Acts::LazyInit<int>::init(pInt,m,3);
  }

  // Option 2: set member variable to pointer returned by some other function
  int* getInt2() const
  {
    const std::function<int*(int)>& f = [](int i) -> int* {return someFunction(i);};
    return Acts::LazyInit<int>::init(pInt,m,f,4);
  }
  
private:
  // lazy-initialized member variable and mutex must be
  // mutable if initialized in a const method
  mutable std::atomic<int*> pInt {0};
  mutable std::mutex m;
};
