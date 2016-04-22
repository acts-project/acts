#ifndef ACTS_DEFAULT_MSG_MACROS_H
#define ACTS_DEFAULT_MSG_MACROS_H 1

// STL include(s)
#include <iostream>

namespace Acts
{
#define MSG_VERBOSE(x) std::cout << "VERBOSE " << (x) << std::endl;
#define MSG_DEBUG(x)   std::cout << "DEBUG   " << (x) << std::endl;
#define MSG_WARNING(x) std::cout << "WARNING " << (x) << std::endl;
#define MSG_INFO(x)    std::cout << "INFO    " << (x) << std::endl;
#define MSG_FATAL(x)   std::cout << "FATAL   " << (x) << std::endl;
}  // end of namespace Acts

#endif // ACTS_DEFAULT_MSG_MACROS_H
