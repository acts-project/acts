#ifndef ACTS_MSG_MACROS_H
#define ACTS_MSG_MACROS_H 1

#ifdef ACTS_MSG_MACROS_PLUGIN
#include ACTS_MSG_MACROS_PLUGIN
#else
  static_assert(false,"no message macros defined");
#endif

#endif // ACTS_MSG_MACROS_H
