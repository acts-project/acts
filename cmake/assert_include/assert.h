#if defined(NDEBUG)
#define _HAD_NDEBUG
#undef NDEBUG
#endif
#include_next <assert.h>
#if defined(_HAD_NDEBUG)
#undef _HAD_NDEBUG
#define NDEBUG
#endif
