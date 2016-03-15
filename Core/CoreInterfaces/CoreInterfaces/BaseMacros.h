///////////////////////////////////////////////////////////////////
// BaseMacros.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_CORE_BASEMACROS_H
#define ACTS_CORE_BASEMACROS_H 1

// these macros simplify/harmonize the retrieval ot TOOLS/SERVICES
#define RETRIEVE_FATAL(handle)                                                                          \
        if (handle.empty()){                                                                            \
            MSG_FATAL("<retrieve> Empty handle for crucial component provided. Arborting.");            \
            return StatusCode::FAILURE;                                                                 \
        } else if (handle.retrieve().isFailure()){                                                  \
            MSG_FATAL("<retrieve> Failed to retrieve InterfaceHandle '" << handle << "'. Arborting.");   \
            return StatusCode::FAILURE;                                                                 \
        } else                                                                                          \
            MSG_DEBUG("<retrieve> Successfully retrieved InterfaceHandle '" << handle <<"'.");                     

// these macros simplify/harmonize the retrieval ot TOOLS/SERVICES
#define RETRIEVE_NONEMPTY_FATAL(handle)                                                                            \
        if (!handle.empty() && handle.retrieve().isFailure()){                                                     \
            MSG_FATAL("<retrieve> Failed to retrieve non-empty InterfaceHandle '" << handle << "'. Arborting.");   \
            return StatusCode::FAILURE;                                                                            \
        } else if (!handle.empty()){                                                                               \
            MSG_DEBUG("<retrieve> Successfully retrieved InterfaceHandle '" << handle <<";.");                     \
        } else                                                                                                     \
            MSG_DEBUG("<retrieve> Empty InterfaceHandle, no retrieval attempted.");                                

#endif
