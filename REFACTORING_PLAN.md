# Detray Conversion Refactoring Plan

## Overview

This document outlines the plan to refactor the current virtual `toDetrayPayload` functions in navigation policies and material classes. The current approach requires forward declarations and creates header dependencies that are problematic when the detray plugin is not available.

## Current Problem

- Functions like `toDetrayPayload` are defined in `Core/src/Navigation/DetrayNavigation.cpp` and `Core/src/Material/DetrayMaterial.cpp`
- These are declared as virtual functions in class headers
- Interface `detray` types are not available in headers, requiring forward declarations in `Core/include/Acts/Geometry/DetrayFwd.hpp`
- Stub implementations throw errors when detray plugin is not available

## Proposed Solution

Replace virtual functions with a generic type-based function dispatch system using function pointers registered at runtime.

## Implementation Plan

### 1. **Generic TypeDispatcher Implementation** ✅

Created `Core/include/Acts/Utilities/TypeDispatcher.hpp` with a reusable template class:

**Key Features:**
- Template parameters: `BaseType` and function signature  
- Register functions by derived type with automatic `dynamic_cast`
- Type-safe dispatch using `std::type_index`
- Support for free functions, lambdas, and function objects

**Usage Example:**
```cpp
TypeDispatcher<BaseClass, ReturnType(Args...)> dispatcher;

// Register function for specific derived type
dispatcher.registerFunction<DerivedType>(myFunction);

// Call with base class reference
ReturnType result = dispatcher(baseObject, args...);
```

Created comprehensive unit tests in `Tests/UnitTests/Core/Utilities/TypeDispatcherTests.cpp`.

### 2. **Create Detray Conversion Dispatchers**

Create global dispatchers for the conversion functions:

```cpp
// In a new header like Acts/Navigation/DetrayConversion.hpp
extern TypeDispatcher<NavigationPolicyBase, 
                     std::unique_ptr<DetraySurfaceGrid>(const SurfaceLookupFunction&, const Logger&)> 
       navigationToDetrayDispatcher;
```

```cpp
// For material conversions
extern TypeDispatcher<ISurfaceMaterial,
                     std::unique_ptr<DetraySurfaceMaterial>()>
       materialToDetrayDispatcher;
```

### 3. **Remove Virtual Functions from Headers**

- Remove `toDetrayPayload` declarations from navigation policy headers
- Remove `toDetrayPayload` declarations from material class headers  
- Keep the forward declarations in `DetrayFwd.hpp` for the return types

### 4. **Register Functions in Plugin**

In the detray plugin initialization (when detray is available):

```cpp
// Register all the conversion functions
navigationToDetrayDispatcher.registerFunction<MultiNavigationPolicy>(convertMultiNavigationPolicy);
navigationToDetrayDispatcher.registerFunction<SurfaceArrayNavigationPolicy>(convertSurfaceArrayNavigationPolicy);
navigationToDetrayDispatcher.registerFunction<MultiLayerNavigationPolicy>(convertMultiLayerNavigationPolicy);
navigationToDetrayDispatcher.registerFunction<TryAllNavigationPolicy>(convertTryAllNavigationPolicy);
navigationToDetrayDispatcher.registerFunction<CylinderNavigationPolicy>(convertCylinderNavigationPolicy);

materialToDetrayDispatcher.registerFunction<BinnedSurfaceMaterial>(convertBinnedSurfaceMaterial);
materialToDetrayDispatcher.registerFunction<HomogeneousSurfaceMaterial>(convertHomogeneousSurfaceMaterial);
// etc.
```

### 5. **Update Call Sites**

Replace current virtual function calls:
```cpp
// Old approach
auto payload = policy->toDetrayPayload(surfaceLookup, logger);
```

With dispatcher calls:
```cpp
// New approach  
auto payload = navigationToDetrayDispatcher(*policy, surfaceLookup, logger);
```

Add proper error handling for unregistered types.

### 6. **Provide Stub Implementation**

When detray plugin is not available:
- The dispatchers remain empty (no functions registered)
- Calling them will throw `std::runtime_error` with a helpful message like: 
  `"No detray conversion function registered for type: <TypeName>. Ensure detray plugin is enabled."`
- This replaces the current stub implementations that throw

## Benefits of This Approach

1. **No virtual functions needed** - eliminates the dependency issue
2. **Type-safe** - compile-time checks ensure correct function signatures and type relationships
3. **Extensible** - easy to add new conversions without modifying base classes  
4. **Clean separation** - conversion logic lives entirely with the detray plugin
5. **Flexible** - supports free functions, lambdas, and function objects
6. **Testable** - each component can be tested independently
7. **Runtime registration** - conversions are only available when detray plugin is loaded

## Implementation Steps

1. ✅ Implement and test `TypeDispatcher` utility class
2. Create detray conversion dispatcher declarations
3. Move existing conversion functions to free functions  
4. Remove virtual `toDetrayPayload` methods from class headers
5. Update detray plugin to register conversion functions
6. Update call sites to use dispatchers
7. Test with and without detray plugin enabled

## Files to Modify

- **New Files:**
  - `Core/include/Acts/Utilities/TypeDispatcher.hpp` ✅
  - `Tests/UnitTests/Core/Utilities/TypeDispatcherTests.cpp` ✅  
  - `Core/include/Acts/Navigation/DetrayConversion.hpp`
  - `Core/include/Acts/Material/DetrayConversion.hpp`

- **Modified Files:**
  - `Core/src/Navigation/DetrayNavigation.cpp` - convert to free functions
  - `Core/src/Material/DetrayMaterial.cpp` - convert to free functions  
  - Navigation policy headers - remove virtual functions
  - Material class headers - remove virtual functions
  - Detray plugin initialization code - register functions
  - Call sites using `toDetrayPayload` - update to use dispatchers

## Migration Path

This refactoring can be done incrementally:
1. Implement the infrastructure (steps 1-2)
2. Add dispatcher calls alongside existing virtual calls
3. Register functions in detray plugin  
4. Switch call sites to use dispatchers
5. Remove virtual function declarations
6. Clean up old implementation code