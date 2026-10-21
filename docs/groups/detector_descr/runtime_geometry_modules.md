@defgroup geometry_module_loading Runtime geometry module loading
@ingroup geometry
@brief Loading tracking geometries from runtime shared libraries.

This topic documents the runtime geometry module loading entry points:

- @ref Acts::loadGeometryModule for modules that do not require user data.
- @ref Acts::loadDD4hepGeometryModule for DD4hep-based modules.

Both loaders validate module ABI compatibility before constructing the
@ref Acts::TrackingGeometry, and keep the shared library loaded for the
lifetime of the returned geometry object.

## Overview

Runtime geometry modules allow tracking geometries to be compiled as
independent shared libraries (`.so`/`.dylib`) and loaded at runtime without
recompiling the host application. This is useful when:

- The geometry changes frequently and recompiling the full framework is costly.
- Different geometry variants need to be swapped at runtime.
- The geometry is developed and distributed independently from the experiment
  framework.

@note Runtime geometry modules rely on `dlopen`/`dlsym` and are
      only supported on Unix-like systems (Linux, macOS).

## The module ABI

Every geometry module must export a single C-linkage entry point:

```c
extern "C" const ActsGeometryModuleV1* acts_geometry_module_v1(void);
```

The returned `ActsGeometryModuleV1` struct (defined in
`Acts/Geometry/GeometryModule.h`) carries four fields:

| Field              | Type                                      | Purpose                                                  |
|--------------------|-------------------------------------------|----------------------------------------------------------|
| `module_abi_tag`   | `const char*`                             | Build-time ABI tag matched against the host library      |
| `user_data_type`   | `const char*`                             | Type name of the extra context the module needs, or null |
| `build`            | `void* (*)(const void* user_data, const void* logger)` | Constructs the geometry; returns a heap-allocated `TrackingGeometry*` or null on failure |
| `destroy`          | `void (*)(void* handle)`                  | Deletes the geometry object returned by `build`          |

You never fill this struct manually. Use the provided helper macros instead
(see @ref geometry_module_writing "Writing a geometry module" below).

## ABI compatibility

Each ACTS build is tagged with an opaque string (`ACTS_GEOMETRY_MODULE_ABI_TAG`)
that encodes the host library ABI. Plain geometry modules use that tag as-is
and are matched against the tag compiled into `ActsCore`.

DD4hep geometry modules use a DD4hep-specific extension of that tag:

```text
${ACTS_GEOMETRY_MODULE_ABI_TAG}|dd4hep-${DD4hep_VERSION}
```

`loadDD4hepGeometryModule` matches `module_abi_tag` against the corresponding
tag compiled into `Acts::PluginDD4hep`, so a DD4hep module must be built
against both the same ACTS build and the same DD4hep version as the host
plugin. Any mismatch causes an immediate `std::runtime_error`.

## Loading a module

### Plain module (no extra context)

Include `Acts/Geometry/GeometryModuleLoader.hpp`, then:

@snippet{trimleft} examples/geometry_module.cpp Load Plain Module

### DD4hep module

Include `ActsPlugins/DD4hep/GeometryModuleLoader.hpp`, then:

@snippet{trimleft} examples/geometry_module.cpp Load DD4hep Module

The loader passes a pointer to `detector` through the opaque `void* user_data`
argument. It validates that the module declares `user_data_type == "dd4hep::Detector"`;
trying to load a plain module with `loadDD4hepGeometryModule` (or vice-versa)
throws a descriptive `std::runtime_error`.

@anchor geometry_module_writing
## Writing a geometry module

### Plain module

1. Create a source file with a build function:

@snippet{trimleft} examples/geometry_module_template.cpp Write Plain Module

2. Register it in CMake using `acts_add_geometry_module`:

@snippet{trimleft} CMakeLists.txt Plain Module CMake

The CMake helper creates a `SHARED` library target, links it against
`Acts::Core`, and injects the `ACTS_GEOMETRY_MODULE_ABI_TAG` compile definition
required by `ACTS_DEFINE_GEOMETRY_MODULE`.

### DD4hep module

1. Create a source file with a build function that takes a `dd4hep::Detector`:

@snippet{trimleft} examples/geometry_module_template_dd4hep.cpp Write DD4hep Module

2. Register it in CMake using `acts_add_dd4hep_geometry_module`:

@snippet{trimleft} CMakeLists.txt DD4hep Module CMake

This links against `Acts::PluginDD4hep` instead of `Acts::Core` and sets
`ACTS_GEOMETRY_MODULE_ABI_TAG` to the DD4hep-specific tag derived from the
installed ACTS tag plus `DD4hep_VERSION`.

## Lifetime management

The `std::shared_ptr<TrackingGeometry>` returned by both loaders uses a custom
deleter that:

1. Calls `descriptor->destroy()` to delete the `TrackingGeometry` object via
   the module's own destructor (important for correct cross-boundary `delete`).
2. Keeps a `shared_ptr<void>` to the `dlopen` handle alive, so the shared
   library is only unloaded **after** the geometry is destroyed — preventing
   use-after-unload of virtual dispatch tables or static data.

## Error handling

All errors are reported as `std::runtime_error`. The loader performs several
checks in order before invoking the build function:

- **File not found.** The path is checked with `std::filesystem::exists` before
  calling `dlopen`. This gives a clearer error than the linker error that
  `dlopen` would produce for a missing file.

- **`dlopen` failure.** If the operating system cannot load the shared library
  (e.g. missing transitive dependencies, wrong architecture), the error string
  from `dlerror` is included in the exception message.

- **Missing entry point.** The loader looks up the symbol
  `acts_geometry_module_v1` via `dlsym`. If it is absent the module was either
  not built with `ACTS_DEFINE_GEOMETRY_MODULE` / `ACTS_DEFINE_DD4HEP_GEOMETRY_MODULE`,
  or the symbol was stripped or hidden by the linker.

- **Null or incomplete descriptor.** The struct returned by the entry point must
  be non-null and have `module_abi_tag` set and both function-pointer fields (`build`, `destroy`)
  non-null. A null or partially initialized descriptor indicates a
  programming error in the module.

- **ABI tag mismatch.** `module_abi_tag` is compared against the tag baked into
  the host loader library at its own build time (`ActsCore` for plain modules,
  `Acts::PluginDD4hep` for DD4hep modules). A mismatch means the module and the
  host were built from different ACTS or DD4hep versions and cannot safely
  interoperate. Rebuild the module against the same ACTS installation and, for
  DD4hep modules, the same DD4hep version.

- **User-data type mismatch.** The `user_data_type` field in the descriptor is
  compared against the type the loader expects. If you call `loadGeometryModule`
  on a module that declares a `user_data_type` (i.e. it needs extra context such
  as a `dd4hep::Detector`), or call a typed loader (e.g. `loadDD4hepGeometryModule`)
  on a plain module, a `std::runtime_error` is thrown with a hint pointing to
  the correct loader to use.

- **`build` returns null.** The module's build function returned a null pointer,
  which means it failed to construct the geometry. Exceptions thrown inside
  `build` are caught by the helper and logged via `ACTS_ERROR` before the null
  is returned, so check the log for the underlying cause.
