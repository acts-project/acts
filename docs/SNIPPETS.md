# Documentation Snippet Pattern

## Overview

ACTS documentation uses **inline code snippets** from working example files rather than duplicating code in markdown files. This ensures documentation examples remain up-to-date and compilable.

## How It Works

### 1. Define Snippets in Source Files

Wrap code sections with special comment markers:

**C++ Example** (`docs/examples/logging.cpp`):

```cpp
//! [Member Logger Pattern]
struct MyClass {
  std::unique_ptr<const Acts::Logger> m_logger;
  const Acts::Logger &logger() const { return *m_logger; }

  MyClass(std::unique_ptr<const Acts::Logger> logger)
      : m_logger(std::move(logger)) {}
};
//! [Member Logger Pattern]
```

**Python Example** (`docs/examples/test_odd.py`):

```python
def test_basic():
    #! [Basic ODD construction]
    import acts.root
    from acts.examples.odd import getOpenDataDetectorDirectory

    odd_dir = getOpenDataDetectorDirectory()
    detector = getOpenDataDetector(materialDecorator=materialDecorator)
    #! [Basic ODD construction]
```

### 2. Reference Snippets in Documentation

Use Doxygen's `@snippet` command in markdown files:

**In `docs/groups/logging.md`:**

```markdown
### Member Logger Pattern

Use this pattern when a class needs persistent logging throughout its lifetime.

@snippet{trimleft} examples/logging.cpp Member Logger Pattern
```

**Options:**

- `@snippet{trimleft}` - Removes leading whitespace (recommended)
- `@snippet` - Preserves original indentation

### 3. Build the Example Code

Ensure snippet files are compiled to catch errors:

**In `docs/CMakeLists.txt`:**

```cmake
# Build examples to validate snippets
add_library(docs-examples SHARED)
target_sources(docs-examples PRIVATE examples/logging.cpp)
target_link_libraries(docs-examples PRIVATE Acts::Core)
```

For Python snippets, use pytest or similar to validate syntax. This is already set up.

## Benefits

- **Always Current**: Documentation reflects actual working code
- **Compile-Time Verification**: Broken examples cause build failures
- **Single Source of Truth**: No duplication between docs and examples
- **Easy Maintenance**: Update code once, documentation updates automatically

## Best Practices

1. **Keep snippets focused**: Extract only the relevant portion of code
2. **Use descriptive names**: Snippet names should clearly indicate their purpose
3. **Test the examples**: Build/run example files to ensure they work
4. **Minimize context**: Snippets should be understandable without reading the entire file
5. **Document around snippets**: Add explanatory text before/after the `@snippet` command

## File Naming Conventions

- Example source files: `docs/examples/*.cpp`, `docs/examples/*.py`
- Documentation files: `docs/groups/*.md`
- Snippet names: Use spaces, describe the pattern/concept (e.g., "Member Logger Pattern")

## Referencing Snippets in Header Files

When documenting classes or components in header files (`.hpp`), you can also use snippets to provide executable examples. This is especially useful for replacing inline code blocks in namespace or class documentation.

**Example in a header file** (`Core/include/Acts/Definitions/Units.hpp`):

```cpp
/// @namespace Acts::UnitConstants
/// @brief Constants and helper literals for physical units.
///
/// Here is how to use unit constants:
///
/// @snippet{trimleft} examples/units.cpp Using Unit Constants
///
/// Or use user-defined literals for more concise code:
///
/// @snippet{trimleft} examples/units.cpp Using Unit Literals
namespace UnitConstants {
  // ... constants ...
}
```

**Benefits for header file documentation:**
- Examples are verified at compile time
- Reduces clutter in header files
- Makes documentation easier to maintain
- Ensures examples use actual working code
