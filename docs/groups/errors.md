@defgroup errors Error definitions

# Error Handling in ACTS

ACTS uses `std::error_code` from the C++ Standard Library for type-safe,
extensible error handling. This approach allows error codes to be returned from
functions without throwing exceptions, while maintaining compatibility with the
standard error handling mechanisms. Error codes are typically used with the
@ref Acts::Result type for ergonomic error propagation.

## std::error_code Overview

`std::error_code` is a platform-dependent error code that consists of:

- An integer error value
- A reference to an error category (error domain)

Various ACTS components define their own strongly-typed error enum class and
register it with the standard library's error code system using
`std::is_error_code_enum` specialization. This allows seamless conversion to
`std::error_code` while preserving type safety.

## Error Enum Pattern

All ACTS error enums follow a consistent pattern:

@snippet{trimleft} examples/errors.cpp Error Enum Pattern

The enums are registered with STL by specializing `std::is_error_code_enum<T>`
to enable implicit conversion to `std::error_code`.

## Usage with Result Type

Error codes are typically used with the @ref Acts::Result type, which
encapsulates either a successful result or an error:

@snippet{trimleft} examples/errors.cpp Usage with Result Type

The @ref Acts::Result type defaults to using `std::error_code` as its error
type, making it compatible with all ACTS error enums.

## Benefits

- **Type Safety**: Each component has its own error enum preventing confusion
  between different error domains
- **No Exceptions**: Errors are returned as values, enabling explicit error
  handling without exception overhead
- **Standard Compatibility**: Integration with `std::error_code` allows use
  with standard library error handling facilities
- **Extensibility**: New error types can be added to any enum without breaking
  binary compatibility
