/* Linux version script for libActsPluginArrow (allowlist).
 *
 * Symbols listed under "global:" are exported from the .so; all other
 * symbols are hidden, regardless of origin (arrow's static archive,
 * template instantiations, system libraries, etc.).
 *
 * Uses an allowlist approach: only symbols matching the patterns below
 * are exported. All others are hidden by default.
 *
 * Pattern syntax: fnmatch-style globs (*, ?, [abc]) against mangled
 * symbol names. On Linux, symbols have a single leading underscore
 * (e.g., _ZN4Acts...).
 */
{
    global:
        /* Anything whose mangled name contains "Acts" — covers all three
         * ACTS namespaces (Acts::, ActsPlugins::, ActsExamples::) plus any
         * template instantiation or vtable / typeinfo / thunk parameterized
         * on ACTS types. C++ mangling is case-sensitive, so this does not
         * match lowercase "acts" or "ACTS". */
        *Acts*;
    local:
        *;
};
