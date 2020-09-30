# Third party software

This contains software that is usually not available through the system
package manager.

**Note** Only include software here that is either a header-only library or
can be built as a static library. Do not use any of this software as part of
the public interface and only include it in pure implementation files to avoid
issues with missing files after installation.

`nlohmann_json` is exempted from this rule as it is handled specially.

## dfelibs

A copy of [dfelibs](https://github.com/msmk0/dfelibs) v20200416, with the
following folders removed:

-   examples
-   unittests

## nlohmann_json

A copy of [nlohmann::json](https://github.com/nlohmann/json) revision 84f19d6, with the
following folders removed:

-   benchmarks
-   doc
-   test
-   include
-   third_party

## autodiff

A copy of [autodiff](https://github.com/autodiff/autodiff), v0.5.11.
In the `CMakeLists.txt` the commands `add_subdirectory(tests)` and 
`add_subdirectory(examples)` have been commented out. All folders/files have been 
removed except the following ones:

-   autodiff (contains the header files)
-   cmake
-   CMakeLists.txt
