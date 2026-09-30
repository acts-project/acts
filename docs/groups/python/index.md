@defgroup python_bindings Python Bindings
@brief Python bindings for the ACTS tracking toolkit, available via PyPI or built from source.

The bindings expose the core library, Fatras, the Examples framework, and (when built from
source) the optional plugins, as a single `acts` Python package.

- @ref python_installation "Installation" — PyPI vs. building from source, what each gives you.
- @ref python_getting_started "Getting started" — the `Sequencer` and the high-level `add*` helpers.
- @ref python_custom_algorithms "Custom algorithms" — writing your own `IAlgorithm`/`IReader` in Python.
- @ref python_data_io "Data reading and writing" — the reader/writer inventory, including the
  ROOT-free path used by the PyPI wheel.
- @ref python_performance_plotting "Performance evaluation and plotting" — extracting efficiency
  and resolution numbers, and plotting them.
- @ref python_caveats "Best practices and caveats" — by-design behavior that is easy to get wrong.
