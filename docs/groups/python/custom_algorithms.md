@defgroup python_custom_algorithms Custom Algorithms
@ingroup python_bindings
@brief Writing your own `IAlgorithm` or `IReader` in Python, and moving data through the whiteboard.

Most readers, algorithms, and writers that come with ACTS are writtend in C++. However, it is also possible to add Python subclasses of `acts.examples.IAlgorithm` and `acts.examples.IReader` directly to the `Sequencer`. This allows to implement custom functionality without touching C++.

# A minimal algorithm

@snippet{trimleft} examples/test_custom_algorithm.py Custom algorithm filtering particles by pT

Points to note:

- The constructor must call the base class constructor (`acts.examples.IAlgorithm.__init__(self,
  name, level)`) before doing anything else.
- `execute(self, context)` is the only method you must override; it returns an
  `acts.examples.ProcessCode` (`SUCCESS`, `SKIP`, `ABORT`, or `END`).
- `self.logger` is available once the base constructor has run.

# Reading and writing the whiteboard

Data is passed between steps via an in-memory `WhiteBoard`, keyed by string. The example above
initializes typed `ReadDataHandle` and `WriteDataHandle` objects with those keys. Inside
`execute`, read with `handle(context.eventStore)` and write with `handle(context, value)`.
`context.eventStore` also gives you `.exists(key)` and `.keys` directly, if you need to probe the
whiteboard rather than go through a typed handle.

# A custom reader

`IReader` follows the same pattern, but implements `availableEvents()` (returning
`(firstEvent, lastEventExclusive)`) instead of taking whiteboard input. `acts.examples.uproot`
ships two ready-to-use examples that read ROOT files without the ROOT plugin — read them for a
complete, working `IReader` implementation, including buffered reads across events.

The [PyPI finding and fitting demo](https://github.com/acts-project/acts/blob/main/Examples/Scripts/Python/pypi_finding_fitting_demo.py)
shows dummy Python track-finding and fitting algorithms in a complete chain.

# Python algorithms and the GIL

Regular CPython uses a Global Interpreter Lock (GIL) that ensures only on thread is executing python code at the same time. `Sequencer.run()` releases the GIL while C++ algorithms execute.

Calls back into Python such as `IAlgorithm.execute` acquire it again. This is not a problem usually, but means, that a custom Python algorithm will run single-threaded and can become a bottleneck, regardless of the `numThreads` configuration of the `Sequencer`.

> [!warning]
> If you define callbacks into C++ in plain python, and C++ code can mutate their states, this can lead to correctness issues.
