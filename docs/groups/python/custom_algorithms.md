@defgroup python_custom_algorithms Custom Algorithms
@ingroup python_bindings
@brief Writing your own `IAlgorithm` or `IReader` in Python, and moving data through the whiteboard.

Every step in a `Sequencer` — reader, algorithm, or writer — implements a small C++ interface.
Python subclasses of `acts.examples.IAlgorithm` and `acts.examples.IReader` are accepted directly
by `Sequencer.addAlgorithm()`/`addReader()`, alongside the built-in bindings, so you can write a
custom step without touching C++.

# A minimal algorithm

@snippet{trimleft} examples/test_custom_algorithm.py Custom algorithm filtering particles by pT

Points to note:

- The constructor must call the base class `__init__` (`acts.examples.IAlgorithm.__init__(self,
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
shows custom Python track-finding and fitting algorithms in a complete chain.

# Python algorithms and the GIL

`Sequencer.run()` releases the Python GIL while C++ algorithms execute. Calls back into Python
(`execute`, `read`, `initialize`, `finalize`, and `name`) acquire it again, so Python steps do not
run Python code in parallel across events. If your algorithm shares state that is not safe to
access concurrently, such as a model or library handle, use `numThreads=1`.
