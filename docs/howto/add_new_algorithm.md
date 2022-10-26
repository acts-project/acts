# Add a new algorithm

## Purpose

The main part of ACTS code is located in the `Core` packages.
To use this code in the ACTS examples framework, an algorithm is needed.
Before doing so, the ideas explored in Examples are typically first developed as algorithms.
In a second step, the essential parts of this code are then moved to the `Core` packages,
which the algorithm then executes.



Code Organisation
------------------


Algorithms reside in `Examples/Algorithms` of the repository.
The code is split into a header with the algorithm class declaration
and a source file containing the implementation.

Assuming that you want to experiment with a new seeding algorithm the files to add would be:
`Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/MySeedingAlgorithm.h`
and
`Examples/Algorithms/TrackFinding/src/MySeedingAlgorithm.cpp`

The `CMakeLists.txt` file needs to be updated as well.

## Algorithm Class

The algorithm class has to inherit from {class}`ActsExamples::BareAlgorithm`
and thus implement this single method:

```cpp
    ProcessCode execute(const AlgorithmContext& ctx) const final override;
```

:::{hint}
There are other types of algorithms for reading inputs from files
or for saving data to files. All these interfaces can be found in
`Examples/Framework/include/ActsExamples/Framework`.
In particular, there are base classes for algorithms for IO operations.
:::

:::{important}
The constructor should ideally follow certain rules. See section on configuration below.
:::

The `execute` method will be called once for each event.
Other methods can also be added to your class.
It is good to remember that the algorithmic code in ACTS is typically stateless,
and if the state needs to be passed between the calls it should be done explicitly.
An algorithm is stateless if does not modify its own attributes while executing.
This way it becomes reentrant and in consequence thread safe.
For instance if the event processing is best organized with such methods:

```cpp

    void prepareBuffers();
    void fillBuffers( const SpacePoint& );
    void extractInfo();
```

that need to pass the data between each other these methods should rather look like this:

```cpp

    struct SeedingBuffers{ ... some buffers ... };
    void prepareBuffers(SeedingBuffers& buffers);
    void fillBuffers(SeedingBuffers& buffers, const SpacePoint& );
    void extractInfo(const SeedingBuffers& buffers);
```

:::{tip}
It is common practice to put the algorithm code in the `ActsExamples` namespace.
:::

## Input and Output

The algorithm would be typically part of some processing sequence
and thus consume and produce some event data.
In the hypothetical example discussed here,
space-points are the input and track seeds an output.
The data can be retrieved in the algorithm using the `get` method of the "store" object of such signature:


```cpp
    template<typename T>
    const T& get(const std::string& name);
    // example
    ctx.eventStore.get<MyType>("A_name");
```
The data is fetched from the "store" that is populated by a preceding algorithm or by the reader.

The data object (or objects) produced by an algorithm can be placed in the store using "store" `add` method:

```cpp

    template<typename T>
    void add(const std::string& name, T&& object);
    // example
    ctx.eventStore.add("A_name", std::move(mydata));
```
The ownership of the object is transferred to the store.
That is, the destruction of this object at the end of event processing is taken care of.

## Configurability

It is customary that an algorithm requires configuration parameters.
For example, in a seeding algorithm these parameters could include which detector layers should be used.
The configuration can be provided to an algorithm through an additional class/structure.
For algorithm it is typically an inner class named `Config`.

That is how the configuration object could look like for `MySeedingAlgorithm`:

```cpp
class MySeedingAlgoritm : ...
    public:
    struct Config {
        std::vector<int> layers; // layers to use by the seeder
        float deltaZ;  // the maximum allowed deviation in r-z plane
    };
    ...
```
:::{tip}
It is customary to put the config structures in the ``Acts`` namespace.
:::

The algorithm constructor would take a `MySeedingAlgorithm::Config` object during
the construction in addition to an argument controlling verbosity of diagnostic messages.

```cpp

    MySeedingAlgorithm::MySeedingAlgorithm( Config cfg, Acts::Logging::Level lvl):
      ActsExamples::BareAlgortihm("MySeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)){...}
```
## Python bindings

In order to use an algorithm in standalone ACTS the algorithm
and the associated config structure need to be accessible from python.
For that, python bindings need to be created using the pybind11 library.
The binding is defined in C++ code in `Examples/Python/src/` directory.
There is one source file per category,
in this particular case the file to edit would be `TrackFinding.cpp`.

The algorithm class and associated config class can be made known to python via such binding definition:

```cpp
    ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::MySeedingAlgorithm,
      "MySeedingAlgorithm",
      layers, deltaZ);
```

The bindings can be tested in a standalone python session:

```python

    from acts .examples import *
    help(MySeedingAlgorithm)
    help(MySeedingConfig)
```
An info about the class and config structure should be printed.

Example empty algorithm
-----------------------
A complete example of an algorithm called `UserAlgorithm` can be found in these two branches/locations:

[Algorithm definition](https://github.com/asalzburger/acts/tree/ws-add-user-algorithm/Examples/Algorithms/Tutorial)

[Python bindings definition](https://github.com/asalzburger/acts/blob/ws-add-user-algorithm-python-bindings/Examples/Python/src/Tutorial.cpp)




