Add a new algorithm
===================

Purpose
-------------

The main part of ACTS code is located in `Core` packages. 
For them to be used in standalone ACTS an algorithm is needed.
Also, new ideas need to be first developed as algorithms.
Typically the essential parts of the code are then moved to the `Core` packages 
and the algorithm only executes it.



Code Organisation
------------------

Algorithms reside in `Examples/Algorithms` are of the repository.
The code is split into header with the algorithm class declaration 
and a source file containing the implementation.

Assuming that you want to experiment with a new seeding algorithm the files to add would be:
`Examples/Algorithms/TrackFinding/include/ActsExamples/TrackFinding/MySeedingAlgorithm.h`
and 
`Examples/Algorithms/TrackFinding/src/MySeedingAlgorithm.h`

Algorithm Class
---------------
The algorithm class has to inherit from {class}`ActsExamples::BareAlgorithm`
and thus implement single method:

.. code-block:: cpp

    ProcessCode execute(const AlgorithmContext& ctx) const final override;

.. hint:: There are other types of algorithms for reading inputs from files
    or for saving data to files. All these interfaces can be found in 
    `Examples/Framework/include/ActsExamples/Framework`.

.. important:: The constructor should ideally follow certain rules. See section on configuration below.

The `execute` method will be called once for each event. 
Obviously other methods can be also added to this class.
It is good to remember that the algorithmic code in ACTS is typically stateless
and if the state needs to be passed between the calls it should be done explicitly.
For instance if the event processing is best organized with a methods:

.. code-block:: cpp

    void prepareBuffers();
    void fillBuffers( const SpacePoint& );
    void extractInfo();

that need to pass the data between each other these methods should rather look like this:

.. code-block:: cpp

    struct SeedingBuffers{ ... some buffers ... };
    void prepareBuffers(SeedingBuffers& buffers);
    void fillBuffers(SeedingBuffers& buffers, const SpacePoint& );
    void extractInfo(const SeedingBuffers& buffers);


..  tip:: It is a common practice to put the algorithm code in ``ActsExamples`` namespace.

Input and Output
------------------

The algorithm would be typically part of some processing sequence and thus consume and produce some event data.
In hypothetical example discussed here space-points cod be the input and track seeds would be an output.
The data can be retrieved in the algorithm using `get` method.

.. code-block:: cpp

    template<typename T>
    const T& get(const std::string& name);

The data is fetched from so called "store" that is populated by preceding algorithm or by the reader.

The data object (or objects) produced by an algorithm can be placed in the store using method `add`.

.. code-block:: cpp

    template<typename T>
    void add(const std::string& name, T&& object);

As the method signature suggests the ownership of the object is transferred to the store.
That is, the destruction of this object at the end of event processing is taken care of.

Configurability
----------------

It is customary that an algorithm requires configuration parameters.
For the sake of an example in hypothetical seeding algorithm it could be parameter steering
which detector layers should be used. 
The configuration can be provided to an algorithm through an additional class/structure.
It can be an inner class of the algorithm or it can be external to it.
One would use definitely external structure if it is shared among several algorithm classes.

For example that is how the configuration object could look like for `MySeedingAlgorithm`:

.. code-block:: cpp

    struct MySeedingConfig {
        std::vector<int> layers; // layers to use by the seeder
        float deltaZ;  // the maximum allowed deviation in r-z plane
    };

.. tip:: It is customary to put the config structures in ``Acts`` namespace.

The algorithm constructor would take `MySeedingConfig` object during the construction in addition 
to the argument controlling verbosity of diagnostic messages.

.. code-block:: cpp

    MySeedingAlgorithm::MySeedingAlgorithm( Acts::MySeedingConfig cfg, Acts::Logging::Level lvl):
      ActsExamples::BareAlgortihm("MySeedingAlgorithm", lvl), 
      m_cfg(std::move(cfg)){...}


Python bindings
---------------
In order to use an algorithm in standalone ACTS 
the algorithm and the associated config structure need to be accessible from python.
For that so called python bindings need to be created using ``pybind11`` library.
The binding is defined in C++ code in `Examples/Python/src/` directory. 
There is one source file per category, in this particular case the file to edit would be `TrackFinding.cpp`.

The configuration structure binding would be defined like this:

.. code-block:: cpp

    using Config = Acts::MySeedingConfig;
    auto c = py::class_<Config>(m, "MySeedingConfig").def(py::init<>()); // this defined the name know in python
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(layers); // makes the layers field accessible in python
    ACTS_PYTHON_MEMBER(deltaZ); // makes the deltaZ accessible
    ACTS_PYTHON_STRUCT_END();
    patchKwargsConstructor(c);

The algorithm class can be made known to python via such binding definition:

.. code-block:: cpp

    auto alg =
        py::class_<ActsExamples::MySeedingAlgorithm, 
                   ActsExamples::BareAlgorithm,
                   std::shared_ptr<ActsExamples::MySeedingAlgorithm>>(
            mex, "MeSeedingAlgorithm")
            .def(py::init<const Acts::MySeedingConfg&, Acts::Logging::Level>(), // makes the constructor callable from python
                 py::arg("config"), py::arg("level")); // defines constructor arguments
        // other methods can be exposed to python (typically config accessor) 

If bindings are defined correctly (and everything compiles) they can be tested in standalone python session (see section on setting up python) by typing:

.. code-block:: python

    from acts .examples import *
    help(MySeedingAlgorithm)
    help(MySeedingConfig)

An info about the class and config structure should be printed.

Example empty algorithm
-----------------------
A complete example of an algorithm called ``UserAlgorithm`` can be found in this two branches/locations:

[Algorithm definition](https://github.com/asalzburger/acts/tree/ws-add-user-algorithm/Examples/Algorithms/Tutorial)

[Python bindings definition](https://github.com/asalzburger/acts/blob/ws-add-user-algorithm-python-bindings/Examples/Python/src/Tutorial.cpp)




