import matplotlib

matplotlib.use("Agg")


def test_plotting(tmp_path):
    #! [Plotting an ACTS histogram]
    import matplotlib.pyplot as plt
    import acts

    # `hist` is a `Histogram1`/`ProfileHistogram1`/`Efficiency1` as returned by
    # `.histograms()["some_key"]` on a performance writer (see
    # Examples/Scripts/Python/pypi_finding_fitting_demo.py for how to build
    # one). `acts._demo_histogram1()` stands in for it here.
    hist = acts._demo_histogram1()

    fig, ax = plt.subplots()
    hist.plot(ax=ax)
    fig.savefig(tmp_path / "demo_histogram.png")
    #! [Plotting an ACTS histogram]

    plt.close(fig)


def test_boost_histogram_conversion():
    #! [Converting to boost-histogram]
    import boost_histogram as bh
    import acts

    hist = acts._demo_histogram1()
    bh_hist = bh.Histogram(
        hist
    )  # any downstream tool that speaks the PlottableHistogram protocol
    #! [Converting to boost-histogram]

    assert bh_hist.ndim == 1
