"""viz package skeleton + matplotlib hard-dependency smoke test (spec §5 c.10)."""


def test_matplotlib_is_a_hard_dependency():
    import matplotlib  # must import without any extra
    assert hasattr(matplotlib, "__version__")
    major, minor = (int(x) for x in matplotlib.__version__.split(".")[:2])
    assert (major, minor) >= (3, 7)   # the pinned floor (pyproject: matplotlib>=3.7)


def test_viz_package_imports():
    import quiverlab.viz as viz
    assert hasattr(viz, "__all__")   # the package exposes its public surface
