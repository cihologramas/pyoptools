from pyoptools.misc.resources import detectCPUs, detectOpenCL


def test_detect_cpus():
    cpus = detectCPUs()
    assert isinstance(cpus, int)
    assert cpus >= 1


def test_detect_opencl():
    has_cl = detectOpenCL()
    assert isinstance(has_cl, bool)
