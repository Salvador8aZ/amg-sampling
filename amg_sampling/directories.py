from pathlib import Path


def qualifyname(directoryname, filename=None):
    directory = Path(directoryname)
    if filename is None:
        return directory
    return directory / filename


def code(filename=None):
    return qualifyname(Path(__file__).parent.resolve(), filename)


def base(filename=None):
    return qualifyname(code().parent, filename)


def data(filename=None):
    return qualifyname(base("data"), filename)


def package_data(filename=None):
    return qualifyname(code("data"), filename)


def tests(filename=None):
    return qualifyname(code("tests"), filename)


def test_data(filename=None):
    return qualifyname(tests("test_data"), filename)
