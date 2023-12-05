# cholla_api

This is a python package to help interact with outputs from [Cholla](https://github.com/cholla-hydro/cholla). 

## Install

To install onto your workstation, you simply need to pip install this package. Currently, this package is on [TestPypi](https://test.pypi.org/project/cholla-api/), with a goal to push onto Pypi. To install, run the following

```
python3 -m pip install --index-url https://test.pypi.org/simple/ --no-deps cholla_api
```

To install from the Github repo, navigate into the `~/cholla_api` directory and run the following

``
python3 -m pip install --user .
``

where we use the `--user` flag to pip install the package onto your home directory, in case you don't have root access.

## Getting Started

After importing `cholla_api` in a python script, we can create a `ChollaSnap` instance that will help in studying a snapshot. To study an entire simulation run, we can also use the `ChollaRun` class. The `ChollaViz` directory holds many helper plotting functions.

## Other Info

Currently, this package works with Cholla simulation runs that use `TYPE=hydro`.

