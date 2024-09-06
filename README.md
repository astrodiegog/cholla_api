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

There are three foundational data type implementations that help manage Cholla outputs: ChollaGrid (information about the domain), ChollaMacroFlags (information on the data simulated), and ChollaBoxHead (information connecting a process domain to the global domain). Under the examples directory, there is a jupyter notebook that will help in getting familiar with typical data structures in this program. Under the ``src/cholla_api`` subdirectory, there are a couple markdown files that describe: the inputs to a Cholla program, the outputs of a Cholla simulation, and our approach in implementing this API.

Under the `data` subdirectory will be information regarding how our API interacts with snapshot files. Under the `snap` subdirectory will be structures describing snapshot-wide information. Under `run` subdirectory will be simulation-wide information regarding the grid domain and the compiler flags used in running a simulation. Under the `OTFanalysis` subdirectory will be data structures that help interact with On-The-Fly analysis outputs from Cholla. Under the `analysis` subdirectory will be how to manipulate data from snapshots in helpful ways. Lastly under the `skewers` subdirectory will be our implementation of what a skewer is within the context of cosmological simulated spectra.

## Distribution

The `dist` files are created when building using PyPA [build](https://github.com/pypa/build).

`tar.gz` serves as the source distribution while `.whl` is a built distribution. 

The source distribution holds metadata and source files to be used by pip when installing, and serves as a hard archive of the package.

The built distribution includes the compiled code, and is used for much quicker installation.

## Features

This python package was developed with cosmological simulation analysis tools. In effect, most features are geared towards that usage

- Compute optical depth along a skewer

## Other Info

Future features in development

- python scripts testing domain decomposition and skewer decomposition
- Compute temperature history 

