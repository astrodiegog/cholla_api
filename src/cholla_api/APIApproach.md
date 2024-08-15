# API Approach

This file is meant to describe the approach and thought process behind the design choices when creating this API for Cholla simulations.

This API was developed with the idea of analyzing cosmological simulations, which means that it is simply not feasible to concatenate datasets from all processes to one large file in order to analyze a cosmological simulation. Cosmological simulations do not necessarily _need_ to be ran on multiple processes, but the API was designed with the idea that a Cholla simulation outputs data from multiple processes, specifically GPU nodes.

We designed this API to be as flexible to whatever simulation was run, relying primarily on information that is specified in the parameter file.


The most basic data structure defined in this API is the ``ChollaGrid`` class:

```bash
class ChollaGrid:
    '''
    Cholla Grid object
        Holds global information about the domain and bounds of a simulation
            run. To complete domain decomposition, also need the number of
            processes simulation was ran with

        Initialized with:
        - nprocs (int): number of MPI processes
        - nx (int): number of cells in x dimension
        - ny (int): number of cells in y dimension
        - nz (int): number of cells in z dimension
        - xmin (float): minimum x coordinate
        - ymin (float): minimum y coordinate
        - zmin (float): minimum z coordinate
        - xmax (float): maximum x coordinate
        - ymax (float): maximum y coordinate
        - zmax (float): maximum z coordinate
    
    '''
```

Given the number of cells in each dimension, as well as the bounding coordinates, we already have enough information to tile processes and perform the domain decomposition given some number of processes because we know that memory is stored in row-major order.

With a ``ChollaGrid`` object defined, we can define the foundational data structure in this API, the ``ChollaBoxHead`` class


```bash
class ChollaBoxHead:
    '''
    Cholla Box Head object
        Holds information regarding the location of the data

        Initialized with:
        - nBox (int): number of the box within snapshot
        - cell_offset_x (int): number of cells offset in x dimension 
        - cell_offset_y (int): number of cells offset in y dimension
        - cell_offset_z (int): number of cells offset in z dimension
        - local_dimx (int): number of cells in the box subvolume in x dimension
        - local_dimy (int): number of cells in the box subvolume in y dimension
        - local_dimz (int): number of cells in the box subvolume in z dimension
    '''
```

This object holds the information that connects a process to the global volume.

That's really cool, but it's all abstract. What if I just want to grab the data held within some file? Well to do that, we tie this header information along with the path to a snapshot to create the ``ChollaBox`` class

```bash
class ChollaBox:
    '''
    Cholla Box object
        Holds information regarding the box for a specific snapshot 

        Initialized with:
            SnapPath (str): path to a snapshot
            ChollaBoxHead (ChollaBoxHead): Cholla Box Head object connecting
                box to rest of volume
    '''
```

We know how each file name is structured as well as the key string used to save the data as Dataset objects, so we define three important methods:

- ``get_XfPath()`` : returns the file path to the hydro data for this box
- ``check_Xkey(key)`` : returns whether a user-provided key is valid to use
- ``get_Xdata(key)`` : returns the key data set, with an optional parameter for the precision to use

where ``X`` is a stand-in for either hydro, gravity, or particle information.

There is some information that is shared between all boxes within a snapshot, so we also have a ``ChollaSnapHead`` class

```bash
class ChollaSnapHead:
    '''
    Cholla Snapshot Head object
        Holds snapshot specific information
        Initialized with:
        - nSnap (int): number of the snapshot within run
    '''   
```

as well as a ``ChollaSnap`` class 

```bash
class ChollaSnap:
    '''
    Cholla Snapshot object
        Holds snapshot specific information and methods to access data for that
            snapshot
            
        Initialized with:
            RunPath (str): path to a directory holding all snapshot directories
            ChollaSnapHead (ChollaSnapHead): Cholla Snap Head object holding
                snapshot specific information
    '''
```

To set some snapshot-wide information to the ``ChollaSnapHead``, we use a specific ``ChollaBox`` to grab that value. For example, the ``set_cosmoinfo(ChollaBox)`` method uses the ``ChollaBox`` to attach the redshift and scale factor to the ``ChollaSnapHead``.

The standard naming scheme in the API is that a ``ChollaXHead`` is an abstraction of what an ``X`` object actually is, down to its barebone essentials. To grab an actual implementation of that ``X`` object, we tie it to some file path with an instance of the abstract object such that ``ChollaX`` is initialized with ``ChollaXHead`` and ``XPath``.

To set global information that holds within all snapshots, we also have a ``ChollaGlobal`` class

```bash
class ChollaGlobal:
    '''
    Cholla Simulation Run object
        Holds information related to an entire simulation run. Initialized with
        - basePath (str): the directory the simulation is held in
        - chGrid (ChollaGrid): ChollaGrid object, holding domain information
        - gravity_flag (bool): whether gravity data was saved
        - particles_flag (bool): whether particle data was saved
        - cosmo_flag (bool): (optional) whether cosmology type was used
        - data_subdir (str): (optional) where data is placed within basePath
    '''
```

This holds enough information to return a ``ChollaSnap`` object given a snapshot number. Given a snapshot number and box id, we have a method to return the appropriate ``ChollaBox`` instance.


Now that is very interesting and cool, but you actually want to _do something_ with the data, not just grab the data. That is what the ``analysis/`` directory is meant to help with. To help in calculating commonly desired values from the raw data outputs, we have developed different calculators.

Firstly, we have created a ``ChollaHydroCalculator`` to deal with common hydro calculations


```bash
class ChollaHydroCalculator:
    '''
    Cholla Calculator object
        Serves as a calculator where the calculated values have some expected
            size and datatype (default is float). Assert that inputs are of same
            shape as dims that was used to initialize this calculator. To 
            complete some analysis, this ChollaCalculator will be the mediator 
            that will act on the primitive saved values. 

        Initialized with:
            dims (tuple): size of data sets to act on
            dtype (np type): (optional) numpy precision to initialize output arrays


        Contains methods to calculate
            hydro velx, vely, velz, velmag
            pressure
            specific internal energy
            temperature
    
    Values are returned in code units unless otherwise specified.
    '''
```

For a higher-order level of abstraction, we can tie a ``ChollaHydroCalculator`` to a ``ChollaBox`` to define the calculation of hydro values specific to a dataset


```bash
class ChollaBoxHydroCalc:
    '''
    Cholla Box Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - ChollaBox (ChollaBox): ChollaBox object
        - dtype (datatype): precision to compute data
    '''
``` 


Likewise, we can define tie a ``ChollaHydroCalculator`` to a ``ChollaSnap`` to create a class that returns calculated values for an entire snapshot

```bash
class ChollaSnapHydroCalc:
    '''
    Cholla Snap Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - ChollaSnap (ChollaSnap): ChollaSnap object
        - ChollaGrid (ChollaGrid): ChollaGrid object
        - dtype (datatype): (optional) precision to calculate data
        - calc_box (bool): (optional) whether to use ChollaBoxHydroCalculator object

        if calc_box is set to True, then we use ChollaBox's place_data method
            to create a global data array
        else, then we instead complete calculation on the global primitive
            data. take care when setting calc_box=False because it can lead
            to loading large numpy arrays into memory !
    '''
``` 

To save some helpful units at cosmological scales, as well as help in converting between proper and comoving coordinates, we have created a cosmological calculator

```bash
class ChollaCosmoGridCalculator:
    '''
    Cholla Cosmology Grid Calculator object
        Serves as a calculator for cosmology-related values given a discretized
            ChollaGrid object. 
        
        Initialized with:
        - ChollaGrid (ChollaGrid): ChollaGrid object
        - OmegaM (float): present-day energy density parameter for matter
        - OmegaR (float): present-day energy density parameter for radiation
        - OmegaK (float): present-day energy density parameter for spatial curvature
        - OmegaL (float): present-day energy density parameter for dark energy
        - w0 (float): constant term in dark energy equation of state
        - wa (float): linear term in dark energy equation of state
        - H0 (float): present-day Hubble parameter in units of [km / s / Mpc]

    Values are returned in code units unless otherwise specified.
    '''
```


This document has described some low-level API to interact with snapshot files. However, Cholla can also output 2 other data files: analysis and skewers !

To describe the barebones of a skewer, we have a ``ChollaOnTheFlySkewerHead`` class


```bash
class ChollaOnTheFlySkewerHead:
    '''
    Cholla On The Fly Skewer Head

    Holds information regarding a specific individual skewer

        Initialized with:
        - skew_id (int): id of the skewer
        - n_i (int): length of the skewer
        - skew_key (str): string to access skewer

    '''
```

where ``skew_key`` refers to which dimension the skewer lies along. An actual implementation of a skewer, would require the path to the file

```bash
class ChollaOnTheFlySkewer:
    '''
    Cholla On The Fly Skewer
    
    Holds skewer specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFSkewerHead (ChollaOnTheFlySkewerHead): header
            information associated with skewer
        - fPath (str): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''
``` 

This object has methods to grab data along that one skewer.

Going one level up, we also have objects to describe the information regarding all of the skewers along one dimension

```bash
class ChollaOnTheFlySkewers_iHead:
    '''
    Cholla On The Fly Skewers_i Head

    Holds information regarding a specific skewer hdf5 group

        Initialized with:
        - n_i (int): length of the skewers
        - n_j (int): length of first dimension spanning cube
        - n_k (int): lenth of second dimension spanning cube
        - n_stride (int): stride cell number between skewers
        - skew_key (str): string to access skewer

    '''
``` 

Tieing this to an actual implementation, we have the ``ChollaOntheFlySkewers_i`` class

```bash
class ChollaOnTheFlySkewers_i:
    '''
    Cholla On The Fly Skewers
    
    Holds skewer specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFSkewersiHead (ChollaOnTheFlySkewers_iHead): header
            information associated with skewer
        - fPath (str): file path to skewers output

    Values are returned in code units unless otherwise specified.
    '''
```

This object has the ability to return a ``ChollaOnTheFlySkewer`` object, as well as all of the skewer data along one dimension.

Lastly, we also have an object that holds all of the information regarding a skewer output

```bash
class ChollaOnTheFlySkewers:
    '''
    Cholla On The Fly Skewers
    
    Holds on-the-fly skewers specific information to an output with methods to 
            access data for that output

        Initialized with:
        - nSkewer (nSkewer): number of the skewer output
        - SkewersPath (str): directory path to skewer output files
        - ChollaGrid (ChollaGrid): grid holding domain information

    Values are returned in code units unless otherwise specified.
    '''
```


To describe the outputs of an analysis file, we have created 2 classes for each of the subgroups in an analysis file. Firstly, to describe the power spectrum, we have a ``ChollaOnTheFlyPowerSpectrumHead`` class initialized with all the information regarding a power spectrum.


```bash
class ChollaOnTheFlyPowerSpectrumHead:
    '''
    Cholla On The Fly Analysis Power Spectrum Head
    
    Holds information regarding the power spectrum calculation

        Initialized with:
        - dlogk (float): differential step in log k-space
        - nfft (int): number of FFT bins
        - nlos (int): number of line-of-sight cells
        - dvHubble (float): differential Hubble flow velocity

    Values are returned in code units unless otherwise specified.
    '''

``` 

An implementation of a power spectrum is that tied to a specific file is specified with a ``ChollaOnTheFlyPowerSpectrum`` class

```bash
class ChollaOnTheFlyPowerSpectrum:
    '''
    Cholla On The Fly Analysis Power Spectrum
    
    Holds power spectrum specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFPowerSpectrumHead (ChollaOnTheFlyPowerSpectrumHead): header
            information for power spectrum
        - fPath (str): file path to analysis output

    Values are returned in code units unless otherwise specified.
    '''
```


Likewise, we have the ``ChollaOnTheFlyPhaseSpaceHead`` class holding information regarding the density and temperature limits as well as the number of histogram for each axis

```bash
class ChollaOnTheFlyPhaseSpaceHead:
    '''
    Cholla On The Fly Analysis Phase Space Head

    Holds information regarding the phase space calculation

        Initialized with:
        - density_min (float): minimum density
        - density_max (float): maximum density
        - T_min (float): minimum temperature
        - T_max (float): maximum temperature
        - n_density (int): number of density bins
        - n_temp (int): number of temperature bins

    Values are returned in code units unless otherwise specified.
    '''
```


The ``ChollaOnTheFlyPhaseSpace`` ties a ``ChollaOnTheFlyPhaseSpaceHead`` with the file path to that object

```bash
class ChollaOnTheFlyPhaseSpace:
    '''
    Cholla On The Fly Analysis Phase Space
    
    Holds phase space specific information to an output with methods to 
            access data for that output

        Initialized with:
        - ChollaOTFPhaseSpaceHead (ChollaOnTheFlyPhaseSpaceHead): header
            information associated with phase space
        - fPath (str): file path to analysis output

    Values are returned in code units unless otherwise specified.
    '''
```

To help interface with an analysis file output, we have created a ``ChollaOnTheFlyAnalysis`` class 


```bash
class ChollaOnTheFlyAnalysis:
    '''
    Cholla On The Fly Analysis
    
    Holds on-the-fly analysis specific information to an output with methods to 
            access data for that output

        Initialized with:
        - nAnalysis (int): number of the analysis output
        - AnalysisPath (str): directory path to analysis output files
        - ChollaGrid (ChollaGrid): grid holding domain information

    Values are returned in code units unless otherwise specified.
    '''
```




When designing these data structures, we attempt to minimize the amount of times that a file is accessed, which is why objects have methods to grab datasets, but don't actually save data onto an attribute of the object.


