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

The last simulation-specific data structure is created to understand what compiler flags were used in making the Cholla executable

```bash
class ChollaMacroFlags:
    '''
    Cholla Macro Flags object
        Holds information related to the compiler flags used for a Cholla
            executable. 
        Initialized with:
        - DualEnergy (bool): whether Gas Internal Energy saved with Dual Energy
            Formalism
        - Gravity (bool): whether gravity was solved on GPU
        - Particles (bool): whether particles were used & tracked on GPU
        - Cosmology (bool): whether equations of moation was solved in 
            comoving frame
        - Cooling_GRACKLE (bool): whether GRACKLE was used for cooling & chemistry
        - Cooling_GPU (bool): whether cooling & chemistry solved onn GPU
        - OTF_Analysis (bool): whether on-the-fly cosmological analysis was 
            performed
    '''
```

This object holding the macro flags helps in understanding what data was actually saved.

That's really cool, but it's all abstract. What if I just want to grab the data held within some file? Well to do that, we tie this header information along with the path to a snapshot to create the ``ChollaBox`` class

```bash
class ChollaBox:
    '''
    Cholla Box object
        Holds information regarding the box for a specific snapshot 

        Initialized with:
            SnapPath (str): path to a snapshot
            chBoxHead (ChollaBoxHead): Cholla Box Head object connecting box to 
                rest of volume
            chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
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
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro compiling information
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

For a higher-order level of abstraction, we can tie a ``ChollaHydroCalculator`` to a ``ChollaBox`` and ``ChollaMacroFlags`` to define the calculation of hydro values specific to a dataset


```bash
class ChollaBoxHydroCalc:
    '''
    Cholla Box Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - chBox (ChollaBox): ChollaBox object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
        - dtype (datatype): (optional) precision to compute data
    '''
``` 


Likewise, we can define tie a ``ChollaHydroCalculator`` to a ``ChollaSnap`` to create a class that returns calculated values for an entire snapshot

```bash
class ChollaSnapHydroCalc:
    '''
    Cholla Snap Hydro Calculator object
        
        This object ties the hydro calculator object to a box's dimensions

        Initialized with:
        - chSnap (ChollaSnap): ChollaSnap object
        - chGrid (ChollaGrid): ChollaGrid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
        - dtype (datatype): (optional) precision to calculate data
        - calc_box (bool): (optional) whether to use ChollaBoxHydroCalculator object

        if calc_box is set to True, then we use ChollaBox's place_data method
            to create a global data array
        else, then we instead complete calculation on the global primitive
            data. take care when setting calc_box=False because it can lead
            to loading large numpy arrays into memory !
    '''
``` 

To save some helpful units at cosmological scales, we have created a cosmology header

```bash
class ChollaCosmologyHead:
    '''
    Cholla Cosmology Head
        Serves as a header object that holds information that helps define a
            specific cosmology
        
        Initialized with:
        - OmegaM (float): present-day energy density parameter for matter
        - OmegaR (float): present-day energy density parameter for radiation
        - OmegaK (float): present-day energy density parameter for spatial curvature
        - OmegaL (float): present-day energy density parameter for dark energy
        - w0 (float): constant term in dark energy equation of state
        - wa (float): linear term in dark energy equation of state
        - H0 (float): present-day Hubble parameter in units of [km / s / Mpc]

    '''
```

To describe cosmological information at a specific scale factor, we have created the ``ChollaSnapCosmologyHead`` header object

```bash
class ChollaSnapCosmologyHead:
    '''
    Cholla Snapshot Cosmology header object
        Serves as a header holding information that combines a ChollaCosmologyHead
            with a specific scale factor with the snapshot header object.
        
        Initialized with:
            snapHead (ChollaSnapHead): provides current redshift
            cosmoHead (ChollaCosmologyHead): provides helpful information of cosmology & units

    Values are returned in code units unless otherwise specified.
    '''
```

Given the cosmology information at a specific scale factor, we can construct a cosmology calculator

```bash
class ChollaCosmoCalculator:
    '''
    Cholla Cosmological Calculator object
        Serves as a calculator for a cosmology at a specific scale factor.
        
        Initialized with:
            snapCosmoHead (ChollaSnapCosmologyHead): provides current redshift
            dims (tuple): size of data sets to act on
            dtype (np type): (optional) numpy precision to initialize output arrays 

    Values are returned in code units unless otherwise specified.
    '''
```

This object can convert between cosmological units and cgs units, as well as from physical and proper values.

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

We have also started to create skewer information that is independent of the On-The-Fly analysis, in the case we want skewers but did not turn on that compiler flags.

To describe a skewer with respect to the entire/global simulation volume, we can instantiate that with information regarding what box the skewer lands on, as well as the local skewer id. To describe a local skewer within a face, we have the following object

```bash
class ChollaSkewerLocalFaceHead:
    '''
    Cholla Skewer Local Face Head

    Holds information regarding a local skewer within face

        Initialized with:
        - localface_id (int): id of the local skewer on face
        - localface_joffset (int): offset along j-axis
        - localface_koffset (int): offset along k-axis
    '''
```

while a face is described as

```bash
class ChollaSkewerFaceHead:
    '''
    Cholla Skewer Face Head

    Holds information regarding a local skewer face

        Initialized with:
        - face_id (int): id of the face
        - face_joffset (int): offset along j-axis
        - face_koffset (int): offset along k-axis
    '''

```

such that we can finally describe an individual skewer with the object

```bash
class ChollaSkewerGlobalHead:
    '''
    Cholla Skewer Global Head

    Holds information regarding a global skewer

        Initialized with:
        - global_id (int): id of the global skewer
        - chFaceHead (ChollaSkewerFaceHead): ChollaSkewerFaceHead object,
            holds info on face within grid
        - chFaceLocalHead (ChollaSkewerLocalFaceHead): ChollaSkewerLocalFaceHead
            object, holds info on skewer within face
        - n_los (int): number of cells along line-of-sight
        - nlos_proc (int): number of processes along line-of-sight
    '''
```

With this header object that describes a skewer in the global simulation volume, we can grab the two $(j,k)$ dimension offsets to describe where the skewer lies on the global simulation, as well as how far the skewer goes in. 

An actual implementation requires simply the path to the snapshots, the ``ChollaSkewerGlobalHead`` for a skewer as well as the simulation grid ``ChollaGrid`` that the skewer lives in. Along the x-direction, a global skewer can be implemented as

```bash
class ChollaSkewerxGlobal:
    '''
    Cholla Skewer x Global object
        Holds information regarding skewer along x direction including methods
            to grab data along skewer

        Initialized with:
        - ChollaSkewerGlobalHead (ChollaSkewerGlobalHead): global head object
        - SnapPath (str): path to snapshot directory
        - ChollaGrid (ChollaGrid): grid object
    '''
```

To describe the outputs of an analysis file, we have created 2 classes for each of the subgroups in an analysis file. Firstly, to describe the power spectrum, we have a ``ChollaOnTheFlyPowerSpectrumHead`` class initialized with all the information regarding a power spectrum.

Whether using the On-The-Fly skewers or the post-simulation skewers, we have an implementation of a skewer cosmological calculator that can compute the optical depth along a skewer. This calculator combines the scale factor from a snapshot, the cosmological framework the simulation lives in, the line-of-sight cells, and the distance between cells to fully describe the necessary ingredients to compute the optical depth.

```bash
class ChollaSkewerCosmoCalculator:
    '''
    Cholla Skewer Calculator object
        Serves as a specific implementaiton of a Cholla Cosmological Calculator
            for a skewer.

        Initialized with:
            snapHead (ChollaSnapHead): provides current redshift
            cosmoHead (ChollaCosmologyHead): provides helpful information of cosmology & units
            n_los (int): number of cells along line-of-sight
            dx (float): comoving distance between cells (kpc)
            dtype (np type): (optional) numpy precision to initialize output arrays
        
        Objects including ghost cells are suffixed with _ghost

    Values are returned in code units unless otherwise specified.
```


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

The abstract description and implementation of a singular skewer is well described, but an entire skewer analysis ties many skewers together. Given the information about cells between skewers, as well as the line-of-sight and the two other dimension number of cells and processes, we can fully describe a skewer analysis object

```bash
class ChollaSkewerAnalysisHead:
    '''
    Cholla Skewer Analysis Head

    Holds information regarding a skewer analysis

        Initialized with:
        - n_stride (int): stride cell number between skewers
        - nlos_global (int): number of line-of-sight global cells
        - nj_global (int): number of global cells along j-dimension
        - nk_global (int): number of global cells along k-dimension
        - nlos_proc (int): number of processes along line-of-sight
        - nj_proc (int): number of processes along j-dimension
        - nk_proc (int): number of processes along k-dimension
    '''
```

With this information, we just need to provide the global id of a skewer, and we can return the ``ChollaSkewerGlobalHead`` of a skewer, which (when tied to a snapshot path and grid) allows us to grab its data.

This is precisely what is done along the x-y-z directions.

```bash
class ChollaSkewerxAnalysis:
    '''
    Cholla Skewer x Analysis object
        Holds information regarding skewer analysis along x direction including
            methods to grab ChollaSkewerxGlobal objects

        Initialized with:
        - n_stride (int): stride cell number between skewers
        - SnapPath (str): path to snapshot directory
        - chGrid (ChollaGrid): grid object
        - chMacroFlags (ChollaMacroFlags): ChollaMacroFlags, holding macro
                compiling information
    '''
```

After creating a ``ChollaSkewerxAnalysis`` object, you can pass through the global id of a skewer to its ``get_skewer`` method, which returns a ``ChollaSkewerxGlobal`` object. This object has enough information to go into the snapshot files and grab the relevant data for it. This can be used to implement skewer analysis post-simulation.



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


