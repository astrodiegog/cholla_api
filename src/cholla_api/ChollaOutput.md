# Cholla Outputs

This file is meant to describe the outputs resulted from Cholla simulations, with an emphasis on cosmology-related information. This document takes lots of inspiration from [Cholla wiki](https://github.com/cholla-hydro/cholla/wiki/Output) on similar information.

## File Type

By default, Cholla will place snapshot information in the 5th version of Hierarchical Data Format, better known as [HDF5](https://github.com/HDFGroup/hdf5). While the official [HDF Group API website](https://docs.hdfgroup.org/hdf5/v1_14/_getting_started.html) is very complete and contains detailed information about the data type, I (Diego G) have personally found the documentation page for [h5py](https://docs.h5py.org/en/stable/), the pythonic interface with HDF5 binary data file format, to be most helpful in getting started with this data format (probably due to how convenient I find Python).

In short, Group objects can hold other Group object or Dataset objects in a hierarchical manner. The File Object acts as the _root group_. A Dataset object functions very similarly to numpy arrays. Group and Dataset objets can have metadata saved to it in the form of Attributes. 

For our cases, Cholla snapshots will (in general) provide global and local information regarding the placement of a box within Attributes header of a file, with the actual raw data being held in Dataset objects.

## Data

This section will describe the files that contain a snapshot's data.

### File Structure

Within the data directory where the parameter text file specified to save snapshot outputs, Cholla will create individual directories for each snapshot. Within each snapshot subdirectory, there is an HDF5 file saved for each process / box. Depending on the compiler flags used in creating the Cholla executable, each process / box may save additional information apart from just the hydro information. For example, the `make.type.cosmology` make file details to run the simulation with gravity and particles, so each process / box will also save the gravity and particles output. 

An example tree structure with ``nprocs`` processes / boxes and ``nsnaps`` number of snapshots detailed to save from the parameter file, is shown below...

```bash
/data/
├── 1
│   ├── 1_gravity.h5.0
│   ├── 1_gravity.h5.1
│   ├── ...
│   ├── 1_gravity.h5.nprocs
│   ├── 1.h5.0
│   ├── 1.h5.1
│   ├── ...
│   ├── 1.h5.nprocs
│   ├── 1_particles.h5.0
│   ├── 1_particles.h5.1
│   ├── ...
│   └── 1_particles.h5.nprocs
├── 2
├── ...
└── nsnaps
```

Files holding the hydro information have the name ``nsnap.h5.proc``, where ``nsnap`` and ``nproc`` specify the snapshot number and process / box number (``4.h5.213`` holds hydro information for process 213 in snapshot 4). Files holding the gravity and particles data suffix the snapshot number with their respective names (``63_gravity.h5.200`` holds gravity information for process 200 in snapshot 63, ``75_particles.h5.3`` holds particle information for process 3 in snapshot 75).


### Hydro Files

Hydro files do not have any subgroups, so their File Objects serve as the root group.

The attributes that are attached are:

- **Macro Flags** : the compiler flags for the Cholla executable
- **bounds** : boundary conditions, further specified in ``src/grid/boundary_conditions.cpp``
- **dims** : number of grid cells in each of the three dimension for entire volume
- **dims_local** : number of grid cells in the box
- **offset** : number of grid cells offset from origin within the entire volume
- **n_fields** : number of fields (conserved and scalars) kept track during simulation, further specified in ``src/grid/grid3D.h``
- **n_step** : number of time steps taken, further specified in ``src/grid/grid3D.h`` 
- **nprocs** : number of processes in each of the three dimensions
- **gamma** : ratio of specific heats 

The following attributes are related to "code units". The reason for "code units" is because Cholla is a versatile program that can run in different units. For a description of the units allowed in Cholla, they are further specified in ``src/global/global.h``.

- **dx** : distance of a cell in each of the three dimensions, in code units
- **dt** : time since last time step, in code units
- **t** : time of the snapshot, in code units
- **density_unit** : conversion from density code units to cgs units [g cm-3]
- **energy_unit** : conversion from energy code units to cgs units [erg]
- **mass_unit** : conversion from mass code units to cgs units [g]
- **length_unit** : conversion from length code units to cgs units [cm]
- **velocity_unit** : conversion from velocity code units to cgs units [cm s-1]
- **time_unit** : conversion from time code units to cgs units [s]

The following conserved fields are saved as datasets in code units

- **density** : mass density in each cell
- **momentum_x** : x-momentum density
- **momentum_y** : y-momentum density
- **momentum_z** : z-momentum density 
- **Energy** : energy density

If the Dual Energy flag was used in compiling the Cholla executable, the thermal energy field will also be saved. To read more on the motivation for saving this value, read the entry on the topic in the [Cholla wiki](https://github.com/cholla-hydro/cholla/wiki/Dual-Energy-Formalism)

- **GasEnergy** : thermal energy density in code units


#### Cosmology Notes

When running with a Cholla executable compiled with the ``COSMOLOGY`` macro flag, these are some important notes regarding the attributes.

The code units used that are attached to the attributes:

- **density_unit** : Solar masses per cubic kiloparsec in grams per cubic centimeter --> $M_{\odot} \rm{kpc}^{-3} = 6.76 \times 10^{-32} \rm{g} \rm{cm}^{-3}$
- **length_unit** : Kiloparsecs in centimeters --> $\rm{kpc} = 3.09 \times 10^{21} \rm{cm}$
- **velocity_unit** : Kiloparsecs per kiloyers --> $\rm{kpc} / \rm{kyr} = 9.78 \times 10^{10} \rm{cm} / \rm{s}$
- **time_unit** : Kiloyears in seconds --> $\rm{kyr} = 3.16 \times 10^{10} \rm{s}$
- **mass_unit** : Solar masses in grams --> $M_{\odot} = 1.99 \times 10^{33} \rm{g}$

**HOWEVER**: this is not what the actual units of the values that are saved. To convert from the density that is saved in snapshots to cgs, the following unit should be used

- **density_unit** : Solar masses per cubic kiloparsecs in grams per cubic centimter, scaled by $h^2$ --> $h^2 M_{\odot} \rm{kpc}^{-3} = 6.76 \times 10^{-32} h^2 \rm{g} \rm{cm}^{-3}$

The velocity units is

- **velocity_unit** : (peculiar velocity) Kilometers per second in centimeters per second --> $\rm{km} \rm{s}^{-1} = 1.00 \times 10^5 \rm{cm} \rm{s}$

To convert from the saved momentum to cgs, the mass density is multiplied by the peculiar velocity

- **momentum_unit** : Product of kilometers per second and solar mass per cubic kiloparsecs in grams per centimeter-squared per second, scaled by $h^2$ --> $h^2 \rm{km} \rm{s}^{-1} M_{\odot} \rm{kpc}^{-3} = (1.00 \times 10^{5} \rm{cm} \rm{s}^{-1})(6.76 \times 10^{-32} h^2 \rm{g} \rm{cm}^{-3}) = 6.76 \times 10^{-27} h^2 \rm{g} \rm{cm}^{-2} \rm{s}^{-1}$

To convert the saved energy to cgs, we scale the momentum unit by another factor of velocity

- **energy_unit** : Product of kilometers per second and kilometers per second and solar mass per cubic kiloparsecs in centiemeters grams per second, scaled by $h^2$ --> $h^2 \rm{km} \rm{s}^{-1} \rm{km} \rm{s}^{-1} M_{\odot} \rm{kpc}^{-3} = 6.76 \times 10^{-22} h^2 \rm{cm} \rm{g} \rm{s}^{-1}$

Likewise, the unit of length and mass should be scaled by $h$

- **length_unit** : Kiloparsecs in centimeters, scaled by $h$ --> $h^{-1}\rm{kpc} = h^{-1} 3.09 \times 10^{21} \rm{cm}$
- **mass_unit** : Solar masses in grams, scaled by $h$ --> $h^{-1} M_{\odot} = h^{-1} 6.76 \times 10^{-32} \rm{g}$

Such that dividing mass unit by length unit-cubed results back in the density unit of $h^2 M_{\odot} \rm{kpc}^{-3}$

The following attributes are also included

- **Current_a** : scale factor at the time of the snapshot
- **Current_z** : redshift at the time of the snapshot
- **H0** : present-day Hubble parameter in units of [km / s / kpc]
- **Omega_L** : present-day energy density parameter for dark energy
- **Omega_M** : present-day energy density parameter for matter

**LASTLY**, the values are ran in super-comoving units and saved in comoving units with their scale factor scaling still attached. We can turn to equations 15-21 of [Martell and Shapiro 1998](https://academic.oup.com/mnras/article/297/2/467/988380) and equations 22-24 of [Teyssier 2002](https://ui.adsabs.harvard.edu/abs/2002A%26A...385..337T/abstract) to notice that we convert from the saved comoving density $\rho_{\rm{comov}}$ and the proper density $\rho_{\rm{proper}}$ with $\rho_{\rm{comov}} = a^3 \rho_{\rm{proper}}$. Likewise, we convert between the saved comoving distance $x_{\rm{comov}}$ (such as the dimensions of the box) and the proper distance $x_{\rm{proper}}$ with $x_{\rm{comov}} = a^{-1}x_{\rm{proper}}$.


### Particle Files

Particle files do not have any subgroups, so their File Objects serve as the root group.

The particle file inherits all the same attributes from the hydro file, with some additions:

- **n_particles_local** : the number of particles in the box
- **particle_mass** : mass of the particles in code units
- **t_particles** : time that the particles use in code units

The following datasets are saved

- **particle_IDs** : id for each particle in the box
- **density** : the density in each cell calculated using the Cloud-in-Cell scheme
- **pos_x** : x-position of particles from origin, in code units
- **pos_y** : y-position of particles from origin, in code units
- **pos_z** : z-position of particles from origin, in code units
- **vel_x** : x-peculiar velocity of particles, in code units
- **vel_y** : y-peculiar velocity of particles, in code units
- **vel_z** : z-peculiar velocity of particles, in code units

The implementation of the Cloud-in-Cell scheme is further specified in ``src/particles/density_CIC.cpp`` To gain an understanding of how the Cloud-in-Cell scheme is computed (going from particle position to discrete density cells), I (Diego G.) have found [Birsdall and Fuss 1968](http://www2.imm.dtu.dk/courses/FortranMPI/Projects/N-body/Birdsall:1969.pdf) and Professor Andrey Kravtsov's [Writing a PM code talk](https://astro.uchicago.edu/~andrey/talks/PM/pm.pdf) particularly helpful.

#### Cosmology Notes

The same discussion on units from hydro also applies to the particles such that the density is in units of $h^2 M_{\odot} \rm{kpc}^{-3}$, velocity is in units of $\rm{km} \rm{s}^{-1}$, position is in units of $h^{-1} \rm{kpc}$, and mass is in units of $h^{-1} M_{\odot}$


### Gravity Files


From line 142 in ``src/grav/grav3D.h``
Device Array containing the gravitational potential of each cell
     * in the grid


## Analysis

This section will describe the files that contain the Lyman $\alpha$ analysis outputs.


### Skewer Files

Files holding the skewer information will have the name ``nSkewerOutput_skewers.h5``, where ``nSkewerOutput`` is the number of the skewer output.

The root Skewer File Object has the following attributes:

- **H0** : present-day Hubble parameter in units of [km / s / Mpc]
- **Omega_L** : present-day energy density parameter for dark energy
- **Omega_M** : present-day energy density parameter for matter
- **Omega_R** : present-day energy density parameter for radiation
- **Omega_K** : present-day energy density parameter for spatial curvature
- **Omega_b** : present-day energy density parameter for  baryonic matter
- **w0** & **wa** : values that specify a time-evolving dark energy equation of state linearly with scale factor, of the form $w(a) = w_0 + w_a (1 - a)$
- **current_a** : scale factor at the time of the skewer output
- **current_z** : redshift at the time of the skewer output
- **Lbox** : global domain length in all three dimensions in units of [kpc / h]

Skewer files have three Groups: ``skewers_x``, ``skewers_y``, and ``skewers_z``.

Each of the three Groups have not attributes.

The following datasets are saved for each skewer Group

- **density** : mass density in each cell, in code units
- **HI_density** : mass density of all Hydrogen species in each cell, in code units
- **HeI_density** : mass density of all Helium species in each cell, in code units
- **los_velocity** : line-of-sight velocity in each cell along the skewer, in code units
- **temperature** : temperature in each cell, in Kelvin
- **vel_Hubble** : the Hubble flow velocity along line-of-sight of skewers, in proper units [km / s]

The density are saved in comoving units, the line-of-sight velocity is the peculiar velocity, and last velocity is the Hubble flow. 

For the group ``skewer_i``, each dataset will have the shape $((n_j / n_{\rm{nstride}})(n_k / n_{\rm{nstride}}), n_i )$, where $n_i$ is the number of grid cells in $i$-dimension.


### Analysis Files

Files holding the analysis information will have the name ``nAnalysisOutput_skewers.h5``, where ``nAnalysisOutput`` is the number of the analysis output.

The root Analysis File Object inherits the same attributes as the Skewer File Object.

Analysis files have two Groups: ``lya_statistics`` and ``phase_diagram``.

The ``lya_statistics`` Group has three attributes:

- **Flux_mean_HI** : mean transmitted flux of Hydrogen in dimensionless units
- **Flux_mean_HeII** : mean transmitted flux of Helium in dimensionless units
- **n_skewers** : total number of skewers used to complete analysis

The ``lya_statistics`` Group holds 1 Group: ``power_spectrum``.

The ``power_spectrum`` Group has no attributes.

The ``power_spectrum`` Group saves the following datasets

- **k_vals** : center $k$ mode values of the $k$ histogram
- **p(k)** : flux power spectrum


The ``phase_diagram`` Group has six attributes:

- **dens_min** : log10 of the minimum density for the phase space diagram
- **dens_max** : log10 of the maximum density for the phase space diagram
- **n_dens** : number of log10-spaced density histogram bins
- **temp_min** : log10 of the minimum temperature for the phase space diagram
- **temp_max** : log10 of the maximum temperature for the phase space diagram
- **n_temp** : number of log10-spaced temperature histogram bins

The ``phase_diagram`` Group saved the following dataset

- **data** : 2D histogram of the phase space diagram

The information held in the ``phase_diagram`` Group let's us reproduce plots similar to [Villasenor et al 2021](https://ui.adsabs.harvard.edu/abs/2021ApJ...912..138V/abstract) Figure 8 and Figure 10. Figure 15 and Figure 16 can be reproduced with information held in ``lya_statistics``. Figure 9 and Figure 13 can be reproduced from Skewer File Objects.



