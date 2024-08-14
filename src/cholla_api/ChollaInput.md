# Cholla Inputs

This file is meant to describe the inputs placed into Cholla simulations, with an emphasis on cosmology-related information. This document takes lots of inspiration from the [Cholla wiki](https://github.com/cholla-hydro/cholla/wiki/Input-File-Parameters) on similar information.


## Parameter File

The parameter file is a text file fed into a Cholla executable that specifies required information to run the desired simulation. Here is an example parameter file for a cosmological simulation:

```
nx=2048
ny=2048
nz=2048
tout=10000
outstep=1000
gamma=1.66666667
#Cosmological Parameters
H0=67.97
Omega_M=0.3069
Omega_L=0.6931
Omega_K=0.0
Omega_R=9.18240850059e-05
Omega_b=0.04853
w0=-1.0
wa=0.0
#analysis
analysis_scale_outputs_file=outputs_DESIRavoux.txt
analysisdir=analysis/
lya_skewers_stride=4
lya_Pk_d_log_k=0.0001
skewersdir=skewers/
#domain properties
xmin=0.0
ymin=0.0
zmin=0.0
xlen=50000.0
ylen=50000.0
zlen=50000.0
xl_bcnd=1
xu_bcnd=1
yl_bcnd=1
yu_bcnd=1
zl_bcnd=1
zu_bcnd=1
init=Read_Grid
# ics directory
indir=ics/
# output snapshot directory
outdir=data/
scale_outputs_file=outputs_DESIRavoux.txt
UVB_rates_file=uvb_rates_V22.txt
# density floor
density_floor=1.0e-05
# temperature floor
temperature_floor=1.0e-03
```

Going line-by-line:

- **ni** : number of grid cells in i-dimension. For 1D simulations, **ny** and **nz** should be set to 1. For 2D simulations, **nz** should be set to 1.
- **tout** & **outstep** : the final output time and time interval for outputs in code units
- **gamma** : ratio of specific heats, feel free to read up on what this value means from its [wikipedia page](https://en.wikipedia.org/wiki/Heat_capacity_ratio)
- **H0** : present-day Hubble parameter in units of [km / s / Mpc]
- **Omega_i** : present-day energy density parameter for i-species. M, L, K, R, relate to matter, dark energy, spacial curvature, and radiation (respectively), feel free to read up on what this value means from the [Friedmann equations wikipedia page](https://en.wikipedia.org/wiki/Friedmann_equations#Density_parameter)
- **w0** & **wa** : values that specify a time-evolving dark energy equation of state linearly with scale factor, of the form $w(a) = w_0 + w_a (1 - a)$
- **analysis_scale_outputs_file** : path to a text file specifying when to complete Lyman $\alpha$ (Ly $\alpha$) statistics
- **analysisdir** : path to a directory to save analysis outputs
- **lya_skewers_stride** : number of cells to jump between skewers
- **lya_Pk_d_log_k** : differential step in log k-space
- **skewersdir** : path to a directory to save skewer outputs
- **imin** : lower boundary in i-dimension direction
- **ilen** : global domain length in i-dimension direction
- **il_bcnd** & **iu_bcnd** : boundary conditions for the lower and upper i-dimension limit. Periodic, reflective, transmissive, custom boundary conditions are specified by 1, 2, 3, 4.
- **init** : name of the initial conditions, further specified in ``src/grid/initial_conditions.cpp``
- **outdir** : path to a directory to save snapshot outputs
- **scale_outputs_file** : path to a text file specifying when to output snapshots
- **UVB_rates_file** : path to a text file specifying background UV photoionization rates of HI, HeI, and HeII
- **density_floor** : minimum density that the simulation won't go under
- **temperature_floor** : minimum temperature that the simulation won't go under 


In this parameter file, we are simulating a $(2048)^3$ cubed box where each side has size $L = 50 \rm{Mpc}$. The cosmological parameter values are from [DESI 2024 VI](https://ui.adsabs.harvard.edu/abs/2024arXiv240403002D/abstract) Table 3 for a Flat $\Lambda$-CDM Universe.


## Domain Decomposition

When running a Cholla executable with a parameter file, we can also specify the number of GPU node processes $N_{\rm{proc}}$ to run the simulation on. With this number and the parameter file, we have enough information to perform a domain decomposition of the entire volume we are simulating. The domain decomposition is performed in ``src/mpi/mpi_routines.cpp`` starting in line 239. Here is a brief overview.

Given the total number of processes ($N_{\rm{proc}}$), we would like to decompose the total grid so that each process works on an individual box, where each individual box has its own much smaller domain. First, we need to tile the processes along each axis. We would like to find the number of processes that span in each dimension. Given the four numbers ($N_{\rm{proc}}$, $n_x$, $n_y$, $n_z$), we want to tile the processes along each axis so that we get three numbers ($n_{\rm{proc},x}$, $n_{\rm{proc},y}$, $n_{\rm{proc},z}$). The number $n_{\rm{proc},i}$ gives the number of processes that can work along the $i$-dimension.

In the case that the simulation is working only in 1 dimension, we place all of the processes along the x-axis. That is ($n_{\rm{proc},x}$, $n_{\rm{proc},y}$, $n_{\rm{proc},z}$) = ($N_{\rm{proc}}$, 1, 1).

In the case that the simulation is working only in 2 dimensions, we find the greatest prime factor of the number of processes and place the greatest number along the x-axis and the other number along the y-axis. That is ($n_{\rm{proc},x}$, $n_{\rm{proc},y}$, $n_{\rm{proc},z}$) = (gpf($N_{\rm{proc}}$), $N_{\rm{proc}}$ / gpf($N_{\rm{proc}}$), 1), where gpf is some greatest prime factor function.

In the case that the simulation is working along all 3 dimensions, we iteratively place the greatest prime factor along the x-y-z dimension. In pseudocode...

```
nproc_tmp = copy(nproc_total)
np_x, np_y, np_z = 1, 1, 1
index = 0
while (nproc_tmp > 1):
    # find gpf of number of processes
    n_gpf = gpf(nproc_tmp)

    # reduce number of processes by gpf
    nproc_tmp = nproc_tmp // n_gpf

    # place process by multiplying, in order of (x,y,z)
    if ((index % 3) == 0):
        np_x = int(np_x * n_gpf)
    elif ((index % 3) == 1):
        np_y = int(np_y * n_gpf)
    else:
        np_z = int(np_z * n_gpf)

    index += 1
```

Now, we have enough information to provide each box with its own subdomain of cells to work with. For each dimension, we can set a local domain of $n_i / n_{\rm{proc},i}$, where $n_i$ is the global number of cells along the $i$-dimesion and $n_{\rm{proc},i}$ is the number of processes along the $i$-dimension. There is some funkiness in the case where the number of processes doesn't evenly divide into the total number of cells, but we will ignore that for now :P For detailed explanation, check out the source code.

Okay great, we know the subdomain of a box, but we still need information about _where_ that box is within the global domain. To do this, we index each box by tiling them along the x-axis, then y-axis, then z-axis. In pseudocode...

```
n = 0 # running index counter
id_x, id_y, id_z = array(nproc_total), array(nproc_total), array(nproc_total)
for k in range(np_z):
    for j in range(np_y):
        for i in range(np_x):
            id_x[n], id_y[n], id_z[n] = i, j, k
            n += 1

# indices of each box will range as : [ (0,0,0), (1,0,0), ... (np_x,0,0), 
#                                     (0,1,0), (1,1,0), ... (np_x,1,0),
#                                     ...
#                                     (0,np_y,0), (1,np_y,0), ... (np_x, np_y, 0),
#                                     ...
#                                     (0,0,1), (1,0,1), ... (np_x,0,1),
#                                     (0,1,1), (1,1,1), ... (np_x,1,1),
#                                     ...
#                                     (0,np_y,1), (1,np_y,1), ... (np_x, np_y, 1),
#                                     ...
#                                     (0,0,np_z), (1,0,np_z), ... (np_x,0,np_z),
#                                     (0,1,np_z), (1,1,np_z), ... (np_x,1,np_z),
#                                     ...
#                                     (0,np_y,np_z), (1,np_y,np_z), ... (np_x, np_y, np_z)]
```

The way that these indices are assigned are in a row-major order, used by Numpy and C/C++. Read more about this style of memory allocation of multidimensional arrays in a linear manner in its [wikipedia page](https://en.wikipedia.org/wiki/Row-_and_column-major_order).

Now that we have the tiling order is determined and local subdomain defined, we can find the offset of a box given the box number by multiplying the local subdomain by its index along each dimension. In pseudocode...

```
offset_x = id_x[box_num] * (nx / np_x)
offset_y = id_y[box_num] * (ny / np_y)
offset_z = id_z[box_num] * (nz / np_z)
```

To summarize,
1. Given the number of processes $N_{\rm{proc}}$ and whether we are simulating in 1, 2, or 3 dimensions, we determine the number of processes along each dimension ($n_{\rm{proc},x}$, $n_{\rm{proc},y}$, $n_{\rm{proc},z}$).
2. For each process, we assign its index along each dimension in a row-major order 
3. To define the offset of a process, we multiply its index along each dimension by the local subdomain dimensions



## Skewer Decomposition

When running a Cholla executable that was compiled with the ``ANALYSIS``, ``PHASE_DIAGRAM``, ``LYA_STATISTICS``, and/or ``OUTPUT_SKEWERS`` macro flags, the simulation will produce skewers of the length of the entire box to compute the one-dimensional flux power spectrum. Computing the Ly $\alpha$ statistics is done across each xyz dimension. Each skewer's Ly $\alpha$ statistics is a pretty time-intensive process, so we do not produce skewers for every cell. Instead, we take an $n_{\rm{stride}}$ number of cell steps between skewers. The Ly $\alpha$ statistics code starts on line 12 of the ``src/analysis/analysis.cpp``, with details in ``src/analysis/lya_statistics.cpp``.


For the analysis along an $i$-dimension, the total number of skewers is $N_{\rm{skewers}} = (n_j / n_{\rm{stride}}) (n_k / n_{\rm{stride}})$, where $n_{j,k}$ is the global number of cells along the $j,k$-dimension.

Within some box, we have $N_{\rm{box, skewers}} = (n_{\rm{box},j} / n_{\rm{stride}}) (n_{\rm{box},k} / n_{\rm{stride}})$ number of skewers, where $n_{\rm{box},j}$ and $n_{\rm{box},k}$ are the number of cells along the $j,k$-dimension within the local domain of the box. These two values are computed on line 1364 of ``src/analysis/lya_statistics.cpp``.


Within some box, the _local_ offset along the $j$-dimension for a skewer with id $\eta_{\rm{local}}$ is calculated as $j_{\rm{local, offset}} = \rm{int}( \eta_{\rm{local}} // ( (n_{\rm{box},j} / n_{\rm{stride}}) ) \cdot n_{\rm{stride}}) )$. The _local_ offset along the $k$-dimension is $k_{\rm{local, offset}} = \eta( \rm{id}_{\rm{local}} \rm{mod} ( (n_{\rm{box},j} / n_{\rm{stride}}) ) \cdot n_{\rm{stride}}) )$. What drives these offsets is how the local skewers are tiled first along the $j$-dimension, then along the $k$-dimension. These values are computed in the function ``Grid3D::Populate_Lya_Skewers_Local`` on line 1226 of ``src/analysis/lya_statistics.cpp``.

Within the global $j-k$ plane, batches of $N_{\rm{box, skewers}}$ are tiled first along the $k$-dimension and $j$-dimension. For some global skewer with $\rm{id}_{\rm{global}}$, the box in which it is found is $\rm{id}_{\rm{box}} = \rm{int} (\rm{id}_{\rm{global}} // N_{\rm{box, skewers}} )$. In effect, the _global_ offsets for some box along the $j$-dimension is $j_{\rm{box, offset}} = \rm{int}( \rm{id}_{\rm{global}} \rm{mod} ( n_{\rm{proc},j} ) \cdot n_{\rm{box},j} )$. Likewise, the offset along the $k$-dimension is $k_{\rm{box, offset}} = \rm{int}( \rm{id}_{\rm{global}} // ( n_{\rm{proc},j} ) \cdot n_{\rm{box},k}) )$. These values are computed in the function ``AnalysisModule::Transfer_Skewers_Data(int axis)`` on line 1037 of ``src/analysis/lya_statistics.cpp``.

For some global skewer with $\rm{id}_{\rm{global}}$, the total $j$-offset is $j_{\rm{global, offset}} = j_{\rm{box, offset}} + j_{\rm{local, offset}}$ and the total $k$-offset is $k_{\rm{global, offset}} = k_{\rm{box, offset}} + k_{\rm{local, offset}}$.

For the analysis along the $x$-dimension, we have $(i,j,k)$ as $(x,y,z)$. For the analysis along the $y$-dimension, we have $(i,j,k)$ as $(y,x,z)$. For the analysis along the $z$-dimension, we have $(i,j,k)$ as $(z,x,y)$. These are chosen in ``Grid3D::Populate_Lya_Skewers_Local`` of ``src/analysis/lya_statistics.cpp``.



## Power Spectrum bins

To create a one-dimensional flux power spectrum, the parameter file let's users specify some **lya_Pk_d_log_k**, which acts as the differential step in log $k$-space. The $k$-bins are specified in the ``Grid3D::Initialize_Power_Spectrum_Measurements(int axis)`` function on line 330 of ``src/analysis/ly_statistics.cpp``. Here is a brief overview of how the bins are created from information within the parameter file.

First we have to define the range the histogram. The minimum is calculated as

$$
\log_{10} k_{\rm{min}} = \log_{10} \left(\frac{2 \pi}{n_{\rm{LOS}} \rm{d}v_{H}} \right)
$$

where $n_{\rm{LOS}}$ is the number of line-of-sight cells and $\rm{d}v_{H}}$ is the differential Hubble flow through one cell. Likewise, the maximum is calculated as

$$
\log_{10} k_{\rm{max}} = \log_{10} \left(\frac{2 \pi (n_{\rm{FFT}} - 1)}{n_{\rm{LOS}} \rm{d}v_{H}} \right)
$$

where $n_{\rm{FFT}} = \rm{int} ( n_x / 2 + 1 )$ is the number of FFT bins.

In the Cholla code, we have 

$$
\log_{10} k_{\rm{start}} = \log_{10} (0.99) + \log_{10} k_{\rm{min}}
$$

The number of $k$-bins is 

$$
n_{k, \rm{bins}} = \rm{int}(1 + \frac{\log_{10}(k_{\rm{max}}) - \log_{10}(k_{\rm{start}}) }{\rm{d} \log_{10} k })
$$


We can recreate the $k$ center values in pseudocode as...

```
iter_arr = range(n_kbins)
kcenters = 10**(log10_kstart + dlogk * (iter_arr + 0.5) )
```





