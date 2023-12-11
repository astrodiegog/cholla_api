# Cholla Visualization

The idea of this `viz` module is that all values have already been loaded and calculated using ChollaRun and ChollaSnap. Now, we would like to visualize our results. There's no more calculation, simply reading files and plotting them.

## ChollaViz

`ChollaViz` - Main parent class that holds the snapshot and number of dimensions. 


## ChollaXDViz

Depending on the number of dimensions used in a simulation run, `ChollaXDViz` use different classes (`Cholla1DViz`, `Cholla2DViz`) that all generally have the same functions:
- pressure
- density
- velocity_x
- velocity_y
- velocity

The ability to compare one of these values within the same simulation run is also possible using the following functions
- density_compare
- pressure_compare
- velocity_compare
- velocityx_compare
- velocityy_compare

## Cholla1DViz, Cholla2DViz, Cholla3DViz

`Cholla2DViz` and `Cholla1DViz` - sub classes created when instantiating `ChollaXDViz`. Each of these have two primary plotting functions: `plot_value` and `plot_value_compare`. These two take in some data, head, plotting format, and plotting arguments.

The plotting format `plt_fmt` is a dictionary that handles the plotting for a specific value and test sample. The general idea is that each simulation run AND each plotting value is going to have different formats. This is meant to be set within its own plotting formatting function that can be saved for tests (`plot_pressure_format` function may be different compared to `plot_density_format` depending on the test sample).
    
The idea is that if you see something you'd like to change after making a plot of a specific test sample, then you create a new key for `plt_fmt` with different values for whatever run or values you're plotting. Within \__init_\_, you can place a default value.

For example, if I wanted to change the title fontsize in density for the 1D version running the sod test, I would add "title_fs" as a key to `plt_fmt`, assign specific font sizes to "sod" within the `density_format` function, and define a default value in \___init__\_.

In the `Cholla1DViz`, `plt_fmt` holds
- value_key: the key used to access a snapshot's data
- y_label OR title: the label on the y-axis OR the title on the plots
- fig_size: the figure size
- yval_fmt: the format that the y value is shown in

In the `Cholla2DViz`, `plt_fmt` holds
- value_key: the key used to access a snapshot's data
- colorbar_fmt: the format on the colorbar
- fig_size: the figure size
- title: the title name

In general, the plotting arguments (`plt_kwargs`) handles what happens with the plot after it has been plotted. Currently, it handles whether the user would like to 1) show the plot and 2) save the plot.


The goal of `ChollaxDViz`, compared to these subclasses, is to handle different dimensions. So `Cholla2DViz` _could_ be used standalone, just need to be careful.


## viz_format

This is supposed to hold any formatting helper functions. Currently holds functions for converting a number to a string that can be used on axes and colorbars.


## Add more visualizations

To add another value to visualize in this framework (like internal energy) that is within the data of the snapshot...

0. Goal is to use `plot_value`, so only special thing is to specify plotting formatting.
1. Figure out some plotting formatting. Add a function for `Cholla1DVizFmt`, `Cholla2DVizFmt`, and `Cholla3DVizFmt` that specifies the formatting. Check `plot_value` function to see the requirements to use `plot_value`.
2. Add a function in `ChollaxDViz` that sets this new plotting function's format. Within function, pass values onto `plot_value`.
3. Add a function in `ChollaViz` that follows the same strategy as the others: define data's key, a string for the key, check the data is in the snapshot, assign a `plt_kwargs` if needed. Set a call to `chollaxDViz.value`

Boom! Now able to easily plot your new value.
