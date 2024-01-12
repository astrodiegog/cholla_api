'''
Produces plots for the hydrodynamic tests using cholla_api
'''

import numpy as np
import h5py
import matplotlib.pyplot as plt

from cholla_api.snap.ChollaSnap import *
from cholla_api.run.ChollaRun import *


plt.style.use("dstyle")
_ = plt.figure()
plt.style.use("dstyle")
_ = plt.figure()

##############
def calc_velocity(data):
    data["vel_x"] = data["momentum_x"] / data["density"]
    data["vel_y"] = data["momentum_y"] / data["density"]
    data["vel_z"] = data["momentum_z"] / data["density"]
    data["vel_mag"] = np.sqrt(data["vel_x"]**2 + data["vel_y"]**2 + data["vel_z"]**2)

def req_velocity(data, head):
    return ["momentum_x", "momentum_y", "momentum_z", "density"]

kwarg_fn = None

velocity_calc = ChollaValueCalc(calc_velocity, req_velocity, kwarg_fn)

##############
def calc_press(data, DE=None, gamma=None):
    if DE:
        data["pressure"] = data["Energy"]*(gamma-1.0)
    else:
        data["pressure"] = (data["Energy"] - 0.5*data["density"]*(data["vel_mag"]**2)) * (gamma - 1.0)

def req_press(data, head):
    return ["vel_mag", "Energy", "density"]

def kwargs_press(kwargs, data, head):
    kwargs["gamma"] = head["gamma"]
    kwargs["DE"] = None

press_calc = ChollaValueCalc(calc_press, req_press, kwargs_press)


tests_1D = [
    "123", "blast", "constant", "sod", "square_wave", 
    "strong_shock", "trac_pen", "two_shocks"]

tests_2D = [
    "KH_disc", "KH_resind", "Rayleigh_Taylor", "sod", 
    "sound_wave", "implosion"]

tests_3D = [
    "KH_resind", "sod_hires"]

test_dir = "/data/users/digarza/cholla_hydro_tests"


plt_kwargs = {"show": False, "save": True}

keys = ["Energy", "density", "momentum_x", "momentum_y", "momentum_z"]
movie_plots = ["density", "velocity", "pressure"]
cholla_calcs = [velocity_calc, press_calc]

imgfbase = "img"


print("Beginning 1D test plots")
for test_name in tests_1D:
    base_dir =  f"{test_dir}/1D_Test/{test_name}"

    ch_run = ChollaRun(base_dir, img_dir=img_dir)
    print(f"\tChollaRun created for {test_name}")
    
    print(f"\tStarting to plot {test_name}")
    ch_run.make_movie(keys, imgfbase, imgftype='png', test_name=test_name,
                      movie_nsnaps=None, valuecalcs=cholla_calcs, 
                      movie_plots=movie_plots, plt_kwargs=plt_kwargs)
    print(f"\tFinished with plotting {test_name}")


print("Beginning 2D test plots")
for test_name in tests_2D:
    base_dir =  f"{test_dir}/2D_Test/{test_name}"

    ch_run = ChollaRun(base_dir, img_dir=img_dir)
    print(f"\tChollaRun created for {test_name}")

    print(f"\tStarting to plot {test_name}")
    ch_run.make_movie(keys, imgfbase, imgftype='png', test_name=test_name,
                      movie_nsnaps=None, valuecalcs=cholla_calcs,
                      movie_plots=movie_plots, plt_kwargs=plt_kwargs)
    print(f"\tFinished with plotting {test_name}")


print("Beginning 3D test plots")
for test_name in tests_3D:
    base_dir =  f"{test_dir}/3D_Test/{test_name}"

    ch_run = ChollaRun(base_dir, img_dir=img_dir)
    print(f"\tChollaRun created for {test_name}")

    print(f"\tStarting to plot {test_name}")
    ch_run.make_movie(keys, imgfbase, imgftype='png', test_name=test_name,
                      movie_nsnaps=None, valuecalcs=cholla_calcs,
                      movie_plots=movie_plots, plt_kwargs=plt_kwargs)
    print(f"\tFinished with plotting {test_name}")


