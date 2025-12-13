import numpy as np
import matplotlib
import math
import timeit
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
import random
from matplotlib.animation import FuncAnimation

GRAVITATIONAL_CONSTANT = 6.67*10**-11
star_x, star_y, star_z = [], [], []
planet_1_x, planet_1_y, planet_1_z = [], [], []
frame_count = 0
time_interval = 60*60*24 # 1 day


def calculate_grav_force(m1, m2, r1, r2):
    grav_force = GRAVITATIONAL_CONSTANT * m1 * m2 * (r1-r2) / abs((r1-r2)**3)
    return grav_force

def calculate_acceleration(grav_force, mass):
    return grav_force / mass

def calculate_new_position(r_k, v_k, a_k, t):
    r_k1 = r_k + v_k * t + a_k * (t**2)
    return r_k1

def calculate_new_velocity(v_k, a_k, a_k1, t):
    v_k1 = v_k + 0.5 * t * (a_k + a_k1)
    return v_k1

"""def calculate_center_of_mass(mass_list, radius_list):
    #Assumes mass and radius are in ordered list
    if len(mass_list) != len(radius_list):
        raise ValueError("Lists are not the same length")
    sum_mass = 0
    sum_r_m = 0
    for i in len(mass_list):
        continue
    return mass_list"""

 

fig, ax = plt.axes(projection="3d") #type of coordinate system
"""x_data = np.arange(0, 50, 0.1)
y_data = np.arange(0, 50, 0.1)
z_data = x_data * y_data
ax.plot(x_data, y_data, z_data)"""


ax.plot(star_x[frame_count], star_y[frame_count], star_z[frame_count])
ax.plot(planet_1_x[frame_count], planet_1_y[frame_count], planet_1_z[frame_count])
def update(frame):
    global graph

    star_x.append()
    star_y.append()
    star_z.append()
    planet_1_x.append()
    planet_1_y.append()
    planet_1_z.appennd()
    # need to append the data

    graph.set_xdata() #put the data in the ()
    graph.set_ydata()
    graph.set_zdata()
anim = FuncAnimation(fig, update, frames=None)
plt.show()
