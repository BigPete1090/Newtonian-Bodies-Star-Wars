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



star_x, star_y, star_z = 0, 0, 0
star_mass = 2e30
star_velocity_x, star_velocity_y, star_velocity_z = 0, 0, 0

planet_1_mass = 5.972e24
planet_1_x, planet_1_y, planet_1_z = 1, 0, 1
planet_1_velocity_x, planet_1_velocity_y, planet_1_velocity_z = 212100, 0, 212100 #m/s

frame_count = 0
time_interval = 60*60*24 # 1 day


def calculate_grav_force(m1, m2, r1, r2, component_1, component_2):
    grav_force = GRAVITATIONAL_CONSTANT * m1 * m2 * (component_1-component_2) / abs((r1-r2)**3)
    return grav_force

def au_to_meter(au):
    m = float(au) * 1.496e11
    return m

def meter_to_au(m):
    au = float(m) / 1.496e11
    return au

def calculate_acceleration(grav_force, mass):
    return grav_force / mass

def calculate_new_position(r_k, v_k, a_k, t):
    r_k1 = r_k + v_k * t + a_k * (t**2)
    return r_k1

def calculate_new_velocity(v_k, a_k, a_k1, t):
    v_k1 = v_k + 0.5 * t * (a_k + a_k1)
    return v_k1

def calculate_distance_from_center(x, y, z):
    distance = np.sqrt(x**2+y**2+z**2)
    return distance

"""def calculate_center_of_mass(mass_list, radius_list):
    #Assumes mass and radius are in ordered list
    if len(mass_list) != len(radius_list):
        raise ValueError("Lists are not the same length")
    sum_mass = 0
    sum_r_m = 0
    for i in len(mass_list):
        continue
    return mass_list"""

for x in range(10):
    print(f"Star location: {star_x}, {star_y}, {star_z}")
    print(f"Planet location: {planet_1_x}, {planet_1_y}, {planet_1_z}")
    print(f"Star: velocities: {star_velocity_x}, {star_velocity_y}, {star_velocity_z}")
    print(f"Planet: velocities: {planet_1_velocity_x}, {planet_1_velocity_y}, {planet_1_velocity_z}")

    star_x, star_y, star_z = au_to_meter(star_x), au_to_meter(star_y), au_to_meter(star_z)

    star_distance = calculate_distance_from_center(star_x, star_y, star_z)
    planet_1_distance = calculate_distance_from_center(planet_1_x, planet_1_y, planet_1_z)
    

    grav_force_x = calculate_grav_force(star_mass, planet_1_mass, star_distance, planet_1_distance, star_x, planet_1_x)
    grav_force_y = calculate_grav_force(star_mass, planet_1_mass, star_distance, planet_1_distance, star_y, planet_1_y)
    grav_force_z = calculate_grav_force(star_mass, planet_1_mass, star_distance, planet_1_distance, star_z, planet_1_z)


    star_acceleration_x = calculate_acceleration(grav_force_x, star_mass)
    planet_1_acceleration_x = calculate_acceleration(grav_force_x, planet_1_mass)
    star_acceleration_y = calculate_acceleration(grav_force_y, star_mass)
    planet_1_acceleration_y = calculate_acceleration(grav_force_y, planet_1_mass)
    star_acceleration_z = calculate_acceleration(grav_force_z, star_mass)
    planet_1_acceleration_z = calculate_acceleration(grav_force_z, planet_1_mass)

    next_star_position_x = calculate_new_position(star_x, star_velocity_x, star_acceleration_x, time_interval)
    next_planet_position_x = calculate_new_position(planet_1_x, planet_1_velocity_x, planet_1_acceleration_x, time_interval)
    next_star_position_y = calculate_new_position(star_y, star_velocity_y, star_acceleration_y, time_interval)
    next_planet_position_y = calculate_new_position(planet_1_y, planet_1_velocity_y, planet_1_acceleration_y, time_interval)
    next_star_position_z = calculate_new_position(star_z, star_velocity_z, star_acceleration_z, time_interval)
    next_planet_position_z = calculate_new_position(planet_1_z, planet_1_velocity_z, planet_1_acceleration_z, time_interval)

    next_star_distance = calculate_distance_from_center(next_star_position_x, next_star_position_y, next_star_position_z)
    next_planet_distance = calculate_distance_from_center(next_planet_position_x, next_planet_position_y, next_planet_position_z)

    next_grav_force_x = calculate_grav_force(star_mass, planet_1_mass, next_star_distance, next_planet_distance, next_star_position_x, next_planet_position_x)
    next_grav_force_y = calculate_grav_force(star_mass, planet_1_mass, next_star_distance, next_planet_distance, next_star_position_y, next_planet_position_y)
    next_grav_force_z = calculate_grav_force(star_mass, planet_1_mass, next_star_distance, next_planet_distance, next_star_position_z, next_planet_position_z)

    next_star_acceleration_x = calculate_acceleration(next_grav_force_x, star_mass)
    next_planet_1_acceleration_x = calculate_acceleration(next_grav_force_x, planet_1_mass)
    next_star_acceleration_y = calculate_acceleration(next_grav_force_y, star_mass)
    next_planet_1_acceleration_y= calculate_acceleration(next_grav_force_y, planet_1_mass)
    next_star_acceleration_z = calculate_acceleration(next_grav_force_z, star_mass)
    next_planet_1_acceleration_z = calculate_acceleration(next_grav_force_z, planet_1_mass)

    star_velocity_x = calculate_new_velocity(star_velocity_x, star_acceleration_x, next_star_acceleration_x, time_interval)
    planet_1_velocity_x = calculate_new_velocity(planet_1_velocity_x, planet_1_acceleration_x, next_planet_1_acceleration_x, time_interval)
    star_velocity_y = calculate_new_velocity(star_velocity_y, star_acceleration_y, next_star_acceleration_y, time_interval)
    planet_1_velocity_y = calculate_new_velocity(planet_1_velocity_y, planet_1_acceleration_y, next_planet_1_acceleration_y, time_interval)
    star_velocity_z= calculate_new_velocity(star_velocity_z, star_acceleration_z, next_star_acceleration_z, time_interval)
    planet_1_velocity_z = calculate_new_velocity(planet_1_velocity_z, planet_1_acceleration_z, next_planet_1_acceleration_z, time_interval)

    

    star_x, star_y, star_z = meter_to_au(next_star_position_x), meter_to_au(next_star_position_y), meter_to_au(next_star_position_z)
    planet_1_x, planet_1_y, planet_1_z = meter_to_au(next_planet_position_x), meter_to_au(next_planet_position_x), meter_to_au(next_planet_position_z)

    print(f"New Star location: {star_x}, {star_y}, {star_z}")
    print(f"New Planet location: {planet_1_x}, {planet_1_y}, {planet_1_z}")
    print(f"Star: New velocities: {star_velocity_x}, {star_velocity_y}, {star_velocity_z}")
    print(f"Planet: New velocities: {planet_1_velocity_x}, {planet_1_velocity_y}, {planet_1_velocity_z}")



# fig, ax = plt.axes(projection="3d") #type of coordinate system
"""x_data = np.arange(0, 50, 0.1)
y_data = np.arange(0, 50, 0.1)
z_data = x_data * y_data
ax.plot(x_data, y_data, z_data)"""


"""ax.plot(star_x[frame_count], star_y[frame_count], star_z[frame_count])
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
plt.show()"""
