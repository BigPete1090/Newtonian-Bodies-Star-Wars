import numpy as np
import matplotlib
import math
import timeit
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d
import random
from matplotlib.animation import FuncAnimation
import csv
from collections import defaultdict


GRAVITATIONAL_CONSTANT = 6.67e-11
all_objects = []

objects = 0
with open('objects.csv', mode='r') as file:
    csv_reader = csv.DictReader(file, delimiter=',')
    for row in csv_reader:
        name_obj = row['Name']
        dict_fun = {"name": row['Name'], 'mass': float(row['Mass (kg)']),
                                'x': float(row['X (AU)']), 'y': float(row['Y (AU)']), 'z': float(row['Z (AU)']), 
                                'velocity_x': float(row['Vx (m/s)']), 'velocity_y': float(row['Vy (m/s)']),  'velocity_z': float(row['Vz (m/s)']), 
                                'color': row['Color']}
        print(f"{dict_fun}")
        all_objects.append(dict_fun)



frame_count = 1
time_interval = (60*60*24) # NEED TO FIX

# fix to accept one argument for r not 2
def calculate_grav_force(m1, m2, r, component_1, component_2):
    grav_force = GRAVITATIONAL_CONSTANT * float(m1) * float(m2) * (float(component_2) - float(component_1)) / (r**3)
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

def calculate_distance_between_objects(x1, y1, z1, x2, y2, z2):
    distance = np.sqrt(((x2-x1)**2)+((y2-y1)**2)+((z2-z1)**2))
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

def run_sim():
    global star_x, star_y, star_z
    global planet_1_x, planet_1_y, planet_1_z
    global star_mass, planet_1_mass
    global star_velocity_x, star_velocity_y, star_velocity_z
    global planet_1_velocity_x, planet_1_velocity_y, planet_1_velocity_z
    global all_objects



    distances_from_center = {}
    grav_forces = {}
    next_distances_from_center = {}
    all_current_accelerations = defaultdict(float)
    all_next_distances_each_direction = defaultdict(float)

    all_next_accelerations = defaultdict(float)
    next_grav_forces = {}



    for i in range(len(all_objects)):
        object_name = all_objects[i]['name']
        object_x_m = au_to_meter(all_objects[i]['x'])
        object_y_m = au_to_meter(all_objects[i]['y'])
        object_z_m = au_to_meter(all_objects[i]['z'])

        object_distance = calculate_distance_from_center(object_x_m, object_y_m, object_z_m)
        distances_from_center[object_name] = object_distance # Should make the object name key


    for i in range(len(all_objects)):
        object_name = all_objects[i]['name']
        object_mass = all_objects[i]['mass']
        object_x = all_objects[i]['x']
        object_y = all_objects[i]['y']
        object_z = all_objects[i]['z']
        object_vx = all_objects[i]["velocity_x"]
        object_vy = all_objects[i]["velocity_y"]
        object_vz = all_objects[i]["velocity_z"]
        

        object_x_m = au_to_meter(object_x)
        object_y_m = au_to_meter(object_y)
        object_z_m = au_to_meter(object_z)

        object_distance = distances_from_center[object_name]
        
        for j in range(len(all_objects)):
            if j != i:
                second_object_name = all_objects[j]['name']
                print(f"Calculating the gravity between object {object_name} and {second_object_name}")
                second_object_mass = all_objects[j]['mass']
                second_object_x = all_objects[j]['x']
                second_object_y = all_objects[j]['y']
                second_object_z = all_objects[j]['z']

                second_object_x_m = au_to_meter(second_object_x)
                second_object_y_m = au_to_meter(second_object_y)
                second_object_z_m = au_to_meter(second_object_z)

                object_distances = calculate_distance_between_objects(object_x_m, object_y_m, object_z_m, 
                                                                      second_object_x_m, second_object_y_m, second_object_z_m)


                # second_object_distance = distances_from_center[second_object_name]

                if (grav_forces.get((object_name, second_object_name, "x"), None) is None or
                    grav_forces.get((object_name, second_object_name, "y")) is None or 
                    grav_forces.get((object_name, second_object_name, "z")) is None): 

                    obj_grav_force_x = calculate_grav_force(object_mass, second_object_mass, object_distances, 
                                                             object_x_m, second_object_x_m)
                    obj_grav_force_y = calculate_grav_force(object_mass, second_object_mass, object_distances, 
                                                             object_y_m, second_object_y_m)
                    obj_grav_force_z = calculate_grav_force(object_mass, second_object_mass, object_distances, 
                                                             object_z_m, second_object_z_m)
                    second_obj_grav_force_x = -obj_grav_force_x
                    second_obj_grav_force_y = -obj_grav_force_y
                    second_obj_grav_force_z = -obj_grav_force_z

                    grav_forces[object_name, second_object_name, "x"] = obj_grav_force_x
                    grav_forces[object_name, second_object_name, "y"] = obj_grav_force_y
                    grav_forces[object_name, second_object_name, "z"] = obj_grav_force_z

                    grav_forces[second_object_name, object_name, 'x'] = second_obj_grav_force_x
                    grav_forces[second_object_name, object_name, 'y'] = second_obj_grav_force_y
                    grav_forces[second_object_name, object_name, 'z'] = second_obj_grav_force_z


                
                    # needs to be inside this loop because if the acceleration due to second object already calculated, it is done
                

                    all_current_accelerations[object_name, 'x'] += calculate_acceleration(obj_grav_force_x, object_mass)
                    all_current_accelerations[object_name, 'y'] += calculate_acceleration(obj_grav_force_y, object_mass)
                    all_current_accelerations[object_name, 'z'] += calculate_acceleration(obj_grav_force_z, object_mass)


        
        #Now that we have the total acceleration based on every object, calculate the position
        next_object_position_x = calculate_new_position(object_x_m, object_vx, all_current_accelerations[object_name, 'x'], time_interval)
        next_object_position_y = calculate_new_position(object_y_m, object_vy, all_current_accelerations[object_name, 'y'], time_interval)
        next_object_position_z = calculate_new_position(object_z_m, object_vz, all_current_accelerations[object_name, 'z'], time_interval)
        next_object_distance = calculate_distance_from_center(next_object_position_x, next_object_position_y, next_object_position_z) #m

        all_next_distances_each_direction[object_name, "x"] = next_object_position_x
        all_next_distances_each_direction[object_name, "y"] = next_object_position_y
        all_next_distances_each_direction[object_name, "z"] = next_object_position_z

        next_distances_from_center[object_name] = next_object_distance



    # Now find next velocity
    for i in range(len(all_objects)):
        
        object_name = all_objects[i]['name']
        object_mass = all_objects[i]['mass']
        object_x = all_objects[i]['x']
        object_y = all_objects[i]['y']
        object_z = all_objects[i]['z']
        object_vx = all_objects[i]["velocity_x"]
        object_vy = all_objects[i]["velocity_y"]
        object_vz = all_objects[i]["velocity_z"]

        next_object_x = all_next_distances_each_direction[object_name, "x"]
        next_object_y = all_next_distances_each_direction[object_name, "y"]
        next_object_z = all_next_distances_each_direction[object_name, "z"]


        print(f"Starting object {object_name}")

        object_ax = all_current_accelerations[(object_name, 'x')]
        object_ay = all_current_accelerations[(object_name, 'y')]
        object_az = all_current_accelerations[(object_name, 'z')]

        object_x_m = au_to_meter(object_x)
        object_y_m = au_to_meter(object_y)
        object_z_m = au_to_meter(object_z)

        next_object_distance = next_distances_from_center[object_name]

        for j in range(len(all_objects)):
            if j != i:

                second_object_name = all_objects[j]['name']
                second_object_mass = all_objects[j]['mass']
                second_object_x = all_objects[j]['x']
                second_object_y = all_objects[j]['y']
                second_object_z = all_objects[j]['z']

                
                next_second_object_x = all_next_distances_each_direction[second_object_name, "x"]
                next_second_object_y = all_next_distances_each_direction[second_object_name, "y"]
                next_second_object_z = all_next_distances_each_direction[second_object_name, "z"]

                second_object_x_m = au_to_meter(second_object_x)
                second_object_y_m = au_to_meter(second_object_y)
                second_object_z_m = au_to_meter(second_object_z)

                object_distances = calculate_distance_between_objects(next_object_x, next_object_y, next_object_z, 
                                                                      next_second_object_x, next_second_object_y, next_second_object_z)
                
                print(f"Calculating the next gravity between object {object_name} and {second_object_name}")
                second_object_distance = distances_from_center[second_object_name]
                if (next_grav_forces.get((object_name, second_object_name, "x"), None) is None or
                    next_grav_forces.get((object_name, second_object_name, "y")) is None or 
                    next_grav_forces.get((object_name, second_object_name, "z")) is None): 

                    obj_grav_force_x = calculate_grav_force(object_mass, second_object_mass, object_distances, 
                                                             next_object_x, next_second_object_x)
                    obj_grav_force_y = calculate_grav_force(object_mass, second_object_mass, object_distance, 
                                                             next_object_y, next_second_object_y)
                    obj_grav_force_z = calculate_grav_force(object_mass, second_object_mass, object_distance, 
                                                             next_object_z, next_second_object_z)
                    second_obj_grav_force_x = -obj_grav_force_x
                    second_obj_grav_force_y = -obj_grav_force_y
                    second_obj_grav_force_z = -obj_grav_force_z

                    next_grav_forces[object_name, second_object_name, "x"] = obj_grav_force_x
                    next_grav_forces[object_name, second_object_name, "y"] = obj_grav_force_y
                    next_grav_forces[object_name, second_object_name, "z"] = obj_grav_force_z

                    next_grav_forces[second_object_name, object_name, 'x'] = second_obj_grav_force_x
                    next_grav_forces[second_object_name, object_name, 'y'] = second_obj_grav_force_y
                    next_grav_forces[second_object_name, object_name, 'z'] = second_obj_grav_force_z

                    all_next_accelerations[object_name, 'x'] += calculate_acceleration(obj_grav_force_x, object_mass)
                    all_next_accelerations[object_name, 'y'] += calculate_acceleration(obj_grav_force_y, object_mass)
                    all_next_accelerations[object_name, 'z'] += calculate_acceleration(obj_grav_force_z, object_mass)



        next_object_acceleration_x = all_next_accelerations[object_name, 'x']
        next_object_acceleration_y = all_next_accelerations[object_name, 'y']
        next_object_acceleration_z = all_next_accelerations[object_name, 'z']

        next_object_velocity_x = calculate_new_velocity(object_vx, object_ax, next_object_acceleration_x, time_interval)
        next_object_velocity_y = calculate_new_velocity(object_vy, object_ay, next_object_acceleration_y, time_interval)
        next_object_velocity_z = calculate_new_velocity(object_vz, object_az, next_object_acceleration_z, time_interval)
        print(f"Found next {object_name} velocity from original position\n ({object_vx}, {object_vy}, {object_vz}) --> ({next_object_velocity_x}, {next_object_velocity_y}, {next_object_velocity_z}))")

        print(f"{object_name} Position all in AU: ({object_x}, {object_y}, {object_z}) --> ({meter_to_au(next_object_position_x)}, {meter_to_au(next_object_position_y)}, {meter_to_au(next_object_position_z)}))\n")
    

        all_objects[i]['x'] = meter_to_au(next_object_position_x)
        all_objects[i]['y'] = meter_to_au(next_object_position_y)
        all_objects[i]['z'] = meter_to_au(next_object_position_z)
        all_objects[i]['velocity_x'] = next_object_velocity_x
        all_objects[i]['velocity_y'] = next_object_velocity_y
        all_objects[i]['velocity_z'] = next_object_velocity_z

        print(f"{all_objects[i]['x']}, {all_objects[i]['y']}, {all_objects[i]['z']}, {all_objects[i]['velocity_x']}, {all_objects[i]['velocity_y']}, {all_objects[i]['velocity_z']}")
    
    return


fig = plt.figure()

ax = fig.add_subplot(projection="3d")

ax.view_init(elev=90, azim=-90)
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_zlim(-10, 10)

traj_lines = []
points = []

for obj in all_objects:
    obj["x_list"] = []
    obj["y_list"] = []
    obj["z_list"] = []

    traj, = ax.plot([], [], [], '-', color=obj["color"], label=obj["name"])
    point, = ax.plot([], [], [], 'o', color=obj["color"])

    traj_lines.append(traj)
    points.append(point)


ax.legend()


def update(frame):
    global frame_count
    print(f"Frame Count")
    run_sim()
    print(f"Simulation Complete for frame: {frame_count}\n\n\n")
    frame_count += 1
    for i, obj in enumerate(all_objects):
        x = obj["x"]
        y = obj["y"]
        z = obj["z"]

        obj["x_list"].append(x)
        obj["y_list"].append(y)
        obj["z_list"].append(z)

        traj_lines[i].set_data(obj["x_list"], obj["y_list"])
        traj_lines[i].set_3d_properties(obj["z_list"])

        points[i].set_data([x], [y])
        points[i].set_3d_properties([z])

    return (*traj_lines, *points)

anim = FuncAnimation(fig, update, frames=2, interval=500)
plt.show()
