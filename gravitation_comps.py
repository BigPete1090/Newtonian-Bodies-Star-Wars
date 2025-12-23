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
initial_energy = 0
showing_energy = 0
current_energy = []

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



frame_count = 0
time_interval = 60*60

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
    r_k1 = r_k + v_k * t + 0.5 * a_k * (t**2)
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

def calculate_kinetic_energy(m, v):
    return 0.5 * m * (v**2)

def calculate_kinetic_energy(m, vx, vy, vz):
    v = np.sqrt((vx**2) + (vy**2) + (vz**2))
    return 0.5 * m * v

def calculate_potential_energy(m1, m2, r):
    u = -1*((GRAVITATIONAL_CONSTANT * m1 * m2)/r)
    return u

def calculate_mechanical_energy(u, k):
    return k + u

"""def calculate_center_of_mass(mass_list, radius_list):
    #Assumes mass and radius are in ordered list
    if len(mass_list) != len(radius_list):
        raise ValueError("Lists are not the same length")
    sum_mass = 0
    sum_r_m = 0
    for i in len(mass_list):
        continue
    return mass_list"""


# Before running sim need to do the calculate center of mass of the whole area and shift everything


def run_sim():
    global star_x, star_y, star_z
    global planet_1_x, planet_1_y, planet_1_z
    global star_mass, planet_1_mass
    global star_velocity_x, star_velocity_y, star_velocity_z
    global planet_1_velocity_x, planet_1_velocity_y, planet_1_velocity_z
    global all_objects
    global initial_energy, current_energy, showing_energy


    k = 0
    u = 0
    distances_from_center = {}
    grav_forces = {}
    next_distances_from_center = {}
    all_current_accelerations = defaultdict(float)
    all_next_distances_each_direction = defaultdict(float)
    new_positions = {}
    new_velocities = {}

    all_next_accelerations = defaultdict(float)
    next_grav_forces = {}
    # Calculate Kinetic Energy
    for i in range(len(all_objects)):
        object_mass = all_objects[i]['mass']
        object_vx = all_objects[i]["velocity_x"]
        object_vy = all_objects[i]["velocity_y"]
        object_vz = all_objects[i]["velocity_z"]
        k += calculate_kinetic_energy(object_mass, object_vx, object_vy, object_vz)

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
        print(f"\nStarting Object: {object_name}")

        

        object_x_m = au_to_meter(object_x)
        object_y_m = au_to_meter(object_y)
        object_z_m = au_to_meter(object_z)

        object_distance = distances_from_center[object_name]
        
        for j in range(len(all_objects)):
            if j > i:
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
                
                u += calculate_potential_energy(object_mass, second_object_mass, object_distances)


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

                    all_current_accelerations[second_object_name, 'x'] += calculate_acceleration(second_obj_grav_force_x, second_object_mass)
                    all_current_accelerations[second_object_name, 'y'] += calculate_acceleration(second_obj_grav_force_y, second_object_mass)
                    all_current_accelerations[second_object_name, 'z'] += calculate_acceleration(second_obj_grav_force_z, second_object_mass)

        
        #Now that we have the total acceleration based on every object, calculate the position
        next_object_position_x = calculate_new_position(object_x_m, object_vx, all_current_accelerations[object_name, 'x'], time_interval)
        next_object_position_y = calculate_new_position(object_y_m, object_vy, all_current_accelerations[object_name, 'y'], time_interval)
        next_object_position_z = calculate_new_position(object_z_m, object_vz, all_current_accelerations[object_name, 'z'], time_interval)

        """print(f"DEBUG: {object_name} calculated next position (meters): x={next_object_position_x}, y={next_object_position_y}, z={next_object_position_z}")
        print(f"DEBUG: {object_name} original position (meters): x={object_x_m}, y={object_y_m}, z={object_z_m}")
        print(f"DEBUG: {object_name} velocity: vx={object_vx}, vy={object_vy}, vz={object_vz}")
        print(f"DEBUG: {object_name} acceleration: ax={all_current_accelerations[object_name, 'x']}, ay={all_current_accelerations[object_name, 'y']}, az={all_current_accelerations[object_name, 'z']}")
        """

        all_next_distances_each_direction[object_name, "x"] = next_object_position_x
        all_next_distances_each_direction[object_name, "y"] = next_object_position_y
        all_next_distances_each_direction[object_name, "z"] = next_object_position_z




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


        object_ax = all_current_accelerations[(object_name, 'x')]
        object_ay = all_current_accelerations[(object_name, 'y')]
        object_az = all_current_accelerations[(object_name, 'z')]

        for j in range(len(all_objects)):
            if j > i:

                second_object_name = all_objects[j]['name']
                second_object_mass = all_objects[j]['mass']

                
                next_second_object_x = all_next_distances_each_direction[second_object_name, "x"]
                next_second_object_y = all_next_distances_each_direction[second_object_name, "y"]
                next_second_object_z = all_next_distances_each_direction[second_object_name, "z"]



                object_distances = calculate_distance_between_objects(next_object_x, next_object_y, next_object_z, 
                                                                      next_second_object_x, next_second_object_y, next_second_object_z)
                
                print(f"Calculating the next gravity between object {object_name} and {second_object_name}")
                if (next_grav_forces.get((object_name, second_object_name, "x"), None) is None or
                    next_grav_forces.get((object_name, second_object_name, "y")) is None or 
                    next_grav_forces.get((object_name, second_object_name, "z")) is None): 

                    obj_grav_force_x = calculate_grav_force(object_mass, second_object_mass, object_distances, 
                                                             next_object_x, next_second_object_x)
                    obj_grav_force_y = calculate_grav_force(object_mass, second_object_mass, object_distances, 
                                                             next_object_y, next_second_object_y)
                    obj_grav_force_z = calculate_grav_force(object_mass, second_object_mass, object_distances, 
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

                    all_next_accelerations[second_object_name, 'x'] += calculate_acceleration(second_obj_grav_force_x, second_object_mass)
                    all_next_accelerations[second_object_name, 'y'] += calculate_acceleration(second_obj_grav_force_y, second_object_mass)
                    all_next_accelerations[second_object_name, 'z'] += calculate_acceleration(second_obj_grav_force_z, second_object_mass)



        next_object_acceleration_x = all_next_accelerations[object_name, 'x']
        next_object_acceleration_y = all_next_accelerations[object_name, 'y']
        next_object_acceleration_z = all_next_accelerations[object_name, 'z']

        

        next_object_velocity_x = calculate_new_velocity(object_vx, object_ax, next_object_acceleration_x, time_interval)
        next_object_velocity_y = calculate_new_velocity(object_vy, object_ay, next_object_acceleration_y, time_interval)
        next_object_velocity_z = calculate_new_velocity(object_vz, object_az, next_object_acceleration_z, time_interval)
        print(f"Found next {object_name} velocity from original position\n ({object_vx:.3f}, {object_vy:.3f}, {object_vz:.3f}) --> ({next_object_velocity_x:.3f}, {next_object_velocity_y:.3f}, {next_object_velocity_z:.3f}))")

        print(f"{object_name} Position all in AU: ({object_x:.3f}, {object_y:.3f}, {object_z:.3f}) --> ({meter_to_au(next_object_x):.3f}, {meter_to_au(next_object_y):.3f}, {meter_to_au(next_object_z):.3f}))")

        new_positions[object_name] = (next_object_x, next_object_y, next_object_z)
        new_velocities[object_name] = (next_object_velocity_x, next_object_velocity_y, next_object_velocity_z)

        """all_objects[i]['x'] = meter_to_au(next_object_position_x)
        all_objects[i]['y'] = meter_to_au(next_object_position_y)
        all_objects[i]['z'] = meter_to_au(next_object_position_z)
        all_objects[i]['velocity_x'] = next_object_velocity_x
        all_objects[i]['velocity_y'] = next_object_velocity_y
        all_objects[i]['velocity_z'] = next_object_velocity_z"""    

    for i in range(len(all_objects)):
        object_name = all_objects[i]['name']

        x_m, y_m, z_m = new_positions[object_name]
        all_objects[i]['x'] = meter_to_au(x_m)
        all_objects[i]['y'] = meter_to_au(y_m)
        all_objects[i]['z'] = meter_to_au(z_m)

        all_objects[i]['velocity_x'], all_objects[i]['velocity_y'], all_objects[i]['velocity_z'] = new_velocities[object_name]

        # find mechanical energy
        mechanical_energy = calculate_mechanical_energy(u, k)
        if frame_count == 0:
            initial_energy = mechanical_energy
        if frame_count // 100 == 0:
            showing_energy == mechanical_energy
        current_energy.append(mechanical_energy)


    return


fig = plt.figure(figsize = (10, 5))


# The arguments 122 mean 1 row, 2 columns, 2nd plot
# Crucially, specify the 'projection' argument
ax1 = fig.add_subplot(1, 2, 1, projection='3d')






ax1.view_init(elev=0, azim=-90)
ax1.set_xlim(-1, 1)
ax1.set_ylim(-3, 3)
ax1.set_zlim(-1, 1)

ax1.set_xlabel('X (AU)')
ax1.set_ylabel('Y (AU)')
ax1.set_zlabel('Z (AU)')

traj_lines = []
points = []

for obj in all_objects:
    obj["x_list"] = []
    obj["y_list"] = []
    obj["z_list"] = []

    traj, = ax1.plot([], [], [], '-', color=obj["color"], label=obj["name"])
    point, = ax1.plot([], [], [], 'o', color=obj["color"])

    traj_lines.append(traj)
    points.append(point)


ax1.legend()

ax2 = fig.add_subplot(1, 2, 2)
x_data_energy = []
y_data_energy = []
line, = ax2.plot(x_data_energy, y_data_energy)

ax2.set_title('Energy vs Frames')
ax2.set_xlabel('Frame count aka days')
ax2.set_ylabel('Mechanical Energy (J)')
ax2.grid(True)


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

    new_x = frame
    new_y = current_energy[-1]
    x_data_energy.append(new_x)
    y_data_energy.append(new_y)
    line.set_data(x_data_energy, y_data_energy)
    ax2.relim()
    ax2.autoscale_view()

    if new_x >= ax2.get_xlim()[1]:
        ax2.set_xlim(new_x - 50, new_x + 50) # Shift the window
        ax2.figure.canvas.draw() # Force redraw of axes
    

    return (*traj_lines, *points, line,)

anim = FuncAnimation(fig, update, frames=30000, interval=0, repeat=False)
plt.tight_layout()
plt.show()
