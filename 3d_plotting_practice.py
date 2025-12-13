import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits import mplot3d



ax = plt.axes(projection="3d") #type of coordinate system

# Single Points
'''
ax = plt.axes(projection="3d") #type of coordinate system
ax.scatter(3, 5, 7) # Scatterplot

'''

# More Points
"""
x_data = np.random.randint(0, 100, (500,))
y_data = np.random.randint(0, 100, (500,))
z_data = np.random.randint(0, 100, (500,))

ax.scatter(x_data, y_data, z_data)
"""


# Plotting a function
"""
x_data = np.arange(0, 50, 0.1)
y_data = np.arange(0, 50, 0.1)
z_data = x_data * y_data
ax.plot(x_data, y_data, z_data)
"""

# Setting a viewing direction
"""
Trying to view from top down for example, use the parameters when dragging the graph that are in the top rihgt
ax.view_init(azim=0, elev=90)
"""







plt.show()





"""
In order to do a .scatter() command, there are other arguments
- Marker = "" is the shape of the things
- Alpha = number between 0-1 is the opacity

Use ax.plot to make a line, ax.scatter to make points

Use ax.set_title("")  To make a title
Use ax.set_xlabel(""), ax.set_ylabel(""), ax.set_zlabel("") to label axis 
Do this all before plt.show()
"""