import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Define the radius of the hypersphere
R = 1

# Define the number of points (the more points, the smoother the sphere)
num_pts = 50

# Define the number of steps in the animation
num_steps = 30

# Create a figure and a 3D subplot
fig = plt.figure(figsize=(8, 8))
ax = fig.add_subplot(111, projection='3d')

# Initialization function: plot the background of each frame
def init():
    # Set the labels and the aspect ratio
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_xlim([-R, R])
    ax.set_ylim([-R, R])
    ax.set_zlim([-R, R])

    return fig,

# Animation function: this is called sequentially
def animate(i):
    ax.clear()
    ax.set_xlim([-R, R])
    ax.set_ylim([-R, R])
    ax.set_zlim([-R, R])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    # The fourth dimension (0 <= psi <= 2*pi)
    psi = 2 * np.pi * i / num_steps

    # 3D slice of the 4D hypersphere
    phi = np.linspace(0, np.pi, num_pts)  # polar angle
    theta = np.linspace(0, 2 * np.pi, num_pts)  # azimuthal angle
    phi, theta = np.meshgrid(phi, theta)

    # Spherical to Cartesian conversion
    x = R * np.sin(phi) * np.cos(theta)
    y = R * np.sin(phi) * np.sin(theta)
    z = R * np.cos(phi) * np.sin(psi)

    # Plot a 3D slice of the 4D hypersphere
    ax.plot_surface(x, y, z, color='b', rstride=1, cstride=1, alpha=0.6)

    return fig,

# Call the animator
anim = FuncAnimation(fig, animate, init_func=init, frames=num_steps, interval=100, blit=False)

# Show the animation
plt.show()
