import pandas as pd
import numpy as np
import pyvista as pv

def write_error_vector_to_vtk(df, filename):
    print("Writing: " + filename )
    x = np.unique(df['X'])
    y = np.unique(df['Y'])
    z = np.array([0.0])

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    height = df['Height_error'].values.reshape(len(y), len(x))
    velocity = np.column_stack((df['X_velocity_error'], df['Y_velocity_error'], np.zeros(len(df))))
    velocity_x = velocity[:, 0].reshape(len(y), len(x))
    velocity_y = velocity[:, 1].reshape(len(y), len(x))
    velocity_z = velocity[:, 2].reshape(len(y), len(x))

    # Create a UniformGrid from the DataFrame
    grid = pv.UniformGrid()
    grid.dimensions = (len(x), len(y), 1)
    grid.origin = (x[0], y[0], z[0])
    grid.spacing = (dx, dy, 1)

    # Add Height and Velocity data to the grid
    grid.point_data["Height_error"] = height.ravel(order='F')
    grid.point_data["Velocity_error_X"] = velocity_x.ravel(order='F')
    grid.point_data["Velocity_error_Y"] = velocity_y.ravel(order='F')
    grid.point_data["Velocity_error_Z"] = velocity_z.ravel(order='F')

    # Save the grid to a VTK file
    grid.save(filename)
    
def write_vector_to_vtk(df, filename):
    print("Writing: " + filename )
    x = np.unique(df['X'])
    y = np.unique(df['Y'])
    z = np.array([0.0])

    dx = x[1] - x[0]
    dy = y[1] - y[0]

    height = df['Height'].values.reshape(len(y), len(x))
    velocity = np.column_stack((df['X_velocity'], df['Y_velocity'], np.zeros(len(df))))
    velocity_x = velocity[:, 0].reshape(len(y), len(x))
    velocity_y = velocity[:, 1].reshape(len(y), len(x))
    velocity_z = velocity[:, 2].reshape(len(y), len(x))

    # Create a UniformGrid from the DataFrame
    grid = pv.UniformGrid()
    grid.dimensions = (len(x), len(y), 1)
    grid.origin = (x[0], y[0], z[0])
    grid.spacing = (dx, dy, 1)

    # Add Height and Velocity data to the grid
    grid.point_data["Height"] = height.ravel(order='F')
    grid.point_data["Velocity_X"] = velocity_x.ravel(order='F')
    grid.point_data["Velocity_Y"] = velocity_y.ravel(order='F')
    grid.point_data["Velocity_Z"] = velocity_z.ravel(order='F')

    # Save the grid to a VTK file
    grid.save(filename)