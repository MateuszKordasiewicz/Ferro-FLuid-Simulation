import numpy as np
import pandas as pd
from scipy.interpolate import griddata

def interpolate_to_common_grid(df, common_grid_shape):
    x_min, x_max = -0.2, 0.2
    y_min, y_max = -0.025, 0.025
    x_grid, y_grid = np.linspace(x_min, x_max, common_grid_shape[0]), np.linspace(y_min, y_max, common_grid_shape[1])
    x_grid, y_grid = np.meshgrid(x_grid, y_grid)

    # Transpose the x and y grids
    x_grid, y_grid = x_grid.T, y_grid.T

    interpolated_dfs = []

    CFL_lens = df['CFL'].unique()

    for CFL_len in CFL_lens:
        print("Interpolating: " + str(CFL_len))

        df_CFL = df[df['CFL'] == CFL_len]
        
        # Swap the order of x and y coordinates when creating the points array
        points = df_CFL[['X', 'Y']].values
        values_height = df_CFL['Height'].values
        values_x_velocity = df_CFL['X_velocity'].values
        values_y_velocity = df_CFL['Y_velocity'].values

        grid_method = "cubic"
        interpolation_method = "nearest"

        height_grid = griddata(points, values_height, (x_grid, y_grid), method=grid_method)
        height_grid_nan_mask = np.isnan(height_grid)
        height_grid[height_grid_nan_mask] = griddata(points, values_height, (x_grid[height_grid_nan_mask], y_grid[height_grid_nan_mask]), method=interpolation_method)

        x_velocity_grid = griddata(points, values_x_velocity, (x_grid, y_grid), method=grid_method)
        x_velocity_grid_nan_mask = np.isnan(x_velocity_grid)
        x_velocity_grid[x_velocity_grid_nan_mask] = griddata(points, values_x_velocity, (x_grid[x_velocity_grid_nan_mask], y_grid[x_velocity_grid_nan_mask]), method=interpolation_method)

        y_velocity_grid = griddata(points, values_y_velocity, (x_grid, y_grid), method=grid_method)
        y_velocity_grid_nan_mask = np.isnan(y_velocity_grid)
        y_velocity_grid[y_velocity_grid_nan_mask] = griddata(points, values_y_velocity, (x_grid[y_velocity_grid_nan_mask], y_grid[y_velocity_grid_nan_mask]), method=interpolation_method)

        df_interpolated = pd.DataFrame({'CFL': np.full(height_grid.size, CFL_len),
                                        'X': x_grid.ravel(),
                                        'Y': y_grid.ravel(),
                                        'Height': height_grid.ravel(),
                                        'X_velocity': x_velocity_grid.ravel(),
                                        'Y_velocity': y_velocity_grid.ravel()})
        interpolated_dfs.append(df_interpolated)

    df_interpolated_all = pd.concat(interpolated_dfs)
    return df_interpolated_all
