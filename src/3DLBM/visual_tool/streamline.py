import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
import argparse

def load_vtk(filename):
    """
    Read the VTK file and return the Structured Grid.
    """
    grid = pv.read(filename)
    return grid

def extract_slice(grid, normal='z'):
    """
    Extract a slice from the grid. By default, the median plane along the z-axis is selected.
    """
    dims = grid.dimensions  # (nx, ny, nz)
    spacing = grid.spacing  # (dx, dy, dz)
    if normal.lower() == 'z':
        z_mid = (dims[2] - 1) / 2.0 * spacing[2]
        slice_plane = grid.slice(normal='z', origin=(0, 0, z_mid))
    elif normal.lower() == 'y':
        y_mid = (dims[1] - 1) / 2.0 * spacing[1]
        slice_plane = grid.slice(normal='y', origin=(0, y_mid, 0))
    elif normal.lower() == 'x':
        x_mid = (dims[0] - 1) / 2.0 * spacing[0]
        slice_plane = grid.slice(normal='x', origin=(x_mid, 0, 0))
    else:
        raise ValueError("Slice direction not supported. Choose from 'x', 'y', or 'z'.")
    return slice_plane

def create_streamline_plot(slice_plane):
    """
    Prepare the streamline plot from the obtained slice.
    Assumes the slice is structured so it can be remapped into a 2D grid for Matplotlib's streamplot.
    """
    # Extract point coordinates and the velocity field
    points = slice_plane.points
    velocity = slice_plane.point_data['velocity']
    x = points[:, 0]
    y = points[:, 1]
    u = velocity[:, 0]
    v = velocity[:, 1]

    unique_x = np.unique(x)
    unique_y = np.unique(y)
    nx = unique_x.size
    ny = unique_y.size

    # Sort points to remap to grid:
    # We create a 2D grid: the reshape function relies on row-major ordering along y, then x.
    try:
        X = x.reshape(ny, nx)
        Y = y.reshape(ny, nx)
        U = u.reshape(ny, nx)
        V = v.reshape(ny, nx)
    except ValueError:
        # If reshape fails, build the grid manually.
        X, Y = np.meshgrid(unique_x, unique_y)
        from scipy.interpolate import griddata
        U = griddata((x, y), u, (X, Y), method='linear')
        V = griddata((x, y), v, (X, Y), method='linear')

    # Creation of the streamline plot
    plt.figure(figsize=(8, 6))
    strm = plt.streamplot(X, Y, U, V, density=3, color='black', linewidth=1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.grid(True)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Display streamlines from a VTK file (Driven Cavity 3D).")
    parser.add_argument("filename", type=str, help="Path to the VTK file to display.")
    parser.add_argument("--normal", type=str, default='z',
                        help="Slice direction ('x', 'y', or 'z'). Default: 'z'")
    args = parser.parse_args()

    # Load the VTK file
    grid = load_vtk(args.filename)
    if 'velocity' not in grid.point_data.keys():
        raise KeyError("The VTK file does not contain the 'velocity' field. Check the input file.")

    # Extract the slice
    slice_plane = extract_slice(grid, normal=args.normal)

    # Create and display the streamline plot
    create_streamline_plot(slice_plane)

if __name__ == '__main__':
    main()
