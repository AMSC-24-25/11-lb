import numpy as np
import matplotlib.pyplot as plt
import imageio
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

# input file name (modifica se il percorso cambia)
input_file = os.path.abspath(os.path.join(__file__, '..', '..', '..', 'output/vel_data.txt'))

iteration_per_frame = 30  # intervallo tra salvataggi nel file

def read_data(file_name):
    with open(file_name, 'r') as f:
        nx = int(f.readline().strip())
        ny = int(f.readline().strip())

        data = [float(line.strip()) for line in f]

    return nx, ny, data

def create_frames(nx, ny, data, num_iterations, vmax):
    frames = []
    for iter in range(num_iterations):
        frame_data = np.array(data[iter * nx * ny:(iter + 1) * nx * ny]).reshape(ny, nx)

        fig, ax = plt.subplots()
        im = ax.imshow(frame_data, cmap='RdBu_r', origin='lower', vmin=0, vmax=vmax)

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)  # puoi regolare size e pad
        plt.colorbar(im, cax=cax, label='Velocity Magnitude')
        ax.set_title(f'Iteration {iter * iteration_per_frame}')

        fig.canvas.draw()
        frame = np.frombuffer(fig.canvas.tostring_rgb(), dtype=np.uint8)
        frame = frame.reshape(fig.canvas.get_width_height()[::-1] + (3,))
        frames.append(frame)

        plt.close(fig)

        percent_complete = (iter + 1) / num_iterations * 100
        print(f"Progress: {percent_complete:.2f}%", end='\r')

    print()
    return frames

def save_gif(frames, output_file):
    imageio.mimsave(output_file, frames, fps=10)  # modifica fps se vuoi

if __name__ == '__main__':
    nx, ny, data = read_data(input_file)
    vmax = max(data)
    num_iterations = len(data) // (nx * ny)

    frames = create_frames(nx, ny, data, num_iterations, vmax)

    output_path = os.path.abspath(os.path.join(__file__, '..', '..', '..', 'output/lbm_simulation.gif'))
    save_gif(frames, output_path)

    print("GIF generated:", os.path.basename(output_path))
    print('\a')
