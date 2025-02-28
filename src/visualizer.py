import numpy as np
import matplotlib.pyplot as plt
import imageio
import sys

# input file name
input_file = 'vel_data.txt'


iteration_per_frame = 20

def read_data(file_name):
    with open(file_name, 'r') as f:
        # reading grid size
        nx = int(f.readline().strip())
        ny = int(f.readline().strip())

        # reading velocities
        num_elements = nx * ny
        data = []
        for line in f:
            data.append(float(line.strip()))

    return nx, ny, data

def create_frames(nx, ny, data, num_iterations):
    frames = []
    for iter in range(num_iterations):
        # legge dati per ogni iterazione e salva in una matrice (ny, nx)
        frame_data = np.array(data[iter * nx * ny:(iter + 1) * nx * ny]).reshape(ny, nx).transpose()
        frame_data = frame_data.T

        plt.imshow(frame_data, cmap='plasma', origin='lower')  # 'origin' è impostato su 'lower' per far partire y da 0 in basso
        plt.colorbar(label='Velocity Magnitude')
        plt.title(f'Iteration {(iter + 1)*iteration_per_frame}')

        # La visualizzazione ora avrà l'asse x da 0 a nx e y da 0 a ny
        plt.pause(0.001)  
        plt.clf()  # Rimuovi il frame precedente per preparare il successivo
        frame = np.frombuffer(plt.gcf().canvas.tostring_rgb(), dtype=np.uint8)
        frame = frame.reshape(plt.gcf().canvas.get_width_height()[::-1] + (3,))
        frames.append(frame)

        # Calcola la percentuale e aggiorna la stessa linea nel terminale
        percent_complete = (iter + 1) / num_iterations * 100
        print(f"Progress: {percent_complete:.2f}%", end='\r')

    print()  # Per assicurarsi che il prompt successivo inizi su una nuova linea
    return frames

def save_video(frames, output_file):
    imageio.mimsave(output_file, frames, fps=24)

if __name__ == '__main__':
    nx, ny, data = read_data(input_file)
    num_iterations = len(data) // (nx * ny)
    frames = create_frames(nx, ny, data, num_iterations)
    save_video(frames, 'lbm_simulation.mp4')

    print("Video generated: lbm_simulation.mp4")
