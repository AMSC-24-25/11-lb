import sys
import numpy as np
import matplotlib.pyplot as plt

def parse_vtk(filename):
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Estrai DIMENSIONS e SPACING
    nx = ny = nz = None
    dx = None
    for line in lines:
        if line.startswith('DIMENSIONS'):
            parts = line.split()
            nx, ny, nz = map(int, parts[1:4])
        elif line.startswith('SPACING'):
            parts = line.split()
            dx = float(parts[1])  # assumiamo passo uguale in tutte le direzioni
        if nx is not None and dx is not None:
            break

    if None in (nx, ny, nz, dx):
        raise ValueError("error.")

    nn = nx * ny * nz

    # Trova il primo LOOKUP_TABLE (fine di SCALARS rho)
    lt_idx = next(i for i, l in enumerate(lines) if l.startswith('LOOKUP_TABLE'))
    i = lt_idx + 1 + nn

    # Trova il blocco VECTORS
    vec_idx = next(j for j in range(i, len(lines)) if lines[j].startswith('VECTORS'))
    i = vec_idx + 1

    # Leggi i vettori di velocità
    vel_vals = []
    count = 0
    while count < nn * 3:
        vals = [float(x) for x in lines[i].split()]
        vel_vals.extend(vals)
        count += len(vals)
        i += 1

    # Rimappa in array 4D
    vel_array = np.array(vel_vals).reshape((nn, 3))
    vel = vel_array.reshape((nz, ny, nx, 3))
    return vel, dx


def plot_profiles(vel, dx, output_prefix, show=False):
    nz, ny, nx, _ = vel.shape
    mid_x = nx // 2
    mid_y = ny // 2
    mid_z = nz // 2

    # Coordinate spaziali
    y = np.arange(ny) * dx
    x = np.arange(nx) * dx

    # Profili centrali
    u_x_vert = vel[mid_z, :, mid_x, 0]
    u_y_horiz = vel[mid_z, mid_y, :, 1]

    def setup_ax():
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(True)
        return fig, ax

    # Profilo verticale
    fig, ax = setup_ax()
    ax.plot(y, u_x_vert, linewidth=0.8, color='black')
    ax.set_xlabel('y')
    ax.set_ylabel(r'$u_y$')
    plt.tight_layout()
    out1 = f'{output_prefix}_vertical_profile.png'
    plt.savefig(out1, dpi=300)
    if show: plt.show()
    plt.close(fig)

    # Profilo orizzontale
    fig, ax = setup_ax()
    ax.plot(x, u_y_horiz, linewidth=0.8, color='black')
    ax.set_xlabel('x')
    ax.set_ylabel(r'$u_x$')
    plt.tight_layout()
    out2 = f'{output_prefix}_horizontal_profile.png'
    plt.savefig(out2, dpi=300)
    if show: plt.show()
    plt.close(fig)

    print('Grafici salvati:', out1, out2)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python velocityprofile.py <file.vtk> [--show]")
        sys.exit(1)
    filename = sys.argv[1]
    show_flag = '--show' in sys.argv
    vel, dx = parse_vtk(filename)
    prefix = filename.rsplit('.', 1)[0]
    plot_profiles(vel, dx, prefix, show=show_flag)
