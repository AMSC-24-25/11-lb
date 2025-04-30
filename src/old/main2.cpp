#include "Lattice.cpp"
#include "Obstacle.cpp"
#include <fstream>
#include <iostream>

int main() {
    // 1. Inizializzazione del reticolo
    Lattice lattice;
    
    // Configurazione parametri
    lattice.nx = 200;
    lattice.ny = 100;
    lattice.x_min = 0.0;
    lattice.x_max = 4.0;
    lattice.y_min = 0.0;
    lattice.y_max = 1.0;
    lattice.tau_lbm = 0.8;
    lattice.u_lbm = 0.1;
    lattice.IBB = true;
    
    // Re-inizializza le strutture dati con i nuovi parametri
    lattice.set_default_lbm();

    // 2. Creazione ostacolo quadrato
    MatrixXd square(4, 2);
    double cx = 2.0, cy = 0.5;  // Centro del dominio in x, centro in y
    double side = 0.2;           // Lato del quadrato
    
    square << cx - side/2, cy - side/2,
              cx + side/2, cy - side/2,
              cx + side/2, cy + side/2,
              cx - side/2, cy + side/2;
    
    Obstacle obstacle("1", square);  // Tag "1"
    lattice.add_obstacle(obstacle);

    // 3. Configurazione condizioni al contorno
    // Velocità in ingresso (parete sinistra)
    lattice.u_left[0].setConstant(lattice.u_lbm);  // Componente x
    lattice.u_left[1].setZero();                  // Componente y

    // Condizioni di uscita (pressione costante a destra)
    lattice.rho_right.setConstant(1.0);

    // 4. Loop temporale
    int it_max = 10000;
    int output_interval = 1000;

    for(int it = 0; it < it_max; ++it) {
        // Calcola campi macroscopici
        lattice.macro();
        
        // Calcola equilibrio
        lattice.equilibrium(lattice.u, lattice.c, lattice.w, lattice.rho, lattice.g_eq);
        
        // Collisione e streaming
        lattice.collision_streaming(
            lattice.g, lattice.g_eq, lattice.g_up,
            lattice.om_p_lbm, lattice.om_m_lbm, lattice.ns,
            lattice.nx, lattice.ny, lattice.lx, lattice.ly
        );
        
        // Applica condizioni al contorno
        // Parete sinistra (Zou-He velocity)
        lattice.zou_he_left_wall_velocity(
            lattice.lx, lattice.ly,
            lattice.u, lattice.u_left,
            lattice.rho, lattice.g
        );
        
        // Parete destra (Zou-He pressure)
        lattice.zou_he_right_wall_pressure(
            lattice.lx, lattice.ly,
            lattice.u, lattice.rho_right, lattice.u_right,
            lattice.rho, lattice.g
        );
        
        // Pareti superiore e inferiore (no-slip)
        lattice.zou_he_top_wall_velocity(
            lattice.lx, lattice.ly,
            lattice.u, lattice.u_top,
            lattice.rho, lattice.g
        );
        lattice.zou_he_bottom_wall_velocity(
            lattice.lx, lattice.ly,
            lattice.u, lattice.u_bot,
            lattice.rho, lattice.g
        );
        
        // Gestione angoli
        lattice.zou_he_bottom_left_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
        lattice.zou_he_top_left_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
        lattice.zou_he_top_right_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
        lattice.zou_he_bottom_right_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);

        // Output progressivo
        if(it % output_interval == 0) {
            std::cout << "Iterazione " << it << "/" << it_max << std::endl;
        }
    }

    // 5. Salvataggio dati
    std::ofstream file("velocity_field.txt");
    file << "x,y,ux,uy\n";
    
    for(int i = 0; i < lattice.nx; ++i) {
        for(int j = 0; j < lattice.ny; ++j) {
            Vector2d coords = lattice.get_coords(i, j);
            double ux = lattice.u[0](i, j);
            double uy = lattice.u[1](i, j);
            
            file << sqrt((ux*ux + uy*uy)) << "\n";
        }
    }
    
    file.close();
    std::cout << "Simulazione completata. Dati salvati in velocity_field.txt" << std::endl;

    return 0;
}