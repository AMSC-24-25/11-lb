#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <Eigen/Dense>
#include "Lattice.cpp"
#include "Obstacle.cpp"

using namespace Eigen;

class StepApp {
public:
    // Parametri
    std::string name = "step";
    double Re_lbm = 500.0;
    double L_lbm = 150;
    double u_lbm = 0.05;
    double rho_lbm = 1.0;
    double t_max = 15.0;
    double x_min = -1.0, x_max = 15.0;
    double y_min = -1.0, y_max = 1.0;
    bool IBB = false;
    std::string stop = "it";
    double obs_cv_ct = 1.0e-3;
    int obs_cv_nb = 1000;
    
    // Parametri di output
    int output_freq = 500;
    int output_it = 0;
    int dpi = 200;
    
    // Parametri calcolati
    double Cs, u_avg, r_cyl, D_lbm, nu_lbm, tau_lbm, dt, dx, dy;
    int nx, ny, it_max, sigma;
    
    // Ostacoli
    std::vector<Obstacle> obstacles;
    Lattice lattice;

    StepApp() {
        compute_lbm_parameters();
        initialize_lattice();
        add_obstacle();
    }

    void compute_lbm_parameters() {
        Cs = 1.0/sqrt(3.0);
        ny = L_lbm;
        u_avg = 2.0*u_lbm/3.0;
        r_cyl = 0.5;
        D_lbm = floor(ny*r_cyl/(y_max - y_min));
        nu_lbm = u_avg*L_lbm/Re_lbm;
        tau_lbm = 0.5 + nu_lbm/(Cs*Cs);
        dt = Re_lbm*nu_lbm/(L_lbm*L_lbm);
        dx = (y_max - y_min)/ny;
        dy = dx;
        nx = floor(ny*(x_max - x_min)/(y_max - y_min));
        it_max = floor(t_max/dt);
        sigma = floor(10*nx);
    }

    void initialize_lattice() {
        lattice.nx = nx;
        lattice.ny = ny;
        lattice.x_min = x_min;
        lattice.x_max = x_max;
        lattice.y_min = y_min;
        lattice.y_max = y_max;
        lattice.tau_lbm = tau_lbm;
        lattice.u_lbm = u_lbm;
        lattice.IBB = IBB;
        lattice.set_default_lbm();
    }

    void add_obstacle() {

    std::string name = "square1";
    int n_pts = 4;
    int n_spts = 100;
    std::string type = "square";
    Eigen::Vector2i size(1, 1);  // Usa un vettore di dimensioni float
    Eigen::Vector2i pos(1, 1);   // Usa un vettore di posizione float

    // Creazione dell'oggetto Obstacle
    Obstacle square1(name, n_pts, n_spts, type, size, pos);
    }

    void set_inlets(int it) {
        double val = it;
        double ret = (1.0 - exp(-val*val/(2.0*sigma*sigma)));
        
        for(int j=0; j<ny; j++) {
            Vector2d pt = lattice.get_coords(0, j);
            Vector2d u = poiseuille(pt);
            lattice.u_left[0](0, j) = ret * u_lbm * u[0];
            lattice.u_left[1](0, j) = 0.0;
        }
        
        lattice.u_top[0].setZero();
        lattice.u_bot[0].setZero();
        lattice.u_right[1].setZero();
        lattice.rho_right.setConstant(rho_lbm);
    }

    Vector2d poiseuille(const Vector2d& pt) {
        double H = y_max - y_min;
        double y = pt[1];
        Vector2d u;
        u[0] = 4.0 * (y_max - y) * (y - y_min) / (H*H);
        u[1] = 0.0;
        return u;
    }

    void run_simulation() {
        for(int it=0; it<it_max; it++) {
            set_inlets(it);
            lattice.macro();
            lattice.equilibrium(lattice.u, lattice.c, lattice.w, lattice.rho, lattice.g_eq);
            lattice.collision_streaming(lattice.g, lattice.g_eq, lattice.g_up, 
                                      lattice.om_p_lbm, lattice.om_m_lbm, lattice.ns, 
                                      lattice.nx, lattice.ny, lattice.lx, lattice.ly);
            set_boundary_conditions();
            
            if(it % output_freq == 0) {
                save_velocity_field(it);
                output_it++;
            }
        }
    }

    void set_boundary_conditions() {
        // Bounce-back per ostacoli
        for(auto& obs : obstacles) {
            lattice.bounce_back_obstacle(lattice.IBB, obs.boundary, lattice.ns, lattice.c, 
                                            obs.ibb, lattice.g_up, lattice.g);
        }
        
        // Condizioni al contorno
        lattice.zou_he_bottom_wall_velocity(lattice.lx, lattice.ly, lattice.u, lattice.u_bot, lattice.rho, lattice.g);
        lattice.zou_he_left_wall_velocity(lattice.lx, lattice.ly, lattice.u, lattice.u_left, lattice.rho, lattice.g);
        lattice.zou_he_top_wall_velocity(lattice.lx, lattice.ly, lattice.u, lattice.u_top, lattice.rho, lattice.g);
        lattice.zou_he_right_wall_pressure(lattice.lx, lattice.ly, lattice.u, lattice.rho_right, lattice.u_right, lattice.rho, lattice.g);
        
        // Angoli
        lattice.zou_he_bottom_left_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
        lattice.zou_he_top_left_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
        lattice.zou_he_top_right_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
        lattice.zou_he_bottom_right_corner_velocity(lattice.lx, lattice.ly, lattice.u, lattice.rho, lattice.g);
    }

    void save_velocity_field(int it) {
        std::ofstream file("velocity_" + std::to_string(it) + ".csv");
        file << "x,y,ux,uy\n";
        
        for(int i=0; i<lattice.nx; i++) {
            for(int j=0; j<lattice.ny; j++) {
                Vector2d pt = lattice.get_coords(i, j);
                file << sqrt(pow(lattice.u[0](i, j),2) + pow(lattice.u[1](i, j),2)) << "\n";
            }
        }
        file.close();
    }
};

int main() {
    StepApp app;
    app.run_simulation();
    return 0;
}