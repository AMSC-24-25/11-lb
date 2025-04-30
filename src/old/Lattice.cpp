
#ifndef Lattice_CPP
#define Lattice_CPP

#include <Eigen/Dense>
#include <vector>
#include <array>
#include <string>
#include <cmath>
#include <iostream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include "Obstacle.cpp"

using namespace Eigen;

// Tipi
using Array2D = ArrayXXd;
using Array3D = std::vector<Array2D>;

class Lattice {
public:
    // Parametri principali
    std::string name = "Lattice";
    double x_min = 0.0, x_max = 1.0;
    double y_min = 0.0, y_max = 1.0;
    int nx = 100, ny = 100;
    double tau_lbm = 1.0;
    double dx = 1.0, dt = 1.0;
    double Cx = dx, Ct = dt, Cr = 1.0;
    double Cn = Cx * Cx / Ct;
    double Cu = Cx / Ct;
    double Cf = Cr * Cx * Cx / Ct;
    int dpi = 100;
    double u_lbm = 0.03;
    double L_lbm = 100;
    double nu_lbm = 0.01;
    double Re_lbm = 100.0;
    double rho_lbm = 1.0;
    bool IBB = false;
    std::string stop = "it";
    double t_max = 1.0;
    int it_max = 1000;
    double obs_cv_ct = 1e-1;
    int obs_cv_nb = 500;

    // Dir output
    std::string results_dir;
    std::string output_dir;
    std::string png_dir;

    // Campi LBM
    int output_it = 0;
    int lx, ly, q = 9;
    double Cs = 1.0 / std::sqrt(3.0);
    double tau_p_lbm, tau_m_lbm, lambda_trt = 0.25;
    double om_p_lbm, om_m_lbm, om_lbm;

    // D2Q9
    Matrix<double, 9, 2> c;
    Vector<double, 9> w;
    std::vector<int> ns;

    // Distribuzioni
    Array3D g, g_eq, g_up;

    // Campi fisici
    Array2D rho;
    std::array<Array2D, 2> u;

    // Condizioni al contorno
    std::array<Array2D, 2> u_left, u_right;
    std::array<Array2D, 2> u_top, u_bot;
    Array2D rho_right;
    Array2D lattice_mask;

    // Costruttore
    Lattice() {
        // Set directory output
        auto now = std::chrono::system_clock::now();
        std::time_t t = std::chrono::system_clock::to_time_t(now);
        std::stringstream ss;
        ss << std::put_time(std::localtime(&t), "%Y-%m-%d_%H_%M_%S");

        results_dir = "./results/";
        output_dir = results_dir + ss.str() + "/";
        png_dir = output_dir + "png/";

        // Parametri e strutture LBM
        set_default_lbm();

        // Info
        print_info();
    }

    // Funzione di inizializzazione LBM
    void set_default_lbm() {
            lx = nx - 1;
            ly = ny - 1;
            Cs = 1.0 / std::sqrt(3.0);

            tau_p_lbm = tau_lbm;
            tau_m_lbm = lambda_trt / (tau_p_lbm - 0.5) + 0.5;
            om_p_lbm = 1.0 / tau_p_lbm;
            om_m_lbm = 1.0 / tau_m_lbm;
            om_lbm   = 1.0 / tau_lbm;

            // D2Q9
            c <<  0,  0,
                1,  0,  -1,  0,
                0,  1,   0, -1,
                1,  1,  -1, -1,
                -1,  1,   1, -1;

            // Pesi
            w = Vector<double, 9>::Ones();
            for (int i = 0; i < 9; ++i) {
                double norm = c.row(i).norm();
                if (norm < 1.1 && i != 0) w(i) = 1.0 / 9.0;
                else if (norm > 1.1)      w(i) = 1.0 / 36.0;
            }
            w(0) = 4.0 / 9.0;

            // Indici bounce-back
            ns = {0, 2, 1, 4, 3, 6, 5, 8, 7};

            // Distribuzioni
            g.resize(9); g_eq.resize(9); g_up.resize(9);
            for (int i = 0; i < 9; ++i) {
                g[i]    = Array2D::Zero(nx, ny);
                g_eq[i] = Array2D::Zero(nx, ny);
                g_up[i] = Array2D::Zero(nx, ny);
            }

            // Campi
            rho = Array2D::Ones(nx, ny);
            for (auto& ui : u) ui = Array2D::Zero(nx, ny);

            for (auto& arr : u_left)  arr = Array2D::Zero(1, ny);
            for (auto& arr : u_right) arr = Array2D::Zero(1, ny);
            for (auto& arr : u_top)   arr = Array2D::Zero(1, nx);
            for (auto& arr : u_bot)   arr = Array2D::Zero(1, nx);

            rho_right = Array2D::Zero(1, ny);
            lattice_mask = Array2D::Zero(nx, ny);
    }

    void print_info() const {
            std::cout << "##################\n";
            std::cout << "### LBM solver ###\n";
            std::cout << "##################\n\n";

            std::cout << "# name       = " << name << "\n";
            std::cout << "# u_lbm      = " << u_lbm << "\n";
            std::cout << "# L_lbm      = " << L_lbm << "\n";
            std::cout << "# nu_lbm     = " << nu_lbm << "\n";
            std::cout << "# Re_lbm     = " << Re_lbm << "\n";
            std::cout << "# tau_p_lbm  = " << tau_p_lbm << "\n";
            std::cout << "# tau_m_lbm  = " << tau_m_lbm << "\n";
            std::cout << "# dt         = " << dt << "\n";
            std::cout << "# dx         = " << dx << "\n";
            std::cout << "# nx         = " << nx << "\n";
            std::cout << "# ny         = " << ny << "\n";
            std::cout << "# IBB        = " << (IBB ? "true" : "false") << "\n\n";
        }

        // Funzione per calcolare i campi macroscopici (densità e velocità)
    void macro() {
            // Calcola densità
            rho = Array2D::Zero(nx, ny);
            for (int q = 0; q < 9; ++q) {
                rho += g[q];  // Somma sui 9 canali per ottenere la densità
            }

            // Calcola velocità
            u[0] = Array2D::Zero(nx, ny);
            u[1] = Array2D::Zero(nx, ny);
            for (int q = 0; q < 9; ++q) {
                u[0] += c(q, 0) * g[q];  // Velocità lungo x
                u[1] += c(q, 1) * g[q];  // Velocità lungo y
            }

            // Normalizza per ottenere la velocità (diviso densità)
            u[0] /= rho;
            u[1] /= rho;
    }

        // Funzione per aggiungere un ostacolo
    void add_obstacle(const Obstacle& obstacle) {
        // Inizializza il tag dell'ostacolo
        int tag = stoi(obstacle.tag);
        std::cout << "# Obstacle " << tag << std::endl;

        // Calcola i confini del poligono
        Array4d poly_bnds = Array4d::Zero();
        poly_bnds(0) = obstacle.polygon.col(0).minCoeff();
        poly_bnds(1) = obstacle.polygon.col(0).maxCoeff();
        poly_bnds(2) = obstacle.polygon.col(1).minCoeff();
        poly_bnds(3) = obstacle.polygon.col(1).maxCoeff();

        // Dichiarazione degli array per l'ostacolo
        std::vector<Array2i> obs;
        std::vector<Array3i> bnd;
        std::vector<double> ibb;

        // Riempie la griglia
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                Array2i pt = get_coords(i, j);

                // Controlla se il punto è dentro il bbox del poligono
                if (pt(0) > poly_bnds(0) && pt(0) < poly_bnds(1) &&
                    pt(1) > poly_bnds(2) && pt(1) < poly_bnds(3)) {
                    if (is_inside(obstacle.polygon, pt)) {
                        lattice_mask(i, j) = tag;
                        obs.push_back(pt);
                    }
                }
            }
        }

        // Stampa il numero di posizioni nell'ostacolo
        std::cout << "# " << obs.size() << " locations in obstacle" << std::endl;

        // Costruisce il bordo dell'ostacolo (cioè il primo strato di fluido)
        for (size_t k = 0; k < obs.size(); ++k) {
            int i = obs[k](0) ;
            int j = obs[k](1) ;

            for (int q = 1; q < 9; ++q) {
                int qb = ns[q];
                int cx = c(q, 0);
                int cy = c(q, 1);
                int ii = i + cx;
                int jj = j + cy;

                if (ii >= nx || jj >= ny) continue;
                if (lattice_mask(ii, jj) == 0) {  // Se non è già un ostacolo
                    bnd.push_back(Array3i(ii, jj, qb));
                }
            }
        }

        // Rimuove i duplicati
        //std::sort(bnd.begin(), bnd.end());
        //bnd.erase(std::unique(bnd.begin(), bnd.end()), bnd.end());

        // Stampa il numero di posizioni sul bordo
        std::cout << "# " << bnd.size() << " locations on boundary" << std::endl;

        // Calcola la distanza al bordo se IBB è True
        if (IBB) {
            for (size_t k = 0; k < bnd.size(); ++k) {
                int i = bnd[k](0) ;
                int j = bnd[k](1) ;
                int q = bnd[k](2) ;
                Array2i pt = get_coords(i, j);
                ArrayXd x = obstacle.polygon.col(0).array() - pt(0);
                ArrayXd y = obstacle.polygon.col(1).array() - pt(1);
                ArrayXd dist = (x.square() + y.square()).sqrt();
                int mpt = dist.minCoeff();
                double mdst = dist(mpt) / (dx * c.row(q).norm());
                ibb.push_back(mdst);
            }
        }

        // Calcola l'area dell'ostacolo
        double area = 0.0;
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                if (lattice_mask(i, j) == tag) {
                    area += dx * dx;
                }
            }
        }

        // Stampa l'area dell'ostacolo
        std::cout << "# Area = " << area << std::endl;
        std::cout << std::endl;
    }

    Eigen::Vector2d get_coords(int i, int j) const {
    double dx = (x_max - x_min) / static_cast<double>(nx - 1);
    double dy = (y_max - y_min) / static_cast<double>(ny - 1);

    double x = x_min + i * dx;
    double y = y_min + j * dy;

    return Eigen::Vector2d(x, y);}

    // Funzione per determinare se un punto è dentro o fuori da un poligono chiuso
    bool is_inside(const Eigen::MatrixXd& poly, const Eigen::ArrayXd& pt) {
        int j = poly.rows() - 1;  // Ultimo punto del poligono
        bool odd_nodes = false;

        // Algoritmo per determinare se il punto è dentro un poligono non convesso
        for (int i = 0; i < poly.rows(); ++i) {
            if (((poly(i, 1) < pt(1) && poly(j, 1) >= pt(1)) || (poly(j, 1) < pt(1) && poly(i, 1) >= pt(1))) &&
                (poly(i, 0) < pt(0) || poly(j, 0) < pt(0))) {

                // Calcola la pendenza
                double slope = (poly(j, 0) - poly(i, 0)) / (poly(j, 1) - poly(i, 1));

                // Verifica il lato
                if ((poly(i, 0) + (pt(1) - poly(i, 1)) * slope) < pt(0)) {
                    odd_nodes = !odd_nodes;
                }
            }
            j = i;  // Passa al prossimo punto
        }

        return odd_nodes;
    }

    // Funzione per calcolo dell'equilibrio
    void equilibrium(
        const std::array<Array2D, 2>& u,             // u[0] e u[1] sono NxM
        const Matrix<double, 9, 2>& c,               // c[q][d]
        const Vector<double, 9>& w,                  // w[q]
        const Array2D& rho,                          // NxM
        Array3D& g_eq                                // g_eq[q] è NxM
    ) {
        int Nx = rho.rows();
        int Ny = rho.cols();

        Array2D v = 1.5 * (u[0].square() + u[1].square());

        for (int q = 0; q < 9; ++q) {
            Array2D t = 3.0 * (u[0] * c(q, 0) + u[1] * c(q, 1));
            g_eq[q] = (1.0 + t + 0.5 * t.square() - v);
            g_eq[q] *= rho * w[q];
        }
    }

    void collision_streaming(
        Array3D& g,                         // distribuzione post-streaming
        const Array3D& g_eq,                // equilibrio
        Array3D& g_up,                      // distribuzione temporanea
        double om_p, double om_m,           // frequenze di rilassamento
        const std::vector<int>& ns,         // mappa q -> q opposto
        int nx, int ny, int lx, int ly      // dimensioni griglia
    ) {
        // Collisione per q = 0
        g_up[0] = (1.0 - om_p) * g[0] + om_p * g_eq[0];
        g[0] = g_up[0];

        // Collisione per q = 1...8
        for (int q = 1; q < 9; ++q) {
            int qb = ns[q];
            g_up[q] =
                (1.0 - 0.5 * (om_p + om_m)) * g[q]
                - 0.5 * (om_p - om_m) * g[qb]
                + 0.5 * (om_p + om_m) * g_eq[q]
                + 0.5 * (om_p - om_m) * g_eq[qb];
        }

        // Streaming
        g[1].block(1, 0, lx, ny) = g_up[1].block(0, 0, lx, ny);             // east
        g[2].block(0, 0, lx, ny) = g_up[2].block(1, 0, lx, ny);             // west
        g[3].block(0, 1, nx, ly) = g_up[3].block(0, 0, nx, ly);             // north
        g[4].block(0, 0, nx, ly) = g_up[4].block(0, 1, nx, ly);             // south
        g[5].block(1, 1, lx, ly) = g_up[5].block(0, 0, lx, ly);             // northeast
        g[6].block(0, 0, lx, ly) = g_up[6].block(1, 1, lx, ly);             // southwest
        g[7].block(0, 1, lx, ly) = g_up[7].block(1, 0, lx, ly);             // northwest
        g[8].block(1, 0, lx, ly) = g_up[8].block(0, 1, lx, ly);             // southeast
    }

    std::pair<double, double> evaluate_drag_lift(
        const MatrixXi& boundary,                  // N x 3
        const std::vector<int>& ns,                // opposites
        const Matrix<int, 9, 2>& c,                // directions
        const Array3D& g_up, const Array3D& g,     // distributions
        double R_ref, double U_ref, double L_ref   // reference values
    ) {
        double fx = 0.0;
        double fy = 0.0;

        int N = boundary.rows();
        for (unsigned int k = 0; k < N; k++) {
            int i  = boundary(k, 0);
            int j  = boundary(k, 1);
            int q  = boundary(k, 2);
            int qb = ns[q];
            int cx = c(q, 0);
            int cy = c(q, 1);

            double g0 = g_up[q](i, j) + g[qb](i, j);
            fx += g0 * cx;
            fy += g0 * cy;
        }

        double Cx = -2.0 * fx / (R_ref * L_ref * U_ref * U_ref);
        double Cy = -2.0 * fy / (R_ref * L_ref * U_ref * U_ref);

        return std::make_pair(Cx, Cy);
    }

    /**
     * Applica le condizioni al contorno bounce-back per un ostacolo (no-slip).
     * Se `IBB == true`, usa interpolated bounce-back (migliore accuratezza).
     * Altrimenti, usa il bounce-back classico (semplice riflessione).
     *
     * @param IBB       Se true, applica interpolated bounce-back
     * @param boundary  Matrice N x 3 con righe (i, j, q): posizione ostacolo + direzione
     * @param ns        Vettore che mappa ciascuna direzione alla sua opposta
     * @param sc        Matrice 9 x 2: vettori velocità per ciascuna direzione
     * @param obs_ibb   Vettore di N double: posizione frazionaria dell’interfaccia (p)
     * @param g_up      Array3D: distribuzioni post-collisione
     * @param g         Array3D: distribuzioni da aggiornare (con bounce-back)
     */
    void bounce_back_obstacle(
        bool IBB,
        const MatrixXi& boundary,           // N x 3 array: (i, j, q)
        const std::vector<int>& ns,        // opposites: q̅ = ns[q]
        const Matrix<int, 9, 2>& sc,       // direzioni discrete D2Q9
        const VectorXd& obs_ibb,           // fattori p per interpolated BB
        const Array3D& g_up,               // distribuzioni post-collisione
        Array3D& g                         // distribuzioni da aggiornare
    ) {
        int N = boundary.rows();  // numero di punti ostacolo

        // Caso 1: Interpolated bounce-back (più accurato)
        if (IBB) {
            for (int k = 0; k < N; ++k) {
                // Estrai la posizione e la direzione
                int i  = boundary(k, 0);
                int j  = boundary(k, 1);
                int q  = boundary(k, 2);
                int qb = ns[q];             // direzione opposta

                // Ottieni la direzione (cx, cy) per qopp (verso l'interno)
                int cx = sc(qb, 0);
                int cy = sc(qb, 1);

                // Coordinate dei vicini secondo la direzione opposta
                int im  = i + cx;
                int jm  = j + cy;
                int imm = i + 2 * cx;
                int jmm = j + 2 * cy;

                // Coefficiente p dell'interfaccia di riflessione (0 < p < 1)
                double p  = obs_ibb(k);
                double pp = 2.0 * p;

                // Formula interpolata dipende da p
                if (p < 0.5) {
                    // Caso p < 0.5 (interfaccia vicina all’ostacolo)
                    g[qb](i, j) =
                        p * (pp + 1.0)          * g_up[q](i, j) +
                        (1.0 + pp) * (1.0 - pp) * g_up[q](im, jm) -
                        p * (1.0 - pp)          * g_up[q](imm, jmm);
                } else {
                    // Caso p >= 0.5 (interfaccia vicina al fluido)
                    g[qb](i, j) =
                        (1.0 / (p * (pp + 1.0)))       * g_up[q](i, j) +
                        ((pp - 1.0) / p)               * g_up[qb](i, j) +
                        ((1.0 - pp) / (1.0 + pp))      * g_up[qb](im, jm);
                }
            }
        }

        // Caso 2: Bounce-back classico (riflessione semplice)
        else {
            for (int k = 0; k < N; ++k) {
                int i  = boundary(k, 0);
                int j  = boundary(k, 1);
                int q  = boundary(k, 2);
                int qb = ns[q];             // direzione opposta

                // Direzione di rimbalzo
                int cx = sc(q, 0);
                int cy = sc(q, 1);
                int ii = i + cx;
                int jj = j + cy;

                // Applica bounce-back classico: riflette g_up[q] in g[qb]
                g[qb](i, j) = g_up[q](i, j);
            }
        }
    }

    // Funzione per calcolo della velocità sulla parete sinistra (Zou-He)
    void zou_he_left_wall_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,         // u[0] e u[1] sono NxM (passed as std::array)
        std::array<Array2D, 2>& u_left,                    // u_left[0,:] è un array NxM
        Array2D& rho,                             // rho è un array NxM
        Array3D& g                                // g[q] è un array NxM per ogni direzione q
    ) {
        double cst1 = 2.0 / 3.0;
        double cst2 = 1.0 / 6.0;
        double cst3 = 1.0 / 2.0;

        // Impostare u[0, 0, :] e u[1, 0, :] a u_left[0, :]
        u[0].row(0) = u_left[0].row(0);
        u[1].row(0) = u_left[1].row(0);

        // Calcolare rho[0, :]
        rho.row(0) = (g[0].row(0) + g[3].row(0) + g[4].row(0) +
                    2.0 * g[2].row(0) + 2.0 * g[6].row(0) +
                    2.0 * g[7].row(0)) / (1.0 - u[0].row(0));

        // Aggiornare g[1]
        g[1].row(0) = g[2].row(0) + cst1 * rho.row(0) * u[0].row(0);

        // Aggiornare g[5]
        g[5].row(0) = g[6].row(0) - cst3 * (g[3].row(0) - g[4].row(0)) +
                    cst2 * rho.row(0) * u[0].row(0) +
                    cst3 * rho.row(0) * u[1].row(0);

        // Aggiornare g[8]
        g[8].row(0) = g[7].row(0) + cst3 * (g[3].row(0) - g[4].row(0)) +
                    cst2 * rho.row(0) * u[0].row(0) -
                    cst3 * rho.row(0) * u[1].row(0);
    }

    void zou_he_right_wall_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,         // u[0] e u[1] sono NxM (passed as std::array)
        const std::array<Array2D, 2>& u_right,                   // u_right[0,:] è un array NxM
        Array2D& rho,                             // rho è un array NxM
        Array3D& g                                // g[q] è un array NxM per ogni direzione q
    ) {
        double cst1 = 2.0 / 3.0;
        double cst2 = 1.0 / 6.0;
        double cst3 = 1.0 / 2.0;

        // Impostare u[0, lx, :] e u[1, lx, :] a u_right[0, :]
        u[0].col(lx) = u_right[0].row(0);
        u[1].col(lx) = u_right[1].row(0);

        // Calcolare rho[lx, :]
        rho.row(lx) = (g[0].row(lx) + g[3].row(lx) + g[4].row(lx) +
                    2.0 * g[1].row(lx) + 2.0 * g[5].row(lx) +
                    2.0 * g[8].row(lx)) / (1.0 + u[0].col(lx));

        // Aggiornare g[2]
        g[2].row(lx) = g[1].row(lx) - cst1 * rho.row(lx) * u[0].col(lx);

        // Aggiornare g[6]
        g[6].row(lx) = g[5].row(lx) + cst3 * (g[3].row(lx) - g[4].row(lx)) -
                    cst2 * rho.row(lx) * u[0].col(lx) -
                    cst3 * rho.row(lx) * u[1].col(lx);

        // Aggiornare g[7]
        g[7].row(lx) = g[8].row(lx) - cst3 * (g[3].row(lx) - g[4].row(lx)) -
                    cst2 * rho.row(lx) * u[0].col(lx) +
                    cst3 * rho.row(lx) * u[1].col(lx);
    }

    // Funzione per calcolo della pressione sulla parete destra (Zou-He)
    void zou_he_right_wall_pressure(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,           // u[0] e u[1] sono NxM (passed as std::array)
        const Array2D& rho_right,                   // rho_right è un array NxM
        const std::array<Array2D, 2>& u_right,                     // u_right[0,:] è un array NxM
        Array2D& rho,                               // rho è un array NxM
        Array3D& g                                  // g[q] è un array NxM per ogni direzione q
    ) {
        double cst1 = 2.0 / 3.0;
        double cst2 = 1.0 / 6.0;
        double cst3 = 1.0 / 2.0;

        // Impostare rho[lx, :]
        rho.row(lx) = rho_right.row(lx);

        // Impostare u[1, lx, :]
        u[1].col(lx) = u_right[1].row(1);

        // Calcolare u[0, lx, :]
        u[0].col(lx) = (g[0].row(lx) + g[3].row(lx) + g[4].row(lx) +
                        2.0 * g[1].row(lx) + 2.0 * g[5].row(lx) +
                        2.0 * g[8].row(lx)) / rho.row(lx) - 1.0;

        // Aggiornare g[2]
        g[2].row(lx) = g[1].row(lx) - cst1 * rho.row(lx) * u[0].col(lx);

        // Aggiornare g[6]
        g[6].row(lx) = g[5].row(lx) + cst3 * (g[3].row(lx) - g[4].row(lx)) -
                    cst2 * rho.row(lx) * u[0].col(lx) -
                    cst3 * rho.row(lx) * u[1].col(lx);

        // Aggiornare g[7]
        g[7].row(lx) = g[8].row(lx) - cst3 * (g[3].row(lx) - g[4].row(lx)) -
                    cst2 * rho.row(lx) * u[0].col(lx) +
                    cst3 * rho.row(lx) * u[1].col(lx);
    }

    // Funzione per il calcolo della velocità sulla parete superiore (Zou-He)
    void zou_he_top_wall_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,           // u[0] e u[1] sono NxM (passed as std::array)
        const std::array<Array2D, 2>& u_top,                       // u_top[0,:] è un array NxM
        Array2D& rho,                               // rho è un array NxM
        Array3D& g                                  // g[q] è un array NxM per ogni direzione q
    ) {
        double cst1 = 2.0 / 3.0;
        double cst2 = 1.0 / 6.0;
        double cst3 = 1.0 / 2.0;

        // Impostare u[0, :, ly] e u[1, :, ly] a u_top[0,:]
        u[0].row(ly) = u_top[0].row(0);
        u[1].row(ly) = u_top[1].row(1);

        // Calcolare rho[:, ly]
        rho.col(ly) = (g[0].col(ly) + g[1].col(ly) + g[2].col(ly) +
                    2.0 * g[3].col(ly) + 2.0 * g[5].col(ly) +
                    2.0 * g[7].col(ly)) / (1.0 + u[1].row(ly));

        // Aggiornare g[4]
        g[4].col(ly) = g[3].col(ly) - cst1 * rho.col(ly) * u[1].row(ly);

        // Aggiornare g[8]
        g[8].col(ly) = g[7].col(ly) - cst3 * (g[1].col(ly) - g[2].col(ly)) +
                    cst3 * rho.col(ly) * u[0].row(ly) -
                    cst2 * rho.col(ly) * u[1].row(ly);

        // Aggiornare g[6]
        g[6].col(ly) = g[5].col(ly) + cst3 * (g[1].col(ly) - g[2].col(ly)) -
                    cst3 * rho.col(ly) * u[0].row(ly) -
                    cst2 * rho.col(ly) * u[1].row(ly);
    }

    // Funzione per il calcolo della velocità sulla parete inferiore (Zou-He)
    void zou_he_bottom_wall_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,           // u[0] e u[1] sono NxM (passed as std::array)
        const std::array<Array2D, 2>& u_bot,                       // u_bot[0,:] è un array NxM
        Array2D& rho,                               // rho è un array NxM
        Array3D& g                                  // g[q] è un array NxM per ogni direzione q
    ) {
        double cst1 = 2.0 / 3.0;
        double cst2 = 1.0 / 6.0;
        double cst3 = 1.0 / 2.0;

        // Impostare u[0, :, 0] e u[1, :, 0] a u_bot[0,:]
        u[0].row(0) = u_bot[0].row(0);
        u[1].row(0) = u_bot[1].row(1);

        // Calcolare rho[:, 0]
        rho.col(0) = (g[0].col(0) + g[1].col(0) + g[2].col(0) +
                    2.0 * g[4].col(0) + 2.0 * g[6].col(0) +
                    2.0 * g[8].col(0)) / (1.0 - u[1].row(0));

        // Aggiornare g[3]
        g[3].col(0) = g[4].col(0) + cst1 * rho.col(0) * u[1].row(0);

        // Aggiornare g[5]
        g[5].col(0) = g[6].col(0) - cst3 * (g[1].col(0) - g[2].col(0)) +
                    cst3 * rho.col(0) * u[0].row(0) +
                    cst2 * rho.col(0) * u[1].row(0);

        // Aggiornare g[7]
        g[7].col(0) = g[8].col(0) + cst3 * (g[1].col(0) - g[2].col(0)) -
                    cst3 * rho.col(0) * u[0].row(0) +
                    cst2 * rho.col(0) * u[1].row(0);
    }

    // Funzione per il calcolo della velocità sulla parete inferiore sinistra (Zou-He)
    void zou_he_bottom_left_corner_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,         // u[0] e u[1] sono NxM (passed as std::array)
        Array2D& rho,                      // rho è un array NxM
        Array3D& g                         // g[q] è un array NxM per ogni direzione q
    ) {
        // Aggiornamenti per u[0,0,0] e u[1,0,0] (velocità angolo in basso a sinistra)
        u[0](0, 0) = u[0](1, 0);  // Impostiamo la velocità x del angolo in basso a sinistra
        u[1](0, 0) = u[1](1, 0);  // Impostiamo la velocità y del angolo in basso a sinistra

        // Aggiornamenti per la densità in basso a sinistra
        rho(0, 0) = rho(1, 0);  // Impostiamo la densità nel angolo in basso a sinistra

        // Aggiornamenti per g[1], g[3], g[5] in basso a sinistra
        g[1](0, 0) = g[2](0, 0) + (2.0 / 3.0) * rho(0, 0) * u[0](0, 0);
        g[3](0, 0) = g[4](0, 0) + (2.0 / 3.0) * rho(0, 0) * u[1](0, 0);
        g[5](0, 0) = g[6](0, 0) + (1.0 / 6.0) * rho(0, 0) * u[0](0, 0) +
                                    (1.0 / 6.0) * rho(0, 0) * u[1](0, 0);

        // Impostiamo g[7] e g[8] a zero
        g[7](0, 0) = 0.0;
        g[8](0, 0) = 0.0;

        // Calcolare g[0] come il residuo per la conservazione della massa
        g[0](0, 0) = rho(0, 0) - g[1](0, 0) - g[2](0, 0) - g[3](0, 0) - g[4](0, 0) -
                    g[5](0, 0) - g[6](0, 0) - g[7](0, 0) - g[8](0, 0);
    }

    // Function to implement Zou-He no-slip top left corner velocity boundary condition
    void zou_he_top_left_corner_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,         // u[0] and u[1] are NxM
        Array2D& rho,                      // rho is an NxM array
        Array3D& g                         // g[q] is an NxM array for each direction q
    ) {
        // Update velocity at the top-left corner (u[0](0, ly) = u[0](1, ly), u[1](0, ly) = u[1](1, ly))
        u[0](0, ly) = u[0](1, ly);  // Set u[0](0, ly) = u[0](1, ly)
        u[1](0, ly) = u[1](1, ly);  // Set u[1](0, ly) = u[1](1, ly)

        // Update density at the top-left corner (rho(0, ly) = rho(1, ly))
        rho(0, ly) = rho(1, ly);  // Set rho(0, ly) = rho(1, ly)

        // Update distribution functions for g at the top-left corner
        g[1](0, ly) = g[2](0, ly) + (2.0 / 3.0) * rho(0, ly) * u[0](0, ly);
        g[4](0, ly) = g[3](0, ly) - (2.0 / 3.0) * rho(0, ly) * u[1](0, ly);
        g[8](0, ly) = g[7](0, ly) + (1.0 / 6.0) * rho(0, ly) * u[0](0, ly) - 
                    (1.0 / 6.0) * rho(0, ly) * u[1](0, ly);

        // Set g[5] and g[6] to zero
        g[5](0, ly) = 0.0;
        g[6](0, ly) = 0.0;

        // Compute mass conservation (g[0](0, ly) = rho(0, ly) - sum(g[1] to g[8]))
        g[0](0, ly) = rho(0, ly) - g[1](0, ly) - g[2](0, ly) - g[3](0, ly) - g[4](0, ly) -
                    g[5](0, ly) - g[6](0, ly) - g[7](0, ly) - g[8](0, ly);
    }

    // Function to implement Zou-He no-slip top right corner velocity boundary condition
    void zou_he_top_right_corner_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,         // u[0] and u[1] are NxM
        Array2D& rho,                      // rho is an NxM array
        Array3D& g                         // g[q] is an NxM array for each direction q
    ) {
        // Update velocity at the top-right corner (u[0](lx, ly) = u[0](lx-1, ly), u[1](lx, ly) = u[1](lx-1, ly))
        u[0](lx, ly) = u[0](lx - 1, ly);  // Set u[0](lx, ly) = u[0](lx-1, ly)
        u[1](lx, ly) = u[1](lx - 1, ly);  // Set u[1](lx, ly) = u[1](lx-1, ly)

        // Update density at the top-right corner (rho(lx, ly) = rho(lx-1, ly))
        rho(lx, ly) = rho(lx - 1, ly);  // Set rho(lx, ly) = rho(lx-1, ly)

        // Update distribution functions for g at the top-right corner
        g[2](lx, ly) = g[1](lx, ly) - (2.0 / 3.0) * rho(lx, ly) * u[0](lx, ly);
        g[4](lx, ly) = g[3](lx, ly) - (2.0 / 3.0) * rho(lx, ly) * u[1](lx, ly);
        g[6](lx, ly) = g[5](lx, ly) - (1.0 / 6.0) * rho(lx, ly) * u[0](lx, ly) - 
                    (1.0 / 6.0) * rho(lx, ly) * u[1](lx, ly);

        // Set g[7] and g[8] to zero
        g[7](lx, ly) = 0.0;
        g[8](lx, ly) = 0.0;

        // Compute mass conservation (g[0](lx, ly) = rho(lx, ly) - sum(g[1] to g[8]))
        g[0](lx, ly) = rho(lx, ly) - g[1](lx, ly) - g[2](lx, ly) - g[3](lx, ly) - g[4](lx, ly) -
                    g[5](lx, ly) - g[6](lx, ly) - g[7](lx, ly) - g[8](lx, ly);
    }

    // Function to implement Zou-He no-slip bottom right corner velocity boundary condition
    void zou_he_bottom_right_corner_velocity(
        int lx, 
        int ly, 
        std::array<Array2D, 2>& u,         // u[0] and u[1] are NxM
        Array2D& rho,                      // rho is an NxM array
        Array3D& g                         // g[q] is an NxM array for each direction q
    ) {
        // Update velocity at the bottom-right corner (u[0](lx, 0) = u[0](lx-1, 0), u[1](lx, 0) = u[1](lx-1, 0))
        u[0](lx, 0) = u[0](lx - 1, 0);  // Set u[0](lx, 0) = u[0](lx-1, 0)
        u[1](lx, 0) = u[1](lx - 1, 0);  // Set u[1](lx, 0) = u[1](lx-1, 0)

        // Update density at the bottom-right corner (rho(lx, 0) = rho(lx-1, 0))
        rho(lx, 0) = rho(lx - 1, 0);  // Set rho(lx, 0) = rho(lx-1, 0)

        // Update distribution functions for g at the bottom-right corner
        g[2](lx, 0) = g[1](lx, 0) - (2.0 / 3.0) * rho(lx, 0) * u[0](lx, 0);
        g[3](lx, 0) = g[4](lx, 0) + (2.0 / 3.0) * rho(lx, 0) * u[1](lx, 0);
        g[7](lx, 0) = g[8](lx, 0) - (1.0 / 6.0) * rho(lx, 0) * u[0](lx, 0) + 
                    (1.0 / 6.0) * rho(lx, 0) * u[1](lx, 0);

        // Set g[5] and g[6] to zero
        g[5](lx, 0) = 0.0;
        g[6](lx, 0) = 0.0;

        // Compute mass conservation (g[0](lx, 0) = rho(lx, 0) - sum(g[1] to g[8]))
        g[0](lx, 0) = rho(lx, 0) - g[1](lx, 0) - g[2](lx, 0) - g[3](lx, 0) - g[4](lx, 0) -
                    g[5](lx, 0) - g[6](lx, 0) - g[7](lx, 0) - g[8](lx, 0);
    }

};
#endif