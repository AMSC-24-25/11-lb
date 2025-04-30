#include <Eigen/Dense>
#include <vector>
#include <array>

using namespace Eigen;

// Alias for clarity
using Array2D = ArrayXXd;
using Array3D = std::vector<Array2D>;  // array 3D: 9 x (NxM)

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
    const Array2D& u_left,                    // u_left[0,:] è un array NxM
    Array2D& rho,                             // rho è un array NxM
    Array3D& g                                // g[q] è un array NxM per ogni direzione q
) {
    double cst1 = 2.0 / 3.0;
    double cst2 = 1.0 / 6.0;
    double cst3 = 1.0 / 2.0;

    // Impostare u[0, 0, :] e u[1, 0, :] a u_left[0, :]
    u[0].row(0) = u_left.row(0);
    u[1].row(0) = u_left.row(0);

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
    const std::array<Array2D, 2>& u,         // u[0] e u[1] sono NxM (passed as std::array)
    const Array2D& u_right,                   // u_right[0,:] è un array NxM
    Array2D& rho,                             // rho è un array NxM
    Array3D& g                                // g[q] è un array NxM per ogni direzione q
) {
    double cst1 = 2.0 / 3.0;
    double cst2 = 1.0 / 6.0;
    double cst3 = 1.0 / 2.0;

    // Impostare u[0, lx, :] e u[1, lx, :] a u_right[0, :]
    u[0].col(lx) = u_right.row(0);
    u[1].col(lx) = u_right.row(0);

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
    const std::array<Array2D, 2>& u,           // u[0] e u[1] sono NxM (passed as std::array)
    const Array2D& rho_right,                   // rho_right è un array NxM
    const Array2D& u_right,                     // u_right[0,:] è un array NxM
    Array2D& rho,                               // rho è un array NxM
    Array3D& g                                  // g[q] è un array NxM per ogni direzione q
) {
    double cst1 = 2.0 / 3.0;
    double cst2 = 1.0 / 6.0;
    double cst3 = 1.0 / 2.0;

    // Impostare rho[lx, :]
    rho.row(lx) = rho_right.row(lx);

    // Impostare u[1, lx, :]
    u[1].col(lx) = u_right.row(1);

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
    const std::array<Array2D, 2>& u,           // u[0] e u[1] sono NxM (passed as std::array)
    const Array2D& u_top,                       // u_top[0,:] è un array NxM
    Array2D& rho,                               // rho è un array NxM
    Array3D& g                                  // g[q] è un array NxM per ogni direzione q
) {
    double cst1 = 2.0 / 3.0;
    double cst2 = 1.0 / 6.0;
    double cst3 = 1.0 / 2.0;

    // Impostare u[0, :, ly] e u[1, :, ly] a u_top[0,:]
    u[0].row(ly) = u_top.row(0);
    u[1].row(ly) = u_top.row(1);

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
    const std::array<Array2D, 2>& u,           // u[0] e u[1] sono NxM (passed as std::array)
    const Array2D& u_bot,                       // u_bot[0,:] è un array NxM
    Array2D& rho,                               // rho è un array NxM
    Array3D& g                                  // g[q] è un array NxM per ogni direzione q
) {
    double cst1 = 2.0 / 3.0;
    double cst2 = 1.0 / 6.0;
    double cst3 = 1.0 / 2.0;

    // Impostare u[0, :, 0] e u[1, :, 0] a u_bot[0,:]
    u[0].row(0) = u_bot.row(0);
    u[1].row(0) = u_bot.row(1);

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