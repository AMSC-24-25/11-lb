#include "3DLBM.hpp"


TDLBM::TDLBM(unsigned int nx_, unsigned int ny_, unsigned int nz_, double U_lid_, double Re_)
    : nx(nx_), ny(ny_), nz(nz_), U_lid(U_lid_), Re(Re_)
{
    dx    = 1.0;
    dt    = 1.0;
    rho0  = 1.0;
    L     = (ny - 1) * dx;
    nu    = (U_lid * L) / Re;            
    tau   = (6.0 * nu + 1.0) / 2.0;      
    ncells = nx * ny * nz;

    e.resize(Q);
    f.resize(ncells * Q, 0.0);
    f_new.resize(ncells * Q, 0.0);
    rho.resize(ncells, 0.0);
    u.resize(ncells, {0.0, 0.0, 0.0});
    u_prev.resize(ncells, {0.0, 0.0, 0.0});

    set_e_values(e);
    init_rho_v(rho, u);
    rho[idx(nx/2, ny/2, nz/2)] *= 1.0001;
    init_f(f, f_new, rho, u, e);
}

double TDLBM::weight(int i) const {
    if (i == 0)            return 1.0/3.0;
    else if (i < 7)        return 1.0/18.0;
    else                   return 1.0/36.0;
}

int TDLBM::get_opposite_direction(int i) const {
    switch(i) {
            case 0: return 0;  
            case 1: return 2;  
            case 2: return 1;  
            case 3: return 4;  
            case 4: return 3;  
            case 5: return 6;  
            case 6: return 5;  
            case 7: return 8;  
            case 8: return 7;
            case 9: return 10;
            case 10: return 9;
            case 11: return 12;
            case 12: return 11;
            case 13: return 14;
            case 14: return 13;
            case 15: return 16;
            case 16: return 15;
            case 17: return 18;
            case 18: return 17;
            default: return 0;
        }
}

int TDLBM::idx(int x, int y, int z) const {
    return x + nx * y + nx * ny * z;
}

int TDLBM::idxi(int x, int y, int z, int i) const {
    return idx(x,y,z) * Q + i;
}

void TDLBM::init_rho_v(std::vector<double>& rho, std::vector<TDouble>& u){
    #pragma omp parallel for collapse(3)
    for (int z=0; z<nz; ++z)
      for (int y=0; y<ny; ++y)
        for (int x=0; x<nx; ++x){
          int id = idx(x,y,z);
          rho[id] = rho0;
          u[id]   = Set_TDouble(0.0,0.0,0.0);
        }
}

void TDLBM::init_f(std::vector<double>& f,
                   std::vector<double>& f_new,
                   const std::vector<double>& rho,
                   const std::vector<TDouble>& u,
                   const std::vector<TDouble>& e)
{
    #pragma omp parallel for collapse(3)
    for (int z=0; z<nz; ++z)
     for (int y=0; y<ny; ++y)
      for (int x=0; x<nx; ++x){
        int cell = idx(x,y,z);
        for (int i=0; i<Q; ++i){
          double fval = feq(rho[cell], weight(i), e[i], u[cell]);
          int pos = idxi(x,y,z,i);
          f[pos]     = fval;
          f_new[pos] = fval;
        }
      }
}

void TDLBM::set_e_values(std::vector<TDouble>& e){
   if (e.size() != Q) e.resize(Q);
    e[0]  = Set_TDouble( 0.0,  0.0,  0.0);  
    e[1]  = Set_TDouble( 1.0,  0.0,  0.0);
    e[2]  = Set_TDouble(-1.0,  0.0,  0.0);
    e[3]  = Set_TDouble( 0.0,  1.0,  0.0);
    e[4]  = Set_TDouble( 0.0, -1.0,  0.0);
    e[5]  = Set_TDouble( 0.0,  0.0,  1.0);
    e[6]  = Set_TDouble( 0.0,  0.0, -1.0);
    e[7]  = Set_TDouble( 1.0,  1.0,  0.0);
    e[8]  = Set_TDouble(-1.0, -1.0,  0.0);
    e[9]  = Set_TDouble(-1.0,  1.0,  0.0);
    e[10] = Set_TDouble( 1.0, -1.0,  0.0);
    e[11] = Set_TDouble( 1.0,  0.0,  1.0);
    e[12] = Set_TDouble(-1.0,  0.0, -1.0);
    e[13] = Set_TDouble( 0.0,  1.0,  1.0);
    e[14] = Set_TDouble( 0.0, -1.0, -1.0);
    e[15] = Set_TDouble(-1.0,  0.0,  1.0);
    e[16] = Set_TDouble( 1.0,  0.0, -1.0);
    e[17] = Set_TDouble( 0.0, -1.0,  1.0);
    e[18] = Set_TDouble( 0.0,  1.0, -1.0);
}

double TDLBM::feq(double rho_, double w,
                  const TDouble& ei, const TDouble& ui) const
{
    double eu = dot(ei, ui);
    double uu = dot(ui, ui);
    return rho_ * w * (1.0 + 3.0*eu + 4.5*eu*eu - 1.5*uu);
}

double TDLBM::compute_rho(const std::vector<double>& f, int x, int y, int z) const {
    double sum = 0.0;
    for (int i=0; i<Q; ++i)
      sum += f[idxi(x,y,z,i)];
    return sum;
}

TDouble TDLBM::compute_u(const std::vector<double>& f, double rho_,
                         const std::vector<TDouble>& e, int x, int y, int z) const
{
    TDouble u_cell = {0.0,0.0,0.0};
    for (int i=0; i<Q; ++i){
      int pos = idxi(x,y,z,i);
      u_cell.x += f[pos] * e[i].x;
      u_cell.y += f[pos] * e[i].y;
      u_cell.z += f[pos] * e[i].z;
    }
    u_cell.x /= rho_;
    u_cell.y /= rho_;
    u_cell.z /= rho_;
    return u_cell;
}

void TDLBM::collide_trt(std::vector<double>& f,
                        std::vector<double>& rho,
                        std::vector<TDouble>& u,
                        const std::vector<TDouble>& e)
{
    double tau_val = get_tau();
    double omega_plus = 1.0 / tau_val;
    double temp = (1.0 / omega_plus - 0.5);
    double omega_minus = 1.0 / ((0.25 / temp) + 0.5);

    #pragma omp parallel for collapse(3)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int cell = idx(x, y, z);
                double rho_cell = compute_rho(f, x, y, z);
                TDouble u_cell = compute_u(f, rho_cell, e, x, y, z);
                rho[cell] = rho_cell;
                u[cell] = u_cell;
                int pos0 = idxi(x, y, z, 0);
                double feq0 = feq(rho_cell, weight(0), e[0], u_cell);
                f[pos0] = f[pos0] - omega_plus * (f[pos0] - feq0);

                for (int i = 1; i < Q; ++i) {
                    int opp = get_opposite_direction(i);
                    if (i < opp) {
                        int pos_i = idxi(x, y, z, i);
                        int pos_opp = idxi(x, y, z, opp);

                        double f_even = 0.5 * (f[pos_i] + f[pos_opp]);
                        double f_odd  = 0.5 * (f[pos_i] - f[pos_opp]);

                        double feq_i   = feq(rho_cell, weight(i), e[i], u_cell);
                        double feq_opp = feq(rho_cell, weight(opp), e[opp], u_cell);

                        double feq_even = 0.5 * (feq_i + feq_opp);
                        double feq_odd  = 0.5 * (feq_i - feq_opp);

                        double f_even_new = f_even - omega_plus * (f_even - feq_even);
                        double f_odd_new  = f_odd  - omega_minus * (f_odd  - feq_odd);

                        f[pos_i]   = f_even_new + f_odd_new;
                        f[pos_opp] = f_even_new - f_odd_new;
                    }
                }
            }
        }
    }
}

void TDLBM::stream(const std::vector<double>& f,
                   std::vector<double>& f_new,
                   const std::vector<TDouble>& e)
{
    std::fill(f_new.begin(), f_new.end(), 0.0);
    #pragma omp parallel for collapse(3)
    for (int z = 0; z < nz; ++z) {
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                for (int i = 0; i < Q; ++i) {
                    int x_dest = x + static_cast<int>(e[i].x);
                    int y_dest = y + static_cast<int>(e[i].y);
                    int z_dest = z + static_cast<int>(e[i].z);
                    
                    if (x_dest >= 0 && x_dest < nx && 
                        y_dest >= 0 && y_dest < ny && 
                        z_dest >= 0 && z_dest < nz) {
                        f_new[idxi(x_dest, y_dest, z_dest, i)] = f[idxi(x, y, z, i)];
                    } else {
                        int i_opposite = get_opposite_direction(i);
                        f_new[idxi(x, y, z, i_opposite)] = f[idxi(x, y, z, i)];
                    }
                }
            }
        }
    }
}

void TDLBM::apply_boundary_conditions(std::vector<double>& f,
                                      const std::vector<TDouble>& e)
{
   #pragma omp parallel
    {
        
        #pragma omp for collapse(2) schedule(static)
        for (int z = 0; z < nz; ++z) {
            for (int x = 0; x < nx; ++x) {
                int y = ny - 1;
                double rho_wall = compute_rho(f, x, y, z);
                TDouble u_wall = {U_lid, 0.0, 0.0};
                for (int i = 0; i < Q; ++i) {
                    if (e[i].y < 0) { 
                        int pos = idxi(x, y, z, i);
                        int i_opposite = get_opposite_direction(i);
                        f[pos] = f[idxi(x, y, z, i_opposite)]
                                 - 6.0 * rho_wall * weight(i) * dot(e[i], u_wall);
                    }
                }
            }
        }

        
        #pragma omp for collapse(2) schedule(static)
        for (int z = 0; z < nz; ++z) {
            for (int x = 0; x < nx; ++x) {
                int y = 0;
                double rho_wall = compute_rho(f, x, y, z);
                TDouble u_wall = {0.0, 0.0, 0.0};
                for (int i = 0; i < Q; ++i) {
                    if (e[i].y > 0) { 
                        int pos = idxi(x, y, z, i);
                        int i_opposite = get_opposite_direction(i);
                        f[pos] = f[idxi(x, y, z, i_opposite)]
                                 - 6.0 * rho_wall * weight(i) * dot(e[i], u_wall);
                    }
                }
            }
        }

        
        #pragma omp for collapse(2) schedule(static)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                int x = 0;
                double rho_wall = compute_rho(f, x, y, z);
                TDouble u_wall = {0.0, 0.0, 0.0};
                for (int i = 0; i < Q; ++i) {
                    if (e[i].x > 0) { 
                        int pos = idxi(x, y, z, i);
                        int i_opposite = get_opposite_direction(i);
                        f[pos] = f[idxi(x, y, z, i_opposite)]
                                 - 6.0 * rho_wall * weight(i) * dot(e[i], u_wall);
                    }
                }
            }
        }

        
        #pragma omp for collapse(2) schedule(static)
        for (int z = 0; z < nz; ++z) {
            for (int y = 0; y < ny; ++y) {
                int x = nx - 1;
                double rho_wall = compute_rho(f, x, y, z);
                TDouble u_wall = {0.0, 0.0, 0.0};
                for (int i = 0; i < Q; ++i) {
                    if (e[i].x < 0) { 
                        int pos = idxi(x, y, z, i);
                        int i_opposite = get_opposite_direction(i);
                        f[pos] = f[idxi(x, y, z, i_opposite)]
                                 - 6.0 * rho_wall * weight(i) * dot(e[i], u_wall);
                    }
                }
            }
        }

        
        #pragma omp for collapse(2) schedule(static)
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int z = 0;
                double rho_wall = compute_rho(f, x, y, z);
                TDouble u_wall = {0.0, 0.0, 0.0};
                for (int i = 0; i < Q; ++i) {
                    if (e[i].z > 0) { 
                        int pos = idxi(x, y, z, i);
                        int i_opposite = get_opposite_direction(i);
                        f[pos] = f[idxi(x, y, z, i_opposite)]
                                 - 6.0 * rho_wall * weight(i) * dot(e[i], u_wall);
                    }
                }
            }
        }

        
        #pragma omp for collapse(2) schedule(static)
        for (int y = 0; y < ny; ++y) {
            for (int x = 0; x < nx; ++x) {
                int z = nz - 1;
                double rho_wall = compute_rho(f, x, y, z);
                TDouble u_wall = {0.0, 0.0, 0.0};
                for (int i = 0; i < Q; ++i) {
                    if (e[i].z < 0) { 
                        int pos = idxi(x, y, z, i);
                        int i_opposite = get_opposite_direction(i);
                        f[pos] = f[idxi(x, y, z, i_opposite)]
                                 - 6.0 * rho_wall * weight(i) * dot(e[i], u_wall);
                    }
                }
            }
        }
    }
}

void TDLBM::evolution(){
    u_prev = u;
    collide_trt(f, rho, u, e);
    stream(f, f_new, e);
    f.swap(f_new);
    apply_boundary_conditions(f, e);
}
