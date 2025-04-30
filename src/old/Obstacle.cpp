#include <Eigen/Dense>
#include <vector>
#include <string>
#ifndef Obstacle_CPP
#define Obstacle_CPP

class Obstacle {
public:
    // Attributi della classe
    std::string name;   // Nome dell'ostacolo
    int n_pts;          // Numero di punti
    int n_spts;         // Numero di sotto-punti
    std::string type;   // Tipo di ostacolo
    Eigen::Vector2i size; // Dimensione dell'ostacolo (due interi per larghezza e altezza)
    Eigen::Vector2i pos;  // Posizione dell'ostacolo

    // Attributi aggiuntivi
    Eigen::MatrixXd polygon; // Poligono che definisce l'ostacolo
    std::string tag;           // Etichetta dell'ostacolo
    Eigen::ArrayXXd area;      // Area dell'ostacolo
    Eigen::ArrayXXd boundary;  // Boundary dell'ostacolo
    Eigen::ArrayXXd ibb;       // Indice di posizione o altro

    // Costruttore della classe
    Obstacle(std::string n, int n_pts, int n_spts, std::string t, Eigen::Vector2i s, Eigen::Vector2i p)
        : name(n), n_pts(n_pts), n_spts(n_spts), type(t), size(s), pos(p) {}

    // Metodi della classe
    void set_polygon(const Eigen::MatrixXd& poly) {
        polygon = poly;
    }

    void set_tag(std::string& t) {
        tag = t;
    }

    void fill(const Eigen::ArrayXXd& a, const Eigen::ArrayXXd& b, const Eigen::ArrayXXd& i) {
        area = a;
        boundary = b;
        ibb = i;
    }
};
#endif