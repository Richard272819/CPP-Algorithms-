#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Structure simple pour représenter un point dans le bulk AdS3 (Plan de Poincaré)
struct PointAdS {
    double x; // Espace de la frontière
    double z; // Direction holographique (z -> 0 c'est la frontière UV)
};

class SimulateurRyuTakayanagi {
private:
    double R_ads;        // Rayon de l'espace Anti-de Sitter
    double G_newton;     // Constante de Newton en 3D
    double cutoff_uv;    // Le fameux cut-off pour éviter la divergence en z = 0
    double L_frontiere;  // Taille du sous-système sur la frontière
    int n_points;        // Nombre de segments pour discrétiser notre courbe

    std::vector<PointAdS> courbe;

public:
    SimulateurRyuTakayanagi(double r, double g, double eps, double l, int pts)
        : R_ads(r), G_newton(g), cutoff_uv(eps), L_frontiere(l), n_points(pts) {
        initialiser_courbe_droite();
    }

    // On commence avec une ligne droite naïve sur le cut-off UV
    void initialiser_courbe_droite() {
        courbe.clear();
        double x_start = -L_frontiere / 2.0;
        double dx = L_frontiere / (n_points - 1);

        for (int i = 0; i < n_points; ++i) {
            courbe.push_back({x_start + i * dx, cutoff_uv});
        }
    }

    // Calcule la longueur propre totale dans la métrique d'AdS3
    double calculer_longueur_ads() {
        double longueur_totale = 0.0;
        for (size_t i = 0; i < courbe.size() - 1; ++i) {
            double dx = courbe[i+1].x - courbe[i].x;
            double dz = courbe[i+1].z - courbe[i].z;
            double z_moy = (courbe[i].z + courbe[i+1].z) / 2.0;

            // Métrique d'AdS3: ds^2 = (R^2 / z^2) * (dx^2 + dz^2)
            double ds = (R_ads / z_moy) * std::sqrt(dx * dx + dz * dz);
            longueur_totale += ds;
        }
        return longueur_totale;
    }

    // Algorithme de relaxation pour minimiser la longueur (trouver la géodésique)
    void optimiser_surface_minimale(int etapes, double taux_apprentissage) {
        // On ne touche pas aux extrémités fixées sur la frontière
        for (int etape = 0; etape < etapes; ++etape) {
            std::vector<PointAdS> nouvelle_courbe = courbe;

            for (int i = 1; i < n_points - 1; ++i) {
                // Évaluation numérique locale des forces géométriques (gradient)
                double original_z = courbe[i].z;
                
                // Petit shift pour calculer la dérivée numérique de la longueur
                courbe[i].z = original_z + 1e-5;
                double len_plus = calculer_longueur_ads();

                courbe[i].z = original_z - 1e-5;
                double len_moins = calculer_longueur_ads();

                courbe[i].z = original_z; // Reset

                double gradient_z = (len_plus - len_moins) / 2e-5;

                // Descente de gradient pour pousser la courbe vers l'intérieur du bulk
                nouvelle_courbe[i].z = original_z - taux_apprentissage * gradient_z;
                
                // Sécurité physique pour pas traverser la frontière
                if (nouvelle_courbe[i].z < cutoff_uv) {
                    nouvelle_courbe[i].z = cutoff_uv;
                }
            }
            courbe = nouvelle_courbe;
        }
    }

    // Application directe de la formule de Ryu-Takayanagi: S = Aire / (4 * G_N)
    double calculer_entropie_holographique() {
        double aire_minimale = calculer_longueur_ads();
        return aire_minimale / (4.0 * G_newton);
    }

    void afficher_profil() {
        std::cout << "\n--- Profil de la surface minimale (Geodésique AdS) ---" << std::endl;
        std::cout << std::setw(12) << "X" << std::setw(15) << "Z (Bulk)" << std::endl;
        // On affiche un échantillon pour pas saturer la console
        int saut = n_points / 10;
        for (int i = 0; i < n_points; i += saut) {
            std::cout << std::setw(12) << courbe[i].x << std::setw(15) << courbe[i].z << std::endl;
        }
    }
};

int main() {
    std::cout << "Lanceur de la simulation de surface minimale Ryu-Takayanagi..." << std::endl;

    // Paramètres physiques
    double R = 1.0;          // Rayon AdS
    double G = 0.05;         // Constante G de Newton
    double epsilon = 0.01;   // Régularisation UV (cutoff)
    double L = 2.0;          // Taille de la région sur la frontière
    int points = 100;        // Résolution de la discrétisation

    SimulateurRyuTakayanagi sim(R, G, epsilon, L, points);

    std::cout << "Longueur initiale (ligne droite): " << sim.calculer_longueur_ads() << " - Clairement pas optimal." << std::endl;

    // Lancement de l'optimisation géométrique
    sim.optimiser_surface_minimale(5000, 0.002);

    std::cout << "Longueur finale (apres relaxation): " << sim.calculer_longueur_ads() << std::endl;
    
    sim.afficher_profil();

    // Résultat final
    double entropie = sim.calculer_entropie_holographique();
    std::cout << "\n=======================================================" << std::endl;
    std::cout << " Entropie d'intrication holographique (S_A): " << entropie << std::endl;
    std::cout << "=======================================================" << std::endl;

    return 0;
}
