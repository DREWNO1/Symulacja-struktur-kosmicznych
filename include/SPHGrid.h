#ifndef SPHGRID_H
#define SPHGRID_H

#include <vector>
#include <cmath>     // Dla std::floor, std::ceil
#include <algorithm> // Dla std::min, std::max
#include <limits>    // Dla std::numeric_limits

#include "Vector3d.h"
#include "Particle.h"  // Zakładamy, że Particle.h definiuje struct Particle i enum ParticleType
#include "SimConfig.h" // Dla SimConfig::H_SPH

// Struktura reprezentująca komórkę siatki SPH
struct GridCell {
    std::vector<size_t> particle_indices; // Indeksy cząstek w tej komórce
};

class SPHGrid {
public:
    SPHGrid();

    // Buduje i wypełnia siatkę cząstkami gazu
    void build(const std::vector<Particle>& particleList, double sph_smoothing_length);

    // Zwraca indeksy cząstek w komórce o podanych współrzędnych (x,y,z) siatki
    // oraz w jej 26 sąsiadujących komórkach.
    // Zwraca pusty wektor, jeśli komórka bazowa jest poza granicami.
    std::vector<size_t> getParticlesInNeighboringCells(int cell_x_base, int cell_y_base, int cell_z_base) const;
    
    // Zwraca indeksy cząstek tylko w DOKŁADNIE tej jednej komórce
    const std::vector<size_t>& getParticlesInCell(int cell_x, int cell_y, int cell_z) const;


    // Oblicza współrzędne komórki siatki dla danej pozycji cząstki
    // Zwraca true jeśli pozycja jest w granicach siatki, false w przeciwnym razie.
    // Współrzędne komórki są zapisywane w out_cell_x, out_cell_y, out_cell_z.
    bool getCellCoordinates(const Vector3d& position, int& out_cell_x, int& out_cell_y, int& out_cell_z) const;

    // Sprawdza, czy siatka została poprawnie zbudowana (ma wymiary > 0)
    bool isBuilt() const;
    
    // Gettery dla wymiarów siatki, jeśli potrzebne na zewnątrz
    int getDimX() const { return grid_dim_x_; }
    int getDimY() const { return grid_dim_y_; }
    int getDimZ() const { return grid_dim_z_; }

private:
    std::vector<GridCell> cells_;         // Wektor komórek siatki
    Vector3d min_bounds_;                 // Minimalne granice przestrzeni siatki
    double cell_size_;                    // Rozmiar boku komórki siatki
    int grid_dim_x_, grid_dim_y_, grid_dim_z_; // Wymiary siatki w liczbie komórek

    // Konwertuje współrzędne komórki (x,y,z) na płaski indeks w wektorze cells_
    size_t getFlatIndex(int x, int y, int z) const;
};

#endif // SPHGRID_H