
#include "Simulation.h"
#include <iostream>


int main(int argc, char* argv[]) {
    Simulation sim_app;

    if (sim_app.initialize("Symulacja Struktur Kosmicznych", 1920, 1080)) {
        std::cout << "Info: Symulacja rozpoczeta." << std::endl;
        std::cout << "Info: Sterowanie: Strzalki (orbita), A/Z (zoom), Lewy Przycisk Myszy+przeciaganie (orbita), Kółko Myszy (zoom), R (reset), P (pauza)" << std::endl;
        sim_app.run();
    } else {
        return 1;
    }
    
    std::cout << "Symulacja zakończona." << std::endl;
    return 0;
}