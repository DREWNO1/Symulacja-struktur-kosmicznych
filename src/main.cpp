// main.cpp
#include "Simulation.h"
#include <iostream>
// #include "SimConfig.h" // Nie jest już bezpośrednio potrzebny tutaj, jeśli Simulation zarządza GUI vars

int main(int argc, char* argv[]) {
    Simulation sim_app;

    // Inicjalizacja wartości GUI w SimConfig jest już niepotrzebna,
    // ponieważ Simulation::initialize() i Simulation::updateSimulationState()
    // będą zarządzać przepływem wartości między polami _gui_ klasy Simulation
    // a właściwymi zmiennymi w SimConfig.

    if (sim_app.initialize("Cosmic Structure Simulation Refactored", 1920, 1080)) {
        std::cout << "Info: Simulation started (Refactored with classes)." << std::endl;
        std::cout << "Info: Controls: Arrow keys (orbit), A/Z (zoom), Mouse LMB+drag (orbit), Mouse Wheel (zoom), R (reset), P (pause)" << std::endl;
        sim_app.run();
    } else {
        // Komunikat o błędzie jest już w Simulation::initialize lub Renderer::initialize
        // std::cerr << "FATAL ERROR: Simulation Initialization failed. Exiting." << std::endl;
        // sim_app.cleanup(); // Destruktor powinien to zrobić
        return 1;
    }

    // sim_app.cleanup(); // Destruktor powinien to zrobić
    std::cout << "Simulation finished." << std::endl;
    return 0;
}