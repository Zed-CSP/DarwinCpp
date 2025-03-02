#include "neat.h"
#include "organism.h"
#include "world.h"
#include <iostream>
#include <chrono>
#include <thread>

int main() {
    std::cout << "Initializing world..." << std::endl;
    
    // Create a world with dimensions 800x600, 50 initial organisms, and 200 food items
    World world(800.0, 600.0, 50, 200);
    
    std::cout << "World initialized. Starting simulation loop..." << std::endl;
    
    // Simulation parameters
    double timeStep = 0.016; // 60 FPS
    int generationLimit = 100;
    
    std::cout << "Starting NEAT organism simulation..." << std::endl;
    std::cout << "Press Ctrl+C to exit" << std::endl;
    
    // Main simulation loop
    while (world.getGeneration() < generationLimit) {
        std::cout << "Updating world (Generation " << world.getGeneration() 
                  << ", Organisms: " << world.getOrganismCount() << ")" << std::endl;
        
        // Update the world
        world.update(timeStep);
        
        // Display stats periodically
        if (static_cast<int>(world.getGeneration()) % 60 == 0) {
            std::cout << "Generation: " << world.getGeneration() 
                      << " | Organisms: " << world.getOrganismCount() 
                      << " | Food: " << world.getFoodCount() << std::endl;
        }
        
        // NOTE: Replace with proper game loop timing in prod
        std::this_thread::sleep_for(std::chrono::milliseconds(static_cast<int>(timeStep * 1000)));
    }
    
    std::cout << "Simulation complete after " << generationLimit << " generations." << std::endl;
    return 0;
}