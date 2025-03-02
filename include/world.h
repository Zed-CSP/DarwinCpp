#pragma once

#include <vector>
#include <memory>
#include <random>
#include "organism.h"

struct Food {
    double x, y;
    double energy;
    bool active;
    
    Food(double x, double y, double energy) 
        : x(x), y(y), energy(energy), active(true) {}
};

class World {
private:
    double width, height;
    std::vector<std::unique_ptr<Organism>> organisms;
    std::vector<Food> foodItems;
    
    int generation;
    int maxOrganisms;
    int initialFoodCount;
    double foodSpawnRate;
    double foodEnergy;
    
    std::mt19937 rng;
    
public:
    World(double width, double height, int initialOrganisms, int initialFood);
    
    void update(double deltaTime);
    void render() const;
    
    // Food management
    void spawnFood(int count);
    bool consumeFood(double x, double y, double radius, double& energyGained);
    
    // Organism management
    void addOrganism(std::unique_ptr<Organism> organism);
    void removeDeadOrganisms();
    
    // Environment queries
    std::vector<Food*> getFoodInRange(double x, double y, double range);
    std::vector<Organism*> getOrganismsInRange(double x, double y, double range);
    
    // Getters
    double getWidth() const { return width; }
    double getHeight() const { return height; }
    int getGeneration() const { return generation; }
    int getOrganismCount() const { return organisms.size(); }
    int getFoodCount() const;
    
    // Random position generator
    double randomX();
    double randomY();
}; 