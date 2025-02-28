#include "world.h"
#include <algorithm>
#include <cmath>

World::World(double width, double height, int initialOrganisms, int initialFood)
    : width(width), height(height), generation(0), 
      maxOrganisms(500), initialFoodCount(initialFood), 
      foodSpawnRate(0.5), foodEnergy(20.0) {
    
    // Initialize random number generator
    std::random_device rd;
    rng = std::mt19937(rd());
    
    // Spawn initial food
    spawnFood(initialFoodCount);
    
    // Create initial organisms
    for (int i = 0; i < initialOrganisms; i++) {
        double x = randomX();
        double y = randomY();
        
        // FIXME: Verify input/output counts match the sense() method implementation
        // Inputs: position(2), angle(1), energy(1), food sensors(8), organism sensors(8), wall sensors(4)
        // Outputs: move forward/back(1), turn(1), eat(1), reproduce(1)
        auto organism = std::make_unique<Organism>(24, 4, x, y);
        organisms.push_back(std::move(organism));
    }
}

void World::update(double deltaTime) {
    // Update all organisms
    for (auto& organism : organisms) {
        if (organism->isAlive()) {
            organism->update(*this, deltaTime);
        }
    }
    
    // Remove dead organisms
    removeDeadOrganisms();
    
    // Spawn new food
    if (std::uniform_real_distribution<double>(0, 1)(rng) < foodSpawnRate * deltaTime) {
        spawnFood(1);
    }
    
    // If all organisms are dead, start a new generation
    if (organisms.empty()) {
        generation++;
        // Create new random organisms
        for (int i = 0; i < initialFoodCount / 10; i++) {
            double x = randomX();
            double y = randomY();
            auto organism = std::make_unique<Organism>(24, 4, x, y);
            organisms.push_back(std::move(organism));
        }
        
        // Replenish food
        spawnFood(initialFoodCount);
    }
    
    // Cap maximum number of organisms
    if (organisms.size() > maxOrganisms) {
        // Sort by energy level (lowest first)
        std::sort(organisms.begin(), organisms.end(), 
            [](const auto& a, const auto& b) {
                return a->getEnergy() < b->getEnergy();
            });
        
        // Remove lowest energy organisms
        organisms.resize(maxOrganisms);
    }
}

void World::render() const {
    // TODO: Implement rendering with graphics library

}

void World::spawnFood(int count) {
    for (int i = 0; i < count; i++) {
        double x = randomX();
        double y = randomY();
        foodItems.push_back(Food(x, y, foodEnergy));
    }
}

bool World::consumeFood(double x, double y, double radius, double& energyGained) {
    for (auto& food : foodItems) {
        if (!food.active) continue;
        
        double dx = food.x - x;
        double dy = food.y - y;
        double distSquared = dx*dx + dy*dy;
        
        if (distSquared <= radius*radius) {
            // Food is within eating radius
            food.active = false;
            energyGained = food.energy;
            return true;
        }
    }
    
    // No food found within radius
    energyGained = 0.0;
    return false;
}

void World::addOrganism(std::unique_ptr<Organism> organism) {
    organisms.push_back(std::move(organism));
}

void World::removeDeadOrganisms() {
    organisms.erase(
        std::remove_if(organisms.begin(), organisms.end(),
            [](const auto& organism) { return !organism->isAlive(); }),
        organisms.end()
    );
}

std::vector<Food*> World::getFoodInRange(double x, double y, double range) {
    std::vector<Food*> result;
    double rangeSquared = range * range;
    
    for (auto& food : foodItems) {
        if (!food.active) continue;
        
        double dx = food.x - x;
        double dy = food.y - y;
        double distSquared = dx*dx + dy*dy;
        
        if (distSquared <= rangeSquared) {
            result.push_back(&food);
        }
    }
    
    return result;
}

std::vector<Organism*> World::getOrganismsInRange(double x, double y, double range) {
    std::vector<Organism*> result;
    double rangeSquared = range * range;
    
    for (auto& organism : organisms) {
        if (!organism->isAlive()) continue;
        
        double dx = organism->getX() - x;
        double dy = organism->getY() - y;
        double distSquared = dx*dx + dy*dy;
        
        if (distSquared <= rangeSquared) {
            result.push_back(organism.get());
        }
    }
    
    return result;
}

int World::getFoodCount() const {
    int count = 0;
    for (const auto& food : foodItems) {
        if (food.active) count++;
    }
    return count;
}

double World::randomX() {
    return std::uniform_real_distribution<double>(0, width)(rng);
}

double World::randomY() {
    return std::uniform_real_distribution<double>(0, height)(rng);
} 