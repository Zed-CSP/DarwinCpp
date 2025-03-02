#include "organism.h"
#include "world.h"
#include <cmath>

Organism::Organism(const Genome& genome, double startX, double startY)
    : brain(genome), x(startX), y(startY), angle(0.0), 
      energy(100.0), speed(30.0), size(5.0), alive(true) {
}

Organism::Organism(int numInputs, int numOutputs, double startX, double startY)
    : x(startX), y(startY), angle(0.0), 
      energy(100.0), speed(30.0), size(5.0), alive(true) {
    brain.createInitialStructure(numInputs, numOutputs);
}

void Organism::update(World& world, double deltaTime) {
    // Lose energy over time (metabolism)
    consumeEnergy(deltaTime * 5.0);
    
    // Die if out of energy
    if (energy <= 0) {
        alive = false;
        return;
    }
    
    // Get sensory input
    std::vector<double> inputs = sense(world);
    
    // Process inputs through neural network
    std::vector<double> outputs = brain.feedForward(inputs);
    
    // Interpret outputs as actions
    if (outputs.size() >= 4) {
        // Movement: forward/backward and turning
        double moveForward = outputs[0] * 2.0 - 1.0;  // Range: -1 to 1
        double turn = outputs[1] * 2.0 - 1.0;         // Range: -1 to 1
        move(moveForward, turn, deltaTime);
        
        // Eating action
        if (outputs[2] > 0.5) {
            eat(world);
        }
        
        // Reproduction action
        if (outputs[3] > 0.5 && energy > 50.0) {
            reproduce(world);
        }
    }
}

std::vector<double> Organism::sense(const World& world) {
    std::vector<double> sensorData;
    
    // Add position and orientation data
    sensorData.push_back(x / world.getWidth());  // Normalized x position
    sensorData.push_back(y / world.getHeight()); // Normalized y position
    sensorData.push_back(angle / (2 * M_PI));    // Normalized angle
    sensorData.push_back(energy / 100.0);        // Normalized energy level
    
    // Food sensors - detect closest food in different directions
    auto foodItems = world.getFoodInRange(x, y, foodSensorRange);
    
    // Initialize food detection in 8 directions (45-degree segments)
    std::vector<double> foodDistances(8, 1.0);  // 1.0 means no food detected (normalized)
    
    for (auto food : foodItems) {
        double dx = food->x - x;
        double dy = food->y - y;
        double distance = std::sqrt(dx*dx + dy*dy);
        
        // Skip if too far
        if (distance > foodSensorRange) continue;
        
        // Calculate angle to food relative to organism's orientation
        double angleToFood = std::atan2(dy, dx) - angle;
        // Normalize to [0, 2π]
        while (angleToFood < 0) angleToFood += 2 * M_PI;
        while (angleToFood >= 2 * M_PI) angleToFood -= 2 * M_PI;
        
        // Determine which of the 8 direction segments this food is in
        int segment = static_cast<int>(angleToFood / (M_PI / 4)) % 8;
        
        // Normalize distance (closer food = higher value)
        double normalizedDistance = 1.0 - (distance / foodSensorRange);
        
        // Update if this food is closer than previously detected food in this segment
        if (normalizedDistance > (1.0 - foodDistances[segment])) {
            foodDistances[segment] = 1.0 - normalizedDistance;
        }
    }
    
    // Add food sensor data
    sensorData.insert(sensorData.end(), foodDistances.begin(), foodDistances.end());
    
    // Organism sensors - detect other organisms
    auto nearbyOrganisms = world.getOrganismsInRange(x, y, visionRange);
    
    // Initialize organism detection in 8 directions
    std::vector<double> organismDistances(8, 1.0);
    
    for (auto organism : nearbyOrganisms) {
        // Skip self
        if (organism == this) continue;
        
        double dx = organism->getX() - x;
        double dy = organism->getY() - y;
        double distance = std::sqrt(dx*dx + dy*dy);
        
        // Skip if too far
        if (distance > visionRange) continue;
        
        // Calculate angle to organism relative to this organism's orientation
        double angleToOrganism = std::atan2(dy, dx) - angle;
        // Normalize to [0, 2π]
        while (angleToOrganism < 0) angleToOrganism += 2 * M_PI;
        while (angleToOrganism >= 2 * M_PI) angleToOrganism -= 2 * M_PI;
        
        // Determine which of the 8 direction segments this organism is in
        int segment = static_cast<int>(angleToOrganism / (M_PI / 4)) % 8;
        
        // Normalize distance (closer organism = higher value)
        double normalizedDistance = 1.0 - (distance / visionRange);
        
        // Update if this organism is closer than previously detected ones in this segment
        if (normalizedDistance > (1.0 - organismDistances[segment])) {
            organismDistances[segment] = 1.0 - normalizedDistance;
        }
    }
    
    // Add organism sensor data
    sensorData.insert(sensorData.end(), organismDistances.begin(), organismDistances.end());
    
    // Wall sensors - detect distance to world boundaries
    double leftDist = x / world.getWidth();
    double rightDist = (world.getWidth() - x) / world.getWidth();
    double topDist = y / world.getHeight();
    double bottomDist = (world.getHeight() - y) / world.getHeight();
    
    sensorData.push_back(leftDist);
    sensorData.push_back(rightDist);
    sensorData.push_back(topDist);
    sensorData.push_back(bottomDist);
    
    return sensorData;
}

void Organism::move(double forward, double turn, double deltaTime) {
    // Update angle based on turn input
    angle += turn * deltaTime * 3.0;
    
    // Normalize angle to [0, 2π]
    while (angle < 0) angle += 2 * M_PI;
    while (angle >= 2 * M_PI) angle -= 2 * M_PI;
    
    // Calculate movement vector
    double dx = std::cos(angle) * forward * speed * deltaTime;
    double dy = std::sin(angle) * forward * speed * deltaTime;
    
    // Update position
    x += dx;
    y += dy;
    
    // Consume energy for movement (proportional to speed)
    double movementCost = std::abs(forward) * deltaTime * 2.0;
    consumeEnergy(movementCost);
}

void Organism::eat(World& world) {
    // Try to consume food within eating radius
    double energyGained = 0.0;
    if (world.consumeFood(x, y, size * 1.5, energyGained)) {
        addEnergy(energyGained);
    }
}

void Organism::reproduce(World& world) {
    // Create a child with a mutated brain
    Genome childBrain = brain;
    childBrain.mutate();
    
    // Calculate spawn position slightly offset from parent
    double spawnDistance = size * 2.0;
    double spawnAngle = angle + M_PI; // Opposite direction
    double childX = x + std::cos(spawnAngle) * spawnDistance;
    double childY = y + std::sin(spawnAngle) * spawnDistance;
    
    // Create the child organism
    auto child = std::make_unique<Organism>(childBrain, childX, childY);
    
    // Transfer some energy to the child
    double energyTransfer = energy * 0.3;
    consumeEnergy(energyTransfer);
    child->addEnergy(energyTransfer * 0.8); // Some energy is lost in reproduction
    
    // Add the child to the world
    world.addOrganism(std::move(child));
}

void Organism::consumeEnergy(double amount) {
    energy -= amount;
    if (energy < 0) energy = 0;
}

void Organism::addEnergy(double amount) {
    energy += amount;
    // Cap maximum energy
    if (energy > 200.0) energy = 200.0;
}

void Organism::render() const {
    // TODO: Implement rendering with graphics library
    // ALPHA SOLUTION: Draw organism as a triangle pointing in direction of angle
    // Size should be proportional to energy level
}