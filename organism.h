#pragma once

#include "neat.h"
#include <vector>
#include <memory>

class World;

class Organism {
private:
    Genome brain;
    double x, y;           // Position
    double angle;          // Orientation
    double energy;         // Energy level
    double speed;
    double size;
    bool alive;
    
    // Sensor ranges
    static constexpr double visionRange = 100.0;
    static constexpr double foodSensorRange = 50.0;
    
public:
    Organism(const Genome& genome, double startX, double startY);
    Organism(int numInputs, int numOutputs, double startX, double startY);
    
    void update(World& world, double deltaTime);
    void render() const;
    
    // Getters
    double getX() const { return x; }
    double getY() const { return y; }
    double getAngle() const { return angle; }
    double getEnergy() const { return energy; }
    double getSize() const { return size; }
    bool isAlive() const { return alive; }
    Genome& getBrain() { return brain; }
    
    // Sensors
    std::vector<double> sense(const World& world);
    
    // Actions
    void move(double forward, double turn, double deltaTime);
    void eat(World& world);
    void reproduce(World& world);
    
    // Energy management
    void consumeEnergy(double amount);
    void addEnergy(double amount);
}; 