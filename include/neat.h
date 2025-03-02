#pragma once

#include <vector>
#include <random>
#include <cmath>
#include <memory>
#include <algorithm>
#include <unordered_map>
#include <string>

// Forward declarations
class Neuron;
class Connection;
class Genome;
class Species;

// Innovation tracking for NEAT
class InnovationHistory {
private:
    static int nextInnovationNumber;
    static std::unordered_map<std::string, int> innovationMap;
    
public:
    static int getInnovation(int fromNode, int toNode);
    static void reset();
};

// Connection between neurons
class Connection {
public:
    int fromNode;
    int toNode;
    double weight;
    bool enabled;
    int innovationNumber;
    
    Connection(int from, int to, double w, bool enabled);
    Connection(const Connection& other);
};

// Neuron in the network
class Neuron {
public:
    int id;
    double value;
    std::vector<int> incomingConnections;
    
    enum Type { INPUT, HIDDEN, OUTPUT };
    Type type;
    
    Neuron(int id, Type t);
    double activate(double x) const;
};

// Neural network genome for NEAT
class Genome {
private:
    std::mt19937 rng;
    
public:
    std::vector<Neuron> neurons;
    std::vector<Connection> connections;
    double fitness;
    
    // Configuration parameters
    double mutateWeightChance = 0.8;
    double weightMutationStep = 0.1;
    double addNodeChance = 0.03;
    double addConnectionChance = 0.05;
    double disableConnectionChance = 0.01;
    
    Genome();
    void createInitialStructure(int numInputs, int numOutputs);
    std::vector<double> feedForward(const std::vector<double>& inputs);
    void mutate();
    
    static Genome crossover(const Genome& parent1, const Genome& parent2);
    static double geneticDistance(const Genome& genome1, const Genome& genome2);
};

// Species for NEAT speciation
class Species {
public:
    std::vector<Genome> members;
    Genome representative;
    double bestFitness;
    int staleness;
    
    Species(const Genome& rep);
    bool isCompatible(const Genome& genome, double threshold);
    void addMember(const Genome& genome);
    void updateStats();
};

// NEAT algorithm implementation
class NEAT {
private:
    std::vector<Species> species;
    std::vector<Genome> population;
    int populationSize;
    int numInputs;
    int numOutputs;
    double speciationThreshold;
    std::mt19937 rng;
    
public:
    NEAT(int popSize, int inputs, int outputs, double threshold = 3.0);
    void speciate();
    void evolve();
    std::vector<Genome>& getPopulation();
}; 