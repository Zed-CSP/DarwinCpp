#include "neat.h"
#include <iostream>
#include <unordered_map>
#include <cmath>

// Initialize static members
int InnovationHistory::nextInnovationNumber = 1;
std::unordered_map<std::string, int> InnovationHistory::innovationMap;

int InnovationHistory::getInnovation(int fromNode, int toNode) {
    std::string key = std::to_string(fromNode) + "-" + std::to_string(toNode);
    if (innovationMap.find(key) == innovationMap.end()) {
        innovationMap[key] = nextInnovationNumber++;
    }
    return innovationMap[key];
}

void InnovationHistory::reset() {
    nextInnovationNumber = 1;
    innovationMap.clear();
}

Connection::Connection() 
    : fromNode(0), toNode(0), weight(0.0), enabled(false), innovationNumber(0) {}

Connection::Connection(int from, int to, double w, bool enabled) 
    : fromNode(from), toNode(to), weight(w), enabled(enabled) {
    innovationNumber = InnovationHistory::getInnovation(from, to);
}

Connection::Connection(const Connection& other) 
    : fromNode(other.fromNode), toNode(other.toNode), 
      weight(other.weight), enabled(other.enabled),
      innovationNumber(other.innovationNumber) {}

Neuron::Neuron(int id, Type t) : id(id), value(0.0), type(t) {}

double Neuron::activate(double x) const {
    // Sigmoid activation function
    return 1.0 / (1.0 + exp(-x));
}

Genome::Genome() : fitness(0.0) {
    std::random_device rd;
    rng = std::mt19937(rd());
}

void Genome::createInitialStructure(int numInputs, int numOutputs) {
    // Create input neurons
    for (int i = 0; i < numInputs; i++) {
        neurons.push_back(Neuron(i, Neuron::INPUT));
    }
    
    // Create output neurons
    for (int i = 0; i < numOutputs; i++) {
        neurons.push_back(Neuron(numInputs + i, Neuron::OUTPUT));
    }
    
    // Connect each input to each output with random weights
    std::uniform_real_distribution<double> weightDist(-1.0, 1.0);
    
    for (int i = 0; i < numInputs; i++) {
        for (int o = 0; o < numOutputs; o++) {
            double weight = weightDist(rng);
            connections.push_back(Connection(i, numInputs + o, weight, true));
        }
    }
}

// Feed forward through the network
std::vector<double> Genome::feedForward(const std::vector<double>& inputs) {
    // Reset all neuron values
    for (auto& neuron : neurons) {
        neuron.value = 0.0;
    }
    
    // Set input values
    for (size_t i = 0; i < inputs.size() && i < neurons.size(); i++) {
        if (neurons[i].type == Neuron::INPUT) {
            neurons[i].value = inputs[i];
        }
    }
    
    // For each neuron, calculate its value based on inputs and connections
    std::unordered_map<int, double> neuronValues;
    
    // Initialize input values in the map
    for (size_t i = 0; i < inputs.size() && i < neurons.size(); i++) {
        if (neurons[i].type == Neuron::INPUT) {
            neuronValues[neurons[i].id] = inputs[i];
        }
    }
    
    // Process all connections to calculate hidden and output neuron values
    for (const auto& connection : connections) {
        if (!connection.enabled) continue;
        
        // Only process if the source neuron has a value
        if (neuronValues.find(connection.fromNode) != neuronValues.end()) {
            // Add the weighted input to the target neuron
            if (neuronValues.find(connection.toNode) == neuronValues.end()) {
                neuronValues[connection.toNode] = 0.0;
            }
            neuronValues[connection.toNode] += neuronValues[connection.fromNode] * connection.weight;
        }
    }
    
    // Apply activation function to hidden and output neurons
    for (auto& neuron : neurons) {
        if (neuron.type != Neuron::INPUT && neuronValues.find(neuron.id) != neuronValues.end()) {
            neuronValues[neuron.id] = neuron.activate(neuronValues[neuron.id]);
        }
    }
    
    // Collect output values
    std::vector<double> outputs;
    for (const auto& neuron : neurons) {
        if (neuron.type == Neuron::OUTPUT && neuronValues.find(neuron.id) != neuronValues.end()) {
            outputs.push_back(neuronValues[neuron.id]);
        }
    }
    
    return outputs;
}

// Mutate the genome
void Genome::mutate() {
    std::uniform_real_distribution<double> rand(0.0, 1.0);
    std::uniform_real_distribution<double> weightDist(-1.0, 1.0);
    
    // Mutate weights
    if (rand(rng) < mutateWeightChance) {
        for (auto& connection : connections) {
            if (rand(rng) < 0.5) {
                // Perturb weight
                connection.weight += (rand(rng) * 2.0 - 1.0) * weightMutationStep;
                // Clamp weight to [-1, 1]
                connection.weight = std::max(-1.0, std::min(1.0, connection.weight));
            } else {
                // Assign new random weight
                connection.weight = weightDist(rng);
            }
        }
    }
    
    // Add a new node by splitting an existing connection
    if (rand(rng) < addNodeChance && !connections.empty()) {
        // Choose a random enabled connection
        std::vector<int> enabledConnections;
        for (size_t i = 0; i < connections.size(); i++) {
            if (connections[i].enabled) {
                enabledConnections.push_back(i);
            }
        }
        
        if (!enabledConnections.empty()) {
            int connIndex = enabledConnections[rand(rng) * enabledConnections.size()];
            Connection& conn = connections[connIndex];
            
            // Disable the original connection
            conn.enabled = false;
            
            // Create a new node
            int newNodeId = neurons.size();
            neurons.push_back(Neuron(newNodeId, Neuron::HIDDEN));
            
            // Create two new connections
            // In -> New with weight 1.0
            connections.push_back(Connection(conn.fromNode, newNodeId, 1.0, true));
            // New -> Out with original weight
            connections.push_back(Connection(newNodeId, conn.toNode, conn.weight, true));
        }
    }
    
    // Add a new connection
    if (rand(rng) < addConnectionChance) {
        // Try to find two unconnected neurons
        int attempts = 0;
        bool connectionAdded = false;
        
        while (!connectionAdded && attempts < 20) {
            attempts++;
            
            // Select source and target neurons
            int sourceIdx = rand(rng) * neurons.size();
            int targetIdx = rand(rng) * neurons.size();
            
            // Skip if connection would create a cycle (input -> input or output -> output)
            if (neurons[sourceIdx].type == Neuron::OUTPUT || neurons[targetIdx].type == Neuron::INPUT) {
                continue;
            }
            
            // Skip if connection already exists
            bool connectionExists = false;
            for (const auto& conn : connections) {
                if (conn.fromNode == neurons[sourceIdx].id && conn.toNode == neurons[targetIdx].id) {
                    connectionExists = true;
                    break;
                }
            }
            
            if (!connectionExists) {
                double weight = weightDist(rng);
                connections.push_back(Connection(neurons[sourceIdx].id, neurons[targetIdx].id, weight, true));
                connectionAdded = true;
            }
        }
    }
    
    // Randomly disable connections
    if (rand(rng) < disableConnectionChance && connections.size() > 1) {
        int connIndex = rand(rng) * connections.size();
        connections[connIndex].enabled = false;
    }
}

// Crossover between two genomes (static method)
Genome Genome::crossover(const Genome& parent1, const Genome& parent2) {
    // Parent with higher fitness passes more of its genes
    const Genome& betterParent = (parent1.fitness >= parent2.fitness) ? parent1 : parent2;
    const Genome& worseParent = (parent1.fitness >= parent2.fitness) ? parent2 : parent1;
    
    Genome child;
    
    // Map innovations to connections for both parents
    std::unordered_map<int, Connection> innovationMap1;
    std::unordered_map<int, Connection> innovationMap2;
    
    for (const auto& conn : betterParent.connections) {
        innovationMap1[conn.innovationNumber] = conn;
    }
    
    for (const auto& conn : worseParent.connections) {
        innovationMap2[conn.innovationNumber] = conn;
    }
    
    // Add all neurons from both parents
    std::unordered_map<int, bool> addedNeurons;
    
    for (const auto& neuron : betterParent.neurons) {
        child.neurons.push_back(neuron);
        addedNeurons[neuron.id] = true;
    }
    
    for (const auto& neuron : worseParent.neurons) {
        if (addedNeurons.find(neuron.id) == addedNeurons.end()) {
            child.neurons.push_back(neuron);
            addedNeurons[neuron.id] = true;
        }
    }
    
    // Inherit connections
    for (const auto& entry : innovationMap1) {
        int innovation = entry.first;
        const Connection& conn = entry.second;
        
        if (innovationMap2.find(innovation) != innovationMap2.end()) {
            // Matching gene - randomly inherit from either parent
            if (std::uniform_real_distribution<double>(0, 1)(child.rng) < 0.5) {
                child.connections.push_back(conn);
            } else {
                child.connections.push_back(innovationMap2[innovation]);
            }
        } else {
            // Disjoint or excess gene from the better parent
            child.connections.push_back(conn);
        }
    }
    
    return child;
}

// Calculate genetic distance between genomes
double Genome::geneticDistance(const Genome& genome1, const Genome& genome2) {
    double c1 = 1.0; // Weight for excess genes
    double c2 = 1.0; // Weight for disjoint genes
    double c3 = 0.4; // Weight for weight differences
    
    // Map innovations to connections for both genomes
    std::unordered_map<int, Connection> innovations1;
    std::unordered_map<int, Connection> innovations2;
    
    for (const auto& conn : genome1.connections) {
        innovations1[conn.innovationNumber] = conn;
    }
    
    for (const auto& conn : genome2.connections) {
        innovations2[conn.innovationNumber] = conn;
    }
    
    // Count disjoint and excess genes
    int disjoint = 0;
    int excess = 0;
    double weightDiff = 0.0;
    int matching = 0;
    
    // Find max innovation number
    int maxInnovation1 = 0;
    int maxInnovation2 = 0;
    
    for (const auto& entry : innovations1) {
        maxInnovation1 = std::max(maxInnovation1, entry.first);
    }
    
    for (const auto& entry : innovations2) {
        maxInnovation2 = std::max(maxInnovation2, entry.first);
    }
    
    // Identify disjoint, excess, and matching genes
    for (const auto& entry : innovations1) {
        int innovation = entry.first;
        if (innovations2.find(innovation) != innovations2.end()) {
            // Matching gene
            matching++;
            weightDiff += std::abs(entry.second.weight - innovations2[innovation].weight);
        } else if (innovation <= maxInnovation2) {
            // Disjoint gene
            disjoint++;
        } else {
            // Excess gene
            excess++;
        }
    }
    
    for (const auto& entry : innovations2) {
        int innovation = entry.first;
        if (innovations1.find(innovation) == innovations1.end()) {
            if (innovation <= maxInnovation1) {
                // Disjoint gene
                disjoint++;
            } else {
                // Excess gene
                excess++;
            }
        }
    }
    
    // Calculate average weight difference for matching genes
    double avgWeightDiff = (matching > 0) ? (weightDiff / matching) : 0;
    
    // Calculate N (normalization factor)
    int N = std::max(innovations1.size(), innovations2.size());
    if (N < 20) N = 1; // Small genomes aren't normalized
    
    // Calculate and return genetic distance
    return (c1 * excess / N) + (c2 * disjoint / N) + (c3 * avgWeightDiff);
}

// Species implementation
Species::Species(const Genome& rep) : representative(rep), bestFitness(0.0), staleness(0) {}

bool Species::isCompatible(const Genome& genome, double threshold) {
    return Genome::geneticDistance(representative, genome) < threshold;
}

void Species::addMember(const Genome& genome) {
    members.push_back(genome);
}

void Species::updateStats() {
    if (members.empty()) return;
    
    // Find best fitness in the species
    double maxFitness = members[0].fitness;
    for (const auto& genome : members) {
        if (genome.fitness > maxFitness) {
            maxFitness = genome.fitness;
        }
    }
    
    // Check if we've improved
    if (maxFitness > bestFitness) {
        bestFitness = maxFitness;
        staleness = 0;
    } else {
        staleness++;
    }
}

// NEAT implementation
NEAT::NEAT(int popSize, int inputs, int outputs, double threshold)
    : populationSize(popSize), numInputs(inputs), numOutputs(outputs),
      speciationThreshold(threshold) {
    
    std::random_device rd;
    rng = std::mt19937(rd());
    
    // Create initial population
    for (int i = 0; i < popSize; i++) {
        Genome genome;
        genome.createInitialStructure(inputs, outputs);
        population.push_back(genome);
    }
    
    // Initial speciation
    speciate();
}

void NEAT::speciate() {
    // Clear existing species
    for (auto& species : species) {
        species.members.clear();
    }
    
    // Assign each genome to a species
    for (auto& genome : population) {
        bool found = false;
        
        for (auto& species : species) {
            if (species.isCompatible(genome, speciationThreshold)) {
                species.addMember(genome);
                found = true;
                break;
            }
        }
        
        if (!found) {
            // Create a new species with this genome as representative
            species.push_back(Species(genome));
            species.back().addMember(genome);
        }
    }
    
    // Remove empty species
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [](const Species& s) { return s.members.empty(); }),
        species.end()
    );
}

void NEAT::evolve() {
    // Calculate fitness for each species
    for (auto& s : species) {
        s.updateStats();
    }
    
    // Remove stale species
    species.erase(
        std::remove_if(species.begin(), species.end(),
            [](const Species& s) { return s.staleness > 15; }),
        species.end()
    );
    
    // Calculate total adjusted fitness
    double totalAdjustedFitness = 0.0;
    for (auto& s : species) {
        double speciesAvgFitness = 0.0;
        for (const auto& genome : s.members) {
            speciesAvgFitness += genome.fitness;
        }
        speciesAvgFitness /= s.members.size();
        totalAdjustedFitness += speciesAvgFitness;
    }
    
    // Create next generation
    std::vector<Genome> newPopulation;
    
    for (auto& s : species) {
        // Calculate number of offspring
        double speciesAvgFitness = 0.0;
        for (const auto& genome : s.members) {
            speciesAvgFitness += genome.fitness;
        }
        speciesAvgFitness /= s.members.size();
        
        int offspring = static_cast<int>((speciesAvgFitness / totalAdjustedFitness) * populationSize);
        
        // Ensure at least one offspring if the species is doing well
        if (offspring == 0 && s.bestFitness > 0) {
            offspring = 1;
        }
        
        // Sort members by fitness
        std::sort(s.members.begin(), s.members.end(),
            [](const Genome& a, const Genome& b) { return a.fitness > b.fitness; });
        
        // Keep the champion unchanged
        if (!s.members.empty() && offspring > 0) {
            newPopulation.push_back(s.members[0]);
            offspring--;
        }
        
        // Create the rest through crossover and mutation
        for (int i = 0; i < offspring; i++) {
            if (s.members.size() == 1) {
                // Only one parent, clone and mutate
                Genome child = s.members[0];
                child.mutate();
                newPopulation.push_back(child);
            } else {
                // Select two parents for crossover
                std::uniform_int_distribution<int> dist(0, std::min(4, (int)s.members.size() - 1));
                int parent1Idx = dist(rng);
                int parent2Idx;
                do {
                    parent2Idx = dist(rng);
                } while (parent2Idx == parent1Idx);
                
                Genome child = Genome::crossover(s.members[parent1Idx], s.members[parent2Idx]);
                child.mutate();
                newPopulation.push_back(child);
            }
        }
    }
    
    // Fill the rest of the population with new random genomes if needed
    while (newPopulation.size() < populationSize) {
        Genome genome;
        genome.createInitialStructure(numInputs, numOutputs);
        newPopulation.push_back(genome);
    }
    
    // Replace the current population
    population = newPopulation;
    
    // Speciate the new population
    speciate();
}

std::vector<Genome>& NEAT::getPopulation() {
    return population;
} 