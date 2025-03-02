#include "neat.h"
#include <iostream>

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
    
    // Process neurons in order (assuming they're topologically sorted)
    for (auto& connection : connections) {
        if (connection.enabled) {
            neurons[connection.toNode].value += 
                neurons[connection.fromNode].value * connection.weight;
        }
    }
    
    // Activate and collect outputs
    std::vector<double> outputs;
    for (auto& neuron : neurons) {
        if (neuron.type == Neuron::OUTPUT) {
            neuron.value = neuron.activate(neuron.value);
            outputs.push_back(neuron.value);
        }
        else if (neuron.type == Neuron::HIDDEN) {
            neuron.value = neuron.activate(neuron.value);
        }
    }
    
    return outputs;
}

void Genome::mutate() {
    std::uniform_real_distribution<double> chance(0.0, 1.0);
    std::uniform_real_distribution<double> weightDist(-1.0, 1.0);
    
    // Mutate weights
    if (chance(rng) < mutateWeightChance) {
        for (auto& connection : connections) {
            if (chance(rng) < 0.9) {
                // Perturb weight
                connection.weight += weightDist(rng) * weightMutationStep;
            } else {
                // Assign new random weight
                connection.weight = weightDist(rng);
            }
        }
    }
    
    // Add a new node
    if (chance(rng) < addNodeChance && !connections.empty()) {
        // Select a random connection to split
        std::uniform_int_distribution<int> connDist(0, connections.size() - 1);
        int connIndex = connDist(rng);
        
        // Disable the selected connection
        connections[connIndex].enabled = false;
        
        // Create a new node
        int newNodeId = neurons.size();
        neurons.push_back(Neuron(newNodeId, Neuron::HIDDEN));
        
        // Add two new connections
        connections.push_back(Connection(
            connections[connIndex].fromNode, 
            newNodeId, 
            1.0, 
            true
        ));
        
        connections.push_back(Connection(
            newNodeId, 
            connections[connIndex].toNode, 
            connections[connIndex].weight, 
            true
        ));
    }
    
    // Add a new connection
    if (chance(rng) < addConnectionChance && neurons.size() > 1) {
        // Try to find two unconnected neurons
        for (int attempts = 0; attempts < 20; attempts++) {
            std::uniform_int_distribution<int> neuronDist(0, neurons.size() - 1);
            int fromIndex = neuronDist(rng);
            int toIndex = neuronDist(rng);
            
            // Ensure we're not connecting output -> input or output -> output
            if (neurons[fromIndex].type == Neuron::OUTPUT) continue;
            
            // Ensure we're not connecting input -> input
            if (neurons[fromIndex].type == Neuron::INPUT && 
                neurons[toIndex].type == Neuron::INPUT) continue;
            
            // Check if connection already exists
            bool connectionExists = false;
            for (const auto& conn : connections) {
                if (conn.fromNode == fromIndex && conn.toNode == toIndex) {
                    connectionExists = true;
                    break;
                }
            }
            
            if (!connectionExists) {
                connections.push_back(Connection(
                    fromIndex, 
                    toIndex, 
                    weightDist(rng), 
                    true
                ));
                break;
            }
        }
    }
    
    // Randomly disable connections
    if (chance(rng) < disableConnectionChance && connections.size() > 1) {
        std::uniform_int_distribution<int> connDist(0, connections.size() - 1);
        int connIndex = connDist(rng);
        connections[connIndex].enabled = !connections[connIndex].enabled;
    }
}

Genome Genome::crossover(const Genome& parent1, const Genome& parent2) {
    // Assume parent1 is the more fit parent
    const Genome& morefit = (parent1.fitness >= parent2.fitness) ? parent1 : parent2;
    const Genome& lessfit = (parent1.fitness >= parent2.fitness) ? parent2 : parent1;
    
    Genome child;
    
    // Copy all neurons from the more fit parent
    child.neurons = morefit.neurons;
    
    // Crossover connections
    std::unordered_map<int, Connection> innovationMap1;
    std::unordered_map<int, Connection> innovationMap2;
    
    for (const auto& conn : morefit.connections) {
        innovationMap1[conn.innovationNumber] = conn;
    }
    
    for (const auto& conn : lessfit.connections) {
        innovationMap2[conn.innovationNumber] = conn;
    }
    
    // Inherit matching genes randomly, disjoint/excess from more fit parent
    for (const auto& pair : innovationMap1) {
        int innovation = pair.first;
        
        if (innovationMap2.find(innovation) != innovationMap2.end()) {
            // Matching gene - randomly choose from either parent
            std::uniform_int_distribution<int> dist(0, 1);
            if (dist(child.rng) == 0) {
                child.connections.push_back(pair.second);
            } else {
                child.connections.push_back(innovationMap2[innovation]);
            }
        } else {
            // Disjoint/excess gene from more fit parent
            child.connections.push_back(pair.second);
        }
    }
    
    return child;
}

double Genome::geneticDistance(const Genome& genome1, const Genome& genome2) {
    const double c1 = 1.0; // Coefficient for excess genes
    const double c2 = 1.0; // Coefficient for disjoint genes
    const double c3 = 0.4; // Coefficient for weight differences
    
    std::unordered_map<int, Connection> innovationMap1;
    std::unordered_map<int, Connection> innovationMap2;
    
    for (const auto& conn : genome1.connections) {
        innovationMap1[conn.innovationNumber] = conn;
    }
    
    for (const auto& conn : genome2.connections) {
        innovationMap2[conn.innovationNumber] = conn;
    }
    
    int disjoint = 0;
    int excess = 0;
    double weightDiff = 0.0;
    int matching = 0;
    
    int maxInnov1 = 0;
    int maxInnov2 = 0;
    
    for (const auto& pair : innovationMap1) {
        maxInnov1 = std::max(maxInnov1, pair.first);
    }
    
    for (const auto& pair : innovationMap2) {
        maxInnov2 = std::max(maxInnov2, pair.first);
    }
    
    // Count matching, disjoint, and excess genes
    for (const auto& pair : innovationMap1) {
        int innovation = pair.first;
        
        if (innovationMap2.find(innovation) != innovationMap2.end()) {
            // Matching gene
            matching++;
            weightDiff += std::abs(pair.second.weight - innovationMap2[innovation].weight);
        } else if (innovation <= maxInnov2) {
            // Disjoint gene
            disjoint++;
        } else {
            // Excess gene
            excess++;
        }
    }
    
    for (const auto& pair : innovationMap2) {
        int innovation = pair.first;
        
        if (innovationMap1.find(innovation) == innovationMap1.end() && 
            innovation <= maxInnov1) {
            // Disjoint gene from genome2
            disjoint++;
        } else if (innovation > maxInnov1) {
            // Excess gene from genome2
            excess++;
        }
    }
    
    // Normalize by size of larger genome
    int N = std::max(innovationMap1.size(), innovationMap2.size());
    if (N < 1) N = 1; // Avoid division by zero
    
    double avgWeightDiff = (matching > 0) ? weightDiff / matching : 0;
    
    return (c1 * excess / N) + (c2 * disjoint / N) + (c3 * avgWeightDiff);
}

Species::Species(const Genome& rep) : representative(rep), bestFitness(0.0), staleness(0) {}

bool Species::isCompatible(const Genome& genome, double threshold) {
    return Genome::geneticDistance(representative, genome) < threshold;
}

void Species::addMember(const Genome& genome) {
    members.push_back(genome);
}

void Species::updateStats() {
    if (members.empty()) return;
    
    // Find the best fitness in this generation
    double maxFitness = members[0].fitness;
    for (const auto& genome : members) {
        maxFitness = std::max(maxFitness, genome.fitness);
    }
    
    // Check if we've improved
    if (maxFitness > bestFitness) {
        bestFitness = maxFitness;
        staleness = 0;
    } else {
        staleness++;
    }
    
    // Select a new random representative
    if (!members.empty()) {
        std::random_device rd;
        std::mt19937 rng(rd());
        std::uniform_int_distribution<int> dist(0, members.size() - 1);
        representative = members[dist(rng)];
    }
}

NEAT::NEAT(int popSize, int inputs, int outputs, double threshold) 
    : populationSize(popSize), numInputs(inputs), numOutputs(outputs),
      speciationThreshold(threshold) {
    
    std::random_device rd;
    rng = std::mt19937(rd());
    
    // Initialize population with basic genomes
    for (int i = 0; i < populationSize; i++) {
        Genome genome;
        genome.createInitialStructure(numInputs, numOutputs);
        population.push_back(genome);
    }
}

void NEAT::speciate() {
    // Clear old species members but keep the species
    for (auto& s : species) {
        s.members.clear();
    }
    
    // Assign each genome to a species
    for (auto& genome : population) {
        bool found = false;
        
        for (auto& s : species) {
            if (s.isCompatible(genome, speciationThreshold)) {
                s.addMember(genome);
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
    
    // Update species stats
    for (auto& s : species) {
        s.updateStats();
    }
}

void NEAT::evolve() {
    // Speciate the population
    speciate();
    
    // Calculate adjusted fitness and total adjusted fitness
    double totalAdjustedFitness = 0.0;
    for (auto& s : species) {
        for (auto& genome : s.members) {
            // Adjust fitness by sharing with species members
            genome.fitness /= s.members.size();
            totalAdjustedFitness += genome.fitness;
        }
    }
    
    // Create new population
    std::vector<Genome> newPopulation;
    
    // Elitism: keep the best genome from each species
    for (auto& s : species) {
        if (s.members.empty()) continue;
        
        // Find best genome in species
        Genome* best = &s.members[0];
        for (auto& genome : s.members) {
            if (genome.fitness > best->fitness) {
                best = &genome;
            }
        }
        
        // Add to new population
        newPopulation.push_back(*best);
    }
    
    // Fill the rest of the population with offspring
    while (newPopulation.size() < populationSize) {
        // Select a species based on adjusted fitness
        Species* selectedSpecies = nullptr;
        double r = std::uniform_real_distribution<double>(0, totalAdjustedFitness)(rng);
        double sum = 0.0;
        
        for (auto& s : species) {
            double speciesFitness = 0.0;
            for (auto& genome : s.members) {
                speciesFitness += genome.fitness;
            }
            
            sum += speciesFitness;
            if (sum >= r) {
                selectedSpecies = &s;
                break;
            }
        }
        
        if (!selectedSpecies || selectedSpecies->members.empty()) {
            // Fallback if no species was selected
            if (!species.empty()) {
                selectedSpecies = &species[0];
            } else {
                // Create a new random genome
                Genome newGenome;
                newGenome.createInitialStructure(numInputs, numOutputs);
                newPopulation.push_back(newGenome);
                continue;
            }
        }
        
        // Create offspring
        Genome child;
        
        if (selectedSpecies->members.size() == 1 || 
            std::uniform_real_distribution<double>(0, 1)(rng) < 0.25) {
            // Asexual reproduction - clone and mutate
            std::uniform_int_distribution<int> dist(0, selectedSpecies->members.size() - 1);
            child = selectedSpecies->members[dist(rng)];
        } else {
            // Sexual reproduction - crossover
            std::uniform_int_distribution<int> dist(0, selectedSpecies->members.size() - 1);
            Genome parent1 = selectedSpecies->members[dist(rng)];
            Genome parent2 = selectedSpecies->members[dist(rng)];
            
            child = Genome::crossover(parent1, parent2);
        }
        
        // Mutate
        child.mutate();
        
        // Add to new population
        newPopulation.push_back(child);
    }
    
    // Replace old population
    population = newPopulation;
}

std::vector<Genome>& NEAT::getPopulation() {
    return population;
} 