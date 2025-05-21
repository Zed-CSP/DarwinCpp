#***WIP***

# NEAT Organism Simulation

A C++ implementation of evolving organisms using NeuroEvolution of Augmenting Topologies (NEAT). Scratch Built Neural Networks in C++.

## Overview

This project simulates a virtual ecosystem where organisms evolve neural networks to control their behavior. Each organism has a brain implemented as a NEAT neural network that evolves over generations through natural selection.

### Key features:
- Neural networks that evolve in both topology and weights
- Organisms that can move, eat, and reproduce
- Food resources that organisms compete for
- Speciation to protect innovation
- Visualization of the simulation (planned)

## Requirements
- C++17 compatible compiler
- Standard Library
- (Optional) SFML or SDL2 for visualization

## Project Structure<br>
├── include/ <br>
│   ├── neat.h           # NEAT algorithm implementation<br>
│   ├── organism.h       # Organism behavior and properties<br>
│   └── world.h          # Simulation environment<br>
├── src/<br>
│   ├── neat.cpp         # NEAT implementation<br>
│   ├── organism.cpp     # Organism implementation<br>
│   ├── world.cpp        # World implementation<br>
│   └── main.cpp         # Entry point<br>
└── README.md<br>

## How It Works

### NEAT Algorithm
The NEAT algorithm evolves neural networks by:
- Adding new neurons and connections
- Adjusting connection weights
- Tracking innovations through historical markers
- Speciation to protect innovation

### Organisms
Each organism has:
- Position and orientation in the world
- Energy level that depletes over time
- Neural network brain that controls behavior
- Sensors to detect food, other organisms, and walls
- Ability to move, eat, and reproduce

### Simulation
The simulation runs in generations:
- Organisms move around seeking food
- Organisms with successful behaviors survive and reproduce
- Neural networks evolve through mutation and crossover
- Over time, more complex and effective behaviors emerge

## Building and Running

```bash
# Clone the repository
git clone https://github.com/yourusername/neat-organism-simulation.git
cd neat-organism-simulation

# Build the project
mkdir build && cd build
cmake ..
make

# Run the simulation
./neat_simulation
```

## Future Improvements

- [ ] Add visualization using SFML/SDL2
- [ ] Implement predator-prey relationships
- [ ] Add different types of food with varying nutritional values
- [ ] Implement more complex environments with obstacles
- [ ] Add detailed statistics tracking and graphing
- [ ] Parallelize simulation for better performance

## References

- Stanley, K. O., & Miikkulainen, R. (2002). Evolving neural networks through augmenting topologies. Evolutionary computation, 10(2), 99-127.
- [NEAT-Python Documentation](https://neat-python.readthedocs.io/en/latest/neat_overview.html)

## License

This project is licensed under the MIT License - see the LICENSE file for details.
