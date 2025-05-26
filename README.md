

# # opto-sim: Quantum Optical Simulator for BB84 Protocol (Physical Abstraction Layer)

This is an early-stage simulator for exploring the **BB84 quantum key distribution (QKD)** protocol at the **physical abstraction layer**, using optical components modeled with mathematical accuracy based on real-world research (within computational and practical limits). .

**Work in Progress**  
This project is in active development and not yet production-ready. However, it serves as a showcase of my interest in quantum technologies and my focus on modeling real-world optical phenomena in Python.

## Project Overview

This simulator aims to mimic the behavior of physical optical components used in QKD experiments—specifically BB84 protocol setups—through mathematical models. It’s intended as a learning tool and experimental platform for simulating realistic noise conditions, signal degradation, and quantum encoding/decoding via optics.

### Current Features

- **Laser Source**  
  - Solid-state lasers modeled after Nd:YAG and Er:YAG systems.

- **Phase Modulators**  
  - Supports X-cut and Y-cut, Z-propagating configurations.

- **Avalanche Photodiodes (APDs)**  
  - Models include shot noise, thermal noise, and random Gaussian/Poisson noise.

- **Optical Components**  
  - Basic tools like beam splitters, couplers, and fiber optic cables.

- **Polarimetry**  
  - Polarimeter with calculation of **Stokes parameters** for polarization state analysis.
  - Poincare Sphere  

### Mathematical Modeling

All models aim to closely reflect published experimental setups and component behaviors found in academic literature. While accuracy is a priority, simplifications are made when necessary to keep computations tractable.

## Dependencies

- `numpy`
- `scipy`
- `matplotlib`

Install dependencies with:

```bash
pip install -r requirements.txt
