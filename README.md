# 2D Electromagnetic Wave Propagation using FDTD

This project simulates the propagation of a 2D electromagnetic wave by solving Maxwell’s time-dependent equations using the Finite Difference Time Domain (FDTD) method. The algorithm is based on Yee’s grid and leapfrog scheme to model electric and magnetic field interactions in space and time.

## 🧪 Overview

- Solves time-dependent Maxwell's equations in 2D
- Uses Yee’s algorithm and leapfrog time-stepping
- Models both electric and magnetic field components (Ez, Hx, Hy)
- Supports two types of boundary conditions:
  - Absorbing (Mur’s method)
  - Dirichlet (perfect electric conductor)

## ⚙️ Parameters

- Grid size: `200 x 200`
- Time step: computed using Courant stability condition
- Spatial resolution: set to 10% of the minimum wavelength
- Source: unit step signal applied at the center of the grid
- Medium: vacuum (ε = ε₀, μ = μ₀)

## 💻 Implementation

- Language: C++/MATLAB
- Central finite differences are used to approximate derivatives
