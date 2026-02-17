# Existence proof of Hopf bifurcation for the Navier-Stokes equations

This repository contains the code corresponding to the proof of a Hopf bifurcation at nu=0.1148 in the Navier-Stokes equations, written by Lindsey van der Aalst and Jan Bouwe van den Berg. The code is heavily inspired by the following code:

J. B. van den Berg, M. Breden, J.-P. Lessard, and L. van Veen. MATLAB code for "Spontaneous periodic orbits in the Navier-Stokes flow", 2019. [https://www.math.vu.nl/~janbouwe/code/navierstokes/](https://www.math.vu.nl/~janbouwe/code/navierstokes/).

---

## How to run the proofs

Run the following file to reproduce the results in the paper:

- [`main_proof.m`](main_proof.m)

Note: The default parameters are currently set to run the "proof" quickly. However, no interval arithmetic is included and the number of modes is too low to yield a successful proof. To run the mathematically rigorous proof, the mode count must be increased, which requires at least 256 GB of RAM. Consequently, the full proof has not yet been executed or verified by the author due to hardware constraints.

---

## Requirements

To obtain rigorous proofs, [INTLAB](https://www.tuhh.de/ti3/rump/) is required.  
You can also run the proof without interval arithmetic, by setting the condition in line 4 of [`main_proof.m`](main_proof.m) to false (as is currently the status)

### Version

- MATLAB: **2024a** 

---

## Contact

For questions about the code, please contact:  
**Lindsey van der Aalst** — l.j.w.van.der.aalst@vu.nl
