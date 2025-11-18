# Boltzmann Equation Solver for Cannibal Dark Sectors

This repository contains a **numerical solver for coupled Boltzmann equations** arising in cannibal dark sector models produced via the freeze-in mechanism, as studied in

> **E. Cervantes**, *Freezing-in cannibal dark sectors* (2024).  
> [arXiv:2407.12104](https://arxiv.org/abs/2407.12104)

The code is written in **Mathematica** and is designed to be fully reproducible.

---

In our work we consider two scenarios: a dark sector composed only by a singlet unstable dark matter candidate $\phi$ (cBE.wl) and a dark sector composed of:

- A **dark matter particle** $S$ with self number changing reactions of the form $3\leftrightarrow 2$,
- A **mediator** $\phi$,

found in the cBE with mediator.wl file. The code contains production from Higgs annihilation/decay as well as production from electroweak gauge bosons. For details see the paper.
The code includes:

- **Freeze-in production** of $S$ and/or $\phi$ from the SM bath
- **Cannibalization processes** in the dark sector (e.g. $phi \phi \phi \leftrightarrow \phi \phi $)
- **Hidden sector temperature evolution and yield**

---

## Features

### Coupled Boltzmann equations

The Boltzmann equation is an integro-partial differential equation, whose solution is the probabilistic phase space distribution function. Solving it in full generality is computationally very expensive and usually not necessary. For instance, when the system is in thermal equilibrium, tracking the temperature and comoving number of particles is sufficient. This code does exactly this, and it is optimized to deal with stiffness during freeze-out. The solver also includes the $3 \leftrightarrow 2$ collision integral tabulated as a function of $m/T$.

### Hidden sector temperature tracking

The solver tracks the evolution of the **dark-sector temperature** $T'$, allowing for:

- Superadiabatic cooling/heating due to cannibalization
- Entropy exchange between the dark matter and the mediator

### Stiff ODE handling

The coupled Boltzmann equations are often **stiff** during cannibal phases.  
We use Mathematica’s ODE solvers with controlled precision and step sizes to obtain stable solutions.

### Relic abundance computation

From the late-time asymptotic value of \( Y_S(x) \), the code computes the dark matter relic abundance:
\[
\Omega_S h^2 = \frac{m_S\, s_0\, Y_S(x \to \infty)}{\rho_{\text{crit}}},
\]
where \( s_0 \) is today’s entropy density and \( \rho_{\text{crit}} \) is the critical density.

---

## Repository Structure

```text
Boltzmann-equations/
│
├── BoltzmannSolver.nb      # Main Mathematica notebook with the coupled ODE solver
├── parameters/             # (Optional) Benchmark parameter files
├── plots/                  # (Optional) Generated figures: Y(x), T'/T, relic density, etc.
└── README.md               # This file
