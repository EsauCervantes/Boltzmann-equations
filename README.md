# Boltzmann Equation Solver for Cannibal Dark Sectors

This repository contains a **numerical solver for coupled Boltzmann equations** arising in cannibal dark sector models produced via the freeze-in mechanism, as studied in

> **E. Cervantes**, *Freezing-in cannibal dark sectors* (2024).  
> [arXiv:2407.12104](https://arxiv.org/abs/2407.12104)

The code is written in **Mathematica** and is designed to be fully reproducible.

---

In our work we consider two scenarios. The first one is a dark sector composed only by a singlet unstable dark matter candidate $\phi$. This one is foun in cBE.wl. The second case is a dark sector composed of:

- A **dark matter particle** $S$ with self number changing reactions of the form $3\leftrightarrow 2$,
- A **mediator** $\phi$,

found in the cBE with mediator.wl file. The code contains production from Higgs annihilation/decay as well as production from electroweak gauge bosons. For details of the evolution of the system is governed by **coupled Boltzmann equations** for the comoving number densities see the paper. The output is given in terms of $x_{DM}$ and $Y$.
The code includes:

- **Freeze-in production** of $S$ and/or $\phi$ from the SM bath
- **Cannibalization processes** in the dark sector (e.g. \( \phi \phi \phi \leftrightarrow \phi \phi \))
- **Mediator decay**  $\phi \to \text{SM}$ and associated entropy injection
- **Hidden sector temperature evolution** \( T'(x) \)
- Proper treatment of **Hubble expansion** and changing relativistic degrees of freedom

---

## Features

### Coupled Boltzmann equations

We solve a system of ODE for the temperature and yield of both S and \(\phi\), where the collision terms \( \mathcal{C}_S \) and \( \mathcal{C}_\phi \) are included in the file collision integrals. The solver also includes:

- Freeze-in production from the SM
- $3 \leftrightarrow 2$ cannibal processes in the dark sector
- Decays and inverse decays of $\phi$

### Hidden sector temperature tracking

The solver also tracks the evolution of the **dark-sector temperature** \( T'(x) \), allowing for:

- Superadiabatic cooling/heating due to cannibalization
- Energy injection from late mediator decays

### Stiff ODE handling

The coupled Boltzmann equations are often **stiff**, especially in freeze-in regimes and during cannibal phases.  
We use Mathematica’s stiffness-aware ODE solvers with controlled precision and step sizes to obtain stable solutions.

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
