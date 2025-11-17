# Boltzmann Equation Solver for Cannibal Dark Sectors

This repository contains a **numerical solver for coupled Boltzmann equations** arising in cannibal dark sector models, as studied in

> **E. Cervantes**, *Freezing-in cannibal dark sectors* (2024).  
> [arXiv:2407.12104](https://arxiv.org/abs/2407.12104)

The code is written in **Mathematica** and is designed to be fully reproducible.

---

We consider two scenarios. The first one is a dark sector composed only by a singlet unstable dark matter candidate $\phi$ in cBE.wl. The second case is a dark sector composed of:

- A **dark matter particle** $S$ with self number changing reactions of the form $3\leftrightarrow 2$,
- A **mediator** $\phi$,

found in the cBE with mediator.wl file. The dark sector is populated via the freeze-in mechanism.

The evolution of the system is governed by **coupled Boltzmann equations** for the comoving number densities
\[
Y_S = \frac{n_S}{s}, \qquad Y_\phi = \frac{n_\phi}{s},
\]
where \( s \) is the SM entropy density and \( x = m_\phi / T \) is used as the time variable.

The code includes:

- **Freeze-in production** of \( S \) and/or \( \phi \) from the SM bath
- **Cannibalization processes** in the dark sector (e.g. \( \phi \phi \phi \to \phi \phi \))
- **Mediator decay** \( \phi \to \text{SM} \) and associated entropy injection
- **Hidden sector temperature evolution** \( T'(x) \)
- Proper treatment of **Hubble expansion** and changing relativistic degrees of freedom

---

## Features

### Coupled Boltzmann equations

We solve a system of the form

```math
\begin{align}
\frac{Y'_S}{Y_S} &= \frac{1}{x\,\tilde H}
\left(
    \langle C_{h\to \phi SS^*} \rangle
  + \langle C_{h\to SS^*} \rangle
  + \langle C_{\phi\phi\leftrightarrow SS^*} \rangle
  + \langle C_{3\leftrightarrow 2} \rangle
\right),
\\[6pt]
-\frac{x_S'}{x_S} &= 
\frac{1}{x\,\tilde H}
\left(
    \langle C_{h\to \phi SS^*} \rangle_2
  + \langle C_{h\to SS^*} \rangle_2
  + \langle C_{\phi S\leftrightarrow \phi S} \rangle_2
  + \langle C_{3\leftrightarrow 2} \rangle_2
\right)
- \frac{Y'_S}{Y_S}
+ \frac{H}{x\,\tilde H} \frac{\langle p^4/E^3 \rangle}{3 T_S}
+ \frac{2 s'}{3 s},
\\[12pt]
\frac{Y'_\phi}{Y_\phi} &= \frac{1}{x\,\tilde H}
\left(
    \langle C_{h\to \phi SS^*} \rangle
  + \langle C_{\text{SM SM}\to \text{SM}\,\phi} \rangle
  + \langle C_{\phi\phi\leftrightarrow SS^*} \rangle
\right),
\\[6pt]
-\frac{x_\phi'}{x_\phi} &= \frac{1}{x\,\tilde H}
\left(
    \langle C_{h\to \phi SS^*} \rangle_2
  + \langle C_{\text{SM SM}\to \text{SM}\,\phi} \rangle_2
  + \langle C_{\phi S\leftrightarrow \phi S} \rangle_2
\right)
- \frac{Y'_\phi}{Y_\phi}
+ \frac{H}{x\,\tilde H} \frac{\langle p^4/E^3 \rangle}{3 T_\phi}
+ \frac{2 s'}{3 s}.
\end{align}


where the collision terms \( \mathcal{C}_S \) and \( \mathcal{C}_\phi \) include:

- Freeze-in production from the SM
- \(2 \leftrightarrow 2\) and \(3 \leftrightarrow 2\) cannibal processes in the dark sector
- Decays and inverse decays of \( \phi \)

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
