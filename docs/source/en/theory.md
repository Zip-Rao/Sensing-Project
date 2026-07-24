# Theory

This page summarises the physics behind the three v1 pipelines. It is a
self-contained overview; for derivations see Gao, Rol, Touzard & Wang (2021,
*PRX Quantum* 2, 040202) and Koch et al. (2007).

## The transmon qubit

A transmon is a Josephson-junction superconducting qubit with Hamiltonian

$$H = 4 E_C\, n^2 - E_J(\Phi)\, \cos\varphi,$$

where $E_C = e^2/2C_\Sigma$ is the charging energy, $n$ and $\varphi$ are the
Cooper-pair-number and phase operators, and the Josephson energy is tuned by an
external flux through a SQUID loop:

$$E_J(\Phi) = E_{J0}\,\lvert\cos(\pi \Phi/\Phi_0)\rvert.$$

In the transmon regime $E_J/E_C \gg 1$, the system is a weakly anharmonic
oscillator:

$$\omega_T(\Phi) = \sqrt{8 E_J(\Phi) E_C} - E_C, \qquad \alpha = -E_C,$$

where $\omega_T$ is the $|0\rangle\!\to\!|1\rangle$ transition frequency and
$\alpha$ the anharmonicity that makes the two-level subspace addressable.

## Flux sensing

An external flux $\Phi(t)$ shifts $E_J$ and hence the qubit frequency
$\omega_T(\Phi)$. A Ramsey sequence $(\pi/2)\!-\!\tau\!-\!(\pi/2)$ accumulates a
coherent phase

$$\phi(\tau) = \int_0^\tau \Delta\omega(t')\,dt' = \int_0^\tau \kappa\,\Phi(t')\,dt',$$

where $\kappa = d\omega_T/d\Phi$ is the flux sensitivity. The excited-state
population then reads

$$p_e(\tau) = \tfrac{1}{2}\bigl(1 - \cos[\phi(\tau)]\bigr),$$

from which $\phi(\tau)$ — and thus $\Phi(t)$ — is reconstructed.

### The sweet spot

Because $E_J(\Phi) \propto \lvert\cos(\pi\Phi/\Phi_0)\rvert$, the frequency is
flat ($d\omega_T/d\Phi = 0$) at integer flux quanta. This *sweet spot* gives
first-order insensitivity to flux noise but zero sensing sensitivity; sensing
protocols bias away from it to a finite $\kappa$.

## Control pulses

Single-qubit gates use a resonant microwave drive. In the rotating frame with
the rotating-wave approximation,

$$H_d \approx \tfrac{\Omega(t)}{2}\bigl(a\,e^{i\phi} + a^\dagger e^{-i\phi}\bigr),$$

with rotation angle $\theta = \int \Omega(t)\,dt$. A $\pi$ pulse flips the
population; a $\pi/2$ pulse prepares a superposition. DRAG shaping (Motzoi 2009)
suppresses leakage into the third level.

## Recommended parameter ranges

| Parameter | Range | Default |
|---|---|---|
| $E_J/h$ | 10–25 GHz | 15 GHz |
| $E_C/h$ | 160–400 MHz | 200 MHz |
| $E_J/E_C$ | ~50 | ~75 |
| $f_{01}$ | 4–8 GHz | ~4.7 GHz |
| $\alpha/h$ | 200–300 MHz | 200 MHz |

These keep $f_{01}$ within commercial AWG/HEMT bandwidth and the anharmonicity
large enough for ~10–20 ns gates without leakage.

## References

- Z. Gao, M. Rol, S. Touzard, C. Wang, *Practical Guide for Building
  Superconducting Quantum Devices*, PRX Quantum **2**, 040202 (2021).
- J. Koch et al., *Charge-insensitive qubit design derived from the Cooper pair
  box*, Phys. Rev. A **76**, 042319 (2007).
- F. Motzoi et al., *Simple Pulses for Elimination of Leakage in Weakly
  Nonlinear Qubits*, Phys. Rev. Lett. **103**, 110501 (2009).
