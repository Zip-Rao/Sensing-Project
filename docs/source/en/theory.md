# Theory

This page summarises the physics behind the three v1 pipelines: from the
transmon device model, through phase accumulation in flux sensing, to the
**short Ramsey-type pulse** and its **response kernel** — the tool that links a
measurement back to the flux waveform being sensed. It is a self-contained
overview keeping only the minimum needed to understand the v1 pipelines; the
full multi-order Volterra expansion and kernel-extraction methods are in the
project report.

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

## Short Ramsey-type pulses

The Ramsey phase formula above assumes the free-evolution time $\tau$ is much
longer than the pulses, so the phase accumulated during the pulses is
negligible. Resolving a fast-varying flux requires a short sequence. A class of
**short Ramsey-type pulses**, motivated by the quantum speed limit, places two
phase-orthogonal control pulses back to back with no separate free-evolution
window (Herb and Degen 2024). This brings the time resolution of a dynamic
signal down to the pulse-duration scale; Herb et al. (2025) reached 1.1 ns on
an NV-center platform. The transmon suits this regime because its microwave
control is programmable at nanosecond timescales.

In this regime the phase during the pulses cannot be neglected, and both the
detuning $\delta\omega(t)$ and the Rabi envelope $\Omega(t)$ vary in time. The
two phase-orthogonal pulses need not be $\pi/2$ rotations: the more general
sequence uses $R_y(\alpha)$ and $R_{\pm x}(\alpha)$. At fixed maximum Rabi rate,
reducing $\alpha$ shortens the control interval and narrows the response kernel,
but also weakens the population response. The angle therefore trades temporal
resolution against sensitivity and acquisition cost. The fixed-axis Bloch
rotation picture no longer applies, so a more general response
description is required.

## The response kernel

Split the rotating-frame Hamiltonian into a **known control term** and an
**unknown perturbation term**:

$$H(t) = H_0(t) + V(t), \qquad H_0(t)=\tfrac{\Omega(t)}{2}\bigl(\cos\phi_d\,\sigma_x+\sin\phi_d\,\sigma_y\bigr), \quad V(t)=\tfrac{\delta\omega(t)}{2}\,\sigma_z,$$

where $H_0$ is the known control sequence and $V(t)$ is driven by the unknown
detuning $\delta\omega(t)$. In the interaction picture defined by $H_0$, the
final measurement expectation $\langle M\rangle_T$ is a functional of
$\delta\omega(t)$. Expanded to first order for a small signal, the change of the
measurement relative to baseline is a **convolution**:

$$\Delta\langle M\rangle \;\approx\; \int_0^T k_1(t)\,\delta\omega(t)\,dt, \qquad k_1(t) = i\,\langle 0|\,[\,W(t),\,M_I(T)\,]\,|0\rangle,$$

with $W(t)=\tfrac12 U_0^\dagger(t)\,\sigma_z\,U_0(t)$ and $U_0$ the propagator of
the control Hamiltonian. The **first-order response kernel** $k_1(t)$ depends
only on the initial state, the control sequence $H_0$, and the measurement
operator $M$. It is the protocol's time-weighting function for the detuning:
its magnitude and sign give the strength and direction with which the detuning
at each instant contributes to the final measurement.

In flux sensing, a small signal makes the detuning approximately linear in the
normalised flux perturbation, $\delta\omega(t)\approx\kappa\,\delta\Phi(t)/\Phi_0$
($\kappa=d\omega_T/d\Phi$), so the measurement and the flux waveform are again a
convolution $\Delta p(t)\approx (k*\Phi)(t)$. This linear relation is the basis
for the inverse problem in
{doc}`waveform reconstruction <examples/waveform_reconstruction>`.

At larger amplitudes the response exceeds the linear approximation, and the
expansion continues to higher-order Volterra kernels $k_n$. The platform
implements these: {py:class}`~sqc.reconstruction.KernelEstimator` estimates
kernels to arbitrary order, including the full off-diagonal
$k_n(t_1,\dots,t_n)$. Transient frequency calibration uses the third-order
kernel $\iiint k_3$ by default to correct the leading nonlinearity. The full
derivation of the Volterra expansion, and the frequency-kernel versus
flux-kernel distinction, are in the project report.

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
- K. Herb, C. L. Degen, *Quantum speed limit in quantum sensing*, Phys. Rev.
  Lett. **133**, 210802 (2024).
- K. Herb et al., *Quantum magnetometry of transient signals with a time
  resolution of 1.1 nanoseconds*, Nat. Commun. **16**, 822 (2025).
