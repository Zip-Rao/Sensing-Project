# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Quantum sensing simulation platform for superconducting Transmon qubits. The platform simulates time-dependent magnetic field sensing via qubit-frequency transduction, supporting multiple sensing protocols. All simulations use natural units (ħ=1) with Fock basis expansion via QuTiP.

## Workflow — before writing any code

1. **Restate your understanding.** What problem are you solving? What is the deliverable? Flag any assumptions you made. If you see a better technical approach, say so directly — the user will decide.
2. **Ask clarifying questions** (at most 3 at a time) until you are 100% confident about:
   - The real goal the user wants to achieve (not just what they said literally).
   - Unspoken constraints or preferences — tech stack, performance requirements, existing code that must remain compatible, parts of the system that must not be touched.
   - Your planned implementation — the core idea and why you chose it.
3. **Do not write code or modify any files** until the user explicitly signals approval (e.g. "可以开始", "go ahead", "proceed").

## Development commands

```bash
# Install dependencies (conda environment, per .vscode/settings.json)
pip install -r requirements.txt

# Core dependency check
python -c "import qutip, numpy, matplotlib; print('ok')"

# Verify project imports
python -c "from src.qubit import TransmonQubit; from src.signal import Signal; print('ok')"

# Run the Gradio web demo
python web_demo.py

# There is no test suite or linting configuration in this project.
```

## Architecture: data flow

```
Signal (signal.py)  →  qubit.qubit_in_mag(Signal)  →  Protocol.evolve(qubit)  →  Analysis methods
                         ↓                                    ↓                        ↓
                  Qubit frequency shifts            mesolve() with              Deconvolution /
                  under external flux               control pulses              numerical inversion
```

## Core modules (`src/`)

### [qubit.py](src/qubit.py) — Qubit physical models
- **`TransmonQubit`** — Central class. Holds EC, EJ, T1, T2, flux, state. Computes frequency and anharmonicity analytically. Key methods:
  - `qubit_in_mag(Phi_signal)` — accepts a `Signal` and pre-computes `self.freq_coeffs` and `self.H_list` (time-dependent Hamiltonian in QuTiP list format). This is the main interface for sensing simulations.
  - `qubit_under_mag(Phi_signal)` — older, slower per-timestep variant; creates a new qubit object at each time point. Prefer `qubit_in_mag()`.
  - `frequency_sensitivity(flux)` — numerical dω/dΦ via central difference.
  - Gate simulation: `simulate_gate()` with DRAG pulses, `ideal_gate()`.
- **`Cavity`** — Multi-mode cavity model with tensor-product operators.
- **`Coupled_System`** — Two qubits coupled via a tunable coupler to a multi-mode cavity. Used for two-qubit gate simulation (iSWAP, CZ).

### [signal.py](src/signal.py) — Magnetic flux signals
- **`Signal`** — Time-domain signal with 8 types (0=zero, 1=constant, 2=sinusoidal, 3=Gaussian, 4=asymmetric impulse, 5=double-peak, 6=basis-expanded, 7=complex wave-packet, 8=custom). Key property: `self.signal` (numpy array), `self.t_list`. Methods: `value_at(t)`, `truncate()`, `update_signal()`, `plot()`.
- **`CompositeSignal`** — Concatenates multiple Signals sequentially.

### [pulse.py](src/pulse.py) — Control pulses
- **`Pulse`** — Single control pulse. Accepts a `Signal` as the Rabi envelope (`Omega`). Builds Hamiltonian in either lab frame (`frame=0`) or rotating frame (`frame=1`), with or without RWA. Returns QuTiP list-format `[op, coeff_array]`.
- **`CompositePulse`** — Sequences multiple Pulses with proper time alignment.
- **Factory functions**: `create_pulse()`, `create_ramsey_pulse()`, `create_echo_pulse()`, `create_cpmg_pulse()`, `create_diff_echo_pulse()`, `create_cryoscope_pulse()`.
- **`get_kernel()`** on `Pulse`/`CompositePulse` computes the control kernel by perturbing with a narrow Gaussian stimulus at each time point.

### [protocal.py](src/protocal.py) — Sensing protocols
- **`Protocal`** (intentional spelling, used project-wide). `type` selects the protocol:
  - `0` Rabi, `1` Ramsey, `2` Differential echo, `3` CPMG (stub), `4` Transient field with sliding measurement, `5` Cryoscope.
- `evolve(qubit)` runs the protocol end-to-end: creates test signals, constructs Hamiltonians, calls `mesolve()`.
- `sliding_measrement()` — slides a control pulse across a signal, used by protocol type 4. Pre-computes qubit time series and Hamiltonians.
- `single_measurement()` — evaluates p_e at one delay time.
- **`Calibration`** — For frequency-flux and phase-flux calibration curves.
- **`IQ_readout()`** — I/Q demodulation by running two Ramsey sequences with π/2 phase offset.

### [analysis.py](src/analysis.py) — Data analysis and field reconstruction
- **`Analysis`** class with methods for extracting physical quantities from mesolve results:
  - `get_signal_from_ramsey_by_iq()` / `by_unwrap()` — Ramsey phase unwrapping and B-field derivation.
  - `get_signal_from_diff_echo()` — Direct formula B = -φ/(2k·κ·t_int).
  - `wiener_deconvolution()` — Linear Wiener deconvolution (FFT-based) of kernel from Δp measurements.
  - `hammerstein_wiener_deconvolution()` — Nonlinear block model: B → ω(B) → φ → Δp.
  - `numerical_inverse()` — Levenberg-Marquardt optimization using full density matrix simulation with basis function parameterization (B-spline, Fourier, Legendre). Includes adjoint-method Jacobian.
  - `get_signal_from_cryoscope()` — Phase-differential method with calibration-curve inversion.
- **Helper functions**: `generate_basis_functions()`, `basis_function_decomposition()`, `forward_simulation()`, `compute_jacobian()` (adjoint), `compute_jacobian_finite_difference()`, `levenberg_marquardt()`.

## Key conventions and gotchas

- **Class name**: `Protocal` is intentionally misspelled everywhere — do not "fix" it.
- **Units**: Frequencies/energies in GHz, time in ns, flux in Φ₀. ħ=1 throughout.
- **Hamiltonian format**: QuTiP list format `[[H0, coeffs], [H1, coeffs_array], ...]` is used pervasively for time-dependent Hamiltonians. `QobjEvo(H, tlist=t_list)` wraps these.
- **`qubit_in_mag()` must be called before simulation** whenever the signal changes — it pre-computes `freq_coeffs` and `H_list` used by `mesolve()`.
- **Kernel computation** uses a hardcoded stimulus amplitude (0.0215) and width (3 ns) in both `Pulse.get_kernel()` and `Analysis.get_kernel()`. These need manual adjustment when qubit parameters change.
- **LM inversion** (numerical_inverse) relies on `qubit_in_mag()` being called before each forward pass. The adjoint Jacobian (`compute_jacobian`) has known convergence issues vs. finite difference — see project TODO 0.3.

## Notebook and results

- [`Simulation.ipynb`](Simulation.ipynb) is the main experimental notebook for protocol testing and visualization.
- [`result/`](result/) contains output data and plotting scripts organized by experiment type (amplitude_scan, lambda_scan, convergence, basis_comparison, etc.).
- [`idea/`](idea/) contains design documents and implementation plans (markdown). The master TODO is [`idea/_TODO_master.md`](idea/_TODO_master.md).
- [`web_demo.py`](web_demo.py) is a Gradio web interface exposing protocols 0, 1, and 4.
