"""sqc.workflows.z_crosstalk — ZCrosstalkWorkflow: dual-qubit crosstalk demo.

End-to-end Z-crosstalk demonstration:
1. Apply flux pulse to qubit A; crosstalk appears at qubit B.
2. Transient-sense on qubit B to reconstruct the parasitic flux.
3. Extract H_BA(omega) by frequency-domain deconvolution.
4. Design compensation and verify crosstalk suppression.

Per phase_5_handbook.md §3.3.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from sqc.control.waveform import Waveform
from sqc.control.flux_signal import FluxSignal
from sqc.workflows.base import Workflow


@dataclass
class ZCrosstalkWorkflow(Workflow):
    """End-to-end Z-crosstalk extraction and compensation demo.

    Setup: 2-qubit chip. Qubit A has its own Z line. Apply a flux pulse
    to A. Qubit B (parked at max-sensitivity point) feels parasitic flux
    via the crosstalk transfer matrix.

    Steps
    -----
    1. Apply transfer matrix: compute phi_B_true from source voltages.
    2. Transient-sense phi_B on qubit B.
    3. Reconstruct phi_B via Wiener deconvolution.
    4. Extract H_BA(omega) = FFT(phi_B_rec) / FFT(V_A).
    5. Design compensation and verify crosstalk suppression.

    Parameters
    ----------
    chip : ChipTopology
        At least 2 qubits (QA, QB) with known transfer_matrix.
    flux_pulse_on_A : Waveform
        Voltage waveform applied to qubit A's Z line.
    true_transfer_matrix : TransferMatrix
        Ground-truth transfer matrix (known in simulation).
    qubit_A_name : str
        Name of the driven qubit. Default "QA".
    qubit_B_name : str
        Name of the sensing qubit. Default "QB".
    t_rabi : np.ndarray or None
        Rabi pulse time axis for transient sensing. If None, uses
        np.linspace(0, 10, 10).
    scan_list : np.ndarray or None
        Pre-computed scan list. If None, defaults to aligning with
        the flux pulse time axis.
    wiener_lambda_reg : float
        Regularization for Wiener deconvolution. Default 1e-3.
    deconv_lambda_reg : float
        Regularization for H_BA extraction. Default 1e-3.
    """

    chip: object  # ChipTopology
    flux_pulse_on_A: Waveform
    true_transfer_matrix: object  # TransferMatrix

    qubit_A_name: str = "QA"
    qubit_B_name: str = "QB"
    t_rabi: Optional[np.ndarray] = None
    scan_list: Optional[np.ndarray] = None
    wiener_lambda_reg: float = 1e-3
    deconv_lambda_reg: float = 1e-3

    # Internal state (populated during run)
    _qubit_A: object = field(default=None, init=False)
    _qubit_B: object = field(default=None, init=False)
    _phi_A: Optional[FluxSignal] = field(default=None, init=False)
    _phi_B_true: Optional[FluxSignal] = field(default=None, init=False)

    def run(self) -> dict:
        """Execute the full Z-crosstalk workflow.

        Returns
        -------
        dict
            Keys:
            - phi_A : FluxSignal
            - phi_B_true : FluxSignal
            - phi_B_reconstructed : FluxSignal
            - H_BA_estimated : np.ndarray (complex)
            - H_BA_true : np.ndarray (complex)
            - fit_error_dB : float
            - parasitic_phase_uncompensated : float
            - parasitic_phase_compensated : float
            - compensation_factor : float
            - phi_B_after_compensation : np.ndarray
        """
        # ---- locate qubits ---------------------------------------------------
        self._qubit_A = self.chip.get_qubit(self.qubit_A_name)
        self._qubit_B = self.chip.get_qubit(self.qubit_B_name)

        # ---- step 1: apply transfer matrix to get on-chip fluxes ------------
        source_voltages = {self.qubit_A_name: self.flux_pulse_on_A}
        on_chip_fluxes = self.true_transfer_matrix.apply(source_voltages)

        self._phi_A = on_chip_fluxes.get(self.qubit_A_name)
        self._phi_B_true = on_chip_fluxes.get(self.qubit_B_name)

        if self._phi_B_true is None:
            raise ValueError(
                f"Transfer matrix has no element for ({self.qubit_B_name}, "
                f"{self.qubit_A_name}). Check true_transfer_matrix.elements."
            )

        # ---- step 2: reconstruct phi_B via transient sensing -----------------
        phi_B_reconstructed = self._reconstruct_phi_B()

        # ---- step 3: extract H_BA by frequency-domain deconvolution ---------
        H_BA_estimated, omega_grid = self._extract_H_BA(phi_B_reconstructed)

        # Get true H_BA for comparison
        try:
            H_BA_true_interp = self._interpolate_H_true(omega_grid)
        except KeyError:
            H_BA_true_interp = np.zeros_like(H_BA_estimated)

        fit_error_dB = self._compute_fit_error(H_BA_estimated, H_BA_true_interp)

        # ---- step 4: design compensation ------------------------------------
        comp_result = self._compensate(phi_B_reconstructed)

        return {
            "phi_A": self._phi_A,
            "phi_B_true": self._phi_B_true,
            "phi_B_reconstructed": phi_B_reconstructed,
            "H_BA_estimated": H_BA_estimated,
            "H_BA_true": H_BA_true_interp,
            "omega_grid": omega_grid,
            "fit_error_dB": fit_error_dB,
            **comp_result,
        }

    # ------------------------------------------------------------------
    # Step 2: reconstruct phi_B via transient sensing + Wiener
    # ------------------------------------------------------------------

    def _reconstruct_phi_B(self) -> FluxSignal:
        """Run transient sensing on qubit B and Wiener-reconstruct phi_B."""
        from sqc.experiments.transient import TransientSensingExperiment
        from sqc.reconstruction.wiener import WienerReconstruction

        # Park qubit B at its optimal work point
        qubit_B_for_exp = self._make_qubit_at_optimal(self._qubit_B)

        # Configure experiment
        t_rabi = self.t_rabi
        if t_rabi is None:
            t_rabi = np.linspace(0, 10, 10)

        # Use the crosstalk flux as the signal to sense
        exp = TransientSensingExperiment(
            qubit=qubit_B_for_exp,
            flux_signal=self._phi_B_true,
            t_rabi=t_rabi,
            scan_list=self.scan_list,
        )
        result = exp.run()

        # Wiener reconstruction
        wiener = WienerReconstruction(lambda_reg=self.wiener_lambda_reg)
        dt = float(self._phi_B_true.t_list[1] - self._phi_B_true.t_list[0])
        kernel = np.asarray(result.data["kernel"])

        phi_rec = wiener.reconstruct(
            measurement=result,
            kernel=kernel,
            dt=dt,
        )
        return phi_rec

    # ------------------------------------------------------------------
    # Step 3: H_BA extraction
    # ------------------------------------------------------------------

    def _extract_H_BA(
        self, phi_B_rec: FluxSignal,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Extract H_BA(omega) = Phi_B(omega) / V_A(omega).

        Uses regularized division to avoid blowup at zeros.
        Resamples phi_B_rec onto V_A's time grid before deconvolution.

        Parameters
        ----------
        phi_B_rec : FluxSignal
            Reconstructed flux at qubit B.

        Returns
        -------
        tuple[np.ndarray, np.ndarray]
            (H_BA_estimated, omega_grid)
        """
        V_A = self.flux_pulse_on_A

        # Resample phi_B_rec onto V_A's time grid
        t_phi = np.asarray(phi_B_rec.t_list)
        phi_arr = np.asarray(phi_B_rec.samples)
        t_V = np.asarray(V_A.t_list)

        phi_B_resampled = np.interp(t_V, t_phi, phi_arr, left=0.0, right=0.0)
        V_A_arr = np.asarray(V_A.samples)

        n_common = len(V_A_arr)
        dt = float(V_A.t_list[1] - V_A.t_list[0]) if len(V_A.t_list) > 1 else 1.0
        omega_grid = 2.0 * np.pi * np.fft.fftfreq(n_common, d=dt)

        Phi_B_omega = np.fft.fft(phi_B_resampled)
        V_A_omega = np.fft.fft(V_A_arr)

        # Regularized division: H = Phi * conj(V) / (|V|^2 + lambda^2)
        lam = self.deconv_lambda_reg
        H_BA = Phi_B_omega * np.conj(V_A_omega) / (np.abs(V_A_omega) ** 2 + lam**2)

        return H_BA, omega_grid

    def _interpolate_H_true(self, omega_grid: np.ndarray) -> np.ndarray:
        """Interpolate H_BA_true onto the same frequency grid as estimated."""
        H_true = self.true_transfer_matrix.H_ji(
            self.qubit_B_name, self.qubit_A_name,
        )
        src_omega = np.asarray(self.true_transfer_matrix.frequency_axis)

        if len(src_omega) <= 1:
            return H_true[0] * np.ones_like(omega_grid) if len(H_true) > 0 else np.zeros_like(omega_grid)

        sort_idx = np.argsort(src_omega)
        src_omega_s = src_omega[sort_idx]
        H_sorted = H_true[sort_idx]

        re = np.interp(omega_grid, src_omega_s, H_sorted.real)
        im = np.interp(omega_grid, src_omega_s, H_sorted.imag)
        return re + 1j * im

    # ------------------------------------------------------------------
    # Step 4: compensation
    # ------------------------------------------------------------------

    def _compensate(self, phi_B_reconstructed: FluxSignal) -> dict:
        """Design and verify crosstalk compensation.

        Computes:
        - parasitic_phase_uncompensated = sum(|phi_B_true|) * dt
        - phi_B_after_compensation = phi_B_true - phi_B_reconstructed
        - parasitic_phase_compensated = sum(|phi_B_after_comp|) * dt
        - compensation_factor = uncom / comp

        Parameters
        ----------
        phi_B_reconstructed : FluxSignal
            Wiener-reconstructed flux at qubit B.

        Returns
        -------
        dict
        """
        phi_true = np.asarray(self._phi_B_true.samples)
        phi_rec = np.asarray(phi_B_reconstructed.samples)

        # Match lengths
        n = min(len(phi_true), len(phi_rec))
        phi_true = phi_true[:n]
        phi_rec = phi_rec[:n]

        dt = float(self._phi_B_true.t_list[1] - self._phi_B_true.t_list[0])

        # Parasitic phase: integral of absolute flux
        phase_uncompensated = float(np.sum(np.abs(phi_true)) * dt)

        # After compensation: residual after subtracting reconstructed
        phi_after_comp = phi_true - phi_rec
        phase_compensated = float(np.sum(np.abs(phi_after_comp)) * dt)

        compensation_factor = phase_uncompensated / max(phase_compensated, 1e-30)

        return {
            "parasitic_phase_uncompensated": phase_uncompensated,
            "parasitic_phase_compensated": phase_compensated,
            "compensation_factor": compensation_factor,
            "phi_B_after_compensation": phi_after_comp,
        }

    # ------------------------------------------------------------------
    # Helper: park qubit at optimal sensitivity point
    # ------------------------------------------------------------------

    @staticmethod
    def _make_qubit_at_optimal(qubit):
        """Return a new qubit object at the optimal flux sensitivity point.

        Works with both src.qubit.TransmonQubit and sqc.devices.transmon.TransmonQubit.
        """
        import numpy as np

        # Determine optimal flux
        # QubitSpec.optimal_work_point() returns float in Phi_0 (e.g., ~0.1512)
        # Legacy TransmonQubit.optimal_work_point() returns arctan(sqrt(2)) in radians (~0.955)
        # We need Phi_0. For legacy, divide by pi to get Phi_0.
        try:
            # Try sqc.devices.transmon approach
            if hasattr(qubit, "spec"):
                spec = qubit.spec()
                optimal = spec.optimal_work_point()
            elif hasattr(qubit, "optimal_work_point"):
                legacy = qubit.optimal_work_point()
                # Legacy returns arctan(sqrt(2)) in rad. Convert to Phi_0.
                optimal = legacy / np.pi
            else:
                optimal = 0.3  # fallback
        except Exception:
            optimal = 0.3

        # Build a new qubit at optimal flux
        n_levels = getattr(qubit, "n_levels", 3)
        T1 = getattr(qubit, "T1", 10000.0)
        T2 = getattr(qubit, "T2", 8000.0)
        EC = getattr(qubit, "EC", 2 * np.pi * 0.2)
        EJ_0 = getattr(qubit, "EJ_0", 2 * np.pi * 15)
        name = getattr(qubit, "name", "Q")

        # Use the sqc TransmonQubit constructor for consistency
        from sqc.devices.transmon import TransmonQubit as SQTransmon

        return SQTransmon(
            EC=EC, EJ=EJ_0, T1=T1, T2=T2,
            flux=optimal, n_levels=n_levels, name=name,
        )

    # ------------------------------------------------------------------
    # Metric
    # ------------------------------------------------------------------

    @staticmethod
    def _compute_fit_error(H_est: np.ndarray, H_true: np.ndarray) -> float:
        """RMSE of estimated vs true transfer function, in dB.

        Parameters
        ----------
        H_est : np.ndarray
            Estimated H_BA(omega).
        H_true : np.ndarray
            True H_BA(omega).

        Returns
        -------
        float
            RMS error in dB.
        """
        err = np.abs(H_est - H_true)
        denom = np.abs(H_true)
        # Avoid division by zero: use mean(abs(H_true)) as floor
        avg_abs = np.mean(denom) if np.mean(denom) > 1e-12 else 1.0
        denom = np.where(denom < 1e-12, avg_abs, denom)
        rms = np.sqrt(np.mean((err / denom) ** 2))
        if rms < 1e-30:
            return -300.0
        return float(20.0 * np.log10(rms))
