"""Unit tests for P5 multi-line transfer matrices."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.control.waveform import Waveform
from sqc.hardware.transfer_matrix import TransferMatrix


pytestmark = pytest.mark.unit


def test_from_dc_matrix_apply_matches_matrix_multiplication():
    transfer = TransferMatrix.from_dc_matrix(
        dc_matrix=np.array([[1.0, 0.05], [0.03, 1.0]]),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )
    t = np.linspace(0.0, 100.0, 1000)
    v_a = Waveform(t_list=t, samples=np.where(t > 10.0, 1.0, 0.0))
    v_b = Waveform(t_list=t, samples=np.zeros_like(t))

    fluxes = transfer.apply({"QA": v_a, "QB": v_b})

    assert fluxes["QA"].samples[-1] == pytest.approx(1.0, abs=1e-12)
    assert fluxes["QB"].samples[-1] == pytest.approx(0.03, abs=1e-12)


def test_transfer_matrix_accessors_split_diagonal_and_crosstalk():
    transfer = TransferMatrix.from_dc_matrix(
        dc_matrix=np.array([[1.0, 0.05], [0.03, 1.0]]),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )

    assert transfer.sources == ["QA", "QB"]
    assert transfer.targets == ["QA", "QB"]
    assert set(transfer.diagonal()) == {"QA", "QB"}
    assert set(transfer.off_diagonal()) == {("QA", "QB"), ("QB", "QA")}
    assert np.allclose(transfer.H_ji("QB", "QA"), 0.03)


def test_transfer_matrix_frequency_domain_delay_matches_fft_convolution():
    n = 64
    dt = 0.5
    t = np.arange(n) * dt
    omega = 2.0 * np.pi * np.fft.fftfreq(n, d=dt)
    transfer = TransferMatrix(
        elements={("QB", "QA"): np.exp(-1j * omega * dt)},
        frequency_axis=omega,
    )
    impulse = np.zeros(n)
    impulse[0] = 1.0

    fluxes = transfer.apply({"QA": Waveform(t_list=t, samples=impulse)})

    expected = np.zeros(n)
    expected[1] = 1.0
    assert np.allclose(fluxes["QB"].samples, expected, atol=1e-12)


def test_transfer_matrix_rejects_mismatched_source_grids():
    transfer = TransferMatrix.from_dc_matrix(
        dc_matrix=np.eye(2),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )
    t = np.linspace(0.0, 10.0, 101)
    shifted = t + 0.01

    with pytest.raises(ValueError):
        transfer.apply(
            {
                "QA": Waveform(t_list=t, samples=np.ones_like(t)),
                "QB": Waveform(t_list=shifted, samples=np.ones_like(shifted)),
            }
        )
