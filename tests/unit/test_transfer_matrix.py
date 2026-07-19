"""tests.unit.test_transfer_matrix — Unit tests for TransferMatrix."""
from __future__ import annotations

import numpy as np
import pytest

from sqc.hardware.transfer_matrix import TransferMatrix
from sqc.control.waveform import Waveform


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def dc_tm_2x2() -> TransferMatrix:
    """2×2 DC transfer matrix with 5% crosstalk."""
    return TransferMatrix.from_dc_matrix(
        dc_matrix=np.array([[1.0, 0.05], [0.03, 1.0]]),
        source_names=["QA", "QB"],
        target_names=["QA", "QB"],
    )


@pytest.fixture
def simple_waveforms() -> dict[str, Waveform]:
    """Two simple waveforms."""
    t = np.linspace(0, 100, 500)
    V_A = Waveform(t_list=t, samples=np.where((t > 20) & (t < 80), 1.0, 0.0))
    V_B = Waveform(t_list=t, samples=np.zeros_like(t))
    return {"QA": V_A, "QB": V_B}


# ---------------------------------------------------------------------------
# from_dc_matrix factory
# ---------------------------------------------------------------------------

class TestFromDcMatrix:
    """Tests for TransferMatrix.from_dc_matrix()."""

    def test_creates_correct_shape(self, dc_tm_2x2):
        assert dc_tm_2x2.sources == ["QA", "QB"]
        assert dc_tm_2x2.targets == ["QA", "QB"]
        assert len(dc_tm_2x2.elements) == 4

    def test_dc_flat_response(self, dc_tm_2x2):
        H = dc_tm_2x2.H_ji("QA", "QA")
        assert np.allclose(H, 1.0)
        H = dc_tm_2x2.H_ji("QB", "QA")
        assert np.allclose(H, 0.03)
        H = dc_tm_2x2.H_ji("QB", "QB")
        assert np.allclose(H, 1.0)
        H = dc_tm_2x2.H_ji("QA", "QB")
        assert np.allclose(H, 0.05)

    def test_shape_mismatch_raises(self):
        with pytest.raises(ValueError, match="does not match"):
            TransferMatrix.from_dc_matrix(
                dc_matrix=np.array([[1.0]]),
                source_names=["QA", "QB"],
                target_names=["QA"],
            )

    def test_n_freq_parameter(self):
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0]]),
            source_names=["S"],
            target_names=["T"],
            n_freq=128,
        )
        assert len(tm.frequency_axis) == 128
        assert len(tm.elements[("T", "S")]) == 128


# ---------------------------------------------------------------------------
# apply consistency
# ---------------------------------------------------------------------------

class TestApply:
    """Tests for TransferMatrix.apply()."""

    def test_apply_zero_crosstalk_returns_self(self, dc_tm_2x2, simple_waveforms):
        """With zero off-diagonal, each target gets only its own source."""
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.0, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )
        fluxes = tm.apply(simple_waveforms)
        assert "QA" in fluxes
        assert "QB" in fluxes
        # QA flux should equal QA voltage
        np.testing.assert_array_almost_equal(
            fluxes["QA"].samples, simple_waveforms["QA"].samples,
        )
        # QB flux should equal QB voltage (all zeros)
        np.testing.assert_array_almost_equal(
            fluxes["QB"].samples, simple_waveforms["QB"].samples,
        )

    def test_apply_with_crosstalk(self, dc_tm_2x2, simple_waveforms):
        """Crosstalk introduces scaled cross-terms."""
        fluxes = dc_tm_2x2.apply(simple_waveforms)

        # QA flux = 1.0 * V_A + 0.05 * V_B = V_A (since V_B=0)
        expected_QA = simple_waveforms["QA"].samples
        np.testing.assert_array_almost_equal(fluxes["QA"].samples, expected_QA)

        # QB flux = 0.03 * V_A + 1.0 * V_B = 0.03 * V_A
        expected_QB = 0.03 * simple_waveforms["QA"].samples
        np.testing.assert_array_almost_equal(fluxes["QB"].samples, expected_QB)

    def test_apply_empty_input(self):
        """Empty source voltage dict returns empty result."""
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0]]),
            source_names=["S"],
            target_names=["T"],
        )
        fluxes = tm.apply({})
        assert fluxes == {}

    def test_apply_missing_element(self, dc_tm_2x2, simple_waveforms):
        """A missing (target, source) pair is treated as zero transfer."""
        # Remove one element
        del dc_tm_2x2.elements[("QB", "QA")]
        fluxes = dc_tm_2x2.apply(simple_waveforms)
        # QB should get zero flux (no transfer from QA, QB is zero)
        assert abs(fluxes["QB"].samples).max() < 1e-12

    def test_apply_preserves_time_axis(self, dc_tm_2x2, simple_waveforms):
        """Output FluxSignal has same t_list as input."""
        fluxes = dc_tm_2x2.apply(simple_waveforms)
        np.testing.assert_array_almost_equal(
            fluxes["QA"].t_list, simple_waveforms["QA"].t_list,
        )

    def test_apply_different_lengths(self):
        """Waveforms of different lengths are zero-padded."""
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[1.0, 0.0], [0.0, 1.0]]),
            source_names=["QA", "QB"],
            target_names=["QA", "QB"],
        )
        t_short = np.linspace(0, 10, 50)
        t_long = np.linspace(0, 20, 200)
        V_A = Waveform(t_list=t_long, samples=np.ones_like(t_long))
        V_B = Waveform(t_list=t_short, samples=np.zeros_like(t_short))
        fluxes = tm.apply({"QA": V_A, "QB": V_B})
        # Result uses the longer t_list
        assert len(fluxes["QA"].t_list) == 200

    def test_apply_single_source(self):
        """Single source, no other lines."""
        tm = TransferMatrix.from_dc_matrix(
            dc_matrix=np.array([[2.5]]),
            source_names=["S0"],
            target_names=["T0"],
        )
        t = np.linspace(0, 50, 100)
        V = Waveform(t_list=t, samples=np.sin(2 * np.pi * 0.1 * t))
        fluxes = tm.apply({"S0": V})
        np.testing.assert_array_almost_equal(fluxes["T0"].samples, 2.5 * V.samples)


# ---------------------------------------------------------------------------
# H_ji / diagonal / off_diagonal access
# ---------------------------------------------------------------------------

class TestAccess:
    """Tests for element access properties."""

    def test_H_ji_access(self, dc_tm_2x2):
        H = dc_tm_2x2.H_ji("QA", "QA")
        assert H.shape == dc_tm_2x2.frequency_axis.shape

    def test_H_ji_missing_raises(self, dc_tm_2x2):
        with pytest.raises(KeyError):
            dc_tm_2x2.H_ji("ZZ", "QA")

    def test_diagonal(self, dc_tm_2x2):
        diag = dc_tm_2x2.diagonal()
        assert set(diag.keys()) == {"QA", "QB"}
        assert np.allclose(diag["QA"], 1.0)
        assert np.allclose(diag["QB"], 1.0)

    def test_off_diagonal(self, dc_tm_2x2):
        off = dc_tm_2x2.off_diagonal()
        assert ("QB", "QA") in off
        assert ("QA", "QB") in off
        assert np.allclose(off[("QB", "QA")], 0.03)
        assert np.allclose(off[("QA", "QB")], 0.05)

    def test_sources_targets_ordering(self, dc_tm_2x2):
        """Sources and targets are returned sorted."""
        assert dc_tm_2x2.sources == ["QA", "QB"]
        assert dc_tm_2x2.targets == ["QA", "QB"]
