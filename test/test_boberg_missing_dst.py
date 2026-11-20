import numpy as np
import pytest

import pymsm


@pytest.mark.parametrize("year,kp", [(2022, 1), (2022, 7), (2025, 7)])
def test_missing_dst_uses_boberg_kp_estimate(monkeypatch, year, kp):
    kwargs = {
        "times": [f"{year}-01-01T12:00:00"],
        "positions": [[450.0, 45.0, 0.0]],
        "kps": np.array([kp]),
    }

    actual = pymsm.PyMSM(**kwargs).getTransmissionFunctions()

    original = pymsm.MagFields.make_lstar

    def with_missing_dst(model, positions, maginput):
        kp_values = maginput["Kp"]
        dst = (
            np.full_like(kp_values, -1e31, dtype=float)
            if isinstance(kp_values, np.ndarray)
            else -1e31
        )
        return original(model, positions, {**maginput, "Dst": dst})

    with monkeypatch.context() as patch:
        patch.setattr(pymsm.MagFields, "make_lstar", with_missing_dst)
        expected = pymsm.PyMSM(**kwargs).getTransmissionFunctions()

    assert 0 < actual[0][0] < 100
    assert actual[3][0] > 0
    for actual_value, expected_value in zip(actual, expected):
        np.testing.assert_allclose(actual_value, expected_value, rtol=1e-10)


def test_mixed_kp_inputs_use_one_missing_dst_per_position():
    times = ["2025-01-01T12:00:00"] * 2
    positions = [[450.0, 45.0, 0.0]] * 2
    kps = np.array([1, 7])
    combined = pymsm.PyMSM(times, positions, kps).getTransmissionFunctions()

    for i in range(len(kps)):
        individual = pymsm.PyMSM(
            times[i:i + 1], positions[i:i + 1], kps[i:i + 1]
        ).getTransmissionFunctions()
        for combined_value, individual_value in zip(combined, individual):
            np.testing.assert_allclose(combined_value[i], individual_value[0])
