import datetime

import numpy as np
import pytest

import pymsm


@pytest.mark.parametrize(
    "date",
    [
        "1954-12-31T12:00:00",
        "1955-01-01T00:00:00",
        "1957-07-01T12:00:00",
        "2020-01-01T00:00:00",
        "2022-07-01T12:00:00",
        "2024-12-31T23:59:59",
        "2025-01-01T00:00:00",
    ],
)
def test_cutoff_interpolates_between_map_epochs(monkeypatch, date):
    model = pymsm.PyMSM([date], [[400.0, 50.0, 0.0]], np.array([1]))
    model.lm[:] = 2.0
    model.radius[:] = 1.0
    monkeypatch.setattr(model.dbMgr, "getMap", lambda *args: None)
    monkeypatch.setattr(
        model, "getRC450km",
        lambda key, lat, lon, lm: float(int(key[:4]) - 1955 + 1),
    )

    def fixed_lstar(position, maginput):
        epoch = int(position["dateTime"][:4])
        model.model.make_lstar_output = {
            "Lm": np.array([2.0 + (epoch - 1955) / 100.0])
        }

    monkeypatch.setattr(model.model, "make_lstar", fixed_lstar)
    result = model.getRc()[0]

    time = datetime.datetime.fromisoformat(date)
    if time.year < 1955:
        interpolated_year = 1955.0
    elif time.year >= 2025:
        interpolated_year = 2025.0
    else:
        start_year = 1955 + 5 * ((time.year - 1955) // 5)
        start = datetime.datetime(start_year, 1, 1)
        end = datetime.datetime(start_year + 5, 1, 1)
        fraction = (time - start) / (end - start)
        interpolated_year = start_year + 5 * fraction
    expected_rc = interpolated_year - 1955 + 1
    expected_lm = 2.0 + (interpolated_year - 1955) / 100.0
    expected = expected_rc * (expected_lm / 2.0) ** 2
    assert result == pytest.approx(expected)
