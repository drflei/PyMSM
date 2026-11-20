"""Calculate and plot vertical cutoff rigidity against Kp."""

import argparse

import matplotlib.pyplot as plt
import numpy as np

import pymsm


def calculate_cutoffs():
    kps = np.arange(10)
    times = ["2020-01-01T00:00:00"] * len(kps)
    positions = [[400.0, 50.0, 0.0]] * len(kps)
    model = pymsm.PyMSM(times=times, positions=positions, kps=kps)
    cutoffs = model.getRc()
    if not np.all(np.isfinite(model.lm) & (model.lm > 0) & (model.lm < 100)):
        raise ValueError("IRBEM returned invalid L values for the Kp sweep")
    if not np.all(np.isfinite(cutoffs) & (cutoffs >= 0)):
        raise ValueError("PyMSM returned invalid cutoff rigidities")
    return kps, cutoffs


def plot_cutoffs(kps, cutoffs, output):
    fig, ax = plt.subplots()
    ax.plot(kps, cutoffs, marker="o")
    ax.set(
        xlabel="Kp",
        ylabel="Vertical cutoff rigidity (GV)",
        title="50° N, 0° E, 400 km; 2020-01-01 00:00 UT",
        xticks=kps,
    )
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(output, dpi=150)
    plt.close(fig)


def test_kp_cutoff_sweep(tmp_path):
    kps, cutoffs = calculate_cutoffs()
    np.testing.assert_array_equal(kps, np.arange(10))
    assert cutoffs.shape == (10,)
    assert np.all(cutoffs > 0)

    for i in (0, 9):
        individual = pymsm.PyMSM(
            ["2020-01-01T00:00:00"], [[400.0, 50.0, 0.0]], kps[i:i + 1]
        ).getRc()
        np.testing.assert_allclose(cutoffs[i], individual[0])

    output = tmp_path / "kp_cutoff_2020.png"
    plot_cutoffs(kps, cutoffs, output)
    assert output.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", default="kp_cutoff_2020.png")
    args = parser.parse_args()
    kps, cutoffs = calculate_cutoffs()
    for kp, cutoff in zip(kps, cutoffs):
        print(f"Kp {kp}: {cutoff:.4f} GV")
    plot_cutoffs(kps, cutoffs, args.output)
    print(f"Plot saved to {args.output}")
