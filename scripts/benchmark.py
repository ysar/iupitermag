import time

import numpy as np
from matplotlib import pyplot as plt

import iupitermag as im

plt.style.use("dark_background")


def bench(label: str, fn, n: int, repeat: int = 5, warmup: int = 1) -> float:
    # warmup
    for _ in range(warmup):
        fn()

    times = []
    for _ in range(repeat):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)

    med = float(np.median(times))
    return med


def bench_for_num(N):
    rng = np.random.default_rng(0)
    positions = np.empty((N, 3), dtype=np.float64)
    positions[:, 0] = rng.uniform(5.0, 35.0, size=N)  # r
    positions[:, 1] = rng.uniform(0.0, np.pi, size=N)  # theta
    positions[:, 2] = rng.uniform(0.0, 2.0 * np.pi, size=N)  # phi

    field = im.InternalField("JRM33")

    def case_calc_field_loop() -> np.ndarray:
        out = np.empty((N, 3), dtype=np.float64)
        for i in range(N):
            r, theta, phi = positions[i]
            out[i, :] = field.calc_field(float(r), float(theta), float(phi))
        return out

    def case_map_calc_field() -> np.ndarray:
        return np.asarray(field.map_calc_field(positions))

    def case_parmap_calc_field() -> np.ndarray:
        return np.asarray(field.parmap_calc_field(positions))

    print("Benchmark: InternalField('JRM33')")
    print(f"N={N:,}")
    print()

    med1 = bench("1) calc_field (loop)", case_calc_field_loop, N)
    med2 = bench("2) map_calc_field", case_map_calc_field, N)
    med3 = bench("3) parmap_calc_field", case_parmap_calc_field, N)
    return (med1, med2, med3)


if __name__ == "__main__":
    len_num = 1000
    num_points = np.geomspace(1, 100_00, len_num, dtype=int)

    medians1 = np.empty(len_num)
    medians2 = np.empty(len_num)
    medians3 = np.empty(len_num)

    for i, n in enumerate(num_points):
        medians1[i], medians2[i], medians3[i] = bench_for_num(n)

    fig, ax = plt.subplots(1, 1, figsize=(6, 3), dpi=200)

    ax.loglog(num_points, medians1, lw=1.0, label="Python loop")
    ax.loglog(num_points, medians2, lw=1.0, label="map_calc_field")
    ax.loglog(num_points, medians3, lw=1.0, label="parmap_calc_field")
    ax.legend(loc="upper left", frameon=False)
    ax.set_xlabel("num points")
    ax.set_ylabel("time [s]")
    plt.show()
    fig.savefig("images/benchmark.png", facecolor="k", bbox_inches="tight")
    plt.show()
    plt.close(fig)
