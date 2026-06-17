import os

import corner
import numpy as np


N_PARM = 8
SKIP_ROWS = 15000000


def read_input_value(name, path="input.ini"):
    prefix = f"{name}:"
    with open(path, "r", encoding="utf-8") as fp:
        for line in fp:
            stripped = line.strip()
            if stripped.startswith(prefix):
                return stripped.split(":", 1)[1].strip().split()[0]
    raise KeyError(f"{name} not found in {path}")


def load_complete_chain(path):
    rows = []
    with open(path, "r", encoding="utf-8") as fp:
        for i, line in enumerate(fp):
            if i < SKIP_ROWS:
                continue
            parts = line.split()
            if len(parts) < N_PARM:
                continue
            rows.append([float(x) for x in parts[:N_PARM]])
    if not rows:
        raise ValueError(f"No complete samples found in {path}")
    return np.asarray(rows, dtype=float)


results_dir = read_input_value("results_dir")
n_beta = int(read_input_value("N_beta"))

# The beta=1 cold chain is the last rank when Beta_Values ends with 1.0.
cold_rank = n_beta - 1
chain_path = os.path.join(results_dir, f"chain2.dat")

chain = load_complete_chain(chain_path)
chain = chain[~np.isnan(chain).any(axis=1)]

plot_chain = chain[:, :N_PARM]

labels = [
    r"$P$ [d]",
    r"$K$ [m s$^{-1}$]",
    r"$\lambda$ [rad]",
    r"$h=\sqrt{e}\sin\omega$",
    r"$k=\sqrt{e}\cos\omega$",
    r"$\dot{\gamma}$ [m s$^{-1}$ d$^{-1}$]",
    r"$\gamma$ [m s$^{-1}$]",
    r"$\sigma_J$ [m s$^{-1}$]",
]
print(f"Reading {chain_path}")
print(f"samples={plot_chain.shape[0]} columns={plot_chain.shape[1]}")

figure = corner.corner(
    plot_chain,
    labels=labels,
    label_kwargs={"fontsize": 13},
    quantiles=[0.16, 0.5, 0.84],
    show_titles=True,
    title_kwargs={"fontsize": 11},
    title_fmt=".5g",
)

out_name = f"pic_corner_{os.path.basename(results_dir)}_chain{cold_rank}.png"
figure.savefig(out_name, dpi=300)
print(f"Saved {out_name}")
