#!/usr/bin/env python3
"""Figure 1 of the paper: how contamination destroys the zero-false-alarm corner.

Left: a pair with negative law Unif(0,1) and positive law Unif(0,2). The largest upper set
with no negative mass is [1, inf); its recall is rho* = P_1([1, inf)) = 1/2.

Right: the contaminated pair of Proposition 1, with negative law (1-eps) P_0 + eps P_1 and
the same positive law. Its total-variation distance from the left pair is at most eps, but
the negative law now reaches the whole support of the positives, so the best
zero-false-alarm rule has recall 0. The exaggerated eps = 0.15 is only for visibility; the
statement is for every eps > 0.

The two recalls are computed here from the piecewise-constant densities, not typed in.
"""
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
EPS = 0.15
# piecewise-constant densities on the cells (0,1) and (1,2)
P1 = (0.5, 0.5)                                     # Unif(0,2)
P0 = (1.0, 0.0)                                     # Unif(0,1)
Q0 = tuple((1 - EPS) * a + EPS * b for a, b in zip(P0, P1))


def recall_of_best_zero_fa_rule(neg, pos):
    """Largest upper set with zero negative mass is the union of the trailing cells where the
    negative density is zero; its recall is the positive mass there (cells have width 1)."""
    mass, cells = 0.0, [1, 0]                       # walk from the right: cell (1,2), then (0,1)
    for idx in cells:
        if neg[idx] > 0:
            break
        mass += pos[idx]
    return mass


rho_p, rho_q = recall_of_best_zero_fa_rule(P0, P1), recall_of_best_zero_fa_rule(Q0, P1)
tv = 0.5 * sum(abs(a - b) for a, b in zip(P0, Q0))
assert abs(rho_p - 0.5) < 1e-12 and rho_q == 0.0 and tv <= EPS + 1e-12


def draw(ax, neg, pos, title, rho, shade):
    xs = [0, 1, 1, 2]
    ax.step([0, 1, 2], [neg[0], neg[1], neg[1]], where="post", color="black", linewidth=1.6,
            label=r"negatives $P_0$" if shade else r"negatives $Q_0^\varepsilon$")
    ax.step([0, 1, 2], [pos[0], pos[1], pos[1]], where="post", color="0.45", linewidth=1.6,
            linestyle="--", label=r"positives $P_1$")
    if shade:
        ax.axvspan(1, 2, color="0.88")
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1.7)
    ax.set_xlabel("score")
    ax.set_title(title)
    ax.text(1.5, 0.2, rf"$\rho^\star={rho:g}$", ha="center")
    return xs


def main():
    plt.rcParams.update({
        "font.size": 10, "axes.labelsize": 10, "axes.titlesize": 10, "legend.fontsize": 10,
        "xtick.labelsize": 10, "ytick.labelsize": 10,
        "font.family": "serif", "mathtext.fontset": "cm", "pdf.fonttype": 42,
    })
    fig, (a, b) = plt.subplots(1, 2, figsize=(6.5, 2.4), sharey=True)
    draw(a, P0, P1, "a pair that attains the corner", rho_p, True)
    draw(b, Q0, P1, rf"contaminated, $\varepsilon={EPS}$, $d\leq\varepsilon$", rho_q, False)
    a.set_ylabel("density")
    a.legend(frameon=False, loc="upper right")
    b.legend(frameon=False, loc="upper right")
    fig.tight_layout()
    fig.savefig(HERE / "contamination.pdf", metadata={"CreationDate": None, "ModDate": None})
    print(f"rho*(P) = {rho_p}, rho*(Q) = {rho_q}, TV between the negative laws = {tv:.3f}")


if __name__ == "__main__":
    main()
