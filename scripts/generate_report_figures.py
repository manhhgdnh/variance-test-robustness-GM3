from pathlib import Path
import csv
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2

ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results" / "tables"
FIGURES = ROOT / "results" / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)


def read_two_column(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    return np.array([float(r["n"]) for r in rows]), np.array([float(r["rejection_percent"]) for r in rows])


# H0 calibration
n, y = read_two_column(TABLES / "null_calibration_reported.csv")
fig, ax = plt.subplots(figsize=(8, 4.8))
ax.plot(n, y, marker="o", label="Reported empirical rejection")
ax.axhline(5, linestyle="--", label="Nominal level = 5%")
ax.set_xscale("log")
ax.set_xticks(n, [str(int(v)) for v in n])
ax.set_xlabel("Sample size n")
ax.set_ylabel("Rejection rate (%)")
ax.set_title("Calibration under H0")
ax.legend()
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(FIGURES / "null_calibration.png", dpi=180)
plt.close(fig)

# Power curves
with open(TABLES / "power_reported.csv", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
ns = np.array([50, 100, 200, 500, 1000], dtype=float)
fig, ax = plt.subplots(figsize=(8, 4.8))
for r in rows:
    vals = [float(r[str(int(v))]) for v in ns]
    ax.plot(ns, vals, marker="o", label=f"epsilon = {r['epsilon']}")
ax.set_xscale("log")
ax.set_xticks(ns, [str(int(v)) for v in ns])
ax.set_ylim(0, 105)
ax.set_xlabel("Sample size n")
ax.set_ylabel("Empirical power (%)")
ax.set_title("Power of the two-sided variance test")
ax.legend()
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(FIGURES / "power_curves.png", dpi=180)
plt.close(fig)

# Exponential robustness
n, y = read_two_column(TABLES / "exponential_reported.csv")
fig, ax = plt.subplots(figsize=(8, 4.8))
ax.plot(n, y, marker="o", label="Exponential data")
ax.axhline(5, linestyle="--", label="Nominal level = 5%")
ax.set_xscale("log")
ax.set_xticks(n, [str(int(v)) for v in n])
ax.set_ylim(0, 40)
ax.set_xlabel("Sample size n")
ax.set_ylabel("Rejection rate (%)")
ax.set_title("Loss of calibration under non-Gaussian data")
ax.legend()
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(FIGURES / "exponential_robustness.png", dpi=180)
plt.close(fig)

# AR(1) robustness
with open(TABLES / "ar1_reported.csv", newline="", encoding="utf-8") as f:
    rows = list(csv.DictReader(f))
fig, ax = plt.subplots(figsize=(8, 4.8))
for r in rows:
    vals = [float(r[str(int(v))]) for v in ns]
    ax.plot(ns, vals, marker="o", label=f"rho = {r['rho']}")
ax.axhline(5, linestyle="--", label="Nominal level = 5%")
ax.set_xscale("log")
ax.set_xticks(ns, [str(int(v)) for v in ns])
ax.set_ylim(0, 45)
ax.set_xlabel("Sample size n")
ax.set_ylabel("Rejection rate (%)")
ax.set_title("Sensitivity to dependence: AR(1)")
ax.legend()
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(FIGURES / "ar1_robustness.png", dpi=180)
plt.close(fig)

# Illustrative chi-square fit with a fixed reproducible seed.
rng = np.random.default_rng(2026)
n_graph = 1000
n_sim = 400
samples = rng.normal(0, 1, size=(n_sim, n_graph))
s2 = samples.var(axis=1, ddof=1)
t = (n_graph - 1) * s2
x = np.linspace(t.min() * 0.98, t.max() * 1.02, 500)
fig, ax = plt.subplots(figsize=(8, 4.8))
ax.hist(t, bins=30, density=True, alpha=0.6, label="Monte Carlo histogram")
ax.plot(x, chi2.pdf(x, df=n_graph - 1), linewidth=2, label=r"$\chi^2_{999}$ density")
ax.set_xlabel(r"$T = (n-1)S^2/\sigma_0^2$")
ax.set_ylabel("Density")
ax.set_title("Distribution of the test statistic under H0")
ax.legend()
fig.tight_layout()
fig.savefig(FIGURES / "chi_square_fit.png", dpi=180)
plt.close(fig)

print(f"Figures written to {FIGURES}")
