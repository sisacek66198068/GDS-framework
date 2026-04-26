import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

G = 1.0

# Reuse daughter material map by recomputing same dynamics
dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

M0 = 1.0
rho_c = 50.0
n = 3.0
A = 4.0
R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

R_A, V_A = 3.0, -0.15
R_B, V_B = 0.05, 0.0
M_A, M_B = M0, 0.0

rows = []

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

for ti in t:
    rhoA = rho_eff(M_A, R_A)
    gamA = Gamma(rhoA)

    grav_A = -G * M_A / (R_A**2 + R_core**2)
    sat = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)
    bounce = k_bounce * sat / (R_A + R_core)
    gate_R = np.exp(-(R_A / (3 * R_core))**2)

    acc_A = grav_A / gamA + bounce

    dM = transfer_strength * sat * gate_R * M_A * dt
    dM = min(dM, M_A)

    M_A -= dM
    M_B += dM

    rhoB = rho_eff(M_B + 1e-12, R_B)

    acc_B = 0.0
    if M_B > 0:
        acc_B += G * M_B / (R_B**2 + R_core**2) * (sat * gate_R)
        acc_B += 0.03 * M_B / (R_B + R_core)

    V_A += acc_A * dt
    R_A += V_A * dt

    if R_A < R_core:
        R_A = R_core
        V_A = abs(V_A) * 0.15

    V_B += acc_B * dt
    R_B += V_B * dt

    H_B = V_B / R_B if R_B > 0 else 0.0

    if M_B > 1e-8 and H_B > 0 and rhoB > 0:
        t_ff = 1.0 / np.sqrt(G * rhoB)
        t_H = 1.0 / H_B
        ratio = t_ff / t_H
    else:
        t_ff = np.inf
        t_H = np.inf
        ratio = np.inf

    rows.append([ti, M_B, R_B, rhoB, H_B, t_ff, t_H, ratio])

df = pd.DataFrame(rows, columns=["t","MB","RB","rhoB","HB","t_ff","t_H","ratio"])
df.to_csv("gds_structure_window.csv", index=False)

valid = df[np.isfinite(df["ratio"])]
allowed = valid[valid["ratio"] < 1.0]

print("=== GDS STRUCTURE WINDOW TEST ===")
print(f"valid samples                    = {len(valid)}")
print(f"samples with t_ff < t_H          = {len(allowed)}")

if len(allowed) > 0:
    print(f"first structure-allowed time     = {allowed['t'].iloc[0]:.6e}")
    print(f"last structure-allowed time      = {allowed['t'].iloc[-1]:.6e}")
    print(f"minimum t_ff/t_H                 = {valid['ratio'].min():.6e}")
    best = valid.loc[valid["ratio"].idxmin()]
    print()
    print("BEST STRUCTURE WINDOW POINT")
    print(best.to_string())
    print()
    print("RESULT:")
    print("There exists an early daughter-region window where local gravitational structure formation is dynamically allowed.")
else:
    print(f"minimum t_ff/t_H                 = {valid['ratio'].min():.6e}")
    print()
    print("RESULT:")
    print("No homogeneous structure-formation window found.")
    print("Next test: local overdensity seeds inside B.")

plt.figure(figsize=(9,5))
plt.plot(valid["t"], valid["ratio"])
plt.axhline(1.0, linestyle="--", label="structure threshold t_ff/t_H = 1")
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("t_ff / t_H")
plt.title("Daughter Region Structure Formation Window")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("structure_window_ratio.png", dpi=200)

print()
print("SAVED:")
print("  gds_structure_window.csv")
print("  structure_window_ratio.png")
