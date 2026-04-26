import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# FIGURE 13 — GDS LIGO ECHO PREDICTION
# FIGURE 14 — GDS JWST STRUCTURE WINDOW
# ============================================================

# -----------------------------
# FIG 13: LIGO echo delay
# -----------------------------
M_solar = 1.98847e30
G = 6.67430e-11
c = 299792458.0

ratio = 25.0  # R_s / R_core from GDS toy model

masses = np.linspace(5, 200, 400)
echo_ms = []

for m in masses:
    M = m * M_solar
    R_s = 2.0 * G * M / c**2
    dt = 2.0 * R_s / c * np.log(ratio)
    echo_ms.append(dt * 1000.0)

echo_ms = np.array(echo_ms)

plt.figure(figsize=(9,5))
plt.plot(masses, echo_ms, linewidth=2)
plt.scatter([10,20,30,50,70,100,150],
            [0.634201,1.268402,1.902603,3.171005,4.439406,6.342009,9.513014],
            s=35,
            label="computed benchmark points")
plt.xlabel("Black-hole remnant mass [$M_\\odot$]")
plt.ylabel("Predicted echo delay [ms]")
plt.title("GDS Finite-Core Prediction: Ringdown Echo Delay")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("MUL/fig13_ligo_echo_prediction.png", dpi=200)

# -----------------------------
# FIG 14: JWST structure window
# -----------------------------
z_start = 11.270862
z_end = 10.820572

plt.figure(figsize=(9,5))
plt.axvspan(z_end, z_start, alpha=0.35, label="GDS structure window")
plt.axvline(8, linestyle="--", label="JWST high-z lower reference")
plt.axvline(15, linestyle="--", label="JWST high-z upper reference")
plt.text((z_start+z_end)/2, 0.55, "GDS window\nz ≈ 10.8–11.3",
         ha="center", va="center", fontsize=12)
plt.xlabel("Redshift z")
plt.ylabel("Structure-formation relevance [normalized]")
plt.title("GDS Early Structure Window in JWST Redshift Range")
plt.xlim(6, 16)
plt.ylim(0, 1)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig("MUL/fig14_jwst_structure_window.png", dpi=200)

print("SAVED:")
print("  MUL/fig13_ligo_echo_prediction.png")
print("  MUL/fig14_jwst_structure_window.png")
