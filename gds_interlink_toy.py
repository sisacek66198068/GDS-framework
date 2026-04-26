import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# GDS INTERLINK TOY MODEL
# Black-hole-like collapse -> saturation core -> white-hole-like expansion
# This is NOT proof of another universe.
# It tests whether a GDS-like saturating field can replace singularity
# by a finite transition channel.
# ============================================================

# Time
dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

# Constants / toy units
G = 1.0
M0 = 1.0

# GDS saturation parameters
rho_c = 50.0
n = 3.0
A = 4.0

# Core regularization
R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

# Initial universe A: collapsing shell
R_A = 3.0
V_A = -0.15

# Universe B: initially tiny expanding seed
R_B = 0.05
V_B = 0.0

# Mass bookkeeping
M_A = M0
M_B = 0.0

hist = []

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    # high density -> saturation / reduced runaway
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

for ti in t:
    rhoA = rho_eff(M_A, R_A)
    gamA = Gamma(rhoA)

    # Collapse acceleration in A, softened by core and saturated by Gamma
    grav_A = -G * M_A / (R_A**2 + R_core**2)

    # Bounce term activates at high density
    sat = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)
    bounce = k_bounce * sat / (R_A + R_core)

    # Net acceleration of collapsing side
    acc_A = grav_A / gamA + bounce

    # Transfer through saturated core
    # only activates when density is high and R_A near core
    gate_R = np.exp(-(R_A / (3*R_core))**2)
    dM = transfer_strength * sat * gate_R * M_A * dt
    dM = min(dM, M_A)

    M_A -= dM
    M_B += dM

    # Universe B receives expansion impulse from transferred mass/field
    rhoB = rho_eff(M_B + 1e-12, R_B)
    gamB = Gamma(rhoB)

    acc_B = 0.0
    if M_B > 0:
        acc_B += +G * M_B / (R_B**2 + R_core**2) * (sat * gate_R)
        acc_B += +0.03 * M_B / (R_B + R_core)

    # Integrate
    V_A += acc_A * dt
    R_A += V_A * dt

    # prevent numerical negative radius; interpret as bounce core
    if R_A < R_core:
        R_A = R_core
        V_A = abs(V_A) * 0.15

    V_B += acc_B * dt
    R_B += V_B * dt

    hist.append([ti, R_A, V_A, M_A, rhoA, gamA, sat, M_B, R_B, V_B])

hist = np.array(hist)
ti, RA, VA, MA, rhoA, gamA, sat, MB, RB, VB = hist.T

# Metrics
print("=== GDS INTERLINK TOY RESULT ===")
print(f"R_A_min              = {RA.min():.6e}")
print(f"rho_A_max            = {rhoA.max():.6e}")
print(f"Gamma_A_min/max      = {gamA.min():.6e} / {gamA.max():.6e}")
print(f"M_A_final            = {MA[-1]:.6e}")
print(f"M_B_final            = {MB[-1]:.6e}")
print(f"R_B_final            = {RB[-1]:.6e}")
print(f"V_B_final            = {VB[-1]:.6e}")
print(f"mass_conservation    = {(MA[-1]+MB[-1]):.12e}")
print()
if RA.min() > 0 and MB[-1] > 1e-6 and RB[-1] > RB[0]:
    print("INTERPRETATION: finite core + transfer + expansion channel appears.")
    print("Toy model supports a regularized bounce/interlink scenario.")
else:
    print("INTERPRETATION: no robust interlink in this parameter set.")

# Plots
plt.figure(figsize=(8,5))
plt.plot(ti, RA, label="Region A radius (collapse side)")
plt.plot(ti, RB, label="Region B radius (expansion side)")
plt.xlabel("time [toy units]")
plt.ylabel("radius [toy units]")
plt.title("GDS Toy Interlink: Collapse Core and Expansion Channel")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("gds_interlink_radius.png", dpi=200)

plt.figure(figsize=(8,5))
plt.plot(ti, MA, label="Mass in region A")
plt.plot(ti, MB, label="Mass in region B")
plt.plot(ti, MA+MB, "--", label="Total mass")
plt.xlabel("time [toy units]")
plt.ylabel("mass [toy units]")
plt.title("Mass Transfer Through Saturated GDS Core")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("gds_interlink_mass.png", dpi=200)

plt.figure(figsize=(8,5))
plt.plot(ti, rhoA, label="Effective density in A")
plt.axhline(rho_c, linestyle="--", label="critical density")
plt.xlabel("time [toy units]")
plt.ylabel("density [toy units]")
plt.yscale("log")
plt.title("Saturation Trigger in the Collapse Region")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("gds_interlink_density.png", dpi=200)

plt.figure(figsize=(8,5))
plt.plot(ti, sat, label="GDS saturation factor")
plt.xlabel("time [toy units]")
plt.ylabel("saturation")
plt.title("Activation of GDS Saturation / Interlink Gate")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("gds_interlink_saturation.png", dpi=200)

print("SAVED:")
print("  gds_interlink_radius.png")
print("  gds_interlink_mass.png")
print("  gds_interlink_density.png")
print("  gds_interlink_saturation.png")
