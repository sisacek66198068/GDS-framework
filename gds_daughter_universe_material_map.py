import numpy as np
import matplotlib.pyplot as plt

# ============================================================
# GDS DAUGHTER UNIVERSE MATERIAL MAP
#
# Goal:
# Track what happens to transferred matter inside region B:
#
# - density evolution
# - dilution during expansion
# - effective Hubble evolution
# - whether structure formation / stabilization is possible
#
# ============================================================

dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

G = 1.0
M0 = 1.0

rho_c = 50.0
n = 3.0
A = 4.0

R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

# Parent region A
R_A = 3.0
V_A = -0.15

# Daughter region B
R_B = 0.05
V_B = 0.0

M_A = M0
M_B = 0.0

hist = []

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

for ti in t:

    # ----------------------------
    # Parent collapse
    # ----------------------------

    rhoA = rho_eff(M_A, R_A)
    gamA = Gamma(rhoA)

    grav_A = -G * M_A / (R_A**2 + R_core**2)

    sat = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)

    bounce = k_bounce * sat / (R_A + R_core)

    acc_A = grav_A / gamA + bounce

    gate_R = np.exp(-(R_A / (3 * R_core))**2)

    dM = transfer_strength * sat * gate_R * M_A * dt
    dM = min(dM, M_A)

    M_A -= dM
    M_B += dM

    # ----------------------------
    # Daughter region B
    # ----------------------------

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

    # FRW-like Hubble parameter
    H_B = V_B / R_B if R_B > 0 else 0.0

    # Toy pressure / temperature proxy
    # expansion should cool / dilute
    pressure_B = rhoB * H_B

    # structure possibility proxy:
    # if density drops but remains finite and H stabilizes
    # matter can potentially reorganize instead of escaping forever
    structure_index = rhoB / (1.0 + abs(H_B))

    hist.append([
        ti,
        M_B,
        R_B,
        V_B,
        H_B,
        rhoB,
        pressure_B,
        structure_index
    ])

hist = np.array(hist)

ti   = hist[:,0]
MB   = hist[:,1]
RB   = hist[:,2]
VB   = hist[:,3]
HB   = hist[:,4]
rhoB = hist[:,5]
PB   = hist[:,6]
SI   = hist[:,7]

late = int(0.8 * len(ti))

rhoB_early = np.max(rhoB)
rhoB_late_mean = np.mean(rhoB[late:])
HB_late_mean = np.mean(HB[late:])
HB_late_std  = np.std(HB[late:])
SI_late_mean = np.mean(SI[late:])

print("=== GDS DAUGHTER UNIVERSE MATERIAL MAP ===")
print()

print("FINAL STATE")
print(f"M_B_final                         = {MB[-1]:.12e}")
print(f"R_B_final                         = {RB[-1]:.12e}")
print(f"V_B_final                         = {VB[-1]:.12e}")
print()

print("DENSITY EVOLUTION")
print(f"rho_B_peak                        = {rhoB_early:.12e}")
print(f"rho_B_late_mean                   = {rhoB_late_mean:.12e}")
print()

print("FRW EXPANSION")
print(f"H_B_late_mean                     = {HB_late_mean:.12e}")
print(f"H_B_late_std                      = {HB_late_std:.12e}")
print()

print("STRUCTURE POSSIBILITY")
print(f"structure_index_late_mean         = {SI_late_mean:.12e}")
print()

if rhoB_late_mean > 0 and HB_late_mean > 0 and HB_late_std < HB_late_mean:
    print("RESULT:")
    print("Transferred matter does not disappear.")
    print("It forms a persistent expanding medium in region B.")
    print("Dynamics allows long-term daughter-universe material evolution.")
else:
    print("RESULT:")
    print("Material evolution is less stable; stronger model needed.")

# ----------------------------
# PLOTS
# ----------------------------

plt.figure(figsize=(9,5))
plt.plot(ti, rhoB)
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("rho_B")
plt.title("Density Evolution of Matter in Daughter Region")
plt.grid(True)
plt.tight_layout()
plt.savefig("daughter_map_01_density.png", dpi=200)

plt.figure(figsize=(9,5))
plt.plot(ti, HB)
plt.xlabel("time")
plt.ylabel("H_B")
plt.title("FRW-like Expansion Rate in Daughter Region")
plt.grid(True)
plt.tight_layout()
plt.savefig("daughter_map_02_HB.png", dpi=200)

plt.figure(figsize=(9,5))
plt.plot(ti, SI)
plt.yscale("log")
plt.xlabel("time")
plt.ylabel("structure index")
plt.title("Material Reorganization Potential in Daughter Region")
plt.grid(True)
plt.tight_layout()
plt.savefig("daughter_map_03_structure.png", dpi=200)

print()
print("SAVED:")
print("  daughter_map_01_density.png")
print("  daughter_map_02_HB.png")
print("  daughter_map_03_structure.png")
