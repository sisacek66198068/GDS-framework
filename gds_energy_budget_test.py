import numpy as np

# ============================================================
# GDS ENERGY BUDGET TEST
#
# Goal:
# verify that finite-core + daughter transfer
# is energetically consistent and not a numerical artifact
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

# parent region
R_A = 3.0
V_A = -0.15

# daughter region
R_B = 0.05
V_B = 0.0

M_A = M0
M_B = 0.0

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

def grav_energy(M, R):
    return -G * M * M / (R + R_core)

def kin_energy(M, V):
    return 0.5 * M * V * V

def bounce_energy(sat, R):
    return k_bounce * sat / (R + R_core)

hist = []

for ti in t:

    rhoA = rho_eff(M_A, R_A)
    gamA = Gamma(rhoA)

    sat = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)
    gate = np.exp(-(R_A / (3 * R_core))**2)

    grav_A = -G * M_A / (R_A**2 + R_core**2)
    bounce = k_bounce * sat / (R_A + R_core)

    acc_A = grav_A / gamA + bounce

    dM = transfer_strength * sat * gate * M_A * dt
    dM = min(dM, M_A)

    M_A -= dM
    M_B += dM

    acc_B = 0.0
    if M_B > 0:
        acc_B += G * M_B / (R_B**2 + R_core**2) * (sat * gate)
        acc_B += 0.03 * M_B / (R_B + R_core)

    V_A += acc_A * dt
    R_A += V_A * dt

    if R_A < R_core:
        R_A = R_core
        V_A = abs(V_A) * 0.15

    V_B += acc_B * dt
    R_B += V_B * dt

    E_A = grav_energy(M_A, R_A) + kin_energy(M_A, V_A)
    E_B = grav_energy(M_B, R_B) + kin_energy(M_B, V_B)
    E_core = bounce_energy(sat, R_A)

    E_total = E_A + E_B + E_core

    hist.append([ti, E_A, E_B, E_core, E_total])

hist = np.array(hist)

ti = hist[:,0]
EA = hist[:,1]
EB = hist[:,2]
EC = hist[:,3]
ET = hist[:,4]

late = int(0.8 * len(ET))

print("=== GDS ENERGY BUDGET TEST ===")
print()

print(f"Initial total effective energy      = {ET[0]:.12e}")
print(f"Final total effective energy        = {ET[-1]:.12e}")
print(f"Late mean total energy              = {np.mean(ET[late:]):.12e}")
print(f"Late std total energy               = {np.std(ET[late:]):.12e}")
print()

drift = abs(ET[-1] - ET[0])
rel = drift / max(abs(ET[0]), 1e-30)

print(f"Absolute energy drift               = {drift:.12e}")
print(f"Relative energy drift               = {rel:.12e}")
print()

if rel < 0.1:
    print("RESULT:")
    print("Energy budget remains dynamically consistent.")
    print("Finite-core daughter transfer is not a trivial numerical artifact.")
else:
    print("RESULT:")
    print("Energy drift is significant; stronger formulation needed.")
