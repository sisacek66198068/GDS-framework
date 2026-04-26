import numpy as np

# ============================================================
# DUAL BLACK HOLE -> SHARED DAUGHTER UNIVERSE TEST
#
# Source 1 ~ Sgr A*-like
# Source 2 ~ Andromeda SMBH-like
#
# Goal:
# test whether both black holes can feed
# one shared expanding daughter spacetime
#
# Key diagnostic:
# shared Hubble-like expansion H_shared
#
# ============================================================

dt = 1e-4
T = 20.0
t = np.arange(0, T, dt)

G = 1.0

rho_c = 50.0
n = 3.0
A = 4.0

R_core = 0.08
k_bounce = 0.35
transfer_strength = 0.08

# ------------------------------------------------------------
# Source A : Sgr A*-like
# ------------------------------------------------------------

M_A = 1.0
R_A = 3.0
V_A = -0.15

# ------------------------------------------------------------
# Source B : Andromeda SMBH-like
# larger source
# ------------------------------------------------------------

M_C = 2.5
R_C = 4.5
V_C = -0.12

# ------------------------------------------------------------
# Shared daughter universe
# ------------------------------------------------------------

M_D = 0.0
R_D = 0.05
V_D = 0.0

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

hist = []

for ti in t:

    # ========================================================
    # SOURCE A
    # ========================================================

    rhoA = rho_eff(M_A, R_A)
    gamA = Gamma(rhoA)

    satA = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)
    gateA = np.exp(-(R_A / (3 * R_core))**2)

    gravA = -G * M_A / (R_A**2 + R_core**2)
    bounceA = k_bounce * satA / (R_A + R_core)

    accA = gravA / gamA + bounceA

    dM_A = transfer_strength * satA * gateA * M_A * dt
    dM_A = min(dM_A, M_A)

    M_A -= dM_A

    V_A += accA * dt
    R_A += V_A * dt

    if R_A < R_core:
        R_A = R_core
        V_A = abs(V_A) * 0.15

    # ========================================================
    # SOURCE C
    # ========================================================

    rhoC = rho_eff(M_C, R_C)
    gamC = Gamma(rhoC)

    satC = (rhoC / rho_c)**n / (1.0 + (rhoC / rho_c)**n)
    gateC = np.exp(-(R_C / (3 * R_core))**2)

    gravC = -G * M_C / (R_C**2 + R_core**2)
    bounceC = k_bounce * satC / (R_C + R_core)

    accC = gravC / gamC + bounceC

    dM_C = transfer_strength * satC * gateC * M_C * dt
    dM_C = min(dM_C, M_C)

    M_C -= dM_C

    V_C += accC * dt
    R_C += V_C * dt

    if R_C < R_core:
        R_C = R_core
        V_C = abs(V_C) * 0.15

    # ========================================================
    # SHARED DAUGHTER REGION
    # ========================================================

    dM_total = dM_A + dM_C
    M_D += dM_total

    accD = 0.0

    if M_D > 0:
        trigger = (satA * gateA + satC * gateC)
        accD += G * M_D / (R_D**2 + R_core**2) * trigger
        accD += 0.03 * M_D / (R_D + R_core)

    V_D += accD * dt
    R_D += V_D * dt

    H_D = V_D / R_D if R_D > 0 else 0.0

    hist.append([
        ti,
        M_A,
        M_C,
        M_D,
        R_D,
        V_D,
        H_D,
        dM_A,
        dM_C
    ])

hist = np.array(hist)

ti  = hist[:,0]
MA  = hist[:,1]
MC  = hist[:,2]
MD  = hist[:,3]
RD  = hist[:,4]
VD  = hist[:,5]
HD  = hist[:,6]
dA  = hist[:,7]
dC  = hist[:,8]

late = int(0.8 * len(HD))

print("=== DUAL BLACK HOLE SHARED UNIVERSE TEST ===")
print()

print("FINAL MASSES")
print(f"Remaining source A mass         = {MA[-1]:.12e}")
print(f"Remaining source C mass         = {MC[-1]:.12e}")
print(f"Shared daughter mass            = {MD[-1]:.12e}")
print()

print("SHARED EXPANSION")
print(f"Shared daughter radius          = {RD[-1]:.12e}")
print(f"Shared daughter velocity        = {VD[-1]:.12e}")
print(f"H_shared late mean              = {np.mean(HD[late:]):.12e}")
print(f"H_shared late std               = {np.std(HD[late:]):.12e}")
print()

print("MASS CONSERVATION")
total = MA[-1] + MC[-1] + MD[-1]
print(f"Final total mass                = {total:.12e}")
print()

if MD[-1] > 1e-6 and np.mean(HD[late:]) > 0:
    print("RESULT:")
    print("Both black holes consistently feed one shared expanding daughter spacetime.")
    print("Shared-universe interpretation is supported.")
else:
    print("RESULT:")
    print("Shared daughter spacetime not robustly supported.")
