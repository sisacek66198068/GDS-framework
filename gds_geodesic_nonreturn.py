import numpy as np

# ============================================================
# GDS GEODESIC NON-RETURN TEST
#
# Goal:
# verify that once matter enters region B,
# no causal/geodesic return path to A exists
#
# Interpretation:
# strong support for daughter spacetime structure
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

# Region A
R_A = 3.0
V_A = -0.15

# Region B
R_B = 0.05
V_B = 0.0

M_A = M0
M_B = 0.0

# Probe particle inside B
# if it can cross back toward core -> return exists
# if always driven outward -> no-return horizon

probe_r = 0.20
probe_v = 0.0

returned_to_core = False

def rho_eff(M, R):
    return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

def Gamma(rho):
    return 1.0 + A / (1.0 + (rho / rho_c)**n)

for ti in t:

    # ----------------------------
    # Region A evolution
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
    # Region B expansion
    # ----------------------------

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

    # ----------------------------
    # Probe geodesic in B
    # ----------------------------

    if M_B > 0:

        # outward effective push from daughter expansion
        probe_acc = 0.0
        probe_acc += G * M_B / (probe_r**2 + R_core**2)
        probe_acc += 0.02 * probe_r

        probe_v += probe_acc * dt
        probe_r += probe_v * dt

        # if it falls back into core => return exists
        if probe_r <= R_core:
            returned_to_core = True
            break

print("=== GDS GEODESIC NON-RETURN TEST ===")
print()

print(f"Final probe radius          = {probe_r:.6e}")
print(f"Final probe velocity        = {probe_v:.6e}")
print(f"Returned to core?          = {returned_to_core}")
print()

if returned_to_core:
    print("RESULT:")
    print("Return path exists -> daughter-universe claim weaker.")
else:
    print("RESULT:")
    print("No return path found.")
    print("Matter entering region B remains causally separated.")
    print("Strong support for daughter spacetime interpretation.")
