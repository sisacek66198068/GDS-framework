import numpy as np

def run(transfer_strength):
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

    R_A, V_A = 3.0, -0.15
    R_B, V_B = 0.05, 0.0
    M_A, M_B = M0, 0.0

    def rho_eff(M, R):
        return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

    def Gamma(rho):
        return 1.0 + A / (1.0 + (rho / rho_c)**n)

    RA_min = 999
    RB_max = R_B
    MB_max = 0.0

    for ti in t:
        rhoA = rho_eff(M_A, R_A)
        gamA = Gamma(rhoA)

        grav_A = -G * M_A / (R_A**2 + R_core**2)
        sat = (rhoA / rho_c)**n / (1.0 + (rhoA / rho_c)**n)
        bounce = k_bounce * sat / (R_A + R_core)
        gate_R = np.exp(-(R_A / (3*R_core))**2)

        acc_A = grav_A / gamA + bounce

        dM = transfer_strength * sat * gate_R * M_A * dt
        dM = min(dM, M_A)

        M_A -= dM
        M_B += dM

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

        RA_min = min(RA_min, R_A)
        RB_max = max(RB_max, R_B)
        MB_max = max(MB_max, M_B)

    return {
        "transfer_strength": transfer_strength,
        "RA_min": RA_min,
        "MA_final": M_A,
        "MB_final": M_B,
        "RB_final": R_B,
        "RB_max": RB_max,
        "VB_final": V_B,
        "mass_total": M_A + M_B,
    }

off = run(0.0)
on  = run(0.08)

print("=== CONTROL TEST: INTERLINK OFF VS ON ===")
print()
print("TRANSFER OFF:")
for k,v in off.items():
    print(f"{k:20s} = {v:.12e}" if isinstance(v,float) else f"{k:20s} = {v}")

print()
print("TRANSFER ON:")
for k,v in on.items():
    print(f"{k:20s} = {v:.12e}" if isinstance(v,float) else f"{k:20s} = {v}")

print()
print("=== DIFFERENCE ===")
print(f"Delta_MB_final = {on['MB_final'] - off['MB_final']:.12e}")
print(f"Delta_RB_final = {on['RB_final'] - off['RB_final']:.12e}")

print()
if off["MB_final"] == 0 and on["MB_final"] > 1e-6 and on["RB_final"] > off["RB_final"]:
    print("RESULT: Expansion of region B requires active interlink transfer.")
    print("Within this model, removed mass from A is quantitatively mapped into B.")
else:
    print("RESULT: Control test is inconclusive.")
