import numpy as np
import pandas as pd

def run_case(rho_c=50.0, n=3.0, A=4.0, R_core=0.08, k_bounce=0.35, transfer_strength=0.08):
    dt = 1e-4
    T = 20.0
    t = np.arange(0, T, dt)

    G = 1.0
    M0 = 1.0

    R_A, V_A = 3.0, -0.15
    R_B, V_B = 0.05, 0.0
    M_A, M_B = M0, 0.0

    RA_min = 1e99
    rhoA_max = 0.0
    gam_min = 1e99
    gam_max = -1e99
    sat_max = 0.0

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

        acc_A = grav_A / gamA + bounce

        gate_R = np.exp(-(R_A / (3*R_core))**2)
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
        rhoA_max = max(rhoA_max, rhoA)
        gam_min = min(gam_min, gamA)
        gam_max = max(gam_max, gamA)
        sat_max = max(sat_max, sat)

    success = (RA_min > 0.0) and (M_B > 1e-6) and (R_B > 0.05) and abs((M_A+M_B)-1.0) < 1e-10

    return {
        "rho_c": rho_c,
        "n": n,
        "A": A,
        "R_core": R_core,
        "k_bounce": k_bounce,
        "transfer_strength": transfer_strength,
        "RA_min": RA_min,
        "rhoA_max": rhoA_max,
        "Gamma_min": gam_min,
        "Gamma_max": gam_max,
        "sat_max": sat_max,
        "MA_final": M_A,
        "MB_final": M_B,
        "RB_final": R_B,
        "VB_final": V_B,
        "mass_total": M_A + M_B,
        "success": success,
    }

rows = []

for rho_c in [20, 50, 100, 200]:
    for n in [2, 3, 4, 6]:
        for A in [1, 2, 4, 8]:
            for k_bounce in [0.15, 0.25, 0.35, 0.50]:
                for transfer_strength in [0.02, 0.05, 0.08, 0.12]:
                    rows.append(run_case(
                        rho_c=rho_c,
                        n=n,
                        A=A,
                        k_bounce=k_bounce,
                        transfer_strength=transfer_strength
                    ))

df = pd.DataFrame(rows)
df.to_csv("gds_interlink_scan_results.csv", index=False)

print("=== GDS INTERLINK PARAMETER SCAN ===")
print("total cases:", len(df))
print("successful:", int(df["success"].sum()))
print("success rate:", float(df["success"].mean()))

print()
print("=== TOP 20 BY MB_final ===")
print(df.sort_values("MB_final", ascending=False).head(20).to_string(index=False))

print()
print("=== TOP 20 BY RB_final ===")
print(df.sort_values("RB_final", ascending=False).head(20).to_string(index=False))

print()
print("SAVED: gds_interlink_scan_results.csv")
