import numpy as np

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

def run_case(seed=0, eps=0.0):
    rng = np.random.default_rng(seed)

    R_A = 3.0 * (1.0 + eps * rng.normal())
    V_A = -0.15 * (1.0 + eps * rng.normal())
    R_B = 0.05 * (1.0 + eps * rng.normal())
    V_B = 0.0
    M_A = M0 * (1.0 + eps * rng.normal())
    M_B = 0.0

    def rho_eff(M, R):
        return M / ((4.0/3.0) * np.pi * (R**3 + R_core**3))

    def Gamma(rho):
        return 1.0 + A / (1.0 + (rho / rho_c)**n)

    RA_min = 1e99
    HB_hist = []

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

        RA_min = min(RA_min, R_A)
        HB_hist.append(V_B / R_B if R_B > 0 else 0.0)

    HB_hist = np.array(HB_hist)
    late = int(0.8 * len(HB_hist))

    return {
        "RA_min": RA_min,
        "MA_final": M_A,
        "MB_final": M_B,
        "RB_final": R_B,
        "HB_late_mean": HB_hist[late:].mean(),
        "mass_total": M_A + M_B,
        "success": (RA_min > 0 and M_B > 1e-6 and R_B > 0.05 and HB_hist[late:].mean() > 0)
    }

eps_values = [0.001, 0.005, 0.01, 0.02, 0.05]
N = 50

rows = []

for eps in eps_values:
    successes = 0
    MBs = []
    RBs = []
    HBs = []
    masses = []

    for seed in range(N):
        r = run_case(seed=seed, eps=eps)
        successes += int(r["success"])
        MBs.append(r["MB_final"])
        RBs.append(r["RB_final"])
        HBs.append(r["HB_late_mean"])
        masses.append(r["mass_total"])

    MBs = np.array(MBs)
    RBs = np.array(RBs)
    HBs = np.array(HBs)
    masses = np.array(masses)

    rows.append([
        eps,
        successes,
        successes / N,
        MBs.mean(),
        MBs.std(),
        RBs.mean(),
        RBs.std(),
        HBs.mean(),
        HBs.std(),
        abs(masses.mean() - 1.0),
        masses.std()
    ])

print("=== GDS PERTURBATION STABILITY TEST ===")
print("eps success/N rate MB_mean MB_std RB_mean RB_std HB_mean HB_std mass_err mass_std")
for r in rows:
    print(
        f"{r[0]:.3e} {int(r[1]):02d}/{N} {r[2]:.3f} "
        f"{r[3]:.6e} {r[4]:.6e} "
        f"{r[5]:.6e} {r[6]:.6e} "
        f"{r[7]:.6e} {r[8]:.6e} "
        f"{r[9]:.6e} {r[10]:.6e}"
    )

print()
if all(r[2] == 1.0 for r in rows):
    print("RESULT:")
    print("Daughter-universe mechanism is stable under tested perturbations.")
    print("It behaves as a robust dynamical attractor, not a fine-tuned trajectory.")
else:
    print("RESULT:")
    print("Some perturbation levels reduce success rate; map stability boundary.")
